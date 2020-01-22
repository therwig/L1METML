import ROOT
from collections import OrderedDict
import numpy as np
from array import array

# open files
f1 = ROOT.TFile.Open("samples/perfNano_DY_reprocessed.root")
t = f1.Get("Events")
t.AddFriend("Events","samples/output_MET_DY.root")

# book histograms
hists=OrderedDict()
def book(h,name,n,a,b,title=""):
    h[name]=ROOT.TH1F("h_"+name,title,n,a,b)
    h[name].Sumw2()
def book2(h,name,nx,ax,bx,ny,ay,by,title=""):
    h[name]=ROOT.TH2F("h2_"+name,title,nx,ax,bx,ny,ay,by)
    h[name].Sumw2()

zpt_n=10
zpt_hi=200
book(hists,"ptz",2*zpt_n,0,zpt_hi)
book2(hists,"ptz_2d"    ,zpt_n,0,zpt_hi,8*zpt_n,0,zpt_hi)
t.Draw("zPt>>"+hists["ptz"].GetName())
t.Draw("zPt:zPt>>"+hists["ptz_2d"].GetName())

met_types=["genMet","L1CaloMet","L1PFMet","L1PuppiMet",
#          "CaloHTMiss","PFHTMiss","PuppiHTMiss",
#           "pred"
]
para_n, para_lo, para_hi = (160,-300,500)
perp_n, perp_lo, perp_hi = (160,-400,400)
for m in met_types:
    book2(hists,m+"_para",zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
    book2(hists,m+"_perp",zpt_n,0,zpt_hi,perp_n,perp_lo,perp_hi)
    t.Draw(m+"_para_z:zPt>>"+hists[m+"_para"].GetName())
    t.Draw(m+"_perp_z:zPt>>"+hists[m+"_perp"].GetName())

# for prediction, must recalculate MVA MET wrt the z boson direction
book2(hists,"pred_para",zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
book2(hists,"pred_perp",zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
ROOT.gROOT.LoadMacro("resultWRTz.C")
t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,1):zPt>>"+hists["pred_para"].GetName())
t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,0):zPt>>"+hists["pred_perp"].GetName())
# xcheck that the genMet results are reproduced
# t.Draw("resultWRTz(genMet_para_puppi,genMet_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,1):zPt>>"+hists["pred_para"].GetName())
# t.Draw("resultWRTz(genMet_para_puppi,genMet_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,0):zPt>>"+hists["pred_perp"].GetName())

if True: #store hists
    fout = ROOT.TFile("tmp/histograms.root","recreate")
    for n in hists:
        hists[n].Write()
    fout.Close()

def DivErrs(n,ne,d,de):
    resp  = np.divide(n, d, out=np.zeros_like(d), where=d!=0)
    fe1=np.divide(ne, n, out=np.zeros_like(n), where=n!=0)
    fe2=np.divide(de, d, out=np.zeros_like(d), where=d!=0)
    respe = resp * np.sqrt(fe1**2+fe2**2)
    return resp,respe

# extract response/resolution graphs from hist
def GraphsFromHist(h,hz,debug=True):
    nx=hz.GetNbinsX()
    z     = np.zeros(nx)
    ze    = np.zeros(nx)
    zeup  = np.zeros(nx)
    zedn  = np.zeros(nx)
    met   = np.zeros(nx) 
    mete  = np.zeros(nx)
    reso  = np.zeros(nx) 
    resoe = np.zeros(nx)

    met_slices=[]
    for i in range(nx):
        bi=i+1 # bin indexing
        lo=hz.GetXaxis().GetBinLowEdge(bi)
        hi=hz.GetXaxis().GetBinUpEdge(bi)
        z_slice=hz.ProjectionY("{}_zpt_{}zpt{}".format(hz.GetName(),int(lo),int(hi)),bi,bi,"e")
        met_slice=h.ProjectionY("{}_{}zpt{}".format(h.GetName(),int(lo),int(hi)),bi,bi,"e")
        met_slices.append(met_slice)

        z[i]    = z_slice.GetMean()
        ze[i]   = z_slice.GetRMS()/np.sqrt(z_slice.GetEntries()) if z_slice.GetEntries() else 0.
        zeup[i] = hi-z_slice.GetMean() # asym errors
        zedn[i] = z_slice.GetMean()-lo # asym errors
        met[i]  = met_slice.GetMean()
        mete[i] = met_slice.GetRMS()/np.sqrt(met_slice.GetEntries()) if met_slice.GetEntries() else 0.
        reso[i] = met_slice.GetRMS()

    # calculate response
    resp, respe = DivErrs(met,mete,z,ze)

    g_resp = ROOT.TGraphAsymmErrors(len(z),array('f',z),array('f',resp),
                                    array('f',zedn), array('f',zeup),
                                    array('f',respe),array('f',respe))
    g_reso = ROOT.TGraphAsymmErrors(len(z),array('f',z),array('f',reso),
                                    array('f',zedn), array('f',zeup),
                                    array('f',resoe),array('f',resoe))
    name="_".join( h.GetName().split('_')[1:] )
    g_resp.SetName("g_{}_resp".format(name))
    g_reso.SetName("g_{}_reso".format(name))

    if debug:
        fdebug=ROOT.TFile("tmp/debug_response_plots.root","recreate")
        ROOT.gStyle.SetOptStat(1111)
        c = ROOT.TCanvas("canv","",750,750)
        suf = ""
        for si,s in enumerate(met_slices):
            s.Write()
            suf="(" if si==1 else (")" if si==nx else "")
            s.Draw()
            c.SaveAs("plots/debug/fits_{}.pdf{}".format(h.GetName(),suf))
        fdebug.Close()
    return g_resp,g_reso

graphs=OrderedDict()
hz=hists["ptz_2d"]

resp, reso = GraphsFromHist(hists["L1PFMet_para"], hz)

for hname in hists:
    h=hists[hname]
    if not (("para" in hname) or ("perp" in hname)):
        continue
    resp, reso = GraphsFromHist(h,hz)
    graphs[hname+"_resp"]=resp
    graphs[hname+"_reso"]=reso

if True: #store graphs
    fout = ROOT.TFile("tmp/graphs.root","recreate")
    for n in graphs:
        graphs[n].Write()
    fout.Close()

# produce summary plots
para_resp = {g:graphs[g] for g in graphs if ("para" in g and "resp" in g)}
para_reso = {g:graphs[g] for g in graphs if ("para" in g and "reso" in g)}
perp_resp = {g:graphs[g] for g in graphs if ("perp" in g and "resp" in g)}
perp_reso = {g:graphs[g] for g in graphs if ("perp" in g and "reso" in g)}

def GetYTitle(x):
    if ("para" in x and "resp" in x): return "<u_{\\parallel}>/<q_{T}>"
    if ("perp" in x and "resp" in x): return "<u_{\\perp}>/<q_{T}>"
    if ("para" in x and "reso" in x): return "\\sigma(u_{\\parallel})"
    if ("perp" in x and "reso" in x): return "\\sigma(u_{\\perp}  )"
    return ""
def Plot(g, outname, 
         ytitle="response", yrange=(0,1)):
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("canv","",750,750)
    pad = ROOT.TPad("pad", "pad", .005, .01, .995, .995)
    pad.Draw()
    pad.cd()

    ROOT.gPad.SetLeftMargin(0.17)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetTopMargin(0.05)

    g.Draw()
    g.SetTitle(";q_{T} [GeV];"+GetYTitle(outname))
    c.SaveAs(outname+".pdf")    
def MultiPlot(gs, outname, 
              ytitle="response", yrange=(0,1)):
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("canv","",750,750)
    pad = ROOT.TPad("pad", "pad", .005, .01, .995, .995)
    pad.Draw()
    # pad2 = ROOT.TPad("pad2", "pad2", .855, .01, .995, .995)
    # pad2.Draw()
    pad.cd()

    ROOT.gPad.SetLeftMargin(0.1)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetTopMargin(0.05)

    # pad2.cd()
    leg = ROOT.TLegend(.15,0.65,0.45,0.95);
    # leg = ROOT.TLegend(.05,0.05,0.88,0.88);
    leg.SetTextFont(42);
    leg.SetHeader("");
    leg.SetNColumns(1);
    # pad.cd()

    mg = ROOT.TMultiGraph()
    for x in gs:
        mg.Add(gs[x])
        gs[x].SetLineWidth(2)
        leg.AddEntry(gs[x],x.split("_")[0])
    mg.Draw("apl PLC PMC")
    mg.SetTitle(";q_{T} [GeV];"+GetYTitle(outname))

    # pad2.cd()
    leg.SetFillStyle(0);
    leg.SetFillColor(0);
    leg.SetBorderSize(0);
    leg.Draw();

    c.SaveAs(outname+".pdf")    

for n in graphs:
    Plot(graphs[n],"plots/test_"+n)

MultiPlot(para_resp,"plots/para_resp")
MultiPlot(para_reso,"plots/para_reso")
MultiPlot(perp_resp,"plots/perp_resp")
MultiPlot(perp_reso,"plots/perp_reso")
