import ROOT
from collections import OrderedDict
import numpy as np
from array import array

def book(h,name,n,a,b,title=""):
    h[name]=ROOT.TH1F("h_"+name,title,n,a,b)
    h[name].Sumw2()
def book2(h,name,nx,ax,bx,ny,ay,by,title=""):
    h[name]=ROOT.TH2F("h2_"+name,title,nx,ax,bx,ny,ay,by)
    h[name].Sumw2()

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
            suf="(" if si==0 else (")" if si==nx-1 else "")
            s.Draw()
            c.SaveAs("plots/debug/fits_{}.pdf{}".format(h.GetName(),suf))
        fdebug.Close()
    return g_resp,g_reso

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

def GetLegEntry(x):
    stem = x.split("_")[0]
    if "pred" in x:
        if "train" in x: stem += "_train"
        if "valid" in x: stem += "_valid"
        if "test" in x: stem += "_test"
        if "all" in x: stem += "_all"
    return stem

def MultiPlot(gs, outname="dummy", 
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
        leg.AddEntry(gs[x], GetLegEntry(x))
        #leg.AddEntry(gs[x],x)
    mg.Draw("apl PLC PMC")
    mg.SetTitle(";q_{T} [GeV];"+GetYTitle(outname))

    # pad2.cd()
    leg.SetFillStyle(0);
    leg.SetFillColor(0);
    leg.SetBorderSize(0);
    leg.Draw();

    c.SaveAs(outname+".pdf")    

