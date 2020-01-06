import ROOT
from collections import OrderedDict

# open files
f1 = ROOT.TFile.Open("samples/perfNano_DY_reprocessed.root")
t = f1.Get("Events")
t.AddFriend("Events","samples/output_MET_DY.root")

# book histograms
hists=OrderedDict()
def book(h,name,n,a,b,title=""):
    h[name]=ROOT.TH1F(name,title,n,a,b)
    h[name].Sumw2()
def book2(h,name,nx,ax,bx,ny,ay,by,title=""):
    h[name]=ROOT.TH2F(name,title,nx,ax,bx,ny,ay,by)
    h[name].Sumw2()

zpt_n=10
zpt_hi=200
book(hists,"ptz",2*zpt_n,0,zpt_hi)
book2(hists,"ptz_2d"    ,zpt_n,0,zpt_hi,8*zpt_n,0,zpt_hi)
t.Draw("zPt>>ptz")
t.Draw("zPt:zPt>>ptz_2d")

met_types=["genMet","L1CaloMet","L1PFMet","L1PuppiMet",
           "CaloHTMiss","PFHTMiss","PuppiHTMiss",
           "pred"]
para_n, para_lo, para_hi = (160,-300,500)
perp_n, perp_lo, perp_hi = (160,-400,400)
for m in met_types:
    book2(hists,m+"_para",zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
    book2(hists,m+"_perp",zpt_n,0,zpt_hi,perp_n,perp_lo,perp_hi)
    t.Draw(m+"_para_puppi:zPt>>"+m+"_para")
    t.Draw(m+"_perp_puppi:zPt>>"+m+"_perp")

if True: #store hists
    fout = ROOT.TFile("tmp/histograms.root","recreate")
    for n in hists:
        hists[n].Write()
    fout.Close()

# proceed to fit the hists from here
