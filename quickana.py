import ROOT
from collections import OrderedDict
import numpy as np
from array import array
from tools import *

# open files
f1 = ROOT.TFile.Open("samples/perfNano_DY_reprocessed.root")
t = f1.Get("Events")
t.AddFriend("Events","samples/output_MET_DY.root")

# book histograms
hists=OrderedDict()
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
    # for ik, kind in enumerate(["train","valid","test"]):
    #     book2(hists,m+"_para_"+kind,zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
    #     book2(hists,m+"_perp_"+kind,zpt_n,0,zpt_hi,perp_n,perp_lo,perp_hi)
    #     t.Draw(m+"_para_z:zPt>>"+hists[m+"_para_"+kind].GetName(),"type={}".format(ik))
    #     t.Draw(m+"_perp_z:zPt>>"+hists[m+"_perp_"+kind].GetName(),"type={}".format(ik))

# for prediction(s), must re-calculate MVA MET wrt the z boson direction by hand
ROOT.gROOT.LoadMacro("resultWRTz.C")
book2(hists,"pred_para_all",zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
book2(hists,"pred_perp_all",zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,1):zPt>>"+hists["pred_para_all"].GetName())
t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,0):zPt>>"+hists["pred_perp_all"].GetName())
for ik, kind in enumerate(["train","valid","test"]):
    book2(hists,"pred_para_"+kind,zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
    book2(hists,"pred_perp_"+kind,zpt_n,0,zpt_hi,para_n,para_lo,para_hi)
    t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,1):zPt>>"+hists["pred_para_"+kind].GetName(),"type=={}".format(ik))
    t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,0):zPt>>"+hists["pred_perp_"+kind].GetName(),"type=={}".format(ik))

if True: #store hists
    fout = ROOT.TFile("tmp/histograms.root","recreate")
    for n in hists:
        hists[n].Write()
    fout.Close()


graphs=OrderedDict()
hz=hists["ptz_2d"]

#resp, reso = GraphsFromHist(hists["L1PFMet_para"], hz)

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
# sel = lambda x : (not "pred" in x or "all" in x) # "all pred" only
# sel = lambda x : (not "pred" in x or "test" in x) # "test pred" only
sel = lambda x : True # "all types"
para_resp = {g:graphs[g] for g in graphs if ("para" in g and "resp" in g and sel(g))}
para_reso = {g:graphs[g] for g in graphs if ("para" in g and "reso" in g and sel(g))}
perp_resp = {g:graphs[g] for g in graphs if ("perp" in g and "resp" in g and sel(g))}
perp_reso = {g:graphs[g] for g in graphs if ("perp" in g and "reso" in g and sel(g))}

for n in graphs:
    Plot(graphs[n],"plots/test_"+n)

MultiPlot(para_resp, outname="plots/para_resp")
MultiPlot(para_reso, outname="plots/para_reso")
MultiPlot(perp_resp, outname="plots/perp_resp")
MultiPlot(perp_reso, outname="plots/perp_reso")



