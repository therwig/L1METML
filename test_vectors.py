import ROOT
from collections import OrderedDict
import numpy as np
from array import array

# open files
f1 = ROOT.TFile.Open("samples/perfNano_DY_reprocessed.root")
t = f1.Get("Events")
t.AddFriend("Events","samples/output_MET_DY.root")

ROOT.gROOT.LoadMacro("resultWRTz.C")
# t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,1):zPt>>"+hists["pred_para"].GetName())
# t.Draw("resultWRTz(pred_para_puppi,pred_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,0):zPt>>"+hists["pred_perp"].GetName())
#xcheck
t.Draw("test_resultWRTz(genMet_para_puppi,genMet_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,1,genMet,genMetPhi,genMet_para_z,genMet_perp_z):zPt>>h1","Entry$<1")
#t.Draw("test_resultWRTz(genMet_para_puppi,genMet_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,0,genMet,genMetPhi):zPt>>h2","Entry$<1")

f1.Close()
