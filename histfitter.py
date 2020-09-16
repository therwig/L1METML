#!/usr/bin/env python
from ROOT import *
from math import sqrt
from collections import OrderedDict
from array import array

import os
os.system("mkdir -p plots")

def plot(g,ytitle="response",yrange=(0,1)):
    gStyle.SetOptStat(0)
    c = TCanvas("canv","",750,750)
    pad = TPad("pad", "pad", .005, .01, .995, .995)
    pad.Draw()
    pad.cd()

    gPad.SetLeftMargin(0.17)
    gPad.SetRightMargin(0.05)
    gPad.SetBottomMargin(0.1)
    gPad.SetTopMargin(0.05)

    #draw this guy first
    b = TH1D("b","",1,0,g.GetX()[g.GetN()-1]+g.GetEXhigh()[g.GetN()-1])
    b.Draw()
    
    b.GetXaxis().SetTitle("q_{T} [GeV]")
    b.GetXaxis().SetTitleOffset(1.2)
    b.GetYaxis().SetTitleOffset(1.5)
    b.GetYaxis().SetLabelSize(0.04)
    b.GetYaxis().SetTitleSize(0.04)
    b.GetYaxis().SetTitle(ytitle)
    b.GetYaxis().SetRangeUser(*yrange)
    b.Draw("axis")

    g.Draw("same")
    c.SaveAs("plots/{}.pdf".format(g.GetName()))

def ArrayDivideErrs(n,ne,d,de):
    x=[]
    xe=[]
    for i in range(len(n)):
        if d[i]==0 or n[i]==0:
            x +=[0.]
            xe+=[0.]
        else:
            x+=[n[i]/d[i]]
            xe+=[n[i]/d[i]*sqrt((ne[i]/n[i])**2+(de[i]/d[i])**2)]
    return x,xe

def RunFits(hz,hpara,hperp,name,fout=None,simple=True):
    #name=h.GetName()
    nx=hz.GetNbinsX()
    z     =[]
    ze    =[]
    zeup  =[]
    zedn  =[]
    para  =[] 
    parae =[]
    resPara =[]
    resParae=[]
    resPerp =[]
    resPerpe=[]

    gStyle.SetOptStat(1111)
    c = TCanvas("canv","",750,750)

    for xi in range(1,nx+1):
        lo=hz.GetXaxis().GetBinLowEdge(xi)
        hi=hz.GetXaxis().GetBinUpEdge(xi)
        z_slice=hz.ProjectionY("{}_zpt_{}zpt{}".format(hz.GetName(),int(lo),int(hi)),xi,xi,"e")
        para_slice=hpara.ProjectionY("{}_para_{}zpt{}".format(hpara.GetName(),int(lo),int(hi)),xi,xi,"e")
        perp_slice=hperp.ProjectionY("{}_perp_{}zpt{}".format(hperp.GetName(),int(lo),int(hi)),xi,xi,"e")

        # fit
        if simple:
            z.append(z_slice.GetMean())
            ze.append(z_slice.GetRMS()/sqrt(z_slice.GetEntries()) if z_slice.GetEntries() else 0.)
            zeup.append(hi-z_slice.GetMean()) # asym errors
            zedn.append(z_slice.GetMean()-lo) # asym errors
            para.append(para_slice.GetMean())
            parae.append(para_slice.GetRMS()/sqrt(para_slice.GetEntries()) if para_slice.GetEntries() else 0.)
            resPara.append(para_slice.GetRMS())
            resPerp.append(perp_slice.GetRMS())
            resParae.append(0.)
            resPerpe.append(0.)
        # else:
        #     x.append((lo+hi)/2) 
        #     xe.append((hi-lo)/2)
        #     fit_lo = hslice.GetBinLowEdge(1)
        #     fit_hi = hslice.GetBinLowEdge(1+hslice.GetNbinsX())
        #     func = TF1("func","gaus",fit_lo,fit_hi)
        #     fit_result = hslice.Fit(func,"QS","",fit_lo,fit_hi)
        #     mean.append(fit_result.Parameter(1)) 
        #     meane.append(fit_result.ParError(1)) 
        #     width.append(fit_result.Parameter(2))
        #     widthe.append(fit_result.ParError(2))
        if fout: 
            fout.cd()
            para_slice.Write()
            perp_slice.Write()
            # write pdfs
            suf = ""
            if xi==1: suf="("
            if xi==nx: suf=")"
            para_slice.Draw()
            c.SaveAs("plots/fits_{}_para.pdf{}".format(name,suf))
            perp_slice.Draw()
            c.SaveAs("plots/fits_{}_perp.pdf{}".format(name,suf))


    resp,respe = ArrayDivideErrs(para,parae,z,ze)

    g_resp = TGraphAsymmErrors(len(z),array('f',z),array('f',resp),
                               array('f',zedn), array('f',zeup),
                               array('f',respe),array('f',respe))
    g_resPara = TGraphAsymmErrors(len(z),array('f',z),array('f',resPara),
                                  array('f',zedn), array('f',zeup),
                                  array('f',resParae),array('f',resParae))
    g_resPerp = TGraphAsymmErrors(len(z),array('f',z),array('f',resPerp),
                                  array('f',zedn), array('f',zeup),
                                  array('f',resPerpe),array('f',resPerpe))
    g_resp.SetName("g_{}_resp".format(name))
    g_resPara.SetName("g_{}_para_res".format(name))
    g_resPerp.SetName("g_{}_perp_res".format(name))

    return g_resp,g_resPara,g_resPerp

def RunAll2(infile, res_file, fit_file, name, ):
    hz    = f_in.Get("ptz_2d")
    hpara = f_in.Get(name+"_para")
    hperp = f_in.Get(name+"_perp")
    gs = RunFits(hz,hpara,hperp,name,fout=fit_file)
    res_file.cd()
    plot(gs[0],"<u_{\\parallel}>/<q_{T}>",(0,2))
    plot(gs[1],"\\sigma(u_{\\parallel})",(0,200))
    plot(gs[2],"\\sigma(u_{\\perp}  )",(0,200))


f_in = TFile("histograms.root","read")
f_out = TFile("results.root","recreate")
f_fits = TFile("fits.root","recreate")

# run analysis
RunAll2(f_in,f_out,f_fits,"pf")
RunAll2(f_in,f_out,f_fits,"puppi")
RunAll2(f_in,f_out,f_fits,"gen")

f_in.Close()
f_fits.Close()
f_out.Close()

