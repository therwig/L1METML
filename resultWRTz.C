float SgnPara(TVector2 x, TVector2 z){
    TVector2 para = x.Proj(z);
    return para.Mod()*(abs(para.DeltaPhi(z))>TMath::Pi()/2 ? 1. : -1.);
    // note here the '<' instead of '>' since pos = parallel here, not antiparallel
    // in other places (z recoil) we use the opposite convention
}
float SgnPerp(TVector2 x, TVector2 z){
    TVector2 perp = x.Norm(z);
    return perp.Mod()*(perp.DeltaPhi(z)>0 ? 1. : -1.);
}

float resultWRTz(float dPuppiPara, float dPuppiPerp, float puppi, float puppiPhi, float z, float zPhi, bool Para){
    // func takes three vectors, Puppi, dPuppi, and Z.
    // Calculates the components of (A+dA) parallel and perpendicular to -Z
    // (note that we choose -Z since MET is antiparallel to the z boson)
    TVector2 dPuppi(-dPuppiPara, dPuppiPerp);
    TVector2 Puppi, Z;
    Puppi.SetMagPhi(puppi,puppiPhi);
    Z.SetMagPhi(z,zPhi);
    dPuppi = dPuppi.Rotate(Puppi.Phi());

    // cout << "Z boson:             "; Z.Print();
    // cout << "L1PuppiMet recoil:   "; Puppi.Print();
    // cout << "correction:          "; dPuppi.Print();
    // cout << "corrected recoil:    "; (dPuppi+Puppi).Print();

    if(Para) return SgnPara(dPuppi+Puppi,-1.*Z);
    else return SgnPerp(dPuppi+Puppi,-1.*Z);
}


// t.Draw("test_resultWRTz(genMet_para_puppi,genMet_perp_puppi,L1PuppiMet,L1PuppiMetPhi,zPt,zPhi,1,genMet,genMetPhi):zPt>>h1","Entry$<1")
float test_resultWRTz(float daPara, float daPerp, float met, float metPhi, float z, float zPhi, 
                      bool Para, float genMet, float genMetPhi, float genMet_para_z, float genMet_perp_z){

    // TVector2 vPuppi;
    // vPuppi.SetMagPhi(met,metPhi);
    // TVector2 vZ;
    // vZ.SetMagPhi(z,zPhi);
    // TVector2 vGen;
    // vGen.SetMagPhi(genMet,genMetPhi);

    // float corr_para = -SgnPara(vGen-vPuppi,vPuppi);
    // float corr_perp = SgnPerp(vGen,vPuppi);
    // TVector2 vCorr(corr_para, corr_perp);
    // vCorr = vCorr.Rotate(vPuppi.Phi());

    // vCorr = vCorr + vPuppi;
    // // vCorr is now PuppiMET _RECOIL_, corrected by some amount, hopefully to make a better MET :)
    // // return the component of this which is aligned (or not with the Z boson)
    // // recall that recoil = -(MET+Z) we we can simply project it onto the z boson

    // float para = SgnPara(vCorr-vZ,vZ);
    // float perp = SgnPerp(vCorr,vZ);

    // if(Para) return para;
    // else return perp;



    TVector2 vPuppi;
    vPuppi.SetMagPhi(met,metPhi);
    TVector2 vZ;
    vZ.SetMagPhi(z,zPhi);
    TVector2 vGen;
    vGen.SetMagPhi(genMet,genMetPhi);

    cout << "Z boson:             "; vZ.Print();
    cout << "genMet recoil:       "; vGen.Print();
    cout << "L1PuppiMet recoil:   "; vPuppi.Print();
    // calc the correction to puppi met needed to match gen met
    cout << " == closure check == " << endl;
    float corr_para = -SgnPara(vGen-vPuppi,vPuppi);
    float corr_perp = SgnPerp(vGen,vPuppi);
    cout << " para =              " << corr_para << endl;
    cout << " perp =              " << corr_perp << endl;
    cout << " whereas precomputed para =   " << -daPara << endl;
    cout << " whereas precomputed perp =   " << daPerp << endl;
    TVector2 vCorr(corr_para, corr_perp);
    cout << "vCorr: "; vCorr.Print();
    vCorr = vCorr.Rotate(vPuppi.Phi());
    cout << "vCorr rotated:       "; vCorr.Print();
    vCorr = vCorr+vPuppi;
    cout << "puppi corr->genMet:  "; vCorr.Print();
    cout << " == FINAL ANSWER == " << endl;
    cout << "genMet:              "; vGen.Print();
    cout << "puppi corr->genMet:  "; vCorr.Print();
    cout << endl;

    float para = SgnPara(vCorr,-1.*vZ);
    float perp = SgnPerp(vCorr,-1.*vZ);
    cout << " == WRT Z == " << endl;
    cout << " genMet_para_z = " << genMet_para_z << endl;
    cout << "         myval = " << para << endl;
    cout << " genMet_perp_z = " << genMet_perp_z << endl;
    cout << "         myval = " << perp << endl;
    

    return 1.;

    // TVector2 v_met;
    // v_met.SetMagPhi(met,metPhi);
    // TVector2 v_z;
    // v_z.SetMagPhi(z,zPhi);
    // TVector2 v_genMet;
    // v_genMet.SetMagPhi(genMet,genMetPhi);

    // v_met = v_met + v_z;
    // v_genMet = v_genMet + v_z;

    // TVector2 v_corr(-daPara, daPerp); //correction to give genMET
    // v_corr = v_corr.Rotate(v_met.Phi());

    // cout << "v_met (plus with z)          "; v_met.Print();
    // cout << "v_genMet (plus with z)       "; v_genMet.Print();
    // cout << "v_corr + v_met  (plus with z)"; (v_corr+v_met).Print();

    // // cout << "v_met (plus with z)" << endl;
    // // v_met.Print();
    // // (v_met+v_z).Print();
    // // cout << "v_genMet (plus with z)" << endl;
    // // v_genMet.Print();
    // // (v_genMet+v_z).Print();
    // // cout << "v_z" << endl;
    // // v_z.Print();
    // // cout << "v_corr  (plus with z)" << endl;
    // // v_corr.Print();
    // // (v_corr+v_z).Print();
    // // cout << "v_corr+v_met  (plus with z)" << endl;
    // // (v_corr+v_met).Print();
    // // (v_corr+v_met+v_z).Print();

    // return 1.;
}
