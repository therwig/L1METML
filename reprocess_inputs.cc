// compile and run with 
// g++ reprocess_inputs.cc `root-config --cflags` `root-config --glibs`
// ./a.out
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using std::cout;
using std::endl;
using std::map;
using std::vector;
using std::pair;

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVector2.h"
#include "TLorentzVector.h"

//
// some helper classes
//

// read values from a tree corresponding to a 
// collection or 4-vectors
class TLVArrayReader {
    TTreeReaderArray<float> *_pt;
    TTreeReaderArray<float> *_eta;
    TTreeReaderArray<float> *_phi;
    TTreeReaderArray<float> *_m;
    TTreeReaderValue<unsigned int> *_n;
    vector<TLorentzVector> _v;
public:
    TLVArrayReader(): _pt(), _eta(), _phi(), _m(), _n(){}
    TLVArrayReader (TTreeReader *r, TString n){
        _pt  = new TTreeReaderArray<float>(*r,n+"_pt");
        _eta = new TTreeReaderArray<float>(*r,n+"_eta");
        _phi = new TTreeReaderArray<float>(*r,n+"_phi");
        _m   = new TTreeReaderArray<float>(*r,n+"_mass");
        _n   = new TTreeReaderValue<unsigned int>(*r,"n"+n);        
    }
    ~TLVArrayReader(){
    }
    inline float pt(unsigned int i){  return (*_pt )[i]; }
    inline float eta(unsigned int i){ return (*_eta)[i]; }
    inline float phi(unsigned int i){ return (*_phi)[i]; }
    inline float m(unsigned int i){   return (*_m  )[i]; }
    //inline int n(){ return *_n; }
    inline int n(){ return (*_pt ).GetSize(); }
    void FillTLVs(){
        _v.clear();
        for(int i=0;i<n();i++){
            TLorentzVector x;
            x.SetPtEtaPhiM((*_pt)[i],(*_eta)[i],(*_phi)[i],(*_m)[i]);
            _v.push_back(x);
        }
    }
    inline TLorentzVector* tlv(unsigned int i){ return &_v[i]; }
};


// helper class to hide the ugly variable reads w/ ttreereader
class MetReader{
    TTreeReaderValue<float> *_pt;
    TTreeReaderValue<float> *_phi;
    TVector2 _v;
public:
    MetReader(): _pt(), _phi() {}
    MetReader(TTreeReader *r, TString n){
        _pt = new TTreeReaderValue<float>(*r,n+"_pt");
        _phi = new TTreeReaderValue<float>(*r,n+"_phi");
    }
    ~MetReader(){
        if(_pt) delete _pt;
        if(_phi) delete _phi;
    }
    inline float pt(){ return **_pt;}
    inline float phi(){ return **_phi;}
    TVector2 GetV(){
        _v.SetMagPhi(pt(),phi());
        return _v;
    }
};

// helper class to hide the ugly variable writes w/ ttree
class MetWriter{
    float _pt;
    float _phi;
    float _para;
    float _perp;
public:
    MetWriter(): _para(), _perp(), _pt(), _phi() {}
    MetWriter(TString n, TTree *t, TString suffix="", bool writePtPhi=false){
        if (suffix.Length()) suffix = "_"+suffix;
        TString paran=n+"_para"+suffix;
        TString perpn=n+"_perp"+suffix;
        if(writePtPhi){
            t->Branch(n,       &_pt, n+"/F");
            t->Branch(n+"Phi", &_phi,n+"Phi/F");
        }
        t->Branch(paran, &_para, paran+"/F");
        t->Branch(perpn, &_perp, perpn+"/F");
    }
    ~MetWriter(){}
    void Set(TVector2 v, float para, float perp){
        _para = para;
        _perp = perp;
        _pt   = v.Mod();
        _phi  = v.Phi();
    }
    void Set(float para, float perp){
        _para = para;
        _perp = perp;
    }
};

//
// some helper functions
//
float SgnPara(const TVector2 recoil, const TVector2 z){
    TVector2 x = recoil.Proj(z);
    return x.Mod()*(fabs(x.DeltaPhi(z))>TMath::Pi()/2 ? 1.: -1.);
}
float SgnPerp(const TVector2 recoil, const TVector2 z){
    TVector2 x = recoil.Norm(z);
    return x.Mod()*(x.DeltaPhi(z)>0 ? 1.: -1.);
}

void BuildMET(TLVArrayReader* a, TVector2& v, float ptcut=0){
    TLorentzVector met;
    a->FillTLVs();
    for(int i=0;i<a->n();i++){
        if(ptcut==0 || a->pt(i) > ptcut)
            met += (*a->tlv(i));
    }
    v.SetMagPhi(met.Pt(), met.Phi());
    // do we need a zero check here?
}


int main(){    
    //
    // input values for precomputed MET
    //
    TFile* fin = new TFile("samples/perfNano_DY.root","read");
    TTree* tin = (TTree*) fin->Get("Events");
    TTreeReader reader(tin);
    vector<TString> met_inputs{"L1CaloMet","L1CaloMetBarrel","L1CaloMetCentral",
            "L1PFMet","L1PFMetBarrel","L1PFMetCentral",
            "L1PuppiMet","L1PuppiMetBarrel","L1PuppiMetCentral",
            "genMet","genMetBarrel","genMetCentral",
            };
    map<TString,MetReader*> input_map;
    for(auto n : met_inputs) input_map[n] = new MetReader(&reader,n);
    // inputs jets for HTmiss computation
    map<TString,TLVArrayReader*> jet_map;
    for(TString s : {"Calo","PF","Puppi"}){
        jet_map[s+"HTMiss"] = new TLVArrayReader(&reader,"L1"+s+"Jets");
    }
    //muons for z reconstruction
    //const bool reqZ=false;
    const bool reqZ=true;
    TLVArrayReader muons(&reader,"PFMu");

    //
    // MET outputs
    //
    TString foutname=reqZ?"samples/perfNano_DY_reprocessed.root":"samples/perfNano_DY_reprocessed_all.root";
    TFile* f = new TFile(foutname,"recreate");
    TTree* t = new TTree("Events","");    
    map<TString,MetWriter*> output_map_z;
    map<TString,MetWriter*> output_map_puppi;
    vector<TString> out_keys;
    for(auto n : met_inputs){
        out_keys.push_back(n);
        if(reqZ) output_map_z[n] = new MetWriter(n,t,"z");
        output_map_puppi[n] = new MetWriter(n,t,"puppi",true); // write also met pt/phi
    }
    for(auto x : jet_map){
        TString n = x.first;
        out_keys.push_back(n);
        if(reqZ) output_map_z[n] = new MetWriter(n,t,"z");
        output_map_puppi[n] = new MetWriter(n,t,"puppi",true);
    }

    // a few other variables to log in the output tree
    float mll; t->Branch("mll", &mll, "mll/F");
    //float ptz; t->Branch("ptz", &ptz, "ptz/F");
    int nmu; t->Branch("nmu", &nmu, "nmu/I");
    int event; t->Branch("event", &event, "event/I");

    float zPt,zPhi; 
    t->Branch("zPt", &zPt, "zPt/F");
    t->Branch("zPhi", &zPhi, "zPhi/F");
    //float genMet,genMetPhi,puppiMet,puppiMetPhi,zPt,zPhi; 
    // t->Branch("genMet", &genMet, "genMet/F");
    // t->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
    // t->Branch("puppiMet", &puppiMet, "puppiMet/F");
    // t->Branch("puppiMetPhi", &puppiMetPhi, "puppiMetPhi/F");


    // eventloop helpers
    TVector2 vGen,vZ,vPuppi,vPuppiRecoil, vin, vrecoil, vPuppiMetProxy, vinMetProxy;
    TLorentzVector tlvZ;

    uint64_t ie=0;
    uint64_t nwritten=0;
    while (reader.Next()){
        ie++;
        //if (ie>10) break;

        event=ie;
        nmu = muons.n();
        mll=-1.;
        zPt=-1.;
        zPhi=-999.;
        if(nmu==2){
            muons.FillTLVs(); // only fill data if needed
            tlvZ = *(muons.tlv(0)) + *(muons.tlv(1));
            mll=tlvZ.M();
            zPt =tlvZ.Pt();
            zPhi=tlvZ.Phi();
            if(reqZ && fabs(tlvZ.M()-91.)<15.){
                vZ.SetMagPhi(tlvZ.Pt(),tlvZ.Phi());
            } else if(reqZ) continue;
        } else if(reqZ) continue;

        vGen = input_map["genMet"]->GetV();
        vPuppi = input_map["L1PuppiMet"]->GetV();
        //vPuppiRecoil = -1.*(vPuppi+vZ);
        vPuppiMetProxy = vPuppi+vZ;
        //be careful about the signs here!

        // precomputed MET
        for(auto name : out_keys){
            // get MET vector (or build from jets)
            if (name.Contains("HTMiss")){
                BuildMET(jet_map[name], vin);
            } else {
                vin = input_map[name]->GetV();
            }
            auto writer_z = output_map_z[name];
            auto writer_puppi = output_map_puppi[name];
            // get MET wrt the PUPPI estimate
            vinMetProxy = (vin+vZ);
            writer_puppi->Set(vinMetProxy,
                              SgnPara(vinMetProxy-vPuppiMetProxy,vPuppiMetProxy),
                              SgnPerp(vinMetProxy,vPuppiMetProxy));
            // get met para, perp components wrt Z
            if(reqZ){
                vrecoil = -1.*(vin+vZ);
                writer_z->Set(SgnPara(vrecoil,vZ), 
                              SgnPerp(vrecoil,vZ));
            }
        }
        // write some extra vars to tree
        // puppiMet    = vPuppiMetProxy.Mag();
        // puppiMetPhi = vPuppiMetProxy.Phi();
        // genMet      = vGen.Mag();
        // genMetPhi   = vGen.Phi();
            
        t->Fill();
        nwritten++;
    }

    //cleanup
    cout << "Wrote " << nwritten << " entries to " << f->GetName() << endl;
    f->Write();
    fin->Close();
    f->Close();

    return 0;
}
