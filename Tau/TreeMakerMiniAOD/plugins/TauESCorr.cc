#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"
#include "TGraphAsymmErrors.h"
//#include <TGraphAsymmErrors.h>
#include <iostream> // std::cerr, std::endl
#include <iomanip>
#include <assert.h> // assert

// Work with DM-dependent energy scale 

TH1* getHistogramTauES () {
    TFile *file = new TFile("/afs/cern.ch/user/a/aoskin/Tau_packeges/CMSSW_10_6_20/src/TauPOG/TauIDSFs/data/TauES_dm_DeepTau2017v2p1VSjet_UL2017.root");
    TH1* hist = dynamic_cast<TH1*>((const_cast<TFile*>(file))->Get("tes"));
    if(!hist) {
        std::cerr << std::endl << "ERROR! Failed to load histogram from input file!" << std::endl;
        assert(0);
    }
    return hist;
}


float getESvsDM(TH1* hist, int dm, const std::string& unc="") {

    Int_t bin   = hist->GetXaxis()->FindBin(dm);
    float tes   = hist->GetBinContent(bin);

    if(unc=="Up")
        tes += hist->GetBinError(bin);
    else if(unc=="Down")
        tes -= hist->GetBinError(bin);
    return tes;
}

// Electrons faking as taus
float getFESvsDMeta(int dm, double eta, const std::string& unc="") {
    // From readme
    TFile *file = new TFile("/afs/cern.ch/user/a/aoskin/Tau_packeges/CMSSW_10_6_20/src/TauPOG/TauIDSFs/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2017ReReco.root");
    // extract TGraph
    //TGraphAsymmErrors* gfes = dynamic_cast<TGraphAsymmErrors*>((const_cast<TFile*>(file))->Get("Graph"));
    TGraphAsymmErrors* gfes = (TGraphAsymmErrors*)file->Get("Graph");
    // or as in git
    // TGraph* gfes_2017=(TGraph*) ffes2017.Get("fes");
    float tes_ele = 1.0;
    //TH1* hist = dynamic_cast<TH1*>((const_cast<TFile*>(file))->Get("tes"));
    if(!gfes) {
        std::cerr << std::endl << "ERROR! Failed to load graph from input file!" << std::endl;
        assert(0);
    }
    if (dm == 0 && abs(eta) < 1.5) {
        tes_ele = gfes->GetY()[0];
        if(unc=="Up")
            tes_ele += gfes->GetErrorYhigh(0);
        else if(unc=="Down")
            tes_ele -= gfes->GetErrorYlow(0);
    }
    if (dm == 1 && abs(eta) < 1.5) {
        tes_ele = gfes->GetY()[1];
        if(unc=="Up")
            tes_ele += gfes->GetErrorYhigh(1);
        else if(unc=="Down")
            tes_ele -= gfes->GetErrorYlow(1);
    }
    if (dm == 0 && abs(eta) >= 1.5) {
        tes_ele = gfes->GetY()[2];
        if(unc=="Up")
            tes_ele += gfes->GetErrorYhigh(2);
        else if(unc=="Down")
            tes_ele -= gfes->GetErrorYlow(2);
    }
    if (dm == 1 && abs(eta) >= 1.5) {
        tes_ele = gfes->GetY()[3];
        if(unc=="Up")
            tes_ele += gfes->GetErrorYhigh(3);
        else if(unc=="Down")
            tes_ele -= gfes->GetErrorYlow(3);
    }
    //
    file->Close();
    return tes_ele;
}