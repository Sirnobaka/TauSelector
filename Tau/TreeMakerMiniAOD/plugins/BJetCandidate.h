//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Definition of class BJetCandidate
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class BJetCandidate {
public :
    //--- variables
    double Pt;
    double E;
    double Eta, Phi;
    math::XYZTLorentzVector FourMomentum;
    double Bprob;
    double BBprob;
    double LepBprob;
    double bProbCSV;
    double bbProbCSV;
    int HadronFlavour;
    bool FromPV;
    double SFLoose;
    double SFMedium;
    double SFTight;
    double SFReshaping;
    double SFReshaping_up;
    double SFReshaping_down;


    //--- constructor & destructor
    BJetCandidate();
    BJetCandidate(const pat::Jet &jet, bool IsFromPV, const BTagCalibration &Bcalibration)
        : Pt(jet.pt()), E(jet.energy()), Eta(jet.eta()), Phi(jet.phi()), FourMomentum(jet.p4()),
        Bprob(jet.bDiscriminator("pfDeepFlavourJetTags:probb")), BBprob(jet.bDiscriminator("pfDeepFlavourJetTags:probbb")),
        LepBprob(jet.bDiscriminator("pfDeepFlavourJetTags:problepb")), bProbCSV(jet.bDiscriminator("pfDeepCSVJetTags:probb")),
        bbProbCSV(jet.bDiscriminator("pfDeepCSVJetTags:probbb")), HadronFlavour(jet.hadronFlavour()),
        FromPV(IsFromPV)
    {

        BTagCalibrationReader reader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
        // calibration instance, btag flavour, measurement type = "comb" (need "ttbar"?)
        reader.load(Bcalibration, BTagEntry::FLAV_B, "comb");
        reader.load(Bcalibration, BTagEntry::FLAV_C, "comb");
        reader.load(Bcalibration, BTagEntry::FLAV_UDSG, "comb");
        // Same but with medium operating point
        BTagCalibrationReader readerMediumOP(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
        readerMediumOP.load(Bcalibration, BTagEntry::FLAV_B, "comb");
        readerMediumOP.load(Bcalibration, BTagEntry::FLAV_C, "comb");
        readerMediumOP.load(Bcalibration, BTagEntry::FLAV_UDSG, "comb");
        // Same but with tight operating point
        BTagCalibrationReader readerTightOP(BTagEntry::OP_TIGHT, "central", {"up", "down"});
        readerTightOP.load(Bcalibration, BTagEntry::FLAV_B, "comb");
        readerTightOP.load(Bcalibration, BTagEntry::FLAV_C, "comb");
        readerTightOP.load(Bcalibration, BTagEntry::FLAV_UDSG, "comb");
        // btag distribution reshaping
        BTagCalibrationReader readerReshaping(BTagEntry::OP_RESHAPING, "central", {"up_hfstats2", "down_hfstats2"});
        readerReshaping.load(Bcalibration, BTagEntry::FLAV_B, "iterativefit");
        readerReshaping.load(Bcalibration, BTagEntry::FLAV_C, "iterativefit");
        readerReshaping.load(Bcalibration, BTagEntry::FLAV_UDSG, "iterativefit");

        if (abs(HadronFlavour) == 5) {
            SFReshaping = readerReshaping.eval_auto_bounds("central", BTagEntry::FLAV_B, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFReshaping_up = readerReshaping.eval_auto_bounds("up_hfstats2", BTagEntry::FLAV_B, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFReshaping_down = readerReshaping.eval_auto_bounds("down_hfstats2", BTagEntry::FLAV_B, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFLoose = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, abs(Eta), Pt, (Bprob + BBprob));
            SFMedium = readerMediumOP.eval_auto_bounds("central", BTagEntry::FLAV_B, abs(Eta), Pt, (Bprob + BBprob));
            SFTight = readerTightOP.eval_auto_bounds("central", BTagEntry::FLAV_B, abs(Eta), Pt, (Bprob + BBprob));
        } else if (abs(HadronFlavour) == 4) {
            SFReshaping = readerReshaping.eval_auto_bounds("central", BTagEntry::FLAV_C, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFReshaping_up = readerReshaping.eval_auto_bounds("up_hfstats2", BTagEntry::FLAV_C, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFReshaping_down = readerReshaping.eval_auto_bounds("down_hfstats2", BTagEntry::FLAV_C, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFLoose = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFMedium = readerMediumOP.eval_auto_bounds("central", BTagEntry::FLAV_C, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFTight = readerTightOP.eval_auto_bounds("central", BTagEntry::FLAV_C, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
        } else {
            SFReshaping = readerReshaping.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFReshaping_up = readerReshaping.eval_auto_bounds("up_hfstats2", BTagEntry::FLAV_UDSG, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFReshaping_down = readerReshaping.eval_auto_bounds("down_hfstats2", BTagEntry::FLAV_UDSG, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFLoose = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFMedium = readerMediumOP.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
            SFTight = readerTightOP.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, abs(Eta), Pt, (Bprob + BBprob + LepBprob));
        }
        if (SFReshaping < 0.001) {
            SFReshaping = 1.;
        }
        if (SFReshaping_up < 0.001) {
            SFReshaping_up = 1.;
        }
        if (SFReshaping_down < 0.001) {
            SFReshaping_down = 1.;
        }
        if (SFLoose < 0.001) {
            SFLoose = 1.;
        }
        if (SFMedium < 0.001) {
            SFMedium = 1.;
        }
        if (SFTight < 0.001) {
            SFTight = 1.;
        }
    }
    BJetCandidate(const pat::Jet &jet, bool IsFromPV)
        : Pt(jet.pt()), E(jet.energy()), Eta(jet.eta()), Phi(jet.phi()), FourMomentum(jet.p4()),
        Bprob(jet.bDiscriminator("pfDeepFlavourJetTags:probb")), BBprob(jet.bDiscriminator("pfDeepFlavourJetTags:probbb")),
        LepBprob(jet.bDiscriminator("pfDeepFlavourJetTags:problepb")), bProbCSV(jet.bDiscriminator("pfDeepCSVJetTags:probb")),
        bbProbCSV(jet.bDiscriminator("pfDeepCSVJetTags:probbb")), HadronFlavour(jet.hadronFlavour()), FromPV(IsFromPV)
    {
        // No btag scales in analyser
        SFReshaping = 1;
        SFReshaping_up = 1;
        SFReshaping_down = 1;
        SFLoose = 1;
        SFMedium = 1;
        SFTight = 1;
    }

    void Monitoring();
    virtual ~BJetCandidate();

    //--- functions

};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Constructor of BJetCandidate object
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
BJetCandidate::BJetCandidate() {

};

BJetCandidate::~BJetCandidate()
{
// Destructor
};

void BJetCandidate::Monitoring() {

    std::cout << "pt      = " << Pt << std::endl
    << "eta       = " << Eta << std::endl
    << "phi       = " << Phi << std::endl
    << "bProb     = " << Bprob << std::endl
    << "bbProb    = " << BBprob << std::endl
    << "LepBprob  = " << LepBprob << std::endl
    << "bProbCSV  = " << bProbCSV << std::endl
    << "bbProbCSV = " << bbProbCSV << std::endl
    << "FromPV    = " << FromPV << std::endl;
    /*
    << "SF reshaping = " << SFReshaping << std::endl
    << "SF reshaping up = " << SFReshaping_up << std::endl
    << "SF reshaping down = " << SFReshaping_down << std::endl
    << "SF loose  = " << SFLoose << std::endl
    << "SF medium = " << SFMedium << std::endl
    << "SF tight  = " << SFTight << std::endl;
    */

};

// Jets in order of decreasing btag

void SortJets(std::vector<BJetCandidate>& items) {
  bool swapped;
  do {
    swapped = false;
    for (unsigned i = 1; i < items.size(); i++) {
      // push worward the jet with larger btag value or larger pt
      if ((items[i-1].Bprob + items[i-1].BBprob + items[i-1].LepBprob) < (items[i].Bprob + items[i].BBprob + items[i].LepBprob)) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      } else if ((items[i-1].Bprob + items[i-1].BBprob + items[i-1].LepBprob) == (items[i].Bprob + items[i].BBprob + items[i].LepBprob)
                  && items[i-1].Pt < items[i].Pt) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      }
    }
  } while (swapped != false);
};