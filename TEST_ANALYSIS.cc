// -*- C++ -*-
#include "lib/SelectPhoton.hh"

namespace Rivet {


  /// Inclusive isolated prompt photon analysis with Atlas data
  class TEST_ANALYSIS : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TEST_ANALYSIS);

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      declare(fs, "FS");

      // Consider the final state jets for the energy density calculation
      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      declare(fj, "KtJetsD05");

      // Consider the leading pt photon with |eta| < 2.37 and pT > 25 GeV
      LeadingParticlesFinalState photonfs(PromptFinalState(FinalState(Cuts::abseta < 2.37 && Cuts::pT > 25*GeV)));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // Book the dsigma/dEt (in eta bins) histograms
      for (size_t i = 0; i < _eta_bins.size() - 1; ++i) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue; // skip this bin
        int offset = i > 2? 0 : 1;
        book(_h_Et_photon[i] ,i + offset, 1, 1);
      }

    }


    /// Return eta bin for either dsigma/dET histogram (area_eta=false) or energy density correction (area_eta=true)
    size_t _getEtaBin(double eta_w, bool area_eta) const {
      const double eta = fabs(eta_w);
      if (!area_eta) {
        return binIndex(eta, _eta_bins);
      } else {
        return binIndex(eta, _eta_bins_areaoffset);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve leading photon
      Particles photons = applyProjection<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();

      if (photons.size() < 1)  vetoEvent;
      const Particle& leadingPhoton = photons[0];
      FourMomentum leadingMomentum = leadingPhoton.momentum();
      // Veto events that wouldn't survived if they were LIV, and look for a missID electron
      if(HasSurvived(leadingMomentum)==0) vetoEvent;

      // Veto events with photon in ECAL crack
      if (inRange(leadingMomentum.abseta(), 1.37, 1.56)) vetoEvent;

      // Compute isolation energy in cone of radius .4 around photon (all particles)
      FourMomentum mom_in_EtCone;
      Particles fs = applyProjection<FinalState>(event, "FS").particles();
      for (const Particle& p : fs) {
        // Check if it's outside the cone of 0.4
        if (deltaR(leadingMomentum, p.momentum()) >= 0.4) continue;
        // Except muons or neutrinos
        if (PID::isNeutrino(p.abspid()) || p.abspid() == PID::MUON) continue;
        // Increment isolation energy
        mom_in_EtCone += p.momentum();
      }
      // Remove the photon energy from the isolation
      mom_in_EtCone -= leadingMomentum;

      // Get the area-filtered jet inputs for computing median energy density, etc.
      vector<double> ptDensity;
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      const FastJets& fast_jets = applyProjection<FastJets>(event, "KtJetsD05");
      const auto clust_seq_area = fast_jets.clusterSeqArea();
      for (const Jet& jet : fast_jets.jets()) {
        const double area = clust_seq_area->area(jet);
        if (area > 1e-3 && jet.abseta() < _eta_bins_areaoffset.back())
          ptDensities.at( _getEtaBin(jet.abseta(), true) ) += jet.pT()/area;
      }
      // Compute the median energy density, etc.
      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; ++b) {
        const int njets = ptDensities[b].size();
        ptDensity += (njets > 0) ? median(ptDensities[b]) : 0.0;
      }
      // Compute the isolation energy correction (cone area*energy density)
      const double etCone_area = PI * sqr(0.4);
      const double correction = ptDensity[_getEtaBin(leadingMomentum.abseta(), true)] * etCone_area;

      // Apply isolation cut on area-corrected value
      // cut is Etiso < 4.8GeV + 4.2E-03 * Et_gamma.
      if (mom_in_EtCone.Et() - correction > 4.8*GeV + 0.0042*leadingMomentum.Et()) vetoEvent;

      // Fill histograms
      const size_t eta_bin = _getEtaBin(leadingMomentum.abseta(), false);
      _h_Et_photon[eta_bin]->fill(leadingMomentum.Et());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double TriggerEfficiency=1.00; //Efficacite de trigger pour un photon
	    double IdentificationEfficiency=0.95; //Efficacite d'identification pour un photon
	    double ReconstructionEfficiency=0.70; //Efficacite de reconstruction pour un photon
      double sf = crossSection() / (picobarn * sumOfWeights()) * TriggerEfficiency * IdentificationEfficiency * ReconstructionEfficiency;
      // double sf = crossSection() / (picobarn * sumOfWeights());
      for (size_t i = 0; i < _eta_bins.size()-1; ++i) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue;
        scale(_h_Et_photon[i], sf);
      }
    }


  private:

    Histo1DPtr _h_Et_photon[5];

    const vector<double> _eta_bins = {0.00, 0.60, 1.37, 1.56, 1.81, 2.37 };
    const vector<double> _eta_bins_areaoffset = {0.0, 1.5, 3.0};

  };


  RIVET_DECLARE_PLUGIN(TEST_ANALYSIS);

}
