/// Skim ntuples using RDataFrame for di-muon analysis
/// Cutflow:
///   0. All events
///   1. L1 Trigger (OR of DoubleMu seeds)
///   2. 2 quality OS (sr) / SS (vr) muons
///   3. + di-pT > 20 GeV
///   4. + dR < 1
///
/// Usage:
///       root -l -b -q 'skim_ntuples.C("sr","/path/to/files/")'
///       root -l -b -q 'skim_ntuples.C("vr", "/path/to/files/")'
///       root -l -b -q 'skim_ntuples.C("sr", "/path/to/files/", true)'  // validate files

void skim_ntuples(TString region = "sr",
                  TString base_dir = "/ceph/cms/store/user/tvami/DiMuonScoutingNtuples/",
                  bool validate = false) {

    ROOT::EnableImplicitMT();

    // Validate region
    if (region != "sr" && region != "vr") {
        std::cerr << "Error: Unknown region '" << region << "'. Use 'sr' or 'vr'.\n";
        return;
    }

    // Output file naming
    TString dataset_name;
    if (base_dir.Contains("ScoutingPFRun3")) {
        // Data path: /ceph/cms/store/data/Run2024E/ScoutingPFRun3/NANOAOD/ScoutNano-v1/2520000
        // Extract: ScoutingPFRun3_Run2024E_2520000
        TObjArray* tokens = base_dir.Tokenize("/");
        TString run_era = "", block = "";
        for (int i = 0; i < tokens->GetEntries(); ++i) {
            TString tok = ((TObjString*)tokens->At(i))->GetString();
            if (tok.BeginsWith("Run20")) run_era = tok;
        }
        TString base_dir_copy = base_dir;
        if (base_dir_copy.EndsWith("/")) base_dir_copy.Remove(base_dir_copy.Length()-1);
        block = base_dir_copy(base_dir_copy.Last('/')+1, base_dir_copy.Length());
        dataset_name = TString::Format("ScoutingPFRun3_%s_%s", run_era.Data(), block.Data());
        delete tokens;
    } else {
        TString base_dir_copy = base_dir;
        if (base_dir_copy.EndsWith("/")) base_dir_copy.Remove(base_dir_copy.Length()-1);
        dataset_name = base_dir_copy(base_dir_copy.Last('/')+1, base_dir_copy.Length());
        dataset_name.ReplaceAll("crab_", "");
    }
    dataset_name.ReplaceAll(".root", "");
    TString output_file = TString::Format("skimmed_%s_%s.root", region.Data(), dataset_name.Data());

    // Build TChain
    TChain chain("Events");

    TString find_cmd = TString::Format("find %s -name '*.root' 2>/dev/null", base_dir.Data());
    FILE* pipe = popen(find_cmd.Data(), "r");
    if (!pipe) {
        std::cerr << "Error: cannot run find command\n";
        return;
    }

    std::cout << "Reading filenames from pipe\n";

    char buffer[2048];
    int nfiles = 0, nzombie = 0;
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        TString filename = buffer;
        filename.ReplaceAll("\n","");
        if (filename.Length() == 0) continue;
        if (validate) {
            TFile* test_file = TFile::Open(filename.Data(), "READ");
            if (test_file && !test_file->IsZombie() && test_file->GetNkeys() > 0) {
                test_file->Close();
                delete test_file;
                chain.Add(filename.Data());
                nfiles++;
            } else {
                std::cerr << "Warning: Skipping zombie: " << filename << "\n";
                nzombie++;
                if (test_file) delete test_file;
            }
        } else {
            chain.Add(filename.Data());
            nfiles++;
        }
        if (nfiles % 10 == 0) std::cout << "Added " << nfiles << " files...\r" << std::flush;
    }
    pclose(pipe);

    if (nzombie > 0) std::cerr << "Skipped " << nzombie << " corrupted files\n";
    if (nfiles == 0) {
        std::cerr << "Error: No ROOT files found.\n";
        return;
    }
    std::cout << "Found " << nfiles << " files\n";

    ROOT::RDataFrame df(chain);

    // --- Define quality object counts (available at all cutflow steps) ---

    // Quality muons: pT > 3 GeV, |eta| < 2.4, trackIso < 0.15, chi2/ndof < 3
    auto df_q = df
        .Define("n_quality_muons", [](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& trackIso,
                const ROOT::VecOps::RVec<float>& normchi2) {
            int n = 0;
            for (size_t i = 0; i < pt.size(); ++i)
                if (pt[i] > 3 && std::abs(eta[i]) < 2.4 && trackIso[i] < 0.15 && normchi2[i] < 3) n++;
            return n;
        }, {"ScoutingMuonVtx_pt", "ScoutingMuonVtx_eta",
            "ScoutingMuonVtx_trackIso", "ScoutingMuonVtx_normchi2"})

        // Quality photons: pT > 5, |eta| < 2.4, H/E < 0.2,
        //   barrel (|eta|<=1.479): sigmaIetaIeta<0.015, ecalIso/E<0.25, hcalIso/E<0.60
        //   endcap (|eta|>1.479):  sigmaIetaIeta<0.045, ecalIso/E<0.10, hcalIso/E<1.00
        .Define("n_quality_photons", [](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& hOverE,
                const ROOT::VecOps::RVec<float>& sieie,
                const ROOT::VecOps::RVec<float>& ecalIso,
                const ROOT::VecOps::RVec<float>& hcalIso) {
            int n = 0;
            for (size_t i = 0; i < pt.size(); ++i) {
                if (pt[i] <= 5 || std::abs(eta[i]) >= 2.4 || hOverE[i] >= 0.2) continue;
                float E = pt[i] * std::cosh(eta[i]);
                if (std::abs(eta[i]) <= 1.479) {
                    if (sieie[i] < 0.015 && ecalIso[i]/E < 0.25 && hcalIso[i]/E < 0.60) n++;
                } else {
                    if (sieie[i] < 0.045 && ecalIso[i]/E < 0.10 && hcalIso[i]/E < 1.00) n++;
                }
            }
            return n;
        }, {"ScoutingPhoton_pt", "ScoutingPhoton_eta", "ScoutingPhoton_hOverE",
            "ScoutingPhoton_sigmaIetaIeta", "ScoutingPhoton_ecalIso", "ScoutingPhoton_hcalIso"})

        // Quality electrons: pT > 5, |eta| < 2.4, H/E < 0.2, |dPhiIn| < 0.06, trackIso/E < 0.001,
        //   barrel: sigmaIetaIeta<0.015, |dEtaIn|<0.008, ecalIso/E<0.25, hcalIso/E<0.40
        //   endcap: sigmaIetaIeta<0.045, |dEtaIn|<0.012, ecalIso/E<0.10, hcalIso/E<0.60
        .Define("n_quality_electrons", [](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& hOverE,
                const ROOT::VecOps::RVec<float>& sieie,
                const ROOT::VecOps::RVec<float>& dEtaIn,
                const ROOT::VecOps::RVec<float>& dPhiIn,
                const ROOT::VecOps::RVec<float>& ecalIso,
                const ROOT::VecOps::RVec<float>& hcalIso,
                const ROOT::VecOps::RVec<float>& trackIso) {
            int n = 0;
            for (size_t i = 0; i < pt.size(); ++i) {
                if (pt[i] <= 5 || std::abs(eta[i]) >= 2.4 || hOverE[i] >= 0.2) continue;
                if (std::abs(dPhiIn[i]) >= 0.06) continue;
                float E = pt[i] * std::cosh(eta[i]);
                if (trackIso[i]/E >= 0.001) continue;
                if (std::abs(eta[i]) <= 1.479) {
                    if (sieie[i] < 0.015 && std::abs(dEtaIn[i]) < 0.008 &&
                        ecalIso[i]/E < 0.25 && hcalIso[i]/E < 0.40) n++;
                } else {
                    if (sieie[i] < 0.045 && std::abs(dEtaIn[i]) < 0.012 &&
                        ecalIso[i]/E < 0.10 && hcalIso[i]/E < 0.60) n++;
                }
            }
            return n;
        }, {"ScoutingElectron_pt", "ScoutingElectron_eta", "ScoutingElectron_hOverE",
            "ScoutingElectron_sigmaIetaIeta", "ScoutingElectron_dEtaIn", "ScoutingElectron_dPhiIn",
            "ScoutingElectron_ecalIso", "ScoutingElectron_hcalIso", "ScoutingElectron_trackIso"})

        // Quality jets (low eta, ScoutingPFJetRecluster): pT > 30, |eta| < 2.5,
        //   neHEF<0.99, neEmEF<0.90, nConstituents>1, muEF<0.80, chHEF>0.01, nCh>0
        .Define("n_quality_jets_loweta", [](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& neHEF,
                const ROOT::VecOps::RVec<float>& neEmEF,
                const ROOT::VecOps::RVec<unsigned char>& nConstituents,
                const ROOT::VecOps::RVec<float>& muEF,
                const ROOT::VecOps::RVec<float>& chHEF,
                const ROOT::VecOps::RVec<int>& nCh) {
            int n = 0;
            for (size_t i = 0; i < pt.size(); ++i)
                if (pt[i] > 30 && std::abs(eta[i]) < 2.5 &&
                    neHEF[i] < 0.99 && neEmEF[i] < 0.90 && nConstituents[i] > 1 &&
                    muEF[i] < 0.80 && chHEF[i] > 0.01 && nCh[i] > 0) n++;
            return n;
        }, {"ScoutingPFJetRecluster_pt", "ScoutingPFJetRecluster_eta",
            "ScoutingPFJetRecluster_neHEF", "ScoutingPFJetRecluster_neEmEF",
            "ScoutingPFJetRecluster_nConstituents", "ScoutingPFJetRecluster_muEF",
            "ScoutingPFJetRecluster_chHEF", "ScoutingPFJetRecluster_nCh"})

        // Quality jets (high eta, ScoutingPFJet): pT > 50, 3.0 < |eta| < 5.0, neutralEMFrac < 0.1
        .Define("n_quality_jets_higheta", [](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& m,
                const ROOT::VecOps::RVec<float>& HFEMEnergy) {
            int n = 0;
            for (size_t i = 0; i < pt.size(); ++i) {
                if (pt[i] <= 50 || std::abs(eta[i]) <= 3.0 || std::abs(eta[i]) >= 5.0) continue;
                float E = std::sqrt(pt[i]*pt[i]*std::cosh(eta[i])*std::cosh(eta[i]) + m[i]*m[i]);
                if (E > 0 && HFEMEnergy[i]/E < 0.1) n++;
            }
            return n;
        }, {"ScoutingPFJet_pt", "ScoutingPFJet_eta", "ScoutingPFJet_m", "ScoutingPFJet_HFEMEnergy"})

        .Define("n_quality_jets", "n_quality_jets_loweta + n_quality_jets_higheta");

    // --- Cutflow ---

    // Step 1: L1 Trigger (OR of 6 DoubleMu seeds)
    auto df_trigger = df_q.Filter(
        "L1_DoubleMu_12_5 || L1_DoubleMu_15_7 || "
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 || L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18 || "
        "L1_DoubleMu4_SQ_OS_dR_Max1p2 || L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        "L1 Trigger");

    // Steps 2-4: Find best qualifying di-muon pair per event
    // Quality muon: pT > 3 GeV, |eta| < 2.4, trackIso < 0.15, chi2/ndof < 3
    // Charge: OS for sr, SS for vr
    // Best pair chosen by highest di-pT
    // Returns RVec<float>{idx1, idx2, dipt, dr, dimass} or all -1 if no pair found
    std::string region_str = region.Data();
    auto df_pairs = df_trigger
        .Define("dimu_pair", [region_str](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& phi,
                const ROOT::VecOps::RVec<float>& mass,
                const ROOT::VecOps::RVec<int>& charge,
                const ROOT::VecOps::RVec<float>& trackIso,
                const ROOT::VecOps::RVec<float>& normchi2) {
            ROOT::VecOps::RVec<float> result = {-1.f, -1.f, -1.f, -1.f, -1.f};

            // Find quality muon indices
            std::vector<int> good;
            for (size_t i = 0; i < pt.size(); ++i) {
                if (pt[i] > 3 && std::abs(eta[i]) < 2.4 && trackIso[i] < 0.15 && normchi2[i] < 3)
                    good.push_back(i);
            }
            if (good.size() < 2) return result;

            float best_dipt = -1;
            for (size_t a = 0; a < good.size(); ++a) {
                for (size_t b = a + 1; b < good.size(); ++b) {
                    int i = good[a], j = good[b];
                    bool os = charge[i] * charge[j] < 0;
                    if (region_str == "sr" && !os) continue;
                    if (region_str == "vr" && os) continue;

                    // Di-muon pT (vector sum)
                    float px = pt[i]*std::cos(phi[i]) + pt[j]*std::cos(phi[j]);
                    float py = pt[i]*std::sin(phi[i]) + pt[j]*std::sin(phi[j]);
                    float dipt = std::sqrt(px*px + py*py);

                    // DeltaR
                    float deta = eta[i] - eta[j];
                    float dphi = phi[i] - phi[j];
                    while (dphi >  M_PI) dphi -= 2*M_PI;
                    while (dphi < -M_PI) dphi += 2*M_PI;
                    float dr = std::sqrt(deta*deta + dphi*dphi);

                    // Invariant mass
                    float px1 = pt[i]*std::cos(phi[i]), py1 = pt[i]*std::sin(phi[i]);
                    float pz1 = pt[i]*std::sinh(eta[i]);
                    float e1  = std::sqrt(pt[i]*pt[i]*std::cosh(eta[i])*std::cosh(eta[i]) + mass[i]*mass[i]);
                    float px2 = pt[j]*std::cos(phi[j]), py2 = pt[j]*std::sin(phi[j]);
                    float pz2 = pt[j]*std::sinh(eta[j]);
                    float e2  = std::sqrt(pt[j]*pt[j]*std::cosh(eta[j])*std::cosh(eta[j]) + mass[j]*mass[j]);
                    float dimass_sq = (e1+e2)*(e1+e2) - (px1+px2)*(px1+px2)
                                    - (py1+py2)*(py1+py2) - (pz1+pz2)*(pz1+pz2);
                    float dimass = dimass_sq > 0 ? std::sqrt(dimass_sq) : 0.f;

                    if (dipt > best_dipt) {
                        best_dipt = dipt;
                        result = {(float)i, (float)j, dipt, dr, dimass};
                    }
                }
            }
            return result;
        }, {"ScoutingMuonVtx_pt", "ScoutingMuonVtx_eta", "ScoutingMuonVtx_phi", "ScoutingMuonVtx_m",
            "ScoutingMuonVtx_charge", "ScoutingMuonVtx_trackIso", "ScoutingMuonVtx_normchi2"})
        .Define("dimu_idx1",  "static_cast<int>(dimu_pair[0])")
        .Define("dimu_idx2",  "static_cast<int>(dimu_pair[1])")
        .Define("dimu_pt",    "dimu_pair[2]")
        .Define("dimu_dr",    "dimu_pair[3]")
        .Define("dimu_mass",  "dimu_pair[4]");

    // Step 2: Has at least one qualifying pair
    TString charge_label = (region == "sr") ? "2 quality OS #mu" : "2 quality SS #mu";
    auto df_dimu = df_pairs.Filter("dimu_idx1 >= 0", charge_label.Data());

    // Step 3: di-pT > 20 GeV
    auto df_dipt = df_dimu.Filter("dimu_pt > 20", "di-p_{T} > 20 GeV");

    // Step 4: dR < 1
    auto df_dr = df_dipt.Filter("dimu_dr < 1", "#DeltaR < 1");

    // Define count actions BEFORE Snapshot so they run in the same event loop
    auto count_all     = df_q.Count();
    auto count_trigger = df_trigger.Count();
    auto count_dimu    = df_dimu.Count();
    auto count_dipt    = df_dipt.Count();
    auto count_dr      = df_dr.Count();

    // --- Book histograms ---

    // 2D: invariant mass at each cutflow step (steps 2-4, where a pair exists)
    // Each event fills at step 2, and also at step 3 if it passes, and also step 4
    auto df_mass2d = df_dimu
        .Define("mass_cf_steps", [](float dipt, float dr) {
            ROOT::VecOps::RVec<double> s = {0.};
            if (dipt > 20) { s.push_back(1.); if (dr < 1) s.push_back(2.); }
            return s;
        }, {"dimu_pt", "dimu_dr"})
        .Define("mass_cf_vals", [](float mass, float dipt, float dr) {
            ROOT::VecOps::RVec<double> v = {(double)mass};
            if (dipt > 20) { v.push_back(mass); if (dr < 1) v.push_back(mass); }
            return v;
        }, {"dimu_mass", "dimu_pt", "dimu_dr"});
    auto h2_mass_cutflow = df_mass2d.Histo2D(
        {"h2_mass_cutflow",
         "Di-muon mass at cutflow steps;Cutflow step;m_{#mu#mu} [GeV]",
         3, -0.5, 2.5, 200, 0, 200},
        "mass_cf_steps", "mass_cf_vals");

    // 1D: invariant mass at each cutflow step
    auto h_mass_dimu = df_dimu.Histo1D({"h_mass_dimu", "m_{#mu#mu} (2 quality #mu);m_{#mu#mu} [GeV];Events", 200, 0, 200}, "dimu_mass");
    auto h_mass_dipt = df_dipt.Histo1D({"h_mass_dipt", "m_{#mu#mu} (di-p_{T} > 20);m_{#mu#mu} [GeV];Events", 200, 0, 200}, "dimu_mass");
    auto h_mass_dr   = df_dr.Histo1D({"h_mass_dr", "m_{#mu#mu} (#DeltaR < 1);m_{#mu#mu} [GeV];Events", 200, 0, 200}, "dimu_mass");

    // Object multiplicities at each cutflow step
    // -- All events --
    auto h_nmu_all  = df_q.Histo1D({"h_nmu_all",  "N quality #mu (all);N_{#mu};Events",  20, -0.5, 19.5}, "n_quality_muons");
    auto h_npho_all = df_q.Histo1D({"h_npho_all", "N quality #gamma (all);N_{#gamma};Events", 20, -0.5, 19.5}, "n_quality_photons");
    auto h_nele_all = df_q.Histo1D({"h_nele_all", "N quality e (all);N_{e};Events",      20, -0.5, 19.5}, "n_quality_electrons");
    auto h_njet_all = df_q.Histo1D({"h_njet_all", "N quality jet (all);N_{jet};Events",  20, -0.5, 19.5}, "n_quality_jets");
    // -- After trigger --
    auto h_nmu_trig  = df_trigger.Histo1D({"h_nmu_trig",  "N quality #mu (trigger);N_{#mu};Events",  20, -0.5, 19.5}, "n_quality_muons");
    auto h_npho_trig = df_trigger.Histo1D({"h_npho_trig", "N quality #gamma (trigger);N_{#gamma};Events", 20, -0.5, 19.5}, "n_quality_photons");
    auto h_nele_trig = df_trigger.Histo1D({"h_nele_trig", "N quality e (trigger);N_{e};Events",      20, -0.5, 19.5}, "n_quality_electrons");
    auto h_njet_trig = df_trigger.Histo1D({"h_njet_trig", "N quality jet (trigger);N_{jet};Events",  20, -0.5, 19.5}, "n_quality_jets");
    // -- After 2 quality muons --
    auto h_nmu_dimu  = df_dimu.Histo1D({"h_nmu_dimu",  "N quality #mu (2#mu);N_{#mu};Events",  20, -0.5, 19.5}, "n_quality_muons");
    auto h_npho_dimu = df_dimu.Histo1D({"h_npho_dimu", "N quality #gamma (2#mu);N_{#gamma};Events", 20, -0.5, 19.5}, "n_quality_photons");
    auto h_nele_dimu = df_dimu.Histo1D({"h_nele_dimu", "N quality e (2#mu);N_{e};Events",      20, -0.5, 19.5}, "n_quality_electrons");
    auto h_njet_dimu = df_dimu.Histo1D({"h_njet_dimu", "N quality jet (2#mu);N_{jet};Events",  20, -0.5, 19.5}, "n_quality_jets");
    // -- After di-pT > 20 --
    auto h_nmu_dipt  = df_dipt.Histo1D({"h_nmu_dipt",  "N quality #mu (di-p_{T}>20);N_{#mu};Events",  20, -0.5, 19.5}, "n_quality_muons");
    auto h_npho_dipt = df_dipt.Histo1D({"h_npho_dipt", "N quality #gamma (di-p_{T}>20);N_{#gamma};Events", 20, -0.5, 19.5}, "n_quality_photons");
    auto h_nele_dipt = df_dipt.Histo1D({"h_nele_dipt", "N quality e (di-p_{T}>20);N_{e};Events",      20, -0.5, 19.5}, "n_quality_electrons");
    auto h_njet_dipt = df_dipt.Histo1D({"h_njet_dipt", "N quality jet (di-p_{T}>20);N_{jet};Events",  20, -0.5, 19.5}, "n_quality_jets");
    // -- After dR < 1 (final) --
    auto h_nmu_dr  = df_dr.Histo1D({"h_nmu_dr",  "N quality #mu (#DeltaR<1);N_{#mu};Events",  20, -0.5, 19.5}, "n_quality_muons");
    auto h_npho_dr = df_dr.Histo1D({"h_npho_dr", "N quality #gamma (#DeltaR<1);N_{#gamma};Events", 20, -0.5, 19.5}, "n_quality_photons");
    auto h_nele_dr = df_dr.Histo1D({"h_nele_dr", "N quality e (#DeltaR<1);N_{e};Events",      20, -0.5, 19.5}, "n_quality_electrons");
    auto h_njet_dr = df_dr.Histo1D({"h_njet_dr", "N quality jet (#DeltaR<1);N_{jet};Events",  20, -0.5, 19.5}, "n_quality_jets");

    // Per-muon distributions at each cutflow step (fills one entry per muon)
    // -- After trigger --
    auto h_mu_eta_trig      = df_trigger.Histo1D({"h_mu_eta_trig",      "#eta (trigger);#eta;Muons",                   100, -3, 3},      "ScoutingMuonVtx_eta");
    auto h_mu_trackIso_trig = df_trigger.Histo1D({"h_mu_trackIso_trig", "trackIso (trigger);trackIso;Muons",            100, 0, 1},       "ScoutingMuonVtx_trackIso");
    auto h_mu_normchi2_trig = df_trigger.Histo1D({"h_mu_normchi2_trig", "#chi^{2}/ndof (trigger);#chi^{2}/ndof;Muons",  100, 0, 20},      "ScoutingMuonVtx_normchi2");
    // -- After 2 quality muons --
    auto h_mu_eta_dimu      = df_dimu.Histo1D({"h_mu_eta_dimu",      "#eta (2#mu);#eta;Muons",                   100, -3, 3},      "ScoutingMuonVtx_eta");
    auto h_mu_trackIso_dimu = df_dimu.Histo1D({"h_mu_trackIso_dimu", "trackIso (2#mu);trackIso;Muons",            100, 0, 1},       "ScoutingMuonVtx_trackIso");
    auto h_mu_normchi2_dimu = df_dimu.Histo1D({"h_mu_normchi2_dimu", "#chi^{2}/ndof (2#mu);#chi^{2}/ndof;Muons",  100, 0, 20},      "ScoutingMuonVtx_normchi2");
    // -- After di-pT > 20 --
    auto h_mu_eta_dipt      = df_dipt.Histo1D({"h_mu_eta_dipt",      "#eta (di-p_{T}>20);#eta;Muons",                   100, -3, 3},      "ScoutingMuonVtx_eta");
    auto h_mu_trackIso_dipt = df_dipt.Histo1D({"h_mu_trackIso_dipt", "trackIso (di-p_{T}>20);trackIso;Muons",            100, 0, 1},       "ScoutingMuonVtx_trackIso");
    auto h_mu_normchi2_dipt = df_dipt.Histo1D({"h_mu_normchi2_dipt", "#chi^{2}/ndof (di-p_{T}>20);#chi^{2}/ndof;Muons",  100, 0, 20},      "ScoutingMuonVtx_normchi2");
    // -- After dR < 1 (final) --
    auto h_mu_eta_dr      = df_dr.Histo1D({"h_mu_eta_dr",      "#eta (#DeltaR<1);#eta;Muons",                   100, -3, 3},      "ScoutingMuonVtx_eta");
    auto h_mu_trackIso_dr = df_dr.Histo1D({"h_mu_trackIso_dr", "trackIso (#DeltaR<1);trackIso;Muons",            100, 0, 1},       "ScoutingMuonVtx_trackIso");
    auto h_mu_normchi2_dr = df_dr.Histo1D({"h_mu_normchi2_dr", "#chi^{2}/ndof (#DeltaR<1);#chi^{2}/ndof;Muons",  100, 0, 20},      "ScoutingMuonVtx_normchi2");

    // Branches to save in the skimmed output
    std::vector<std::string> branches_to_keep = {
        // Event-level
        "event", 
        "run", 
        "luminosityBlock",
        // L1 Triggers
        "L1_DoubleMu_12_5", 
        "L1_DoubleMu_15_7",
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", 
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",
        "L1_DoubleMu4_SQ_OS_dR_Max1p2", 
        "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        // Muons
        "nScoutingMuonVtx",
        "ScoutingMuonVtx_pt", 
        "ScoutingMuonVtx_eta",
        "ScoutingMuonVtx_phi",
        "ScoutingMuonVtx_m",
        "ScoutingMuonVtx_charge",
        "ScoutingMuonVtx_trackIso",
        "ScoutingMuonVtx_nValidPixelHits",
        "ScoutingMuonVtx_nTrackerLayersWithMeasurement",
        "ScoutingMuonVtx_trk_chi2",
        "ScoutingMuonVtx_trk_ndof",
        "ScoutingMuonVtx_normchi2",
        // MET
        "ScoutingMET_pt",
        "ScoutingMET_phi",
        // Electrons
        "nScoutingElectron",
        "ScoutingElectron_pt",
        "ScoutingElectron_eta",
        "ScoutingElectron_phi",
        "ScoutingElectron_m",
        "ScoutingElectron_hOverE",
        "ScoutingElectron_sigmaIetaIeta",
        "ScoutingElectron_dEtaIn",
        "ScoutingElectron_dPhiIn",
        "ScoutingElectron_hcalIso",
        "ScoutingElectron_ecalIso",
        "ScoutingElectron_trackIso",
        "ScoutingElectron_bestTrack_chi2overndf",
        // Photons
        "nScoutingPhoton",
        "ScoutingPhoton_pt",
        "ScoutingPhoton_eta",
        "ScoutingPhoton_phi",
        "ScoutingPhoton_m",
        "ScoutingPhoton_hOverE",
        "ScoutingPhoton_sigmaIetaIeta",
        "ScoutingPhoton_ecalIso",
        "ScoutingPhoton_hcalIso",
        // Reclustered PF Jets
        "nScoutingPFJetRecluster",
        "ScoutingPFJetRecluster_pt",
        "ScoutingPFJetRecluster_eta",
        "ScoutingPFJetRecluster_phi",
        "ScoutingPFJetRecluster_mass",
        "ScoutingPFJetRecluster_nCh",
        "ScoutingPFJetRecluster_nConstituents",
        "ScoutingPFJetRecluster_nNh",
        "ScoutingPFJetRecluster_nPhotons",
        "ScoutingPFJetRecluster_neHEF",
        "ScoutingPFJetRecluster_neEmEF",
        "ScoutingPFJetRecluster_chHEF",
        "ScoutingPFJetRecluster_chEmEF",
        "ScoutingPFJetRecluster_muEF",
        "ScoutingPFJetRecluster_particleNet_prob_b",
        "ScoutingPFJetRecluster_particleNet_prob_g",
        "ScoutingPFJetRecluster_particleNet_prob_uds",
        "ScoutingPFJetRecluster_particleNet_prob_undef",
        // PF Jets
        "nScoutingPFJet",
        "ScoutingPFJet_pt",
        "ScoutingPFJet_eta",
        "ScoutingPFJet_phi",
        "ScoutingPFJet_m",
        "ScoutingPFJet_chargedHadronEnergy",
        "ScoutingPFJet_neutralHadronEnergy",
        "ScoutingPFJet_photonEnergy",
        "ScoutingPFJet_electronEnergy",
        "ScoutingPFJet_muonEnergy",
        "ScoutingPFJet_HFEMEnergy",
        "ScoutingPFJet_HFEMMultiplicity",
        // Computed columns
        "dimu_pt",
        "dimu_dr",
        "dimu_mass",
        "n_quality_muons",
        "n_quality_photons",
        "n_quality_electrons",
        "n_quality_jets",
        "n_quality_jets_loweta",
        "n_quality_jets_higheta",
    };

    // Snapshot triggers the event loop for all booked actions
    df_dr.Snapshot("Events", output_file.Data(), branches_to_keep);

    auto report = df_dr.Report();
    report->Print();

    // Save cutflow histogram
    TFile* f = TFile::Open(output_file.Data(), "UPDATE");

    TH1F* h_cutflow = new TH1F("h_cutflow", "Cutflow;Cut;Events", 5, 0, 5);
    h_cutflow->GetXaxis()->SetBinLabel(1, "All events");
    h_cutflow->GetXaxis()->SetBinLabel(2, "Trigger");
    h_cutflow->GetXaxis()->SetBinLabel(3, charge_label.Data());
    h_cutflow->GetXaxis()->SetBinLabel(4, "di-p_{T} > 20 GeV");
    h_cutflow->GetXaxis()->SetBinLabel(5, "#DeltaR < 1");

    double n_total = *count_all;
    h_cutflow->SetBinContent(1, n_total);
    h_cutflow->SetBinContent(2, *count_trigger);
    h_cutflow->SetBinContent(3, *count_dimu);
    h_cutflow->SetBinContent(4, *count_dipt);
    h_cutflow->SetBinContent(5, *count_dr);

    h_cutflow->SetMinimum(0);
    h_cutflow->SetMaximum(1.3 * n_total);
    h_cutflow->Write();

    // Write 2D mass vs cutflow step histogram with bin labels
    auto* h2_ptr = h2_mass_cutflow.GetPtr();
    h2_ptr->GetXaxis()->SetBinLabel(1, charge_label.Data());
    h2_ptr->GetXaxis()->SetBinLabel(2, "di-p_{T} > 20 GeV");
    h2_ptr->GetXaxis()->SetBinLabel(3, "#DeltaR < 1");
    h2_ptr->Write();

    // Write 1D mass histograms
    h_mass_dimu->Write();
    h_mass_dipt->Write();
    h_mass_dr->Write();

    // Write object multiplicity histograms
    std::vector<ROOT::RDF::RResultPtr<TH1D>> obj_histos = {
        h_nmu_all, h_npho_all, h_nele_all, h_njet_all,
        h_nmu_trig, h_npho_trig, h_nele_trig, h_njet_trig,
        h_nmu_dimu, h_npho_dimu, h_nele_dimu, h_njet_dimu,
        h_nmu_dipt, h_npho_dipt, h_nele_dipt, h_njet_dipt,
        h_nmu_dr, h_npho_dr, h_nele_dr, h_njet_dr,
    };
    for (auto& h : obj_histos) h->Write();

    // Write per-muon distribution histograms
    std::vector<ROOT::RDF::RResultPtr<TH1D>> mu_histos = {
        h_mu_eta_trig, h_mu_trackIso_trig, h_mu_normchi2_trig,
        h_mu_eta_dimu, h_mu_trackIso_dimu, h_mu_normchi2_dimu,
        h_mu_eta_dipt, h_mu_trackIso_dipt, h_mu_normchi2_dipt,
        h_mu_eta_dr,   h_mu_trackIso_dr,   h_mu_normchi2_dr,
    };
    for (auto& h : mu_histos) h->Write();

    f->Close();
    delete f;

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Output written to: " << output_file << "\n";
    std::cout << std::string(60, '=') << "\n";
}


/* Documentation of original branches (for reference):

\begin{table}[htbp]
\centering
\caption{ROOT branches used in the analysis, grouped by object type.}
\label{tab:root_leaves}
\begin{tabular}{ll}
\toprule
\textbf{Category} & \textbf{Branch Name} \\
\midrule
\multicolumn{2}{l}{\textit{Event-level}} \\
& {event} \\
& {run} \\
& {luminosityBlock} \\
\midrule
\multicolumn{2}{l}{\textit{L1 Triggers}} \\
& {L1_DoubleMu_12_5} \\
& {L1_DoubleMu_15_7} \\
& {L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7} \\
& {L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18} \\
& {L1_DoubleMu4_SQ_OS_dR_Max1p2} \\
& {L1_DoubleMu4p5_SQ_OS_dR_Max1p2} \\
\midrule
\multicolumn{2}{l}{\textit{Muons ({ScoutingMuonVtx})}} \\
& {nScoutingMuonVtx} \\
& {ScoutingMuonVtx_pt}, {ScoutingMuonVtx_eta}, {ScoutingMuonVtx_phi}, {ScoutingMuonVtx_m}, {ScoutingMuonVtx_charge} \\
& {ScoutingMuonVtx_trackIso} \\
& {ScoutingMuonVtx_nValidPixelHits} \\
& {ScoutingMuonVtx_nTrackerLayersWithMeasurement} \\
& {ScoutingMuonVtx_trk_chi2}, {ScoutingMuonVtx_trk_ndof}, {ScoutingMuonVtx_normchi2} \\
\midrule
\multicolumn{2}{l}{\textit{Missing Transverse Energy ({ScoutingMET})}} \\
& {ScoutingMET_pt}, {ScoutingMET_phi} \\
\midrule
\multicolumn{2}{l}{\textit{Electrons ({ScoutingElectron})}} \\
& {nScoutingElectron} \\
& {ScoutingElectron_pt}, {ScoutingElectron_eta}, {ScoutingElectron_phi}, {ScoutingElectron_m} \\
& {ScoutingElectron_hOverE}, {ScoutingElectron_sigmaIetaIeta} \\
& {ScoutingElectron_dEtaIn}, {ScoutingElectron_dPhiIn} \\
& {ScoutingElectron_hcalIso}, {ScoutingElectron_ecalIso}, {ScoutingElectron_trackIso} \\
& {ScoutingElectron_bestTrack_chi2overndf} \\
\midrule
\multicolumn{2}{l}{\textit{Photons ({ScoutingPhoton})}} \\
& {nScoutingPhoton} \\
& {ScoutingPhoton_pt}, {ScoutingPhoton_eta}, {ScoutingPhoton_phi}, {ScoutingPhoton_m} \\
& {ScoutingPhoton_hOverE}, {ScoutingPhoton_sigmaIetaIeta} \\
& {ScoutingPhoton_ecalIso}, {ScoutingPhoton_hcalIso} \\
\midrule
\multicolumn{2}{l}{\textit{Reclustered PF Jets ({ScoutingPFJetRecluster})}} \\
& {nScoutingPFJetRecluster} \\
& {ScoutingPFJetRecluster_pt}, {ScoutingPFJetRecluster_eta}, {ScoutingPFJetRecluster_phi}, {ScoutingPFJetRecluster_mass} \\
& {ScoutingPFJetRecluster_nCh}, {ScoutingPFJetRecluster_nConstituents}, {ScoutingPFJetRecluster_nNh}, {ScoutingPFJetRecluster_nPhotons} \\
& {ScoutingPFJetRecluster_neHEF}, {ScoutingPFJetRecluster_neEmEF}, {ScoutingPFJetRecluster_chHEF}, {ScoutingPFJetRecluster_chEmEF}, {ScoutingPFJetRecluster_muEF} \\
& {ScoutingPFJetRecluster_particleNet_prob_b}, {ScoutingPFJetRecluster_particleNet_prob_g} \\
& {ScoutingPFJetRecluster_particleNet_prob_uds}, {ScoutingPFJetRecluster_particleNet_prob_undef} \\
\midrule
\multicolumn{2}{l}{\textit{PF Jets ({ScoutingPFJet})}} \\
& {nScoutingPFJet} \\
& {ScoutingPFJet_pt}, {ScoutingPFJet_eta}, {ScoutingPFJet_phi}, {ScoutingPFJet_m} \\
& {ScoutingPFJet_chargedHadronEnergy}, {ScoutingPFJet_neutralHadronEnergy} \\
& {ScoutingPFJet_photonEnergy}, {ScoutingPFJet_electronEnergy}, {ScoutingPFJet_muonEnergy} \\
& {ScoutingPFJet_HFEMEnergy}, {ScoutingPFJet_HFEMMultiplicity} \\
\bottomrule
\end{tabular}
\end{table}

*/
