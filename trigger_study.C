/// Trigger study using RDataFrame for di-muon analysis
/// Based on skim_ntuples.C but only saves histograms (no event output)
///
/// Cutflow:
///   0. All events
///   1. L1 Trigger (OR of DoubleMu seeds)
///   2. 2 quality OS (sr) / SS (vr) muons
///   3. + di-pT > 20 GeV
///
/// Usage:
///       root -l -b -q 'trigger_study.C("sr","/path/to/files/")'
///       root -l -b -q 'trigger_study.C("vr", "/path/to/files/")'
///       root -l -b -q 'trigger_study.C("sr", "/path/to/files/", true)'  // validate files
///       root -l -b -q 'trigger_study.C("sr", "/path/to/files/", false, 0, 100)'  // first 100 files
///       root -l -b -q 'trigger_study.C("sr", "/path/to/files/", false, 2, 100)'  // files 200-299

void trigger_study(TString region = "sr",
                  TString base_dir = "/ceph/cms/store/user/tvami/DiMuonPlusX/",
                  bool validate = true,
                  int job_index = 0,
                  int files_per_job = 0) {  // 0 = process all files

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
    if (files_per_job > 0)
        dataset_name += TString::Format("_job%d", job_index);
    TString output_file = TString::Format("trigger_study_%s_%s.root", region.Data(), dataset_name.Data());

    // Build TChain
    TChain chain("Events");

    // Collect all filenames and sort for deterministic ordering
    TString find_cmd = TString::Format("find %s -name '*.root' 2>/dev/null | sort", base_dir.Data());
    FILE* pipe = popen(find_cmd.Data(), "r");
    if (!pipe) {
        std::cerr << "Error: cannot run find command\n";
        return;
    }

    std::cout << "Reading filenames from pipe\n";

    std::vector<std::string> all_files;
    char buffer[2048];
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        TString filename = buffer;
        filename.ReplaceAll("\n","");
        if (filename.Length() > 0) all_files.push_back(filename.Data());
    }
    pclose(pipe);

    std::cout << "Found " << all_files.size() << " total files in directory\n";

    // Select file slice for this job
    size_t start = 0, end = all_files.size();
    if (files_per_job > 0) {
        start = job_index * files_per_job;
        end = std::min(start + (size_t)files_per_job, all_files.size());
        if (start >= all_files.size()) {
            std::cerr << "Error: job_index " << job_index << " out of range (only "
                      << all_files.size() << " files, " << files_per_job << " per job).\n";
            return;
        }
        std::cout << "Processing files " << start << " to " << end - 1
                  << " (job " << job_index << ", " << files_per_job << " per job)\n";
    }

    int nfiles = 0, nzombie = 0;
    for (size_t i = start; i < end; ++i) {
        TString filename = all_files[i].c_str();
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

    if (nzombie > 0) std::cerr << "Skipped " << nzombie << " corrupted files\n";
    if (nfiles == 0) {
        std::cerr << "Error: No ROOT files found.\n";
        return;
    }
    std::cout << "Processing " << nfiles << " files\n";

    ROOT::RDataFrame df(chain);

    // --- Define quality object counts (available at all cutflow steps) ---

    // Quality muons: pT > 3 GeV, |eta| < 2.4, trackIso < 0.15, chi2/ndof < 3
    auto df_q = df
        .Define("n_qual_muons", [](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& trackIso,
                const ROOT::VecOps::RVec<float>& normchi2) {
            int n = 0;
            for (size_t i = 0; i < pt.size(); ++i)
                if (pt[i] > 3 && std::abs(eta[i]) < 2.4 && trackIso[i] < 0.15 && normchi2[i] < 3) n++;
            return n;
        }, {"ScoutingMuonVtx_pt", "ScoutingMuonVtx_eta",
            "ScoutingMuonVtx_trackIso", "ScoutingMuonVtx_normchi2"});

    // --- Cutflow ---

    // Step 1: L1 Trigger (OR of 6 DoubleMu seeds)
    auto df_trigger = df_q.Filter(
        "L1_DoubleMu_12_5 || L1_DoubleMu_15_7 || "
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 || L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18 || "
        "L1_DoubleMu4_SQ_OS_dR_Max1p2 || L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        "L1 Trigger");

    // Steps 2-3: Find best qualifying di-muon pair per event
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

    // Define count actions
    auto count_all     = df_q.Count();
    auto count_trigger = df_trigger.Count();
    auto count_dimu    = df_dimu.Count();
    auto count_dipt    = df_dipt.Count();

    // --- Book histograms ---

    // 2D: invariant mass at each cutflow step (steps 2-3, where a pair exists)
    auto df_mass2d = df_dimu
        .Define("mass_cf_steps", [](float dipt) {
            ROOT::VecOps::RVec<double> s = {0.};
            if (dipt > 20) s.push_back(1.);
            return s;
        }, {"dimu_pt"})
        .Define("mass_cf_vals", [](float mass, float dipt) {
            ROOT::VecOps::RVec<double> v = {(double)mass};
            if (dipt > 20) v.push_back(mass);
            return v;
        }, {"dimu_mass", "dimu_pt"});
    auto h2_mass_cutflow = df_mass2d.Histo2D(
        {"h2_mass_cutflow",
         "Di-muon mass at cutflow steps;Cutflow step;m_{#mu#mu} [GeV]",
         2, -0.5, 1.5, 200, 0, 200},
        "mass_cf_steps", "mass_cf_vals");

    // 1D: invariant mass at each cutflow step
    auto h_mass_dimu = df_dimu.Histo1D({"h_mass_dimu", "m_{#mu#mu} (2 quality #mu);m_{#mu#mu} [GeV];Events", 200, 0, 200}, "dimu_mass");
    auto h_mass_dipt = df_dipt.Histo1D({"h_mass_dipt", "m_{#mu#mu} (di-p_{T} > 20);m_{#mu#mu} [GeV];Events", 200, 0, 200}, "dimu_mass");

    // --- Trigger efficiency histograms ---
    // Measured in events passing orthogonal HTT triggers
    auto df_trig_eff = df_q
        .Filter(
            "L1_HTT120_SingleLLPJet40 || "
            "L1_HTT120er || "
            "L1_HTT160_SingleLLPJet50 || "
            "L1_HTT160er || "
            "L1_HTT200_SingleLLPJet60 || "
            "L1_HTT200er || "
            "L1_HTT240_SingleLLPJet70 || "
            "L1_HTT255er || "
            "L1_HTT280er || L1_HTT280er_QuadJet_70_55_40_35_er2p5 || "
            "L1_HTT320er || L1_HTT320er_QuadJet_70_55_40_40_er2p5 || "
            "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 || "
            "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 || "
            "L1_HTT360er || "
            "L1_HTT400er || "
            "L1_HTT450er",
            "Orthogonal HTT trigger")
        .Define("trigger_or",
            "(int)(L1_DoubleMu_12_5 || L1_DoubleMu_15_7 || "
            "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 || L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18 || "
            "L1_DoubleMu4_SQ_OS_dR_Max1p2 || L1_DoubleMu4p5_SQ_OS_dR_Max1p2)")
        .Define("pretrig_dimu", [region_str](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& eta,
                const ROOT::VecOps::RVec<float>& phi,
                const ROOT::VecOps::RVec<float>& mass,
                const ROOT::VecOps::RVec<int>& charge,
                const ROOT::VecOps::RVec<float>& trackIso,
                const ROOT::VecOps::RVec<float>& normchi2) {
            ROOT::VecOps::RVec<float> result = {-1.f, -1.f, -1.f};
            std::vector<int> good;
            for (size_t i = 0; i < pt.size(); ++i) {
                if (pt[i] > 3 && std::abs(eta[i]) < 2.4 && trackIso[i] < 0.15 && normchi2[i] < 3)
                    good.push_back(i);
            }
            if (good.size() < 2) return result;
            float best_dipt = -1;
            int best_i = -1, best_j = -1;
            for (size_t a = 0; a < good.size(); ++a) {
                for (size_t b = a + 1; b < good.size(); ++b) {
                    int i = good[a], j = good[b];
                    bool os = charge[i] * charge[j] < 0;
                    if (region_str == "sr" && !os) continue;
                    if (region_str == "vr" && os) continue;
                    float px = pt[i]*std::cos(phi[i]) + pt[j]*std::cos(phi[j]);
                    float py = pt[i]*std::sin(phi[i]) + pt[j]*std::sin(phi[j]);
                    float dipt = std::sqrt(px*px + py*py);
                    if (dipt > best_dipt) {
                        best_dipt = dipt;
                        best_i = i; best_j = j;
                    }
                }
            }
            if (best_i < 0) return result;
            float lead_pt = std::max(pt[best_i], pt[best_j]);
            float sublead_pt = std::min(pt[best_i], pt[best_j]);
            return ROOT::VecOps::RVec<float>{lead_pt, sublead_pt, best_dipt};
        }, {"ScoutingMuonVtx_pt", "ScoutingMuonVtx_eta", "ScoutingMuonVtx_phi", "ScoutingMuonVtx_m",
            "ScoutingMuonVtx_charge", "ScoutingMuonVtx_trackIso", "ScoutingMuonVtx_normchi2"})
        .Define("any_muon_info", [](
                const ROOT::VecOps::RVec<float>& pt,
                const ROOT::VecOps::RVec<float>& phi) {
            ROOT::VecOps::RVec<float> result = {-1.f, -1.f, -1.f};
            if (pt.empty()) return result;
            auto idx = ROOT::VecOps::Argsort(pt, [](float a, float b){ return a > b; });
            result[0] = pt[idx[0]];
            if (pt.size() >= 2) {
                result[1] = pt[idx[1]];
                float px = pt[idx[0]]*std::cos(phi[idx[0]]) + pt[idx[1]]*std::cos(phi[idx[1]]);
                float py = pt[idx[0]]*std::sin(phi[idx[0]]) + pt[idx[1]]*std::sin(phi[idx[1]]);
                result[2] = std::sqrt(px*px + py*py);
            }
            return result;
        }, {"ScoutingMuonVtx_pt", "ScoutingMuonVtx_phi"})
        .Define("any_leading_mu_pt",    "any_muon_info[0]")
        .Define("any_subleading_mu_pt", "any_muon_info[1]")
        .Define("any_dimu_pt",          "any_muon_info[2]")
        .Define("leading_mu_pt",    "pretrig_dimu[0]")
        .Define("subleading_mu_pt", "pretrig_dimu[1]")
        .Define("pretrig_dipt",     "pretrig_dimu[2]")
        // Individual trigger flags
        .Define("trig_DoubleMu_12_5",                  "(int)(L1_DoubleMu_12_5)")
        .Define("trig_DoubleMu_15_7",                  "(int)(L1_DoubleMu_15_7)")
        .Define("trig_DoubleMu4p5er2p0_SQ_OS_Mass_Min7",  "(int)(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7)")
        .Define("trig_DoubleMu4p5er2p0_SQ_OS_Mass_7to18", "(int)(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18)")
        .Define("trig_DoubleMu4_SQ_OS_dR_Max1p2",         "(int)(L1_DoubleMu4_SQ_OS_dR_Max1p2)")
        .Define("trig_DoubleMu4p5_SQ_OS_dR_Max1p2",       "(int)(L1_DoubleMu4p5_SQ_OS_dR_Max1p2)");

    // Trigger names and column names for looping
    struct TrigInfo { std::string col; std::string label; };
    std::vector<TrigInfo> triggers = {
        {"trigger_or",                                "DoubleMu Triggers OR"},
        {"trig_DoubleMu_12_5",                        "DoubleMu_12_5"},
        {"trig_DoubleMu_15_7",                        "DoubleMu_15_7"},
        {"trig_DoubleMu4p5er2p0_SQ_OS_Mass_Min7",     "DoubleMu4p5er2p0_SQ_OS_Mass_Min7"},
        {"trig_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",    "DoubleMu4p5er2p0_SQ_OS_Mass_7to18"},
        {"trig_DoubleMu4_SQ_OS_dR_Max1p2",            "DoubleMu4_SQ_OS_dR_Max1p2"},
        {"trig_DoubleMu4p5_SQ_OS_dR_Max1p2",          "DoubleMu4p5_SQ_OS_dR_Max1p2"},
    };

    // 1D histogram: which individual triggers fire (one entry per passing trigger per event)
    auto df_trig_bins = df_trig_eff
        .Define("trig_fired_bins", [](int t0, int t1, int t2, int t3, int t4, int t5) {
            ROOT::VecOps::RVec<double> bins = {0.};  // "All" bin always filled
            if (t0) bins.push_back(1.);
            if (t1) bins.push_back(2.);
            if (t2) bins.push_back(3.);
            if (t3) bins.push_back(4.);
            if (t4) bins.push_back(5.);
            if (t5) bins.push_back(6.);
            return bins;
        }, {"trig_DoubleMu_12_5", "trig_DoubleMu_15_7",
            "trig_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", "trig_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",
            "trig_DoubleMu4_SQ_OS_dR_Max1p2", "trig_DoubleMu4p5_SQ_OS_dR_Max1p2"});
    auto h_trig_fired = df_trig_bins.Histo1D(
        {"h_trig_fired", "Individual L1 trigger counts;Trigger;Events", 7, -0.5, 6.5},
        "trig_fired_bins");

    // Book 2D trigger efficiency histograms for each trigger (leading pT, subleading pT, di-pT)
    // -- Inclusive: any muons, no quality/charge cuts --
    auto df_trig_eff_1mu = df_trig_eff.Filter("any_leading_mu_pt >= 0", ">=1 muon");
    auto df_trig_eff_2mu_any = df_trig_eff.Filter("any_subleading_mu_pt >= 0", ">=2 muons");

    std::vector<ROOT::RDF::RResultPtr<TH2D>> trigeff_histos;
    for (auto& t : triggers) {
        trigeff_histos.push_back(df_trig_eff_1mu.Histo2D(
            {("h2_trigeff_leading_" + t.label).c_str(),
             (t.label + " vs leading #mu p_{T};p_{T}^{lead} [GeV];" + t.label).c_str(),
             100, 0, 100, 2, -0.5, 1.5},
            "any_leading_mu_pt", t.col));
        trigeff_histos.push_back(df_trig_eff_2mu_any.Histo2D(
            {("h2_trigeff_subleading_" + t.label).c_str(),
             (t.label + " vs subleading #mu p_{T};p_{T}^{sublead} [GeV];" + t.label).c_str(),
             100, 0, 100, 2, -0.5, 1.5},
            "any_subleading_mu_pt", t.col));
        trigeff_histos.push_back(df_trig_eff_2mu_any.Histo2D(
            {("h2_trigeff_dipt_" + t.label).c_str(),
             (t.label + " vs di-#mu p_{T};p_{T}^{#mu#mu} [GeV];" + t.label).c_str(),
             100, 0, 200, 2, -0.5, 1.5},
            "any_dimu_pt", t.col));
    }

    // -- With quality + charge: 2 quality muons with OS (sr) / SS (vr) --
    auto df_trig_eff_qual = df_trig_eff.Filter("pretrig_dimu[0] >= 0", "Quality OS/SS pair");

    std::vector<ROOT::RDF::RResultPtr<TH2D>> trigeff_histos_2mu;
    for (auto& t : triggers) {
        trigeff_histos_2mu.push_back(df_trig_eff_qual.Histo2D(
            {("h2_trigeff_leading_2mu_" + t.label).c_str(),
             (t.label + " vs leading #mu p_{T} (qual);p_{T}^{lead} [GeV];" + t.label).c_str(),
             100, 0, 100, 2, -0.5, 1.5},
            "leading_mu_pt", t.col));
        trigeff_histos_2mu.push_back(df_trig_eff_qual.Histo2D(
            {("h2_trigeff_subleading_2mu_" + t.label).c_str(),
             (t.label + " vs subleading #mu p_{T} (qual);p_{T}^{sublead} [GeV];" + t.label).c_str(),
             100, 0, 100, 2, -0.5, 1.5},
            "subleading_mu_pt", t.col));
        trigeff_histos_2mu.push_back(df_trig_eff_qual.Histo2D(
            {("h2_trigeff_dipt_2mu_" + t.label).c_str(),
             (t.label + " vs di-#mu p_{T} (qual);p_{T}^{#mu#mu} [GeV];" + t.label).c_str(),
             100, 0, 200, 2, -0.5, 1.5},
            "pretrig_dipt", t.col));
    }

    // --- Write histograms only (no event Snapshot) ---

    TFile* f = TFile::Open(output_file.Data(), "RECREATE");

    // Cutflow histogram
    TH1F* h_cutflow = new TH1F("h_cutflow", "Cutflow;Cut;Events", 4, 0, 4);
    h_cutflow->GetXaxis()->SetBinLabel(1, "All events");
    h_cutflow->GetXaxis()->SetBinLabel(2, "Trigger");
    h_cutflow->GetXaxis()->SetBinLabel(3, charge_label.Data());
    h_cutflow->GetXaxis()->SetBinLabel(4, "di-p_{T} > 20 GeV");

    double n_total = *count_all;
    h_cutflow->SetBinContent(1, n_total);
    h_cutflow->SetBinContent(2, *count_trigger);
    h_cutflow->SetBinContent(3, *count_dimu);
    h_cutflow->SetBinContent(4, *count_dipt);

    h_cutflow->SetMinimum(0);
    h_cutflow->SetMaximum(1.3 * n_total);
    h_cutflow->Write();

    // Write 2D mass vs cutflow step histogram with bin labels
    auto* h2_ptr = h2_mass_cutflow.GetPtr();
    h2_ptr->GetXaxis()->SetBinLabel(1, charge_label.Data());
    h2_ptr->GetXaxis()->SetBinLabel(2, "di-p_{T} > 20 GeV");
    h2_ptr->Write();

    // Write 1D mass histograms
    h_mass_dimu->Write();
    h_mass_dipt->Write();

    // Write individual trigger fired histogram
    auto* h_tf = h_trig_fired.GetPtr();
    h_tf->GetXaxis()->SetBinLabel(1, "All");
    h_tf->GetXaxis()->SetBinLabel(2, "DoubleMu_12_5");
    h_tf->GetXaxis()->SetBinLabel(3, "DoubleMu_15_7");
    h_tf->GetXaxis()->SetBinLabel(4, "DoubleMu4p5er2p0_SQ_OS_Mass_Min7");
    h_tf->GetXaxis()->SetBinLabel(5, "DoubleMu4p5er2p0_SQ_OS_Mass_7to18");
    h_tf->GetXaxis()->SetBinLabel(6, "DoubleMu4_SQ_OS_dR_Max1p2");
    h_tf->GetXaxis()->SetBinLabel(7, "DoubleMu4p5_SQ_OS_dR_Max1p2");
    h_tf->Write();

    // Write trigger efficiency histograms (inclusive)
    for (auto& h : trigeff_histos) h->Write();

    // Write trigger efficiency histograms (>= 2 quality muons)
    for (auto& h : trigeff_histos_2mu) h->Write();

    f->Close();
    delete f;

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Output written to: " << output_file << "\n";
    std::cout << std::string(60, '=') << "\n";
}
