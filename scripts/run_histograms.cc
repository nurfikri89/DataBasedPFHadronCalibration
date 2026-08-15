// ROOT RDataFrame -kirjasto nopeaan, lazy datan käsittelyyn
#include "ROOT/RDataFrame.hxx"
// ROOTin histogrammi- ja tiedostokirjastot
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
// C++ standardikirjastoja tiedostonkäsittelyyn, tulostukseen, datarakenteisiin ja ajastukseen
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <tuple>
#include <chrono>
#include <cmath>        // std::abs
#include "TH1D.h"       // käytät TH1D:tä
#include "TSystem.h"    // gSystem->Exec
#include "TString.h"    // Form(...)
#include "ROOT/RVec.hxx"

#include <nlohmann/json.hpp>
#include <map>
#include <vector>
#include <fstream>
#include "vetomap_utils.h"
// #include "./scripts/vetomap_utils.h"
#include "ROOT/RDFHelpers.hxx"

// golden JSON -apufunktiot
std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> LoadGoldenJSON(const std::string &filename) {
    std::ifstream f(filename);
    nlohmann::json j;
    f >> j;
    std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> goodLumis;
    for (auto &[runStr, ranges] : j.items()) {
        unsigned int run = std::stoul(runStr);
        for (auto &range : ranges)
            goodLumis[run].push_back({range[0], range[1]});
    }
    return goodLumis;
}

bool PassesGoldenJSON(unsigned int run, unsigned int lumi,
                      const std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> &goodLumis) {
    auto it = goodLumis.find(run);
    if (it == goodLumis.end()) return false;
    for (auto &r : it->second)
        if (lumi >= r.first && lumi <= r.second) return true;
    return false;
}

// MÄÄRITELLÄÄN BINITYS pionien momentille (p = pT * cosh(η))
// tarkka binitys pienillä p-arvioilla auttaa havaitsemaan responsen muutoksia
const double trkPBins[] = {
  3.0, 3.3, 3.6, 3.9, 4.2, 4.6, 5.0, 5.5, 6.0, 6.6, 7.2, 7.9, 8.6,
  9.3, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,16.5, 18.0, 20.0, 22.0, 25.0,
  30.0, 35.0, 40.0, 45.0, 55.0, 65.0, 80.0, 100.0
};
const int nTrkPBins = sizeof(trkPBins)/sizeof(double) - 1;


/* 100 bins of width 1 from 0 to 100 -> 101 edges */
static const double trkP100Bins[] = {
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
  30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
  40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
  50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
  60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
  70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
  80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
  90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
  100
};
const int nTrkP100Bins = sizeof(trkP100Bins)/sizeof(double) - 1;

/* 50 bins of width 0.1 from -2.5 to 2.5 -> 51 edges */
static const double trkEta50Bins[] = {
  -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6,
  -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6,
  -0.5, -0.4, -0.3, -0.2, -0.1,  0.0,  0.1,  0.2,  0.3,  0.4,
   0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,
   1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4,
   2.5
};
const int nTrkEta50Bins = sizeof(trkEta50Bins)/sizeof(double) - 1;

/* 50 bins of width 0.05 from 0 to 2.5 -> 51 edges */
static const double trkEP50Bins[51] = {
  0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
  0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
  1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45,
  1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95,
  2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45,
  2.50
};
const int nTrkEP50Bins = sizeof(trkEP50Bins)/sizeof(double) - 1;

const int nphi = 72;
const double vphi[nphi+1] = {
  -3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443,
  -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658,
  -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873,
  -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
  0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785,
  0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658,
  1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531,
  2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142
};

double etabins[] = {
  -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305,
  -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
  0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
  1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5
};
const int netabins = sizeof(etabins)/sizeof(etabins[0])-1;

double absetabins[] = {
  0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
  1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5
};
const int nabsetabins = sizeof(absetabins)/sizeof(absetabins[0])-1;



// helper: käännetään /store/... XRootD-URL:ksi (HIP SE). Jätä valmiit root://-URLit sellaisenaan
std::string NormalizeInputPath(const std::string &s) {
  if (s.rfind("root://", 0) == 0) return s;                // jo valmis
  if (s.rfind("/store/user/nbinnorj", 0) == 0) return "root://hip-cms-se.csc.fi/" + s;  // HUOM tupla-slash
  if (s.rfind("/eos/cms/store/", 0) == 0) return "root://eoscms.cern.ch/" + s;  // HUOM tupla-slash
  return s;
}


// PÄÄFUNKTIO: tekee histogrammit ja kirjoittaa ne tiedostoon
// void run_histograms(const char* listfile, const char* tag, const char* outpath) {
void run_histograms(const char* listfile, const char* tag, const char* outpath = "") {


  // Lataa veto map (aja kerran per prosessi)
  // try {
  //   veto::Load("/eos/user/m/mmarjama/my_pion_analysis/scripts/JetVetoMap.root", "h_jetVetoMap");
  //   std::cout << "[VETO] Loaded veto map." << std::endl;
  // } catch (const std::exception& e) {
  //   std::cerr << e.what() << "\n[VETO] Proceeding without veto (map missing)." << std::endl;
  // }

  // aloitetaan ajanotto ajon kestolle
  auto start = std::chrono::high_resolution_clock::now();

  // monisäikeisyys
  ROOT::EnableImplicitMT(4);

  // LUETAAN ROOT-tiedostojen polut listasta (robusti)
  std::ifstream infile(listfile);
  if (!infile.is_open()) {
    std::cerr << "[run_histograms] ERROR: cannot open list file: " << listfile << "\n";
    return;
  }
  std::string line;
  std::vector<std::string> files;
  files.reserve(1000);
  while (std::getline(infile, line)) {
    if (!line.empty() && line[0] == '#') continue;
    auto norm = NormalizeInputPath(line);
    if (!norm.empty())
      files.push_back(norm);
  }

  // std::cout << "Processing " << files.size() << " files..." << std::endl;

  // // Luo RDataFrame
  // ROOT::RDataFrame df("Events", files);


  std::cout << "Checking files for validity..." << std::endl;

  // Luo uusi lista vain ehjistä tiedostoista
  std::vector<std::string> goodFiles;
  std::vector<std::string> badFiles;

  for (const auto& f : files) {
    goodFiles.push_back(f);
  }

  // for (const auto& f : files) {
  //   try {
  //     // Yritetään avata tiedosto vain testimielessä
  //     std::unique_ptr<TFile> testFile(TFile::Open(f.c_str()));
  //     if (!testFile || testFile->IsZombie() || testFile->GetNkeys() == 0) {
  //       badFiles.push_back(f);
  //       continue;
  //     }
  //     goodFiles.push_back(f);
  //   } catch (...) {
  //     badFiles.push_back(f);
  //   }
  // }

  // Tulostetaan yhteenveto
  std::cout << "Good files: " << goodFiles.size() << " / " << files.size() << std::endl;
  if (!badFiles.empty()) {
    std::cerr << "[Warning] Skipping " << badFiles.size() << " bad (zombie) files:\n";
    for (const auto& bf : badFiles) std::cerr << "  " << bf << std::endl;
  }

  if (goodFiles.empty()) {
    std::cerr << "[ERROR] No valid input files found, aborting.\n";
    return;
  }

  std::cout << "Creating RDataFrame with " << goodFiles.size() << " valid files..." << std::endl;
  TChain treeEvents("Events");
  for (const auto& f : goodFiles) {
    treeEvents.Add(f.c_str());
  }

  // Luo RDataFrame vain ehjistä tiedostoista
  // ROOT::RDataFrame df("Events", goodFiles);

  // keräysvektorit
  std::vector<ROOT::RDF::RResultPtr<TH1D>>     h1_histos;
  std::vector<ROOT::RDF::RResultPtr<TH2D>>     h2_histos;
  std::vector<ROOT::RDF::RResultPtr<TProfile>> profiles;
  std::vector<ROOT::RDF::RResultPtr<TH3D>>     h3_histos;

  ROOT::RDataFrame df(treeEvents);
  ROOT::RDF::Experimental::AddProgressBar(df);

  // // --- Golden JSON (vain datalle) ---
  // auto golden = LoadGoldenJSON("/eos/user/c/cmsdqm/www/CAF/certification/Collisions25/Cert_Collisions2025_391658_397595_Golden.json");

  // auto df_json = df.Filter([&golden](unsigned int run, unsigned int lumi) {
  //     return PassesGoldenJSON(run, lumi, golden);
  // }, {"run", "luminosityBlock"}, "Golden JSON filter");

  // --- Päätä automaattisesti, tarvitaanko Golden JSON -filtteriä ---
  // Logiikka:
  //  - Data = polussa "ZeroBias"
  //  - 2025 = polussa "2025"
  //  - JSON filtteriä käytetään VAIN, jos (Data && 2025)
  //  - MC:lle ja 2024-datalle EI Golden JSONia

  bool applyGoldenJSON = false;

  if (!goodFiles.empty()) {
      const std::string& sample = goodFiles.front();

      bool isData = (sample.find("ZeroBias") != std::string::npos);
      bool is2025 = (sample.find("2025")    != std::string::npos);

      if (isData && is2025) {
          // applyGoldenJSON = true;
          applyGoldenJSON = false;
      }
  }

  ROOT::RDF::RNode df_json(df);

  if (applyGoldenJSON) {
    std::cout << "[INFO] Applying 2025 Golden JSON filter..." << std::endl;

    auto golden = LoadGoldenJSON(
      "/eos/user/c/cmsdqm/www/CAF/certification/Collisions25/"
      "Cert_Collisions2025_391658_397595_Golden.json"
    );

    df_json = df.Filter(
      [&golden](unsigned int run, unsigned int lumi) {
          return PassesGoldenJSON(run, lumi, golden);
      },
      {"run", "luminosityBlock"},
      "Golden JSON filter"
    );
  } else {
    std::cout << "[INFO] Skipping Golden JSON filter (MC or non-2025 data)" << std::endl;
    // identiteettifiltteri, selkeäksi lokiin
    df_json = df.Filter("1", "No JSON filter");
  }

  // --- Event cleaning flags (standard CMS “tight” set) ---
  auto df_clean = df_json.Filter(
    "Flag_goodVertices && "
    "Flag_globalSuperTightHalo2016Filter && "
    "Flag_EcalDeadCellTriggerPrimitiveFilter && "
    "Flag_BadPFMuonFilter && "
    "Flag_BadPFMuonDzFilter && "
    "Flag_hfNoisyHitsFilter && "
    "Flag_eeBadScFilter && "
    "Flag_ecalBadCalibFilter",
    "Event cleaning flags"
  );

  auto df_eta_precut = df_clean
    .Define("etaMaskBarrel", "abs(IsoChHadPFCand_eta) < 1.3");
    .Define("etaMaskEndcap", "not(etaMaskBarrel) && abs(IsoChHadPFCand_eta) < 2.5");
    .Define("etaMask",       "etaMaskBarrel || etaMaskEndcap")

  // --- Tästä eteenpäin jatkuu normaali analyysiketju ---
  auto df_eta = df_eta_precut
    // Isoloidut pionit
    .Filter("Sum(isoMask) > 0", "Has isolated charged pions")

  auto df_iso = df_eta
    // Isoloitujen kandidaattien η ja φ (tarvitaan veto-maskiin)
    .Define("etaIso",     "IsoChHadPFCand_eta[etaMask]")
    .Define("phiIso",     "IsoChHadPFCand_phi[etaMask]")
    .Define("ptIso",      "IsoChHadPFCand_pt[etaMask]")
    .Define("pIso",       "ptIso * cosh(etaIso)")
    .Define("hcalIso",    "IsoChHadPFCand_hcalEnergy[etaMask]")
    .Define("ecalIso",    "IsoChHadPFCand_ecalEnergy[etaMask]")
    .Define("rawEcalIso", "IsoChHadPFCand_rawEcalEnergy[etaMask]")
    .Define("rawHcalIso", "IsoChHadPFCand_rawHcalEnergy[etaMask]")
    // --- Summat ja E/p (vartioi nollajaot) ---
    .Define("detIsoE",     "ecalIso + hcalIso")
    .Define("rawDetIsoE",  "rawEcalIso + rawHcalIso")
    .Define("ep_def",      "ROOT::VecOps::Where(pIso > 0, (ecalIso + hcalIso)/pIso, -1.0)")
    .Define("ep_raw",      "ROOT::VecOps::Where(pIso > 0, (rawEcalIso + rawHcalIso)/pIso, -1.0)")
    // HCAL E/E_raw ja leikkausmaski
    .Define("hcalCorrFactor", "ROOT::VecOps::Where(rawHcalIso > 0, hcalIso/rawHcalIso, -1.0)")
    .Define("passHCAL",       "hcalCorrFactor > 0.9")
    //
    .Define("fracEcal", "ecalIso / pIso")
    .Define("fracHcal", "hcalIso / pIso")
    // Hadronityypit: JIT tulkitsee tämän elementtikohtaiseksi operaatioksi
    .Define("isHadH",    "(fracEcal <= 0)  && (fracHcal > 0.f)")
    .Define("isHadEH",   "(fracEcal > 0.f) && (fracEcal * pIso >= 1.f) && (fracHcal > 0)")
    .Define("isHadMIP",  "(fracEcal > 0.f) && (fracEcal * pIso < 1.f)  && (fracHcal > 0)")
    .Define("isHadE",    "(fracEcal > 0.f) && (fracHcal <= 0)")
    //
    .Define("HadHTrk_ptIso",       "ptIso[isHadH]")
    .Define("HadHTrk_pIso",        "pIso[isHadH]")
    .Define("HadHTrk_etaIso",      "etaIso[isHadH]")
    .Define("HadEHTrk_ptIso",      "ptIso[isHadEH]")
    .Define("HadEHTrk_pIso",       "pIso[isHadEH]")
    .Define("HadEHTrk_etaIso",     "etaIso[isHadEH]")
    .Define("HadMIPTrk_ptIso",     "ptIso[isHadMIP]")
    .Define("HadMIPTrk_pIso",      "pIso[isHadMIP]")
    .Define("HadMIPTrk_etaIso",    "etaIso[isHadMIP]")
    .Define("HadETrk_ptIso",       "ptIso[isHadE]")
    .Define("HadETrk_pIso",        "pIso[isHadE]")
    .Define("HadETrk_etaIso",      "etaIso[isHadE]")
    //
    .Define("BarrelTrk_etaIso",          "etaIso[etaMaskBarrel]")
    .Define("BarrelTrk_phiIso",          "phiIso[etaMaskBarrel]")
    .Define("BarrelTrk_ptIso",           "ptIso[etaMaskBarrel]")
    .Define("BarrelTrk_pIso",            "pIso[etaMaskBarrel]")
    .Define("BarrelTrk_ep_def",          "ep_def[etaMaskBarrel]")
    .Define("BarrelTrk_ep_raw",          "ep_raw[etaMaskBarrel]")
    .Define("BarrelTrk_hcalCorrFactor",  "hcalCorrFactor[etaMaskBarrel]")
    .Define("BarrelTrk_passHCAL",        "passHCAL[etaMaskBarrel]")
    .Define("BarrelTrk_isHadH",          "isHadH[etaMaskBarrel]")
    .Define("BarrelTrk_isHadEH",         "isHadEH[etaMaskBarrel]")
    .Define("BarrelTrk_isHadMIP",        "isHadMIP[etaMaskBarrel]")
    .Define("BarrelTrk_isHadE",          "isHadE[etaMaskBarrel]")
    //
    .Define("EndcapTrk_etaIso",          "etaIso[etaMaskEndcap]")
    .Define("EndcapTrk_phiIso",          "phiIso[etaMaskEndcap]")
    .Define("EndcapTrk_ptIso",           "ptIso[etaMaskEndcap]")
    .Define("EndcapTrk_pIso",            "pIso[etaMaskEndcap]")
    .Define("EndcapTrk_ep_def",          "ep_def[etaMaskEndcap]")
    .Define("EndcapTrk_ep_raw",          "ep_raw[etaMaskEndcap]")
    .Define("EndcapTrk_hcalCorrFactor",  "hcalCorrFactor[etaMaskEndcap]")
    .Define("EndcapTrk_passHCAL",        "passHCAL[etaMaskEndcap]")
    .Define("EndcapTrk_isHadH",          "isHadH[etaMaskEndcap]")
    .Define("EndcapTrk_isHadEH",         "isHadEH[etaMaskEndcap]")
    .Define("EndcapTrk_isHadMIP",        "isHadMIP[etaMaskEndcap]")
    .Define("EndcapTrk_isHadE",          "isHadE[etaMaskEndcap]")

    // Raw E/p maskattuna hadronityypeittäin (1D-jakaumia varten)
    .Define("rawEpIso_H",   "ROOT::VecOps::Where(isHadH,   ep_raw_all, -1.0)")
    .Define("rawEpIso_E",   "ROOT::VecOps::Where(isHadE,   ep_raw_all, -1.0)")
    .Define("rawEpIso_MIP", "ROOT::VecOps::Where(isHadMIP, ep_raw_all, -1.0)")
    .Define("rawEpIso_EH",  "ROOT::VecOps::Where(isHadEH,  ep_raw_all, -1.0)");

    // .Define("vetoMask", "veto::Mask(etaIso, phiIso)")
    // Vetomaski (1 => hylätään)
    // Käytä tästä eteenpäin vain hyväksyttyjä: isoMask && !vetoMask
    // .Define("keepMask", R"(
    //   ROOT::VecOps::RVec<char> keep(isoMask.size(), false);
    //   size_t j = 0;
    //   for (size_t i = 0; i < isoMask.size(); ++i) {
    //     if (isoMask[i]) {
    //       // Jos tämä iso-kandidaatti ei ole veto-alueella, pidä se
    //       keep[i] = (j < vetoMask.size()) ? !vetoMask[j] : true;
    //       ++j;
    //     }
    //   }
    //   return keep;
    // )")


    // Käytetään jatkossa vain hyväksyttyjä kandidaatteja




    // 2) Hadronimaskit + niistä johdetut vektorit (yhtenä ketjuna)


// 2b) Per-tyyppiset maskit HCAL-cutilla + viipaleet (X=E/p default, Y=HCALcorr, Z=p)
// auto df_iso2 = df_iso
//   // yhdistelmämaskit: hadronityyppi AND passHCAL
//   .Define("mask_H_cut",  "isHadH   && passHCAL")
//   .Define("mask_E_cut",  "isHadE   && passHCAL")
//   .Define("mask_MIP_cut","isHadMIP && passHCAL")
//   .Define("mask_EH_cut", "isHadEH  && passHCAL")

//   // viipaleet: X=E/p (default), Y=HCALcorr, Z=p  (HCAL cut)
//   .Define("ep_def_H_cut",   "ep_def_all[mask_H_cut]")
//   .Define("ep_def_E_cut",   "ep_def_all[mask_E_cut]")
//   .Define("ep_def_MIP_cut", "ep_def_all[mask_MIP_cut]")
//   .Define("ep_def_EH_cut",  "ep_def_all[mask_EH_cut]")

//   .Define("hcal_H_cut",     "hcalCorrFactor[mask_H_cut]")
//   .Define("hcal_E_cut",     "hcalCorrFactor[mask_E_cut]")
//   .Define("hcal_MIP_cut",   "hcalCorrFactor[mask_MIP_cut]")
//   .Define("hcal_EH_cut",    "hcalCorrFactor[mask_EH_cut]")

//   .Define("p_H_cut",        "pIso[mask_H_cut]")
//   .Define("p_E_cut",        "pIso[mask_E_cut]")
//   .Define("p_MIP_cut",      "pIso[mask_MIP_cut]")
//   .Define("p_EH_cut",       "pIso[mask_EH_cut]");

// Per-tyyppiset maskit ILMAN HCAL-cuttia + viipaleet (X=E/p default, Y=HCALcorr, Z=p)
// auto df_iso_nocut = df_iso
//   // yhdistelmämaskit: vain hadronityyppi (ei passHCAL)
//   .Define("mask_H_all",   "isHadH")
//   .Define("mask_E_all",   "isHadE")
//   .Define("mask_MIP_all", "isHadMIP")
//   .Define("mask_EH_all",  "isHadEH")

//   // viipaleet: X=E/p (default), Y=HCALcorr, Z=p  (ilman leikkausta)
//   .Define("ep_def_H_all",   "ep_def_all[mask_H_all]")
//   .Define("ep_def_E_all",   "ep_def_all[mask_E_all]")
//   .Define("ep_def_MIP_all", "ep_def_all[mask_MIP_all]")
//   .Define("ep_def_EH_all",  "ep_def_all[mask_EH_all]")

//   .Define("hcal_H_all",     "hcalCorrFactor[mask_H_all]")
//   .Define("hcal_E_all",     "hcalCorrFactor[mask_E_all]")
//   .Define("hcal_MIP_all",   "hcalCorrFactor[mask_MIP_all]")
//   .Define("hcal_EH_all",    "hcalCorrFactor[mask_EH_all]")

//   .Define("p_H_all",        "pIso[mask_H_all]")
//   .Define("p_E_all",        "pIso[mask_E_all]")
//   .Define("p_MIP_all",      "pIso[mask_MIP_all]")
//   .Define("p_EH_all",       "pIso[mask_EH_all]");

// // RAW-versiot (no cut): vaihdetaan X-akseleiksi ep_raw_all[..]
// auto df_iso_raw_all = df_iso_nocut
//   .Define("ep_raw_H_all",   "ep_raw_all[mask_H_all]")
//   .Define("ep_raw_E_all",   "ep_raw_all[mask_E_all]")
//   .Define("ep_raw_MIP_all", "ep_raw_all[mask_MIP_all]")
//   .Define("ep_raw_EH_all",  "ep_raw_all[mask_EH_all]");

// // RAW + HCAL-cut versiot: ketjuta df_iso2:sta (mask_*_cut on määritelty siellä)
// auto df_iso_raw_cut = df_iso2
//   .Define("ep_raw_H_cut",   "ep_raw_all[mask_H_cut]")
//   .Define("ep_raw_E_cut",   "ep_raw_all[mask_E_cut]")
//   .Define("ep_raw_MIP_cut", "ep_raw_all[mask_MIP_cut]")
//   .Define("ep_raw_EH_cut",  "ep_raw_all[mask_EH_cut]");


// 1D: perusjakaumat, jotta plot_histograms.cc:n "histNames" löytyy
h1_histos.push_back(df_iso.Histo1D({"h_pt_iso",       "Isolated track p_{T};p_{T} [GeV];Events",  nTrkPBins,  trkP100Bins},      "ptIso"));
h1_histos.push_back(df_iso.Histo1D({"h_p_iso",        "Isolated track p;p [GeV];Events",          nTrkPBins,  trkP100Bins},      "pIso"));
h1_histos.push_back(df_iso.Histo1D({"h_eta_iso",      "Isolated track #eta;#eta;Events",          nTrkEta50Bins, trkEta50Bins},  "etaIso"));
h1_histos.push_back(df_iso.Histo1D({"h_ep_iso",       "Default E/p;E/p;Events",                   nTrkEP50Bins,  trkEP50Bins},   "ep_def_all"));
h1_histos.push_back(df_iso.Histo1D({"h_raw_ep_iso",   "Raw E/p;E/p;Events",                       nTrkEP50Bins,  trkEP50Bins},   "ep_raw_all"));

h2_histos.push_back(df_iso.Histo2D({
  "h2_trkpt_vs_trketa",";track pT [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "ptIso", "etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadH_trkpt_vs_trketa","HadH;track pT [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadHTrk_ptIso", "HadHTrk_etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadEH_trkpt_vs_trketa","HadEH;track pT [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadEHTrk_ptIso", "HadEHTrk_etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadMIP_trkpt_vs_trketa","HadMIP;track pT [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadMIPTrk_ptIso", "HadMIPTrk_etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadE_trkpt_vs_trketa","HadE;track pT [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadETrk_ptIso", "HadETrk_etaIso")
);

// h3_histos.push_back(df_iso.Histo3D({
//   "h3_trkpt_vs_trketa_vs_trkphi",";track pT [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "ptIso", "etaIso", "phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadH_trkpt_vs_trketa_vs_trkphi","HadH;track pT [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadHTrk_ptIso", "HadHTrk_etaIso", "HadHTrk_phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadEH_trkpt_vs_trketa_vs_trkphi","HadEH;track pT [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadEHTrk_ptIso", "HadEHTrk_etaIso", "HadEHTrk_phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadMIP_trkpt_vs_trketa_vs_trkphi","HadMIP;track pT [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadMIPTrk_ptIso", "HadMIPTrk_etaIso", "HadMIPTrk_phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadE_trkpt_vs_trketa_vs_trkphi","HadE;track pT [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadETrk_ptIso", "HadETrk_etaIso", "HadETrk_phiIso")
// );

h2_histos.push_back(df_iso.Histo2D({
  "h2_trkp_vs_trketa",";track p [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "pIso", "etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadH_trkp_vs_trketa","HadH;track p [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadHTrk_pIso", "HadHTrk_etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadEH_trkp_vs_trketa","HadEH;track p [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadEHTrk_pIso", "HadEHTrk_etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadMIP_trkp_vs_trketa","HadMIP;track p [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadMIPTrk_pIso", "HadMIPTrk_etaIso")
);
h2_histos.push_back(df_iso.Histo2D({
  "h2_HadE_trkp_vs_trketa","HadE;track p [GeV];track #eta",
  nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins},
  "HadETrk_pIso", "HadETrk_etaIso")
);


// h3_histos.push_back(df_iso.Histo3D({
//   "h3_trkp_vs_trketa_vs_trkphi",";track p [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "pIso", "etaIso", "phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadH_trkp_vs_trketa_vs_trkphi","HadH;track p [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadHTrk_pIso", "HadHTrk_etaIso", "HadHTrk_phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadEH_trkp_vs_trketa_vs_trkphi","HadEH;track p [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadEHTrk_pIso", "HadEHTrk_etaIso", "HadEHTrk_phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadMIP_trkp_vs_trketa_vs_trkphi","HadMIP;track p [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadMIPTrk_pIso", "HadMIPTrk_etaIso", "HadMIPTrk_phiIso")
// );
// h3_histos.push_back(df_iso.Histo2D({
//   "h3_HadE_trkp_vs_trketa_vs_trkphi","HadE;track p [GeV];track #eta;track #phi",
//   nTrkPBins, trkP100Bins, nTrkEta50Bins, trkEta50Bins, nphi, vphi},
//   "HadETrk_pIso", "HadETrk_etaIso", "HadETrk_phiIso")
// );

// 1D: E/p per hadronityyppi (sekä default että raw)
// df_iso = df_iso
//   .Define("defEpIso_H",   "ROOT::VecOps::Where(isHadH,   ep_def_all, -1.0)")
//   .Define("defEpIso_E",   "ROOT::VecOps::Where(isHadE,   ep_def_all, -1.0)")
//   .Define("defEpIso_MIP", "ROOT::VecOps::Where(isHadMIP, ep_def_all, -1.0)")
//   .Define("defEpIso_EH",  "ROOT::VecOps::Where(isHadEH,  ep_def_all, -1.0)");

// h1_histos.push_back(df_iso.Histo1D({"h_ep_isHadH",   "Default E/p: HCAL only;E/p;Events", 40, 0.05, 2.5}, "defEpIso_H"));
// h1_histos.push_back(df_iso.Histo1D({"h_ep_isHadE",   "Default E/p: ECAL only;E/p;Events", 40, 0.05, 2.5}, "defEpIso_E"));
// h1_histos.push_back(df_iso.Histo1D({"h_ep_isHadMIP", "Default E/p: MIP;E/p;Events",       40, 0.05, 2.5}, "defEpIso_MIP"));
// h1_histos.push_back(df_iso.Histo1D({"h_ep_isHadEH",  "Default E/p: ECAL+HCAL;E/p;Events", 40, 0.05, 2.5}, "defEpIso_EH"));

// h1_histos.push_back(df_iso.Histo1D({"h_raw_ep_isHadH",   "Raw E/p: HCAL only;E/p;Events", 40, 0.05, 2.5}, "rawEpIso_H"));
// h1_histos.push_back(df_iso.Histo1D({"h_raw_ep_isHadE",   "Raw E/p: ECAL only;E/p;Events", 40, 0.05, 2.5}, "rawEpIso_E"));
// h1_histos.push_back(df_iso.Histo1D({"h_raw_ep_isHadMIP", "Raw E/p: MIP;E/p;Events",       40, 0.05, 2.5}, "rawEpIso_MIP"));
// h1_histos.push_back(df_iso.Histo1D({"h_raw_ep_isHadEH",  "Raw E/p: ECAL+HCAL;E/p;Events", 40, 0.05, 2.5}, "rawEpIso_EH"));

// // Fraktioprofiilit: no cut
// profiles.push_back(df_iso.Profile1D({"h_frac_H_all",   "Fraction (no cut): HCAL only;p (GeV);fraction",    nTrkPBins, trkPBins}, "pIso", "isHadH_float"));
// profiles.push_back(df_iso.Profile1D({"h_frac_E_all",   "Fraction (no cut): ECAL only;p (GeV);fraction",    nTrkPBins, trkPBins}, "pIso", "isHadE_float"));
// profiles.push_back(df_iso.Profile1D({"h_frac_MIP_all", "Fraction (no cut): MIP;p (GeV);fraction",          nTrkPBins, trkPBins}, "pIso", "isHadMIP_float"));
// profiles.push_back(df_iso.Profile1D({"h_frac_EH_all",  "Fraction (no cut): ECAL+HCAL;p (GeV);fraction",    nTrkPBins, trkPBins}, "pIso", "isHadEH_float"));

// // Fraktioprofiilit: HCAL E/E_raw > 0.9 cut
// profiles.push_back(df_iso.Profile1D({"h_frac_H_cut",   "Fraction (HCAL E/E_{raw}>0.9): HCAL only;p (GeV);fraction",    nTrkPBins, trkPBins}, "pIso_cut", "isHadH_float_cut"));
// profiles.push_back(df_iso.Profile1D({"h_frac_E_cut",   "Fraction (HCAL E/E_{raw}>0.9): ECAL only;p (GeV);fraction",    nTrkPBins, trkPBins}, "pIso_cut", "isHadE_float_cut"));
// profiles.push_back(df_iso.Profile1D({"h_frac_MIP_cut", "Fraction (HCAL E/E_{raw}>0.9): MIP;p (GeV);fraction",          nTrkPBins, trkPBins}, "pIso_cut", "isHadMIP_float_cut"));
// profiles.push_back(df_iso.Profile1D({"h_frac_EH_cut",  "Fraction (HCAL E/E_{raw}>0.9): ECAL+HCAL;p (GeV);fraction",    nTrkPBins, trkPBins}, "pIso_cut", "isHadEH_float_cut"));


// 2D & profiilit: S1–S4 (nimet yhtenevät myöhempien plottien kanssa)
// S1: Raw/p, HCAL cut
// h2_histos.push_back(df_iso.Histo2D({"h2_ep_vs_p_S1_raw_cut", "S1: Raw E/p vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//                                     nTrkPBins, trkPBins, 40, 0.05, 2.25}, "pIso", "ep_raw_cut"));
// profiles .push_back(df_iso.Profile1D({"prof_ep_vs_p_S1_raw_cut", "S1: <Raw E/p> vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//                                       nTrkPBins, trkPBins}, "pIso", "ep_raw_cut"));
// h3_histos.push_back(df_iso.Histo3D({"h3_S1_raw_cut_ep_vs_hcal_p",
//                                     "S1: Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (cut);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//                                     50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_raw_cut", "hcalCorrFactor", "pIso"));

// S2: Default/p, HCAL cut
// h2_histos.push_back(df_iso.Histo2D({"h2_ep_vs_p_S2_def_cut", "S2: Default E/p vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//                                     nTrkPBins, trkPBins, 40, 0.05, 2.25}, "pIso", "ep_def_cut"));
// profiles .push_back(df_iso.Profile1D({"prof_ep_vs_p_S2_def_cut", "S2: <Default E/p> vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//                                       nTrkPBins, trkPBins}, "pIso", "ep_def_cut"));
// h3_histos.push_back(df_iso.Histo3D({"h3_S2_def_cut_ep_vs_hcal_p",
//                                     "S2: Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (cut);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//                                     50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_def_cut", "hcalCorrFactor", "pIso"));

// S3: Raw/p, no cut
// h2_histos.push_back(df_iso.Histo2D({"h2_ep_vs_p_S3_raw_all", "S3: Raw E/p vs p (no cut);p (GeV);E/p",
//                                     nTrkPBins, trkPBins, 40, 0.05, 2.25}, "pIso", "ep_raw_all"));
// profiles .push_back(df_iso.Profile1D({"prof_ep_vs_p_S3_raw_all", "S3: <Raw E/p> vs p (no cut);p (GeV);E/p",
//                                       nTrkPBins, trkPBins}, "pIso", "ep_raw_all"));
// h3_histos.push_back(df_iso.Histo3D({"h3_S3_raw_ep_vs_hcal_p",
//                                     "S3: Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p;E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//                                     50, 0.0, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_raw_all", "hcalCorrFactor", "pIso"));

// S4: Default/p, no cut
// h2_histos.push_back(df_iso.Histo2D({"h2_ep_vs_p_S4_def_all", "S4: Default E/p vs p (no cut);p (GeV);E/p",
//                                     nTrkPBins, trkPBins, 40, 0.05, 2.25}, "pIso", "ep_def_all"));
// profiles .push_back(df_iso.Profile1D({"prof_ep_vs_p_S4_def_all", "S4: <Default E/p> vs p (no cut);p (GeV);E/p",
//                                       nTrkPBins, trkPBins}, "pIso", "ep_def_all"));
// h3_histos.push_back(df_iso.Histo3D({"h3_S4_def_ep_vs_hcal_p",
//                                     "S4: Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p;E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//                                     50, 0.0, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_def_all", "hcalCorrFactor", "pIso"));

// Per-tyyppiset 3D:t: Default E/p + HCAL cut (plot_histograms.cc käyttää näitä)
// h3_histos.push_back(df_iso2.Histo3D({"h3_resp_corr_p_isHadH",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, HCAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_def_H_cut", "hcal_H_cut", "p_H_cut"));
// h3_histos.push_back(df_iso2.Histo3D({"h3_resp_corr_p_isHadE",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, ECAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_def_E_cut", "hcal_E_cut", "p_E_cut"));
// h3_histos.push_back(df_iso2.Histo3D({"h3_resp_corr_p_isHadMIP",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, MIP);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_def_MIP_cut", "hcal_MIP_cut", "p_MIP_cut"));
// h3_histos.push_back(df_iso2.Histo3D({"h3_resp_corr_p_isHadEH",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, ECAL+HCAL);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_def_EH_cut", "hcal_EH_cut", "p_EH_cut"));

// Per-tyyppiset 3D:t: Default E/p, no cut (Y 0.0–2.5)
// h3_histos.push_back(df_iso_nocut.Histo3D({"h3_resp_def_p_isHadH_all",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, HCAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_def_H_all", "hcal_H_all", "p_H_all"));
// h3_histos.push_back(df_iso_nocut.Histo3D({"h3_resp_def_p_isHadE_all",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_def_E_all", "hcal_E_all", "p_E_all"));
// h3_histos.push_back(df_iso_nocut.Histo3D({"h3_resp_def_p_isHadMIP_all",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, MIP);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_def_MIP_all", "hcal_MIP_all", "p_MIP_all"));
// h3_histos.push_back(df_iso_nocut.Histo3D({"h3_resp_def_p_isHadEH_all",
//   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL+HCAL);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_def_EH_all", "hcal_EH_all", "p_EH_all"));

// Per-tyyppiset 3D:t: Raw E/p, no cut (Y 0.0–2.5)
// h3_histos.push_back(df_iso_raw_all.Histo3D({"h3_resp_raw_p_isHadH_all",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, HCAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_raw_H_all", "hcal_H_all", "p_H_all"));
// h3_histos.push_back(df_iso_raw_all.Histo3D({"h3_resp_raw_p_isHadE_all",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_raw_E_all", "hcal_E_all", "p_E_all"));
// h3_histos.push_back(df_iso_raw_all.Histo3D({"h3_resp_raw_p_isHadMIP_all",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, MIP);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_raw_MIP_all", "hcal_MIP_all", "p_MIP_all"));
// h3_histos.push_back(df_iso_raw_all.Histo3D({"h3_resp_raw_p_isHadEH_all",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL+HCAL);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.0, 2.5,  30, 0.0, 100.0}, "ep_raw_EH_all", "hcal_EH_all", "p_EH_all"));

// Per-tyyppiset 3D:t: Raw E/p + HCAL cut (pyysit nämä nimillä "h3_resp_raw_p_isHad*_cut")
// h3_histos.push_back(df_iso_raw_cut.Histo3D({"h3_resp_raw_p_isHadH_cut",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, HCAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_raw_H_cut", "hcal_H_cut", "p_H_cut"));
// h3_histos.push_back(df_iso_raw_cut.Histo3D({"h3_resp_raw_p_isHadE_cut",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, ECAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_raw_E_cut", "hcal_E_cut", "p_E_cut"));
// h3_histos.push_back(df_iso_raw_cut.Histo3D({"h3_resp_raw_p_isHadMIP_cut",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, MIP);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_raw_MIP_cut", "hcal_MIP_cut", "p_MIP_cut"));
// h3_histos.push_back(df_iso_raw_cut.Histo3D({"h3_resp_raw_p_isHadEH_cut",
//   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (HCAL cut, ECAL+HCAL);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//   50, 0.05, 2.5,  50, 0.9, 2.5,  30, 0.0, 100.0}, "ep_raw_EH_cut", "hcal_EH_cut", "p_EH_cut"));


// 2D & profiilit: S1–S4 (nimet yhtenevät myöhempien plottien kanssa)

// ------------------------------------------------------------
// S1: Raw/p, HCAL cut — inklusiivinen (kaikki hadronityypit)
// ------------------------------------------------------------
// h2_histos.push_back(df_iso.Histo2D(
//   {"h2_ep_vs_p_S1_raw_cut",
//    "S1: Raw E/p vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "pIso", "ep_raw_cut"));

// profiles.push_back(df_iso.Profile1D(
//   {"prof_ep_vs_p_S1_raw_cut",
//    "S1: <Raw E/p> vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "pIso", "ep_raw_cut"));

// h3_histos.push_back(df_iso.Histo3D(
//   {"h3_S1_raw_cut_ep_vs_hcal_p",
//    "S1: Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (cut);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//    50, 0.0, 2.5,
//    50, 0.0, 2.5,
//    150,0.0, 150.0},
//   "ep_raw_cut", "hcalCorrFactor", "pIso"));

// ------------------------------------------------------------
// S1: Raw/p, HCAL cut — per hadron type (H, E, MIP, EH)
// (raw E/p + HCAL cut: df_iso_raw_cut)
// ------------------------------------------------------------
// h2_histos.push_back(df_iso_raw_cut.Histo2D(
//   {"h2_ep_vs_p_S1_raw_cut_isHadH",
//    "S1 (HCAL cut, HCAL only): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_H_cut", "ep_raw_H_cut"));
// profiles.push_back(df_iso_raw_cut.Profile1D(
//   {"prof_ep_vs_p_S1_raw_cut_isHadH",
//    "S1 (HCAL cut, HCAL only): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_H_cut", "ep_raw_H_cut"));

// h2_histos.push_back(df_iso_raw_cut.Histo2D(
//   {"h2_ep_vs_p_S1_raw_cut_isHadE",
//    "S1 (HCAL cut, ECAL only): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_E_cut", "ep_raw_E_cut"));
// profiles.push_back(df_iso_raw_cut.Profile1D(
//   {"prof_ep_vs_p_S1_raw_cut_isHadE",
//    "S1 (HCAL cut, ECAL only): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_E_cut", "ep_raw_E_cut"));

// h2_histos.push_back(df_iso_raw_cut.Histo2D(
//   {"h2_ep_vs_p_S1_raw_cut_isHadMIP",
//    "S1 (HCAL cut, MIP): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_MIP_cut", "ep_raw_MIP_cut"));
// profiles.push_back(df_iso_raw_cut.Profile1D(
//   {"prof_ep_vs_p_S1_raw_cut_isHadMIP",
//    "S1 (HCAL cut, MIP): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_MIP_cut", "ep_raw_MIP_cut"));

// h2_histos.push_back(df_iso_raw_cut.Histo2D(
//   {"h2_ep_vs_p_S1_raw_cut_isHadEH",
//    "S1 (HCAL cut, ECAL+HCAL): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_EH_cut", "ep_raw_EH_cut"));
// profiles.push_back(df_iso_raw_cut.Profile1D(
//   {"prof_ep_vs_p_S1_raw_cut_isHadEH",
//    "S1 (HCAL cut, ECAL+HCAL): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_EH_cut", "ep_raw_EH_cut"));


// ------------------------------------------------------------
// S2: Default/p, HCAL cut — inklusiivinen
// ------------------------------------------------------------
// h2_histos.push_back(df_iso.Histo2D(
//   {"h2_ep_vs_p_S2_def_cut",
//    "S2: Default E/p vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "pIso", "ep_def_cut"));

// profiles.push_back(df_iso.Profile1D(
//   {"prof_ep_vs_p_S2_def_cut",
//    "S2: <Default E/p> vs p (HCAL E/E_{raw}>0.9);p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "pIso", "ep_def_cut"));

// h3_histos.push_back(df_iso.Histo3D(
//   {"h3_S2_def_cut_ep_vs_hcal_p",
//    "S2: Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (cut);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//    50, 0.0, 2.5,
//    50, 0., 2.5,
//    150,0.0, 150.0},
//   "ep_def_cut", "hcalCorrFactor", "pIso"));

// ------------------------------------------------------------
// S2: Default/p, HCAL cut — per hadron type (H, E, MIP, EH)
// (default E/p + HCAL cut: df_iso2)
// ------------------------------------------------------------
// h2_histos.push_back(df_iso2.Histo2D(
//   {"h2_ep_vs_p_S2_def_cut_isHadH",
//    "S2 (HCAL cut, HCAL only): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.05, 2.25},
//   "p_H_cut", "ep_def_H_cut"));
// profiles.push_back(df_iso2.Profile1D(
//   {"prof_ep_vs_p_S2_def_cut_isHadH",
//    "S2 (HCAL cut, HCAL only): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_H_cut", "ep_def_H_cut"));

// h2_histos.push_back(df_iso2.Histo2D(
//   {"h2_ep_vs_p_S2_def_cut_isHadE",
//    "S2 (HCAL cut, ECAL only): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_E_cut", "ep_def_E_cut"));
// profiles.push_back(df_iso2.Profile1D(
//   {"prof_ep_vs_p_S2_def_cut_isHadE",
//    "S2 (HCAL cut, ECAL only): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_E_cut", "ep_def_E_cut"));

// h2_histos.push_back(df_iso2.Histo2D(
//   {"h2_ep_vs_p_S2_def_cut_isHadMIP",
//    "S2 (HCAL cut, MIP): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_MIP_cut", "ep_def_MIP_cut"));
// profiles.push_back(df_iso2.Profile1D(
//   {"prof_ep_vs_p_S2_def_cut_isHadMIP",
//    "S2 (HCAL cut, MIP): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_MIP_cut", "ep_def_MIP_cut"));

// h2_histos.push_back(df_iso2.Histo2D(
//   {"h2_ep_vs_p_S2_def_cut_isHadEH",
//    "S2 (HCAL cut, ECAL+HCAL): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_EH_cut", "ep_def_EH_cut"));
// profiles.push_back(df_iso2.Profile1D(
//   {"prof_ep_vs_p_S2_def_cut_isHadEH",
//    "S2 (HCAL cut, ECAL+HCAL): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_EH_cut", "ep_def_EH_cut"));


// ------------------------------------------------------------
// S3: Raw/p, no cut — inklusiivinen
// ------------------------------------------------------------
// h2_histos.push_back(df_iso.Histo2D(
//   {"h2_ep_vs_p_S3_raw_all",
//    "S3: Raw E/p vs p (no cut);p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "pIso", "ep_raw_all"));

// profiles.push_back(df_iso.Profile1D(
//   {"prof_ep_vs_p_S3_raw_all",
//    "S3: <Raw E/p> vs p (no cut);p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "pIso", "ep_raw_all"));

// h3_histos.push_back(df_iso.Histo3D(
//   {"h3_S3_raw_ep_vs_hcal_p",
//    "S3: Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p;E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//    50, 0.0, 2.5,
//    50, 0.0, 2.5,
//    150,0.0, 150.0},
//   "ep_raw_all", "hcalCorrFactor", "pIso"));

// ------------------------------------------------------------
// S3: Raw/p, no cut — per hadron type (df_iso_raw_all)
// ------------------------------------------------------------
// h2_histos.push_back(df_iso_raw_all.Histo2D(
//   {"h2_ep_vs_p_S3_raw_all_isHadH",
//    "S3 (no cut, HCAL only): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_H_all", "ep_raw_H_all"));
// profiles.push_back(df_iso_raw_all.Profile1D(
//   {"prof_ep_vs_p_S3_raw_all_isHadH",
//    "S3 (no cut, HCAL only): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_H_all", "ep_raw_H_all"));

// h2_histos.push_back(df_iso_raw_all.Histo2D(
//   {"h2_ep_vs_p_S3_raw_all_isHadE",
//    "S3 (no cut, ECAL only): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_E_all", "ep_raw_E_all"));
// profiles.push_back(df_iso_raw_all.Profile1D(
//   {"prof_ep_vs_p_S3_raw_all_isHadE",
//    "S3 (no cut, ECAL only): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_E_all", "ep_raw_E_all"));

// h2_histos.push_back(df_iso_raw_all.Histo2D(
//   {"h2_ep_vs_p_S3_raw_all_isHadMIP",
//    "S3 (no cut, MIP): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_MIP_all", "ep_raw_MIP_all"));
// profiles.push_back(df_iso_raw_all.Profile1D(
//   {"prof_ep_vs_p_S3_raw_all_isHadMIP",
//    "S3 (no cut, MIP): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_MIP_all", "ep_raw_MIP_all"));

// h2_histos.push_back(df_iso_raw_all.Histo2D(
//   {"h2_ep_vs_p_S3_raw_all_isHadEH",
//    "S3 (no cut, ECAL+HCAL): Raw E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.50},
//   "p_EH_all", "ep_raw_EH_all"));
// profiles.push_back(df_iso_raw_all.Profile1D(
//   {"prof_ep_vs_p_S3_raw_all_isHadEH",
//    "S3 (no cut, ECAL+HCAL): <Raw E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_EH_all", "ep_raw_EH_all"));


// ------------------------------------------------------------
// S4: Default/p, no cut — inklusiivinen
// ------------------------------------------------------------
// h2_histos.push_back(df_iso.Histo2D(
//   {"h2_ep_vs_p_S4_def_all",
//    "S4: Default E/p vs p (no cut);p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.50},
//   "pIso", "ep_def_all"));

// profiles.push_back(df_iso.Profile1D(
//   {"prof_ep_vs_p_S4_def_all",
//    "S4: <Default E/p> vs p (no cut);p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "pIso", "ep_def_all"));

// h3_histos.push_back(df_iso.Histo3D(
//   {"h3_S4_def_ep_vs_hcal_p",
//    "S4: Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p;E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
//    50, 0.0, 2.5,
//    50, 0.0, 2.5,
//    150,0.0, 150.0},
//   "ep_def_all", "hcalCorrFactor", "pIso"));

// ------------------------------------------------------------
// S4: Default/p, no cut — per hadron type (df_iso_nocut)
// ------------------------------------------------------------
// h2_histos.push_back(df_iso_nocut.Histo2D(
//   {"h2_ep_vs_p_S4_def_all_isHadH",
//    "S4 (no cut, HCAL only): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_H_all", "ep_def_H_all"));
// profiles.push_back(df_iso_nocut.Profile1D(
//   {"prof_ep_vs_p_S4_def_all_isHadH",
//    "S4 (no cut, HCAL only): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_H_all", "ep_def_H_all"));

// h2_histos.push_back(df_iso_nocut.Histo2D(
//   {"h2_ep_vs_p_S4_def_all_isHadE",
//    "S4 (no cut, ECAL only): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_E_all", "ep_def_E_all"));
// profiles.push_back(df_iso_nocut.Profile1D(
//   {"prof_ep_vs_p_S4_def_all_isHadE",
//    "S4 (no cut, ECAL only): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_E_all", "ep_def_E_all"));

// h2_histos.push_back(df_iso_nocut.Histo2D(
//   {"h2_ep_vs_p_S4_def_all_isHadMIP",
//    "S4 (no cut, MIP): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_MIP_all", "ep_def_MIP_all"));
// profiles.push_back(df_iso_nocut.Profile1D(
//   {"prof_ep_vs_p_S4_def_all_isHadMIP",
//    "S4 (no cut, MIP): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_MIP_all", "ep_def_MIP_all"));

// h2_histos.push_back(df_iso_nocut.Histo2D(
//   {"h2_ep_vs_p_S4_def_all_isHadEH",
//    "S4 (no cut, ECAL+HCAL): Default E/p vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins,
//    40, 0.0, 2.5},
//   "p_EH_all", "ep_def_EH_all"));
// profiles.push_back(df_iso_nocut.Profile1D(
//   {"prof_ep_vs_p_S4_def_all_isHadEH",
//    "S4 (no cut, ECAL+HCAL): <Default E/p> vs p;p (GeV);E/p",
//    nTrkPBins, trkPBins},
//   "p_EH_all", "ep_def_EH_all"));


// ------------------------------------------------------------
// Per-tyyppiset 3D:t: Default E/p, no cut (Y 0.0–2.5)
// ------------------------------------------------------------
h3_histos.push_back(df_iso.Histo3D(
  {"h3_barrel_resp_def_p_isHadH_all",
   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, HCAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_def_H_all", "hcal_H_all", "p_H_all"));
h3_histos.push_back(df_iso.Histo3D(
  {"h3_barrel_resp_def_p_isHadE_all",
   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_def_E_all", "hcal_E_all", "p_E_all"));
h3_histos.push_back(df_iso.Histo3D(
  {"h3_barrel_resp_def_p_isHadMIP_all",
   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, MIP);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_def_MIP_all", "hcal_MIP_all", "p_MIP_all"));
h3_histos.push_back(df_iso.Histo3D(
  {"h3_barrel_resp_def_p_isHadEH_all",
   "Default E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL+HCAL);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_def_EH_all", "hcal_EH_all", "p_EH_all"));


// ------------------------------------------------------------
// Per-tyyppiset 3D:t: Raw E/p, no cut (Y 0.0–2.5)
// ------------------------------------------------------------
h3_histos.push_back(df_iso_raw_all.Histo3D(
  {"h3_resp_raw_p_isHadH_all",
   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, HCAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_raw_H_all", "hcal_H_all", "p_H_all"));
h3_histos.push_back(df_iso_raw_all.Histo3D(
  {"h3_resp_raw_p_isHadE_all",
   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL only);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_raw_E_all", "hcal_E_all", "p_E_all"));
h3_histos.push_back(df_iso_raw_all.Histo3D(
  {"h3_resp_raw_p_isHadMIP_all",
   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, MIP);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_raw_MIP_all", "hcal_MIP_all", "p_MIP_all"));
h3_histos.push_back(df_iso_raw_all.Histo3D(
  {"h3_resp_raw_p_isHadEH_all",
   "Raw E/p vs E_{HCAL}/E_{HCAL}^{raw} vs p (no cut, ECAL+HCAL);E/p;E_{HCAL}/E_{HCAL}^{raw};p (GeV)",
   50, 0.0, 2.5,
   50, 0.0, 2.5,
   150,0.0, 150.0},
  "ep_raw_EH_all", "hcal_EH_all", "p_EH_all"));

  auto cutReport = df_iso.Report();
  cutReport->Print();
//  KIRJOITA TIEDOSTOON (kirjoitus laukaisee myös evaluoinnin)
// root tiedoston avaus
std::string finalPath;
if (outpath && std::string(outpath).size()>0) {
  finalPath = outpath;
  // tee hakemisto, jos ei olemassa
  std::string dir = gSystem->DirName(finalPath.c_str());
  if (dir != "" && dir != "." && dir != "/") {
    gSystem->mkdir(dir.c_str(), true);
  }
} else {
  // std::string outdir = "histograms";
  // gSystem->mkdir(outdir.c_str(), true);
  // finalPath = outdir + "/" + std::string(tag) + ".root";
  // std::string outdir = "../histograms";
  // std::string outdir = "/eos/user/m/mmarjama/my_pion_analysis/histograms";
  std::string outdir = "./histograms";
  gSystem->mkdir(outdir.c_str(), true);
  finalPath = outdir + "/" + std::string(tag) + ".root";
}



TFile out(finalPath.c_str(), "RECREATE");
out.cd();
// kirjoitetaan root tiedostoon (kirjoitus aloittaa myös evaluoinnin) ---
for (auto& h : h1_histos) h->Write();
for (auto& h : h2_histos) h->Write();
for (auto& h : h3_histos) h->Write();
for (auto& p : profiles)  p->Write();

// for (auto& h : h1_histos) h->GetValue();
// for (auto& h : h2_histos) h->GetValue();
// for (auto& h : h3_histos) h->GetValue();
// for (auto& p : profiles)  p->GetValue();


out.Close();

// tulostetaan kokonaisajankesto
auto end = std::chrono::high_resolution_clock::now();
std::cout << "Wrote " << finalPath << " | Elapsed time: "
          << std::chrono::duration<double>(end - start).count() << " s\n";
}