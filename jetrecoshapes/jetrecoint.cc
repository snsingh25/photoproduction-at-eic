// =============================================================================
//            Integrated Jet Shape Analysis (psi at fixed r)
// =============================================================================
// Description: Calculates the integrated jet shape psi(r=0.3) for dijet events.
//              Classifies jets as Thin (Quark-like) or Thick (Gluon-like).
//              Includes strict event validation and formatted output.
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <map>

using namespace std;
using namespace fastjet;

// =============================================================================
// CONFIGURATION
// =============================================================================
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera300_pT7/hera300_pT7.root";
const string output_pdf = "integrated_jet_shapes_clean.pdf";

// Jet Clustering
const double R = 1.0;
const double etaMin = -1.0;
const double etaMax = 4.0;

// Dijet Selection
const double leading_jet_etMin = 10.0;
const double subleading_jet_etMin = 7.0;

// Shape Analysis
const double JET_SHAPE_RADIUS = 0.3; // The fixed r value for psi(r)

// Physics Classification Cuts
const double THIN_JET_CUT = 0.8;   // psi > 0.8  -> Thin (Quark)
const double THICK_JET_CUT = 0.6;  // psi < 0.6  -> Thick (Gluon)
const double SOFT_JET_CUT = 0.01;  // psi < 0.01 -> Artifact/Reject

// Process Type Enum (based on your file structure)
enum ProcessType {
    UNKNOWN = 0,
    QQ = 1,
    GG = 2,
    GQ = 3
};

// =============================================================================
// HELPER: CALCULATE PSI(r)
// =============================================================================
double calculateIntegratedJetShape(const PseudoJet& jet, double r_cone) {
    if (jet.Et() <= 0.0) return 0.0;

    double energy_within_r = 0.0;
    vector<PseudoJet> constituents = jet.constituents();
    
    for (const auto& constituent : constituents) {
        // Use FastJet's built-in delta_R to handle phi wrapping automatically and correctly
        if (constituent.delta_R(jet) <= r_cone) {
            energy_within_r += constituent.Et();
        }
    }
    
    return energy_within_r / jet.Et();
}

// =============================================================================
// MAIN
// =============================================================================
int main() {
    
    // 1. Setup Style
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    cout << "=============================================================================\n";
    cout << "STARTING INTEGRATED JET SHAPE ANALYSIS (psi at r=" << JET_SHAPE_RADIUS << ")\n";
    cout << "=============================================================================\n";

    // 2. Initialize Histograms
    map<int, TH1F*> histos;
    histos[QQ] = new TH1F("h_QQ", "QQ Jets", 50, 0, 1);
    histos[GG] = new TH1F("h_GG", "GG Jets", 50, 0, 1);
    histos[GQ] = new TH1F("h_GQ", "GQ Jets", 50, 0, 1);

    // Style: QQ (Red, Hatch), GG (Blue, Hatch), GQ (Green, Hatch)
    histos[QQ]->SetLineColor(kRed);     histos[QQ]->SetFillColor(kRed);     histos[QQ]->SetFillStyle(3004); histos[QQ]->SetLineWidth(2);
    histos[GG]->SetLineColor(kBlue);    histos[GG]->SetFillColor(kBlue);    histos[GG]->SetFillStyle(3005); histos[GG]->SetLineWidth(2);
    histos[GQ]->SetLineColor(kGreen+2); histos[GQ]->SetFillColor(kGreen+2); histos[GQ]->SetFillStyle(3006); histos[GQ]->SetLineWidth(2);

    // 3. Initialize Statistics Counters
    struct Stats {
        int total_clean = 0;
        int thin = 0;
        int thick = 0;
    };
    map<int, Stats> jet_stats;
    
    struct EventStats {
        int total = 0;
        int tt_events = 0; // Both jets are thin
    };
    map<int, EventStats> evt_stats;

    int total_events_read = 0;
    int dijet_events_selected = 0;
    int rejected_soft_events = 0;

    // 4. Open File
    TFile* file = TFile::Open(input_filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "ERROR: Cannot open file: " << input_filename << endl;
        return 1;
    }

    // Access specific Tree path from your code
    TDirectoryFile* dir = (TDirectoryFile*)file->Get("Combined_Events");
    if (!dir) { cerr << "ERROR: Directory Combined_Events not found." << endl; return 1; }
    
    TTree* tree = (TTree*)dir->Get("Combined_Events");
    if (!tree) { cerr << "ERROR: Tree Combined_Events not found." << endl; return 1; }

    // 5. Branch Setup
    vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
    int processType = 0;

    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("eta", &eta); // Used for pre-check if needed, though we recalc eta from PxPyPz
    tree->SetBranchAddress("processType", &processType);

    Long64_t n_entries = tree->GetEntries();
    // cout << "Found " << n_entries << " entries." << endl;

    // 6. Event Loop
    for (Long64_t i = 0; i < n_entries; ++i) {
        if (i % 100000 == 0 && i > 0) cout << "Processing event " << i << "..." << endl;
        
        tree->GetEntry(i);
        total_events_read++;

        if (!px || px->empty()) continue;

        // A. Reconstruct Particles
        vector<PseudoJet> particles;
        for (size_t j = 0; j < px->size(); ++j) {
            if ((*energy)[j] > 0) {
                particles.push_back(PseudoJet((*px)[j], (*py)[j], (*pz)[j], (*energy)[j]));
            }
        }
        if (particles.size() < 2) continue;

        // B. Cluster Jets
        ClusterSequence cs(particles, JetDefinition(antikt_algorithm, R));
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

        // C. Select Jets (Eta & Subleading Et cut)
        vector<PseudoJet> selected_jets;
        for (const auto& jet : jets) {
            if (jet.Et() > subleading_jet_etMin && jet.eta() > etaMin && jet.eta() < etaMax) {
                selected_jets.push_back(jet);
            }
        }

        // D. Dijet Event Selection
        // Strict Requirement: Exactly 2 jets, leading jet must pass higher threshold
        if (selected_jets.size() != 2) continue;
        if (selected_jets[0].Et() < leading_jet_etMin) continue;

        // E. Calculate Shapes & Validate Event
        // CRITICAL FIX: Calculate both psis first. If ANY jet is soft, reject the WHOLE event 
        // to maintain consistency between event counters and histograms.
        vector<double> psis;
        bool clean_event = true;

        for (const auto& jet : selected_jets) {
            double psi = calculateIntegratedJetShape(jet, JET_SHAPE_RADIUS);
            
            // Validity Checks
            if (psi < 0.0 || psi > 1.0) { clean_event = false; break; }
            if (psi < SOFT_JET_CUT)     { clean_event = false; break; } // Reject soft artifacts

            psis.push_back(psi);
        }

        if (!clean_event) {
            rejected_soft_events++;
            continue;
        }

        // F. Fill Histograms & Statistics
        // If we reached here, the event is valid, clean, and dijet
        dijet_events_selected++;
        
        // Event Level Classification
        bool both_thin = (psis[0] > THIN_JET_CUT && psis[1] > THIN_JET_CUT);
        evt_stats[processType].total++;
        if (both_thin) evt_stats[processType].tt_events++;

        // Jet Level Filling
        for (double psi : psis) {
            if (histos.find(processType) != histos.end()) {
                histos[processType]->Fill(psi);
                
                jet_stats[processType].total_clean++;
                if (psi > THIN_JET_CUT) jet_stats[processType].thin++;
                if (psi < THICK_JET_CUT) jet_stats[processType].thick++;
            }
        }
    }

    file->Close();

    // 7. Visualization
    TCanvas* c1 = new TCanvas("c1", "Integrated Jet Shapes", 800, 800);
    
    // Auto-scale Y axis
    double max_val = max({histos[QQ]->GetMaximum(), histos[GG]->GetMaximum(), histos[GQ]->GetMaximum()});
    histos[QQ]->SetMaximum(max_val * 1.2);
    histos[QQ]->SetTitle(Form("Integrated Jet Shape #psi(r=%.1f); #psi(r); Entries", JET_SHAPE_RADIUS));
    histos[QQ]->GetYaxis()->SetTitleOffset(1.3);

    histos[QQ]->Draw("HIST");
    histos[GG]->Draw("HIST SAME");
    histos[GQ]->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.15, 0.70, 0.45, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(histos[QQ], Form("QQ (%d jets)", jet_stats[QQ].total_clean), "f");
    legend->AddEntry(histos[GG], Form("GG (%d jets)", jet_stats[GG].total_clean), "f");
    legend->AddEntry(histos[GQ], Form("GQ (%d jets)", jet_stats[GQ].total_clean), "f");
    legend->Draw();

    c1->Print(output_pdf.c_str());

    // 8. Print Results
    auto print_separator = []() { cout << "-----------------------------------------------------------------------------" << endl; };
    
    cout << "\n";
    cout << "=============================================================================\n";
    cout << "ANALYSIS STATISTICS\n";
    cout << "=============================================================================\n";
    cout << "Total events read:       " << total_events_read << endl;
    cout << "Rejected (Soft Jets):    " << rejected_soft_events << " events (containing artifact jets)" << endl;
    cout << "Final Selected Dijets:   " << dijet_events_selected << endl;
    
    cout << "\n>>> EVENT CLASSIFICATION (TT = Both jets Thin)" << endl;
    print_separator();
    cout << left << setw(10) << "Type" << setw(15) << "Total Events" << setw(15) << "TT Events" << setw(15) << "Fraction (%)" << endl;
    print_separator();
    
    int types[] = {QQ, GG, GQ};
    string names[] = {"QQ", "GG", "GQ"};
    
    for(int i=0; i<3; ++i) {
        int t = types[i];
        double frac = (evt_stats[t].total > 0) ? 100.0 * evt_stats[t].tt_events / evt_stats[t].total : 0.0;
        cout << left << setw(10) << names[i] 
             << setw(15) << evt_stats[t].total 
             << setw(15) << evt_stats[t].tt_events 
             << setw(15) << fixed << setprecision(1) << frac << endl;
    }
    
    cout << "\n>>> JET CLASSIFICATION (Thin: psi > " << THIN_JET_CUT << ", Thick: psi < " << THICK_JET_CUT << ")" << endl;
    print_separator();
    cout << left << setw(10) << "Type" << setw(15) << "Total Jets" << setw(15) << "Thin Jets" << setw(15) << "Thick Jets" << endl;
    print_separator();

    int tot_thin = 0, tot_thick = 0;
    for(int i=0; i<3; ++i) {
        int t = types[i];
        cout << left << setw(10) << names[i] 
             << setw(15) << jet_stats[t].total_clean
             << setw(15) << jet_stats[t].thin 
             << setw(15) << jet_stats[t].thick << endl;
        
        tot_thin += jet_stats[t].thin;
        tot_thick += jet_stats[t].thick;
    }
    print_separator();
    
    // Efficiency/Purity
    if (tot_thin > 0 && tot_thick > 0) {
        cout << "\n>>> PERFORMANCE METRICS" << endl;
        cout << "Thin Jet Purity (is QQ?):  " << fixed << setprecision(2) << (100.0 * jet_stats[QQ].thin / tot_thin) << " %" << endl;
        cout << "Thick Jet Purity (is GG?): " << fixed << setprecision(2) << (100.0 * jet_stats[GG].thick / tot_thick) << " %" << endl;
    }

    cout << "\nPlot saved to: " << output_pdf << endl;
    cout << "=============================================================================\n";

    // Cleanup
    delete c1;
    delete legend;
    for(auto& h : histos) delete h.second;

    return 0;
}