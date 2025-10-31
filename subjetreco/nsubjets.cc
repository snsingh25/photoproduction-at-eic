#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TDirectoryFile.h"
#include "TLine.h"
#include "TError.h"
#include "TBranch.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;

// Analysis parameters
const double R = 1.0;
const double ycut = 0.0005;
const double leading_jet_etMin = 10.0;
const double subleading_jet_etMin = 7.0;
// const double jet_etMin = 10.0;
const double etaMin = -1.0;
const double etaMax = 3.0;

int main() {
    gROOT->SetBatch(kTRUE);
    
    cout << "\n========================================" << endl;
    cout << "  Subjet Multiplicity Analysis" << endl;
    cout << "========================================" << endl;
    cout << "Jet Algorithm: anti-kT, R = " << R << endl;
    cout << "Subjet Algorithm: kT, ycut = " << ycut << endl;
    cout << "Leading Jet ET > " << leading_jet_etMin << " GeV" << endl;
    cout << "Sub-Leading Jet ET > " << subleading_jet_etMin << " GeV" << endl;
    // cout << "All Jets ET > " << jet_etMin << " GeV" << endl;
    cout << "Jet Eta Range: " << etaMin << " < eta < " << etaMax << endl;
    cout << "========================================\n" << endl;

    // Open ROOT file
    const string inputFile = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt7GeV/hera300_pT7/hera300_pT7.root";
    TFile* file = TFile::Open(inputFile.c_str(), "READ");
    
    if (!file || file->IsZombie()) {
        cerr << "ERROR: Cannot open file: " << inputFile << endl;
        return 1;
    }
    
    cout << "Successfully opened: " << inputFile << "\n" << endl;
    
    // Subprocess names and colors
    vector<string> subprocess_names = {"QQ_Events", "GG_Events", "GQ_Events"};
    vector<string> subprocess_labels = {"Quarks", "Gluons", "Quark-Gluon"};
    vector<int> colors = {kRed, kBlue, kGreen+2};
    vector<int> fill_styles = {3345, 3354, 3395}; // Hatched fill patterns
    
    // Create histograms
    TH1D* h_qq = new TH1D("h_qq", "Subjet Multiplicity", 10, 0, 10);
    TH1D* h_gg = new TH1D("h_gg", "Subjet Multiplicity", 10, 0, 10);
    TH1D* h_gq = new TH1D("h_gq", "Subjet Multiplicity", 10, 0, 10);
    
    vector<TH1D*> histograms = {h_qq, h_gg, h_gq};
    
    // Set histogram styles
    for (int i = 0; i < 3; i++) {
        histograms[i]->SetLineColor(colors[i]);
        histograms[i]->SetFillColor(colors[i]);
        histograms[i]->SetFillStyle(fill_styles[i]);
        histograms[i]->SetLineWidth(2);
        histograms[i]->GetXaxis()->SetTitle("n_{subjets} (y_{cut} = 0.0005)");
        histograms[i]->GetYaxis()->SetTitle("Normalised Number of Jets");
        histograms[i]->GetXaxis()->SetTitleSize(0.045);
        histograms[i]->GetYaxis()->SetTitleSize(0.045);
        histograms[i]->GetXaxis()->SetLabelSize(0.04);
        histograms[i]->GetYaxis()->SetLabelSize(0.04);
    }
    
    // Loop over three subprocesses
    for (int iProc = 0; iProc < 3; iProc++) {
        
        cout << "\n----------------------------------------" << endl;
        cout << "Processing: " << subprocess_names[iProc] << endl;
        cout << "----------------------------------------" << endl;
        
        // Access directory and tree
        TDirectoryFile* dir = (TDirectoryFile*)file->Get(subprocess_names[iProc].c_str());
        if (!dir) {
            cerr << "ERROR: Cannot find directory: " << subprocess_names[iProc] << endl;
            continue;
        }
        
        TTree* tree = (TTree*)dir->Get(subprocess_names[iProc].c_str());
        if (!tree) {
            cerr << "ERROR: Cannot find tree in directory: " << subprocess_names[iProc] << endl;
            continue;
        }
        
        // Set up branch addresses
        vector<float>* px = nullptr;
        vector<float>* py = nullptr;
        vector<float>* pz = nullptr;
        vector<float>* energy = nullptr;
        vector<float>* eta = nullptr;
        int processType = 0;
        
        tree->SetBranchAddress("px", &px);
        tree->SetBranchAddress("py", &py);
        tree->SetBranchAddress("pz", &pz);
        tree->SetBranchAddress("energy", &energy);
        tree->SetBranchAddress("eta", &eta);
        tree->SetBranchAddress("processType", &processType);
        
        int total_events = 0;
        int selected_events = 0;
        int total_jets = 0;
        
        Long64_t nEntries = tree->GetEntries();
        cout << "Total events in tree: " << nEntries << endl;
        
        // Event loop
        for (Long64_t i = 0; i < nEntries; ++i) {

            // if (i >= 10000) break;
            
            if (i % 50000 == 0 && i > 0) {
                cout << "  Processed " << i << " events..." << endl;
            }
            
            tree->GetEntry(i);
            total_events++;
            
            // Check if particles exist
            if (!px || px->empty()) continue;
            
            // Create particles for FastJet
            vector<PseudoJet> particles;
            for (size_t j = 0; j < px->size(); ++j) {
                if ((*energy)[j] <= 0) continue;
                particles.push_back(PseudoJet((*px)[j], (*py)[j], (*pz)[j], (*energy)[j]));
            }
            
            if (particles.empty()) continue;
            
            // Cluster jets with anti-kT algorithm
            JetDefinition jet_def(antikt_algorithm, R);
            ClusterSequence cs(particles, jet_def);
            vector<PseudoJet> jets = cs.inclusive_jets();
            
            // Select Jets on eta cuts
            vector<PseudoJet> selected_jets;
            for (const auto& jet : jets) {
                if (jet.eta() > etaMin && jet.eta() < etaMax) {
                    selected_jets.push_back(jet);
                }
            }
            
            // Sort jets by ET (descending)
            sort(selected_jets.begin(), selected_jets.end(), 
                 [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
            
            // Require at least one jet
            // if (selected_jets.empty()) continue;
            // Require Dijets
            if (selected_jets.size() != 2) continue;
            
            // Apply leading and subleading ET cuts
            if (selected_jets[0].Et() < leading_jet_etMin) continue;  // Leading > 10 GeV
            if (selected_jets[1].Et() < subleading_jet_etMin) continue;  // Subleading > 7 GeV
            
            selected_events++;
            
            // Loop over ALL selected jets in this event
            for (const auto& jet : selected_jets) {
                
                // Get jet constituents
                vector<PseudoJet> constituents = jet.constituents();
                
                if (constituents.size() < 2) continue; // Need at least 2 constituents
                
                // Recluster constituents with kT algorithm
                JetDefinition kt_def(kt_algorithm, R);
                ClusterSequence cs_kt(constituents, kt_def);
                
                // Get exclusive subjets with ycut
                vector<PseudoJet> subjets = cs_kt.exclusive_jets_ycut(ycut);
                int n_subjets = subjets.size();
                
                // Fill histogram
                histograms[iProc]->Fill(n_subjets);
                total_jets++;
            }
            
        } // End of event loop
        
        cout << "----------------------------------------" << endl;
        cout << "Summary for " << subprocess_labels[iProc] << ":" << endl;
        cout << "  Total events processed: " << total_events << endl;
        cout << "  Events passing cuts: " << selected_events << endl;
        cout << "  Total jets analyzed: " << total_jets << endl;
        cout << "  Jets in histogram: " << histograms[iProc]->GetEntries() << endl;
        cout << "----------------------------------------" << endl;
        
    } // End of subprocess loop
    
    // Normalize histograms
    cout << "\nNormalizing histograms..." << endl;
    for (int i = 0; i < 3; i++) {
        double integral = histograms[i]->Integral();
        if (integral > 0) {
            histograms[i]->Scale(1.0 / integral);
            cout << subprocess_labels[i] << " normalized (integral was " << integral << ")" << endl;
        }
    }

    cout << "\n========================================" << endl;
    cout << "  Creating Plot" << endl;
    cout << "========================================\n" << endl;

    // SAFETY CHECKS - Add these
    cout << "Checking histograms..." << endl;
    if (!h_qq) { cout << "ERROR: h_qq is null!" << endl; return 1; }
    if (!h_gg) { cout << "ERROR: h_gg is null!" << endl; return 1; }
    if (!h_gq) { cout << "ERROR: h_gq is null!" << endl; return 1; }

    cout << "QQ entries: " << h_qq->GetEntries() << endl;
    cout << "GG entries: " << h_gg->GetEntries() << endl;
    cout << "GQ entries: " << h_gq->GetEntries() << endl;

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kError;
    
    // Create canvas
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "Subjet Multiplicity", 900, 700);
    c1->cd(); 
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.08);
    c1->SetBottomMargin(0.12);
    
    // Find maximum for y-axis
    // double max_val = 0;
    // for (int i = 0; i < 3; i++) {
    //     double this_max = histograms[i]->GetMaximum();
    //     if (this_max > max_val) max_val = this_max;
    // }
    double max_val = TMath::Max(h_qq->GetMaximum(), TMath::Max(h_gg->GetMaximum(), h_gq->GetMaximum()));
    
    // Draw histograms
    h_qq->SetMaximum(max_val * 1.3);
    h_qq->SetMinimum(0);
    h_qq->Draw("HIST");
    h_gg->Draw("HIST SAME");
    h_gq->Draw("HIST SAME");
    
    // Create legend
    TLegend* leg = new TLegend(0.65, 0.65, 0.90, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.04);
    leg->AddEntry(h_qq, "Quarks", "f");
    leg->AddEntry(h_gg, "Gluons", "f");
    leg->AddEntry(h_gq, "Quark-Gluon", "f");
    leg->Draw();
    
    // Add text with cuts
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    // latex.DrawLatex(0.15, 0.85, "E_{T}^{jet} > 10 GeV");
    // latex.DrawLatex(0.15, 0.80, "-1 < #eta_{jet} < 3");
    latex.DrawLatex(0.15, 0.85, "Dijet Events");
    latex.DrawLatex(0.15, 0.80, Form("E_{T}^{leading} > %.0f GeV", leading_jet_etMin));
    latex.DrawLatex(0.15, 0.75, Form("E_{T}^{subleading} > %.0f GeV", subleading_jet_etMin));
    latex.DrawLatex(0.15, 0.70, Form("%.0f < #eta_{jet} < %.0f", etaMin, etaMax));
    
    // Save canvas
    // c1->SaveAs("subjet_multiplicity.pdf");
    // TString filename = Form("sJM_Et%.0fGeV_eta%.0fto%.0f.pdf", jet_etMin, etaMin, etaMax);
    TString filename = Form("sJM_Dijet_Et_L%.0f_SL_%.0fGeV_eta%.0fto%.0f.pdf", leading_jet_etMin, subleading_jet_etMin, etaMin, etaMax);
    c1->SaveAs(filename);
    cout << "Plot saved as: " << filename << endl;
    
    // Print histogram statistics
    cout << "\n========================================" << endl;
    cout << "  Histogram Statistics" << endl;
    cout << "========================================" << endl;
    for (int i = 0; i < 3; i++) {
        cout << subprocess_labels[i] << ":" << endl;
        cout << "  Mean n_subjets: " << histograms[i]->GetMean() << endl;
        cout << "  RMS: " << histograms[i]->GetRMS() << endl;
        cout << "  Total entries: " << histograms[i]->GetEntries() << endl;
        cout << endl;
    }
    
    // Clean up
    delete c1;
    for (auto h : histograms) delete h;
    
    cout << "========================================" << endl;
    cout << "  Analysis Complete!" << endl;
    cout << "========================================\n" << endl;
    
    return 0;
    file->Close();
}