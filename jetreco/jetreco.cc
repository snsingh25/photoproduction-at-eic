// =============================================================================
// Jet Reconstruction Analysis for HERA Photoproduction Events
// =============================================================================
// Description: Reconstructs jets from different subprocess categories using 
//              anti-kT algorithm and saves basic jet properties with dijet cuts
//
// Input: ROOT file with photoproduction events from improved event generator
// Output: ROOT file with reconstructed dijet events for each process category
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TObjArray.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;
using namespace fastjet;

// =============================================================================
// DUAL OUTPUT CLASS FOR CONSOLE AND FILE LOGGING
// =============================================================================
class DualOutput {
private:
    std::ofstream logFile;
    bool fileOpen;
    
public:
    DualOutput(const std::string& filename) : fileOpen(false) {
        logFile.open(filename.c_str());
        fileOpen = logFile.is_open();
        if (!fileOpen) {
            std::cerr << "Warning: Could not open log file: " << filename << std::endl;
        }
    }
    
    ~DualOutput() {
        if (fileOpen) {
            logFile.close();
        }
    }
    
    template<typename T>
    DualOutput& operator<<(const T& data) {
        std::cout << data;  // Output to console
        if (fileOpen) {
            logFile << data;  // Output to file
            logFile.flush();  // Ensure immediate write
        }
        return *this;
    }
    
    // Handle special cases like std::endl, std::flush
    DualOutput& operator<<(std::ostream& (*manip)(std::ostream&)) {
        std::cout << manip;
        if (fileOpen) {
            logFile << manip;
            logFile.flush();
        }
        return *this;
    }
};

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================

// Analysis mode selection
const bool DIJET_ONLY = true;    // Set to true for dijet analysis, false for all jets

// Jet clustering parameters 
const double R = 1.0;             // Jet radius
const double etMin = 0.0;         // Minimum ET for any jet (HERA = 17, EIC = 10), ALL jets must pass this Base Cut
const double etaMin = -4.0;        // Minimum eta for any jet
const double etaMax = 4.0;         // Maximum eta for any jet

// Dijet-specific cuts (only used when DIJET_ONLY = true)
const double leading_jet_et_min = 10.0;      // Minimum ET for leading jet (GeV)
const double subleading_jet_et_min = 7.0;    // Minimum ET for subleading jet (GeV)
const double leading_jet_eta_max = 4.5;      // Maximum |eta| for leading jet
const double subleading_jet_eta_max = 4.5;   // Maximum |eta| for subleading jet

// Additional dijet cuts
// const double dijet_mass_min = 0.0;          // Minimum dijet invariant mass (GeV)
// const double delta_phi_min = 2.5;            // Minimum azimuthal separation (radians)
// const double delta_eta_max = 4.0;            // Maximum rapidity separation

// Particle selection cuts
const double particle_pt_min = 0.1;     // Minimum particle pT (GeV)
const double particle_eta_max = 5.0;    // Maximum particle |eta|

// Event progress reporting
const int report_every = 100000;    // Report progress every N events

// Input Events
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt5GeV/eic141_pt5GeV/eic141.root";

// =============================================================================
// OUTPUT FILENAME GENERATION
// =============================================================================
string generateOutputFilename(const string& input_filename, double R, double etMin) {
    string baseFileName = input_filename.substr(input_filename.find_last_of("/") + 1);
    baseFileName = baseFileName.substr(0, baseFileName.find_last_of("."));  // Remove .root extension
    string mode_str, leadjet_str, subleadjet_str;
    if (DIJET_ONLY) {
        mode_str = "dijets";
        leadjet_str = to_string((int)leading_jet_et_min);
        subleadjet_str = to_string((int)subleading_jet_et_min);
        return mode_str + "_" + baseFileName + "_R" + to_string((int)(R * 10)) + "_EtMin" + to_string((int)etMin) + "_" + leadjet_str + "_" + subleadjet_str + ".root";
    } else {
        mode_str = "alljets";
        return mode_str + "_" + baseFileName + "_R" + to_string((int)(R * 10)) + "_EtMin" + to_string((int)etMin) + ".root";
    }
}

string generateLogFilename(const string& root_filename) {
    string log_filename = root_filename;
    // Replace .root with .log
    size_t pos = log_filename.find_last_of(".");
    if (pos != string::npos) {
        log_filename = log_filename.substr(0, pos) + ".log";
    } else {
        log_filename += ".log";
    }
    return log_filename;
}

const string output_filename = generateOutputFilename(input_filename, R, etMin);
const string log_filename = generateLogFilename(output_filename);

// =============================================================================
// JET RECONSTRUCTION FUNCTION
// =============================================================================
vector<PseudoJet> reconstructJets(const vector<float>& px, const vector<float>& py, 
                                  const vector<float>& pz, const vector<float>& energy,
                                  const vector<float>& eta) {
    
    // Check input vector sizes
    if (px.empty() || px.size() != py.size() || px.size() != pz.size() || 
        px.size() != energy.size() || px.size() != eta.size()) {
        return vector<PseudoJet>();
    }
    
    // Create input particles for FastJet
    vector<PseudoJet> input_particles;
    
    for (size_t i = 0; i < px.size(); ++i) {
        // Apply particle-level cuts
        double pt = sqrt(px[i]*px[i] + py[i]*py[i]);
        if (pt < particle_pt_min) continue;
        if (abs(eta[i]) > particle_eta_max) continue;
        
        // Check for valid energy and momentum
        if (energy[i] <= 0) continue;
        
        // Create PseudoJet
        PseudoJet particle(px[i], py[i], pz[i], energy[i]);
        particle.set_user_index(i);  // Store original index
        input_particles.push_back(particle);
    }
    
    if (input_particles.size() < 2) {
        return vector<PseudoJet>();  
    }
    
    // Define anti-kT algorithm
    JetDefinition jet_def(antikt_algorithm, R);
    
    // Run clustering
    ClusterSequence cs(input_particles, jet_def);
    
    // Get jets and apply cuts
    vector<PseudoJet> jets = cs.inclusive_jets();
    vector<PseudoJet> selected_jets;
    
    for (const auto& jet : jets) {
        double et = jet.Et();
        double jet_eta = jet.eta();
        
        // Apply basic jet-level cuts
        if (et > etMin && jet_eta > etaMin && jet_eta < etaMax) {
            selected_jets.push_back(jet);
        }
    }
    
    // Sort jets by ET (descending)
    sort(selected_jets.begin(), selected_jets.end(), 
         [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
    
    return selected_jets;
}

// =============================================================================
// DIJET SELECTION FUNCTION
// =============================================================================
bool passesDijetCuts(const vector<PseudoJet>& jets) {
    // Need at least 2 jets
    if (jets.size() < 2) return false;
    
    PseudoJet leading_jet = jets[0];
    PseudoJet subleading_jet = jets[1];
    
    // Leading jet cuts
    if (leading_jet.Et() < leading_jet_et_min) return false;
    if (abs(leading_jet.eta()) > leading_jet_eta_max) return false;
    
    // Subleading jet cuts
    if (subleading_jet.Et() < subleading_jet_et_min) return false;
    if (abs(subleading_jet.eta()) > subleading_jet_eta_max) return false;
    
    // Calculate dijet properties
    TLorentzVector p4_leading, p4_subleading;
    p4_leading.SetPxPyPzE(leading_jet.px(), leading_jet.py(), leading_jet.pz(), leading_jet.E());
    p4_subleading.SetPxPyPzE(subleading_jet.px(), subleading_jet.py(), subleading_jet.pz(), subleading_jet.E());
    
    // // Dijet invariant mass cut
    // double dijet_mass = (p4_leading + p4_subleading).M();
    // if (dijet_mass < dijet_mass_min) return false;
    
    // // Azimuthal separation cut (back-to-back topology)
    // double delta_phi = abs(leading_jet.phi() - subleading_jet.phi());
    // if (delta_phi > M_PI) delta_phi = 2*M_PI - delta_phi;  // Wrap to [0, π]
    // if (delta_phi < delta_phi_min) return false;
    
    // // Rapidity separation cut
    // double delta_eta = abs(leading_jet.eta() - subleading_jet.eta());
    // if (delta_eta > delta_eta_max) return false;
    
    return true;
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================
int main() {

    // Create dual output for console and file logging
    DualOutput dout(log_filename);

    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    dout << "=============================================================================\n";
    dout << "JET RECONSTRUCTION ANALYSIS FOR HERA PHOTOPRODUCTION EVENTS\n";
    dout << "=============================================================================\n";
    dout << "Analysis Parameters:\n";
    dout << "  Analysis Mode:    " << (DIJET_ONLY ? "DIJET ONLY" : "ALL JETS") << "\n";
    dout << "  Jet Algorithm:    Anti-kT\n";
    dout << "  Jet Radius:       " << R << "\n";
    dout << "  Min Jet ET:       " << etMin << " GeV\n";
    dout << "  Jet Eta Range:    [" << etaMin << ", " << etaMax << "]\n";
    dout << "  Min Particle pT:  " << particle_pt_min << " GeV\n";
    
    if (DIJET_ONLY) {
        dout << "\nDijet Selection Cuts:\n";
        dout << "  Leading Jet ET:   > " << leading_jet_et_min << " GeV\n";
        dout << "  Subleading Jet ET:> " << subleading_jet_et_min << " GeV\n";
        dout << "  Leading Jet |η|:  < " << leading_jet_eta_max << "\n";
        dout << "  Subleading |η|:   < " << subleading_jet_eta_max << "\n";
        // dout << "  Dijet Mass:       > " << dijet_mass_min << " GeV\n";
        // dout << "  Δφ(j1,j2):        > " << delta_phi_min << " rad\n";
        // dout << "  Δη(j1,j2):        < " << delta_eta_max << "\n";
    }
    dout << "=============================================================================\n\n";
    
    dout << "Successfully opened input file: " << input_filename << "\n";
    dout << "Log file: " << log_filename << "\n\n";
    
    // Create output file
    TFile* output_file = new TFile(output_filename.c_str(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        cerr << "ERROR: Cannot create output file: " << output_filename << endl;
        return 1;
    }
    
    // Define event categories to process
    vector<string> categories = {"QQ_Events", "GG_Events", "GQ_Events", 
                                "Combined_Events", "Resolved_Events", "Direct_Events"};
    
    // Overall statistics
    int total_events_processed = 0;
    int total_events_selected = 0;
    
    // Process each category
    for (const string& category : categories) {
        
        dout << "Processing " << category << "...\n";
        
        // Get tree from input file
        TDirectory* input_dir = (TDirectory*)input_file->Get(category.c_str());
        if (!input_dir) {
            dout << "  Warning: Directory " << category << " not found, skipping...\n";
            continue;
        }
        
        TTree* input_tree = (TTree*)input_dir->Get(category.c_str());
        if (!input_tree) {
            dout << "  Warning: Tree " << category << " not found, skipping...\n";
            continue;
        }
        
        Long64_t n_entries = input_tree->GetEntries();
        dout << "  Found " << n_entries << " events\n";
        
        if (n_entries == 0) {
            dout << "  No events found, skipping...\n";
            continue;
        }
        
        // Set up input branches - check if branches exist first
        vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
        
        if (!input_tree->GetBranch("px") || !input_tree->GetBranch("py") || 
            !input_tree->GetBranch("pz") || !input_tree->GetBranch("energy") || 
            !input_tree->GetBranch("eta")) {
            dout << "  Error: Required branches not found in tree, skipping...\n";
            continue;
        }
        
        input_tree->SetBranchAddress("px", &px);
        input_tree->SetBranchAddress("py", &py);
        input_tree->SetBranchAddress("pz", &pz);
        input_tree->SetBranchAddress("energy", &energy);
        input_tree->SetBranchAddress("eta", &eta);
        
        // Create output directory and tree
        output_file->cd();
        TDirectory* output_dir = output_file->mkdir(category.c_str());
        output_dir->cd();
        
        TTree* output_tree = new TTree(("jets_" + category).c_str(), 
                                      ("Reconstructed Jets for " + category).c_str());
        
        // Output variables - basic jet info
        Int_t eventID, n_jets;
        vector<float> jet_et, jet_px, jet_py, jet_pz, jet_eta, jet_phi, jet_mass;
        
        // Dijet-specific variables
        Bool_t passes_dijet_cuts;
        Float_t leading_jet_et, subleading_jet_et;
        Float_t leading_jet_eta, subleading_jet_eta;
        Float_t leading_jet_phi, subleading_jet_phi;
        Float_t dijet_mass, dijet_delta_phi, dijet_delta_eta;
        Float_t dijet_pt, dijet_y;
        
        // Set up output branches - basic jet info
        output_tree->Branch("eventID", &eventID);
        output_tree->Branch("n_jets", &n_jets);
        output_tree->Branch("jet_et", &jet_et);
        output_tree->Branch("jet_px", &jet_px);
        output_tree->Branch("jet_py", &jet_py);
        output_tree->Branch("jet_pz", &jet_pz);
        output_tree->Branch("jet_eta", &jet_eta);
        output_tree->Branch("jet_phi", &jet_phi);
        output_tree->Branch("jet_mass", &jet_mass);
        
        // Dijet-specific branches
        output_tree->Branch("passes_dijet_cuts", &passes_dijet_cuts);
        output_tree->Branch("leading_jet_et", &leading_jet_et);
        output_tree->Branch("subleading_jet_et", &subleading_jet_et);
        output_tree->Branch("leading_jet_eta", &leading_jet_eta);
        output_tree->Branch("subleading_jet_eta", &subleading_jet_eta);
        output_tree->Branch("leading_jet_phi", &leading_jet_phi);
        output_tree->Branch("subleading_jet_phi", &subleading_jet_phi);
        output_tree->Branch("dijet_mass", &dijet_mass);
        output_tree->Branch("dijet_delta_phi", &dijet_delta_phi);
        output_tree->Branch("dijet_delta_eta", &dijet_delta_eta);
        output_tree->Branch("dijet_pt", &dijet_pt);
        output_tree->Branch("dijet_y", &dijet_y);
        
        // Event counters
        int events_processed = 0;
        int events_with_jets = 0;
        int events_selected = 0;
        int total_jets = 0;
        
        // Process events
        for (Long64_t i = 0; i < n_entries; ++i) {
            
            if (i % report_every == 0) {
                dout << "    Processing event " << i << "/" << n_entries 
                     << " (" << (100*i/n_entries) << "%)\r" << flush;
            }
            
            input_tree->GetEntry(i);
            
            // Check if vectors are valid and not empty
            if (!px || !py || !pz || !energy || !eta) {
                continue;
            }
            
            // Check if all vectors have the same size
            if (px->size() != py->size() || px->size() != pz->size() || 
                px->size() != energy->size() || px->size() != eta->size()) {
                continue;
            }
            
            events_processed++;
            
            // Clear output vectors
            jet_et.clear();
            jet_px.clear();
            jet_py.clear();
            jet_pz.clear();
            jet_eta.clear();
            jet_phi.clear();
            jet_mass.clear();
            
            // Reconstruct jets
            vector<PseudoJet> jets = reconstructJets(*px, *py, *pz, *energy, *eta);
            
            // Basic event info
            eventID = i;
            n_jets = jets.size();
            
            // Initialize dijet variables
            passes_dijet_cuts = false;
            leading_jet_et = subleading_jet_et = -999.0;
            leading_jet_eta = subleading_jet_eta = -999.0;
            leading_jet_phi = subleading_jet_phi = -999.0;
            dijet_mass = dijet_delta_phi = dijet_delta_eta = -999.0;
            dijet_pt = dijet_y = -999.0;
            
            // Fill jet info and check dijet criteria
            bool fill_event = false;
            
            if (n_jets >= 2) {
                events_with_jets++;
                
                // Check dijet cuts
                passes_dijet_cuts = passesDijetCuts(jets);
                
                // Fill leading and subleading jet info
                leading_jet_et = jets[0].Et();
                subleading_jet_et = jets[1].Et();
                leading_jet_eta = jets[0].eta();
                subleading_jet_eta = jets[1].eta();
                leading_jet_phi = jets[0].phi();
                subleading_jet_phi = jets[1].phi();
                
                // Calculate dijet system properties
                TLorentzVector p4_leading, p4_subleading;
                p4_leading.SetPxPyPzE(jets[0].px(), jets[0].py(), jets[0].pz(), jets[0].E());
                p4_subleading.SetPxPyPzE(jets[1].px(), jets[1].py(), jets[1].pz(), jets[1].E());
                TLorentzVector dijet_system = p4_leading + p4_subleading;
                
                dijet_mass = dijet_system.M();
                dijet_pt = dijet_system.Pt();
                dijet_y = dijet_system.Rapidity();
                
                dijet_delta_phi = abs(jets[0].phi() - jets[1].phi());
                if (dijet_delta_phi > M_PI) dijet_delta_phi = 2*M_PI - dijet_delta_phi;
                dijet_delta_eta = abs(jets[0].eta() - jets[1].eta());
                
                // Decide whether to fill this event
                if (DIJET_ONLY) {
                    fill_event = passes_dijet_cuts;
                } else {
                    fill_event = true;  // Fill all events with ≥2 jets for "all jets" mode
                }
            }
            
            // Fill jet vectors based on mode
            if (fill_event) {
                events_selected++;
                
                if (DIJET_ONLY) {
                    // Only save leading and subleading jets
                    total_jets += 2;
                    for (int j = 0; j < 2; ++j) {
                        jet_et.push_back(jets[j].Et());
                        jet_px.push_back(jets[j].px());
                        jet_py.push_back(jets[j].py());
                        jet_pz.push_back(jets[j].pz());
                        jet_eta.push_back(jets[j].eta());
                        jet_phi.push_back(jets[j].phi());
                        jet_mass.push_back(jets[j].m());
                    }
                } else {
                    // Save all jets
                    total_jets += n_jets;
                    for (const auto& jet : jets) {
                        jet_et.push_back(jet.Et());
                        jet_px.push_back(jet.px());
                        jet_py.push_back(jet.py());
                        jet_pz.push_back(jet.pz());
                        jet_eta.push_back(jet.eta());
                        jet_phi.push_back(jet.phi());
                        jet_mass.push_back(jet.m());
                    }
                }
                
                output_tree->Fill();
            }
        }
        
        dout << "\n  Events processed: " << events_processed << "\n";
        dout << "  Events with ≥2 jets: " << events_with_jets 
             << " (" << (100.0*events_with_jets/events_processed) << "%)\n";
        
        if (DIJET_ONLY) {
            dout << "  Events passing dijet cuts: " << events_selected 
                 << " (" << (100.0*events_selected/events_processed) << "%)\n";
        } else {
            dout << "  Events selected: " << events_selected 
                 << " (" << (100.0*events_selected/events_processed) << "%)\n";
        }
        
        dout << "  Total jets saved: " << total_jets << "\n";
        dout << "  Average jets/selected event: " 
             << (events_selected > 0 ? (double)total_jets/events_selected : 0) << "\n\n";
        
        // Update global counters
        total_events_processed += events_processed;
        total_events_selected += events_selected;
        
        // Write tree
        output_tree->Write();
        
        // Create basic histograms
        TH1F* h_n_jets = new TH1F(("h_n_jets_" + category).c_str(), 
                                 ("Number of Jets per Event - " + category + ";N_{jets};Events").c_str(), 
                                 10, 0, 10);
        TH1F* h_jet_et = new TH1F(("h_jet_et_" + category).c_str(), 
                                 ("Jet E_{T} Distribution - " + category + ";E_{T} [GeV];Jets").c_str(), 
                                 100, 0, 200);
        TH1F* h_jet_eta = new TH1F(("h_jet_eta_" + category).c_str(), 
                                  ("Jet #eta Distribution - " + category + ";#eta;Jets").c_str(), 
                                  80, -4, 4);
        
        // Dijet-specific histograms
        TH1F* h_dijet_mass = new TH1F(("h_dijet_mass_" + category).c_str(),
                                     ("Dijet Mass - " + category + ";M_{jj} [GeV];Events").c_str(),
                                     100, 0, 300);
        TH1F* h_dijet_delta_phi = new TH1F(("h_dijet_delta_phi_" + category).c_str(),
                                          ("Dijet #Delta#phi - " + category + ";#Delta#phi [rad];Events").c_str(),
                                          50, 0, M_PI);
        
        // Fill histograms by reading back the tree
        vector<float> *hist_jet_et = nullptr, *hist_jet_eta = nullptr;
        Int_t hist_n_jets;
        Bool_t hist_passes_cuts;
        Float_t hist_dijet_mass, hist_dijet_dphi;
        
        output_tree->SetBranchAddress("n_jets", &hist_n_jets);
        output_tree->SetBranchAddress("jet_et", &hist_jet_et);
        output_tree->SetBranchAddress("jet_eta", &hist_jet_eta);
        output_tree->SetBranchAddress("passes_dijet_cuts", &hist_passes_cuts);
        output_tree->SetBranchAddress("dijet_mass", &hist_dijet_mass);
        output_tree->SetBranchAddress("dijet_delta_phi", &hist_dijet_dphi);
        
        for (Long64_t i = 0; i < output_tree->GetEntries(); ++i) {
            output_tree->GetEntry(i);
            h_n_jets->Fill(hist_n_jets);
            
            if (hist_jet_et && hist_jet_eta) {
                for (size_t j = 0; j < hist_jet_et->size(); ++j) {
                    h_jet_et->Fill((*hist_jet_et)[j]);
                    h_jet_eta->Fill((*hist_jet_eta)[j]);
                }
            }
            
            if (hist_passes_cuts && hist_dijet_mass > 0) {
                h_dijet_mass->Fill(hist_dijet_mass);
                h_dijet_delta_phi->Fill(hist_dijet_dphi);
            }
        }
        
        h_n_jets->Write();
        h_jet_et->Write();
        h_jet_eta->Write();
        h_dijet_mass->Write();
        h_dijet_delta_phi->Write();
    }
    
    // Write and close files
    output_file->Write();
    output_file->Close();
    input_file->Close();
    
    dout << "=============================================================================\n";
    dout << "ANALYSIS COMPLETED SUCCESSFULLY!\n";
    dout << "=============================================================================\n";
    dout << "Output file: " << output_filename << "\n";
    dout << "Log file: " << log_filename << "\n";
    dout << "\nOverall Statistics:\n";
    dout << "  Total events processed: " << total_events_processed << "\n";
    dout << "  Total events selected: " << total_events_selected 
         << " (" << (100.0*total_events_selected/total_events_processed) << "%)\n";
    
    dout << "\nFor each category, the following are saved:\n";
    dout << "  - Jet trees with basic properties (ET, px, py, pz, eta, phi, mass)\n";
    dout << "  - Dijet-specific variables (masses, angles, cuts flags)\n";
    dout << "  - Basic distribution histograms\n";
    dout << "\nJets are reconstructed using anti-kT algorithm with:\n";
    dout << "  - R = " << R << "\n";
    dout << "  - ET > " << etMin << " GeV\n";
    dout << "  - " << etaMin << " < eta < " << etaMax << "\n";
    
    if (DIJET_ONLY) {
        dout << "\nDijet selection applied with:\n";
        dout << "  - Leading jet ET > " << leading_jet_et_min << " GeV\n";
        dout << "  - Subleading jet ET > " << subleading_jet_et_min << " GeV\n";
        // dout << "  - Dijet mass > " << dijet_mass_min << " GeV\n";
    }
    dout << "=============================================================================\n";
    
    return 0;
}