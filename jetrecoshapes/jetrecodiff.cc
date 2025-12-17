// =============================================================================
//            Differential Jet Shapes Analysis for Photoproduction Events
// =============================================================================
// Description: Calculates differential jet shapes for r = 0.1 to 1.0
//              for different subprocess categories (QQ, GG, GQ) and 
//              photoproduction types (Combined, Direct, Resolved)
//
// Input: ROOT file with photoproduction events
// Output: Console and file printout of differential jet shape values
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>
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
// Input file
const string input_filename = "/Users/siddharthsingh/Analysis/ph-new/evt/allevents_pt5GeV/eic141_pt5/eic141_pt5.root";

// Analysis mode selection
const bool DIJET_ONLY = false;    // Set to true for dijet analysis, false for all jets

// Jet clustering parameters 
const double R = 1.0;             // Jet radius
const double jet_etMin = 10.0;    // Minimum ET for ALL jets (base cut)
const double etaMin = 1.5;       // Minimum eta 
const double etaMax = 2.0;        // Maximum eta

// Inclusive jet cut (only used when DIJET_ONLY = false)
const double inclusive_leading_jet_et_min = 0.0;  // Leading jet ET cut for inclusive
const bool flag_inclusiveleadjet = false;

// Kinematic cuts (following ZEUS analysis)
// const double y_min = 0.2;           // Minimum inelasticity
// const double y_max = 0.85;          // Maximum inelasticity
// const double Q2_max = 1.0;          // Maximum Q² (GeV²) for quasi-real photons

// Dijet-specific cuts (only used when DIJET_ONLY = true)
const double leading_jet_et_min = 17.0;      // Minimum ET for leading jet (GeV)
const double subleading_jet_et_min = 12.0;   // Minimum ET for subleading jet (GeV)

// Differential jet shape parameters
const double delta_r = 0.10;      // Annulus width (±0.05 around r) it will be divided by 2 in the code
const double r_center_min = 0.1;         // Minimum r value
const double r_center_max = 1.0;         // Maximum r value
const double r_step = 0.1;        // Step size for r

// Particle selection cuts (minimal cuts - include soft particles)
// const double particle_eta_max = 5.0;    // Maximum particle |eta|

// Progress reporting
const int report_every = 100000;   // Report progress every N events

// =============================================================================
// LOG FILENAME GENERATION
// =============================================================================
string generateLogFilename(const string& input_filename, double etaMin, double etaMax, 
                          double leading_et, double subleading_et, double etMin) {
    // Extract dataset name from input filename
    string baseFileName = input_filename.substr(input_filename.find_last_of("/") + 1);
    baseFileName = baseFileName.substr(0, baseFileName.find_last_of("."));  // Remove .root extension
    
    // File naming scheme, for 2 decimal digits after eta
    auto formatFloat = [](double value) -> string {
    stringstream ss;
    ss << fixed << setprecision(2) << value;
    return ss.str();
    };

    // Generate log filename with parameters
    string eta_str = "eta_" + formatFloat(etaMin) + "_" + formatFloat(etaMax);
    string j1_str = "j1_" + to_string((int)leading_et);
    string j2_str = "j2_" + to_string((int)subleading_et);
    string etMin_str = "etMin_" + to_string((int)etMin);
    if (DIJET_ONLY == true){
        return "diff_dijets" + baseFileName + "_" + eta_str + "_" + j1_str + "_" + j2_str + "_" + etMin_str + ".log";
    }
    if (DIJET_ONLY == false){
        return "diff_alljets_" + baseFileName + "_" + eta_str + "_" + etMin_str + ".log";
    }
}
    

const string log_filename = generateLogFilename(input_filename, etaMin, etaMax, 
                                               leading_jet_et_min, subleading_jet_et_min, jet_etMin);

// =============================================================================
// DIFFERENTIAL JET SHAPE CALCULATION (FOR DIJET OR ALL-JET EVENTS)
// =============================================================================

double calculateDifferentialJetShape(const vector<PseudoJet>& jets, double r) {
    
    // Check jet requirements based on analysis mode
    if (DIJET_ONLY && jets.size() != 2) return 0.0;
    if (!DIJET_ONLY && jets.empty()) return 0.0;
    
    double R1 = r - (delta_r / 2.0);  // Inner radius
    double R2 = r + (delta_r / 2.0);  // Outer radius
    double rho_total = 0.0;
    int n_jets = jets.size();
    
    // Loop over all jets (2 for dijet mode, all for inclusive mode)
    for (int i = 0; i < n_jets; i++) {
        
        vector<PseudoJet> constituents = jets[i].constituents();
        double select_const_Et = 0.0;
        
        // Loop over constituents and select ones in annulus [R1, R2]
        for (unsigned int j = 0; j < constituents.size(); j++) {
            
            double dphi = jets[i].phi() - constituents[j].phi();
            double drap = jets[i].rapidity() - constituents[j].rapidity();
            double const_r = sqrt(dphi*dphi + drap*drap);
            
            // Add ET of constituents in the annulus
            if (const_r >= R1 && const_r <= R2) {
                select_const_Et += constituents[j].Et();
            }
        }
        
        // Calculate contribution from this jet
        double jet_contribution = select_const_Et / jets[i].Et();
        rho_total += (jet_contribution / delta_r);
    }
    
    return rho_total;
}

// =============================================================================
// DIJET SELECTION FUNCTION
// =============================================================================
bool passesDijetCuts(const vector<PseudoJet>& jets) {
    // Need exactly 2 jets for dijet mode
    if (jets.size() != 2) return false;
    
    PseudoJet leading_jet = jets[0];
    PseudoJet subleading_jet = jets[1];
    
    // Leading jet cuts
    if (leading_jet.Et() < leading_jet_et_min) return false;
    if (leading_jet.eta() < etaMin || leading_jet.eta() > etaMax) return false;
    
    // Subleading jet cuts
    if (subleading_jet.Et() < subleading_jet_et_min) return false;
    if (subleading_jet.eta() < etaMin || subleading_jet.eta() > etaMax) return false;
    
    return true;
}

// =============================================================================
// JET RECONSTRUCTION AND SHAPE ANALYSIS FUNCTION
// =============================================================================

vector<double> reconstructJetsAndCalculateShapes(const vector<float>& px, const vector<float>& py, 
                                                 const vector<float>& pz, const vector<float>& energy,
                                                 const vector<float>& eta) {
    
    // Check input vector sizes
    if (px.empty() || px.size() != py.size() || px.size() != pz.size() || 
        px.size() != energy.size() || px.size() != eta.size()) {
        return vector<double>();
    }
    
    // Create input particles for FastJet
    vector<PseudoJet> input_particles;
    
    for (size_t i = 0; i < px.size(); ++i) {
        // Apply minimal particle-level cuts
        // if (abs(eta[i]) > particle_eta_max) continue;
        if (energy[i] <= 0) continue;
        
        // Create PseudoJet
        PseudoJet particle(px[i], py[i], pz[i], energy[i]);
        particle.set_user_index(i);
        input_particles.push_back(particle);
    }
    
    if (input_particles.size() < 2) {
        return vector<double>();
    }
    
    // Define anti-kT algorithm
    JetDefinition jet_def(antikt_algorithm, R);
    
    // Run clustering - ClusterSequence stays in scope for constituent access
    ClusterSequence cs(input_particles, jet_def);
    
    // Get all jets
    vector<PseudoJet> jets = cs.inclusive_jets();
    
    // Sort jets by ET (descending) to identify leading jet
    sort(jets.begin(), jets.end(), 
         [](const PseudoJet& a, const PseudoJet& b) { return a.Et() > b.Et(); });
    
    // Apply cuts based on analysis mode
    vector<PseudoJet> selected_jets;
    
    for (const auto& jet : jets) {
        double et = jet.Et();
        double jet_eta = jet.eta();
        
        // Apply basic cuts: ET > jet_etMin and eta in range
        if (et > jet_etMin && jet_eta > etaMin && jet_eta < etaMax) {
            selected_jets.push_back(jet);
        }
    }

    // Additional check: leading jet must pass higher threshold
    if (!selected_jets.empty() && flag_inclusiveleadjet) {
        if (selected_jets[0].Et() < inclusive_leading_jet_et_min) {
            return vector<double>();  // Reject event if leading jet too soft
        }
    }

    // Check event selection based on analysis mode
    if (DIJET_ONLY) {
        // For dijet mode: require exactly 2 jets passing dijet cuts
        if (!passesDijetCuts(selected_jets)) {
            return vector<double>();
        }
    } else {
        // For inclusive mode: require at least 1 jet
        if (selected_jets.empty()) {
            return vector<double>();
        }
    }
    
    // Calculate differential jet shapes for all r center values
    vector<double> shapes;
    int n_r_points = 10; //int((r_center_max - r_center_min) / r_step) + 1;

    for (int r_idx = 0; r_idx < n_r_points; r_idx++) {
        double r_center = r_center_min + r_idx * r_step;  // Center of annulus
        double diff_shape = calculateDifferentialJetShape(selected_jets, r_center);
        shapes.push_back(diff_shape);
    }
    
    return shapes;
}

// =============================================================================
// MAIN ANALYSIS FUNCTION
// =============================================================================

int main() {
    
    // Create dual output for console and file logging
    DualOutput dout(log_filename);
    
    dout << "=============================================================================\n";
    dout << "DIFFERENTIAL JET SHAPES ANALYSIS FOR PHOTOPRODUCTION EVENTS\n";
    dout << "=============================================================================\n";
    dout << "Analysis Parameters:\n";
    dout << "  Analysis Mode:      " << (DIJET_ONLY ? "DIJET ONLY (exactly 2 jets required)" : "INCLUSIVE (all jets)") << "\n";
    dout << "  Jet Algorithm:      Anti-kT\n";
    dout << "  Jet Radius:         " << R << "\n";
    dout << "  Base Jet ET:        > " << jet_etMin << " GeV (all jets)\n";
    // dout << "  Inelasticity (y):   [" << y_min << ", " << y_max << "]\n";
    
    if (DIJET_ONLY) {
        dout << "  Leading Jet ET:     > " << leading_jet_et_min << " GeV\n";
        dout << "  Subleading Jet ET:  > " << subleading_jet_et_min << " GeV\n";
    }
    else {
        dout << "  Leading Jet ET:     > " << inclusive_leading_jet_et_min << " GeV (inclusive)\n";
    }
    
    dout << "  Jet Eta Range:      [" << etaMin << ", " << etaMax << "]\n";
    // dout << "  Max Particle |η|:   " << particle_eta_max << "\n";
    dout << "  r Range:             [" << r_center_min << ", " << r_center_max << "] (annulus centers) with step " << r_step << "\n";
    dout << "  Annulus Width:      ±" << (delta_r/2.0) << " around r\n";
    dout << "=============================================================================\n\n";
    
    // Open input file
    TFile* input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_filename << endl;
        return 1;
    }
    
    dout << "Successfully opened input file: " << input_filename << "\n";
    dout << "Log file: " << log_filename << "\n\n";
    
    // Define event categories to process
    vector<string> categories = {"QQ_Events", "GG_Events", "GQ_Events", 
                                "Combined_Events", "Resolved_Events", "Direct_Events"};
    
    // Storage for differential jet shapes
    map<string, vector<vector<double>>> diff_shapes; // [category][r_index][event_values]
    
    // Initialize storage
    int n_r_points = 10; //int((r_center_max - r_center_min) / r_step) + 1;
    for (const string& category : categories) {
        diff_shapes[category].resize(n_r_points);
    }
    
    // Overall statistics
    int total_events_processed = 0;
    int total_dijets_selected = 0;
    
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
        // dout << "  Found " << n_entries << " events\n";
        
        if (n_entries == 0) {
            dout << "  No events found, skipping...\n";
            continue;
        }
        
        // Set up input branches
        vector<float> *px = nullptr, *py = nullptr, *pz = nullptr, *energy = nullptr, *eta = nullptr;
        // float inelasticity, scatteredElectronEnergy, photonEnergy; 
        
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
        // input_tree->SetBranchAddress("inelasticity", &inelasticity);
        // input_tree->SetBranchAddress("scatteredElectronEnergy", &scatteredElectronEnergy);
        // input_tree->SetBranchAddress("photonEnergy", &photonEnergy);    
        
        // Event counters
        int events_processed = 0;
        int events_with_jets = 0;
        int events_selected = 0;
        
        // Process events
        for (Long64_t i = 0; i < n_entries; ++i) {
            
            // if (i % report_every == 0) {
            //     dout << "    Processing event " << i << "/" << n_entries 
            //          << " (" << (100*i/n_entries) << "%)\r" << flush;
            // }
            
            input_tree->GetEntry(i);

            // Apply ZEUS kinematic cuts (before any particle-level processing)
            // if (inelasticity < y_min || inelasticity > y_max) continue;
            
            // Check if vectors are valid
            if (!px || !py || !pz || !energy || !eta) continue;
            if (px->size() != py->size() || px->size() != pz->size() || 
                px->size() != energy->size() || px->size() != eta->size()) continue;
            
            events_processed++;
            
            // Reconstruct jets and calculate all differential shapes in one go
            vector<double> shapes = reconstructJetsAndCalculateShapes(*px, *py, *pz, *energy, *eta);
            
            if (shapes.empty()) continue;
            
            events_with_jets++;
            events_selected++;  // Event passed all selection criteria
            
            // Store the shapes for each r value
            for (int r_idx = 0; r_idx < (int)shapes.size() && r_idx < n_r_points; r_idx++) {
                diff_shapes[category][r_idx].push_back(shapes[r_idx]);
            }
        }
        
        dout << "\n  Events processed: " << events_processed << "\n";
        dout << "  Events with jets: " << events_with_jets 
             << " (" << (100.0*events_with_jets/events_processed) << "%)\n";
        
        if (DIJET_ONLY) {
            dout << "  Events with exactly 2 jets passing dijet cuts: " << events_selected 
                 << " (" << (100.0*events_selected/events_processed) << "%)\n\n";
        } else {
            dout << "  Events selected (≥1 jet): " << events_selected 
                 << " (" << (100.0*events_selected/events_processed) << "%)\n\n";
        }
             
        // Update global counters
        total_events_processed += events_processed;
        total_dijets_selected += events_selected;
    }
    
    // Calculate and print results
    dout << "=============================================================================\n";
    dout << "DIFFERENTIAL JET SHAPE RESULTS\n";
    dout << "=============================================================================\n";
    
    if (DIJET_ONLY) {
        dout << "Requirement: Exactly 2 jets passing all dijet cuts\n";
        dout << "Leading Jet: ET > " << leading_jet_et_min << " GeV\n";
        dout << "Subleading Jet: ET > " << subleading_jet_et_min << " GeV\n";
        dout << "Both Jets: " << etaMin << " < η < " << etaMax << "\n";
    } else {
        dout << "Requirement: At least 1 jet passing cuts\n";
        dout << "All Jets: ET > " << jet_etMin << " GeV, " << etaMin << " < η < " << etaMax << "\n";
    }
    
    // dout << "Particles: |η| < " << particle_eta_max << "\n";
    dout << "=============================================================================\n\n";
    
    // Print header
    dout << setw(6) << "r";
    for (const string& category : categories) {
        dout << setw(15) << category;
    }
    dout << "\n";
    
    // Print separator
    dout << setw(6) << "-----";
    for (const string& category : categories) {
        dout << setw(15) << "-------------";
    }
    dout << "\n";
    
    // Print results for each r value
    for (int r_idx = 0; r_idx < n_r_points; r_idx++) {
        double r_center = r_center_min + r_idx * r_step;  // Center of annulus
        dout << setw(6) << setprecision(2) << fixed << r_center;  // Show 2 decimal places

        for (const string& category : categories) {
            if (!diff_shapes[category][r_idx].empty()) {
                // Calculate average differential jet shape for this r
                double sum = 0.0;
                for (double val : diff_shapes[category][r_idx]) {
                    sum += val;
                }
                double average = sum / diff_shapes[category][r_idx].size();
                dout << setw(15) << setprecision(4) << fixed << average;
            } else {
                dout << setw(15) << "N/A";
            }
        }
        dout << "\n";
    }
    
    dout << "\n=============================================================================\n";
    dout << "Number of events contributing to each category:\n";
    dout << "=============================================================================\n";
    
    for (const string& category : categories) {
        if (!diff_shapes[category][0].empty()) {
            dout << setw(20) << category << ": " << diff_shapes[category][0].size() << " events\n";
        }
    }
    
    dout << "\n=============================================================================\n";
    dout << "OVERALL STATISTICS\n";
    dout << "=============================================================================\n";
    dout << "Total events processed: " << total_events_processed << "\n";
    
    if (DIJET_ONLY) {
        dout << "Total dijet events selected: " << total_dijets_selected 
             << " (" << (100.0*total_dijets_selected/total_events_processed) << "%)\n";
    } else {
        dout << "Total events selected: " << total_dijets_selected 
             << " (" << (100.0*total_dijets_selected/total_events_processed) << "%)\n";
    }
    
    dout << "\n=============================================================================\n";
    dout << "ANALYSIS COMPLETED SUCCESSFULLY!\n";
    dout << "Log file: " << log_filename << "\n";
    dout << "=============================================================================\n";
    
    input_file->Close();
    return 0;
}