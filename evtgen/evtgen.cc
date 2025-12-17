// =============================================================================
// Pure Photoproduction Event Generation Script
// =============================================================================
// Description: Generates photoproduction events using PYTHIA8 and stores
//              final state particles categorized by process type
//
// Output: ROOT file with particle data for different process types
// Author: Siddharth Singh
// =============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "Pythia8/Pythia.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <iomanip>

using namespace Pythia8;
using namespace std;

// =============================================================================
// ANALYSIS CONFIGURATION
// =============================================================================

struct EventConfig {
    // Event generation
    int nEvents = 1000000;           // Number of events to generate
    
    // Beam configuration
    double protonEnergy = 100;     // GeV (100, 275, 275, 820)
    double electronEnergy = 10;    // GeV (10, 10, 18, 27.5)
    
    // Physics cuts
    double pTHatMin = 5.0;          // Minimum partonic pT
    double Q2max = 1.0;             // Maximum Q2 for photoproduction
    double pT0Ref = 3.0;            // Tuned for photon-hadron interactions
    
    void print() const {
        cout << "Event Generation Configuration:\n";
        cout << "  Number of Events: " << nEvents << "\n";
        cout << "  Proton Energy:    " << protonEnergy << " GeV\n";
        cout << "  Electron Energy:  " << electronEnergy << " GeV\n";
        cout << "  pT Hat Min:       " << pTHatMin << " GeV\n";
        cout << "  Q2 Max:           " << Q2max << " GeV^2\n";
    }
};

// =============================================================================
// SUBPROCESS CLASSIFICATION
// =============================================================================

int classifySubprocess(int processCode) {
    // QQ: Quark-initiated processes
    if (processCode == 112 || processCode == 114 || processCode == 116 || 
        processCode == 121 || processCode == 122 || processCode == 123 || 
        processCode == 124 || processCode == 271 || processCode == 272 || 
        processCode == 273 || processCode == 281 || processCode == 282 || 
        processCode == 283 || processCode == 284 || processCode == 285) {
        return 1; // QQ
    }
    // GG: Gluon-gluon processes
    else if (processCode == 111 || processCode == 115) {
        return 2; // GG
    }
    // GQ: Gluon-quark processes  
    else if (processCode == 113 || processCode == 274) {
        return 3; // GQ
    }
    else {
        return 0; // Unclassified
    }
}

// =============================================================================
// MAIN FUNCTION
// =============================================================================

int main() {
    
    try {
        
        // =====================================================================
        // CONFIGURATION
        // =====================================================================
        
        EventConfig config;
        
        const string outputFileName = "eic64.root";
        
        cout << "=============================================================================\n";
        cout << "PHOTOPRODUCTION EVENT GENERATION\n";
        cout << "=============================================================================\n";
        cout << "Output File: " << outputFileName << "\n";
        config.print();
        cout << "=============================================================================\n\n";

        // =====================================================================
        // PYTHIA INITIALIZATION
        // =====================================================================
        
        Pythia pythia;
        Settings& settings = pythia.settings;
        const Info& info = pythia.info;

        // Reduce output verbosity
        pythia.readString("Init:showMultipartonInteractions = off");
        pythia.readString("Init:showChangedSettings = off");
        pythia.readString("Init:showChangedParticleData = off");
        pythia.readString("Next:numberCount = 100000");
        pythia.readString("Next:numberShowInfo = 0");
        pythia.readString("Next:numberShowProcess = 0");
        pythia.readString("Next:numberShowEvent = 0");

        // Beam configuration
        pythia.readString("Beams:frameType = 2");
        pythia.readString("Beams:idA = 2212");          // Proton
        pythia.readString("Beams:idB = 11");            // Electron
        pythia.readString("Beams:eA = " + to_string(config.protonEnergy));
        pythia.readString("Beams:eB = " + to_string(config.electronEnergy));
        pythia.readString("PDF:beamB2gamma = on");      // Enable photon beam

        // Event settings for photoproduction
        settings.mode("Photon:ProcessType", 0);         // Automatic mix of direct/resolved
        pythia.readString("Photon:Q2max = " + to_string(config.Q2max));
        // pythia.readString("MultipartonInteractions:pT0Ref = " + to_string(config.pT0Ref));
        pythia.readString("PhaseSpace:pTHatMin = " + to_string(config.pTHatMin));
        
        // Enable relevant processes
        pythia.readString("HardQCD:all = on");          // Resolved processes
        pythia.readString("PhotonParton:all = on");     // Direct processes

        // Initialize generator
        if (!pythia.init()) {
            throw runtime_error("Failed to initialize PYTHIA");
        }

        cout << "PYTHIA initialized successfully.\n\n";

        // =====================================================================
        // CREATE OUTPUT FILE
        // =====================================================================
        
        unique_ptr<TFile> outputFile(new TFile(outputFileName.c_str(), "RECREATE"));
        if (!outputFile || outputFile->IsZombie()) {
            throw runtime_error("Failed to create output file: " + outputFileName);
        }
        
        // Create directory structure
        TDirectory* dirQQ = outputFile->mkdir("QQ_Events");
        TDirectory* dirGG = outputFile->mkdir("GG_Events");
        TDirectory* dirGQ = outputFile->mkdir("GQ_Events");
        TDirectory* dirCombined = outputFile->mkdir("Combined_Events");
        TDirectory* dirResolved = outputFile->mkdir("Resolved_Events");
        TDirectory* dirDirect = outputFile->mkdir("Direct_Events");
        TDirectory* dirConfig = outputFile->mkdir("Analysis_Config");

        // =====================================================================
        // CREATE TREES FOR EVENT DATA
        // =====================================================================
        
        struct EventData {
            // Event information
            Int_t eventID, processCode, processType, photonMode;
            Bool_t isResolved, isDirect;
            Float_t crossSection, inelasticity;
            Float_t photonEnergy, scatteredElectronEnergy;
            
            // Final state particles
            vector<Float_t> px, py, pz, energy;
            vector<Float_t> pT, eta, phi, mass;
            vector<Int_t> pdgId, status;
            vector<Bool_t> isCharged, isHadron;
            
            // Scattered electron
            Float_t scElectron_px, scElectron_py, scElectron_pz, scElectron_E;
            Float_t scElectron_pT, scElectron_eta, scElectron_phi;
            
            void clear() {
                px.clear(); py.clear(); pz.clear(); energy.clear();
                pT.clear(); eta.clear(); phi.clear(); mass.clear();
                pdgId.clear(); status.clear();
                isCharged.clear(); isHadron.clear();
            }
        };
        
        auto createEventTree = [](TDirectory* dir, const string& name) -> pair<TTree*, EventData*> {
            dir->cd();
            TTree* tree = new TTree(name.c_str(), ("Events for " + name).c_str());
            EventData* data = new EventData();
            
            // Event information
            tree->Branch("eventID", &data->eventID);
            tree->Branch("processCode", &data->processCode);
            tree->Branch("processType", &data->processType);
            tree->Branch("photonMode", &data->photonMode);
            tree->Branch("isResolved", &data->isResolved);
            tree->Branch("isDirect", &data->isDirect);
            tree->Branch("crossSection", &data->crossSection);
            tree->Branch("inelasticity", &data->inelasticity);
            tree->Branch("photonEnergy", &data->photonEnergy);
            tree->Branch("scatteredElectronEnergy", &data->scatteredElectronEnergy);
            
            // Final state particles
            tree->Branch("px", &data->px);
            tree->Branch("py", &data->py);
            tree->Branch("pz", &data->pz);
            tree->Branch("energy", &data->energy);
            tree->Branch("pT", &data->pT);
            tree->Branch("eta", &data->eta);
            tree->Branch("phi", &data->phi);
            tree->Branch("mass", &data->mass);
            tree->Branch("pdgId", &data->pdgId);
            tree->Branch("status", &data->status);
            tree->Branch("isCharged", &data->isCharged);
            tree->Branch("isHadron", &data->isHadron);
            
            // Scattered electron
            tree->Branch("scElectron_px", &data->scElectron_px);
            tree->Branch("scElectron_py", &data->scElectron_py);
            tree->Branch("scElectron_pz", &data->scElectron_pz);
            tree->Branch("scElectron_E", &data->scElectron_E);
            tree->Branch("scElectron_pT", &data->scElectron_pT);
            tree->Branch("scElectron_eta", &data->scElectron_eta);
            tree->Branch("scElectron_phi", &data->scElectron_phi);
            
            return make_pair(tree, data);
        };
        
        // Create trees
        pair<TTree*, EventData*> qqPair = createEventTree(dirQQ, "QQ_Events");
        TTree* treeQQ = qqPair.first;
        EventData* dataQQ = qqPair.second;
        
        pair<TTree*, EventData*> ggPair = createEventTree(dirGG, "GG_Events");
        TTree* treeGG = ggPair.first;
        EventData* dataGG = ggPair.second;
        
        pair<TTree*, EventData*> gqPair = createEventTree(dirGQ, "GQ_Events");
        TTree* treeGQ = gqPair.first;
        EventData* dataGQ = gqPair.second;
        
        pair<TTree*, EventData*> combinedPair = createEventTree(dirCombined, "Combined_Events");
        TTree* treeCombined = combinedPair.first;
        EventData* dataCombined = combinedPair.second;
        
        pair<TTree*, EventData*> resolvedPair = createEventTree(dirResolved, "Resolved_Events");
        TTree* treeResolved = resolvedPair.first;
        EventData* dataResolved = resolvedPair.second;
        
        pair<TTree*, EventData*> directPair = createEventTree(dirDirect, "Direct_Events");
        TTree* treeDirect = directPair.first;
        EventData* dataDirect = directPair.second;

        // =====================================================================
        // EVENT COUNTERS
        // =====================================================================
        
        int totalEvents = 0;
        int validEvents = 0;
        int qqEvents = 0, ggEvents = 0, gqEvents = 0, unclassifiedEvents = 0;
        int resolvedEvents = 0, directEvents = 0;
        
        double totalCrossSection = 0.0;
        double qqCrossSection = 0.0, ggCrossSection = 0.0, gqCrossSection = 0.0;
        
        cout << "Starting event generation...\n";
        cout << "Progress: ";

        // =====================================================================
        // EVENT GENERATION LOOP
        // =====================================================================
        
        for (int iEvent = 0; iEvent < config.nEvents; ++iEvent) {
            
            // Progress indicator
            if (iEvent % (config.nEvents/10) == 0) {
                cout << (100 * iEvent / config.nEvents) << "% ";
                cout.flush();
            }
            
            totalEvents++;
            
            // Generate next event
            if (!pythia.next()) continue;
            
            // Debug first event
            if (iEvent == 0) {
                cout << "\n--- First Event Debug ---\n";
                pythia.process.list();
                cout << "------------------------\n";
            }
            
            validEvents++;
            
            // Extract event information
            int processCode = pythia.info.code();
            int processType = classifySubprocess(processCode);
            int photonMode = info.photonMode();
            double crossSection = pythia.info.weight();
            totalCrossSection += crossSection;
            
            bool isResolved = (photonMode == 1);
            bool isDirect = (photonMode == 2);
            
            // Update counters
            if (processType == 1) {
                qqEvents++;
                qqCrossSection += crossSection;
            } else if (processType == 2) {
                ggEvents++;
                ggCrossSection += crossSection;
            } else if (processType == 3) {
                gqEvents++;
                gqCrossSection += crossSection;
            } else {
                unclassifiedEvents++;
            }
            
            if (isResolved) resolvedEvents++;
            if (isDirect) directEvents++;
            
            // Helper function to fill event data
            auto fillEventData = [&](EventData* data) {
                data->clear();
                
                // Event information
                data->eventID = iEvent;
                data->processCode = processCode;
                data->processType = processType;
                data->photonMode = photonMode;
                data->isResolved = isResolved;
                data->isDirect = isDirect;
                data->crossSection = crossSection;
                
                // Initialize electron and photon energies
                data->photonEnergy = 0.0;
                data->scatteredElectronEnergy = 0.0;
                data->inelasticity = 0.0;
                
                // Initialize scattered electron variables
                data->scElectron_px = data->scElectron_py = data->scElectron_pz = 0.0;
                data->scElectron_E = data->scElectron_pT = data->scElectron_eta = data->scElectron_phi = 0.0;
                
                // Extract particle information
                for (int i = 0; i < pythia.event.size(); ++i) {
                    const Particle& particle = pythia.event[i];
                    
                    // Find incoming beam photon energy
                    if (particle.id() == 22 && particle.status() == -13) {
                        data->photonEnergy = particle.e();
                    }
                    
                    // Find scattered beam electron
                    if (particle.id() == 11 && particle.mother1() == 0 && particle.mother2() == 0) {
                        data->scElectron_px = particle.px();
                        data->scElectron_py = particle.py();
                        data->scElectron_pz = particle.pz();
                        data->scElectron_E = particle.e();
                        data->scElectron_pT = particle.pT();
                        data->scElectron_eta = particle.eta();
                        data->scElectron_phi = particle.phi();
                        data->scatteredElectronEnergy = particle.e();
                        
                        // Calculate inelasticity
                        if (config.electronEnergy > 0) {
                            data->inelasticity = (config.electronEnergy - data->scatteredElectronEnergy) / config.electronEnergy;
                        }
                        continue;
                    }
                    
                    // Store final state particles
                    if (particle.isFinal()) {
                        data->px.push_back(particle.px());
                        data->py.push_back(particle.py());
                        data->pz.push_back(particle.pz());
                        data->energy.push_back(particle.e());
                        data->pT.push_back(particle.pT());
                        data->eta.push_back(particle.eta());
                        data->phi.push_back(particle.phi());
                        data->mass.push_back(particle.m());
                        data->pdgId.push_back(particle.id());
                        data->status.push_back(particle.status());
                        data->isCharged.push_back(particle.isCharged());
                        data->isHadron.push_back(particle.isHadron());
                    }
                }
            };
            
            // Fill combined tree (all events)
            fillEventData(dataCombined);
            treeCombined->Fill();
            
            // Fill subprocess-specific trees
            if (processType == 1) {
                fillEventData(dataQQ);
                treeQQ->Fill();
            } else if (processType == 2) {
                fillEventData(dataGG);
                treeGG->Fill();
            } else if (processType == 3) {
                fillEventData(dataGQ);
                treeGQ->Fill();
            }
            
            // Fill photoproduction-type trees
            if (isResolved) {
                fillEventData(dataResolved);
                treeResolved->Fill();
            }
            if (isDirect) {
                fillEventData(dataDirect);
                treeDirect->Fill();
            }
        }
        
        cout << "100%\n\n";

        // =====================================================================
        // SAVE CONFIGURATION AND STATISTICS
        // =====================================================================
        
        dirConfig->cd();
        TTree* configTree = new TTree("GenerationConfig", "Event Generation Configuration");
        
        Double_t protonEnergy = config.protonEnergy;
        Double_t electronEnergy = config.electronEnergy;
        Double_t pTHatMin = config.pTHatMin;
        Double_t Q2max = config.Q2max;
        // Double_t pT0Ref = config.pT0Ref;
        Int_t nEventsGenerated = config.nEvents;
        Int_t validEventsGenerated = validEvents;
        
        configTree->Branch("protonEnergy", &protonEnergy);
        configTree->Branch("electronEnergy", &electronEnergy);
        configTree->Branch("pTHatMin", &pTHatMin);
        configTree->Branch("Q2max", &Q2max);
        // configTree->Branch("pT0Ref", &pT0Ref);
        configTree->Branch("nEventsGenerated", &nEventsGenerated);
        configTree->Branch("validEventsGenerated", &validEventsGenerated);
        
        configTree->Fill();
        
        // Create summary histogram
        TH1F* hEventCounts = new TH1F("hEventCounts", 
            "Event Counts by Process Type;Process Type;Number of Events", 6, 0, 6);
        hEventCounts->GetXaxis()->SetBinLabel(1, "QQ");
        hEventCounts->GetXaxis()->SetBinLabel(2, "GG");
        hEventCounts->GetXaxis()->SetBinLabel(3, "GQ");
        hEventCounts->GetXaxis()->SetBinLabel(4, "Unclassified");
        hEventCounts->GetXaxis()->SetBinLabel(5, "Resolved");
        hEventCounts->GetXaxis()->SetBinLabel(6, "Direct");
        
        hEventCounts->SetBinContent(1, qqEvents);
        hEventCounts->SetBinContent(2, ggEvents);
        hEventCounts->SetBinContent(3, gqEvents);
        hEventCounts->SetBinContent(4, unclassifiedEvents);
        hEventCounts->SetBinContent(5, resolvedEvents);
        hEventCounts->SetBinContent(6, directEvents);

        // =====================================================================
        // FINAL STATISTICS
        // =====================================================================
        
        pythia.stat();
        
        cout << "\n=============================================================================\n";
        cout << "EVENT GENERATION SUMMARY\n";
        cout << "=============================================================================\n";
        cout << "Total Events Generated:    " << totalEvents << "\n";
        cout << "Valid Events:              " << validEvents << "\n";
        cout << "\nEvent Classification:\n";
        cout << "  QQ Events:               " << qqEvents 
             << " (" << setprecision(1) << (100.0 * qqEvents / validEvents) << "%)\n";
        cout << "  GG Events:               " << ggEvents 
             << " (" << (100.0 * ggEvents / validEvents) << "%)\n";
        cout << "  GQ Events:               " << gqEvents 
             << " (" << (100.0 * gqEvents / validEvents) << "%)\n";
        cout << "  Unclassified:            " << unclassifiedEvents 
             << " (" << (100.0 * unclassifiedEvents / validEvents) << "%)\n";
        cout << "\nPhotoproduction Type:\n";
        cout << "  Resolved Events:         " << resolvedEvents 
             << " (" << (100.0 * resolvedEvents / validEvents) << "%)\n";
        cout << "  Direct Events:           " << directEvents 
             << " (" << (100.0 * directEvents / validEvents) << "%)\n";
        
        // Cross section validation
        cout << scientific << setprecision(3);
        cout << "\nCross Section Information:\n";
        cout << "  Total Cross Section:     " << (totalCrossSection / validEvents) << " mb\n";
        cout << "  QQ Cross Section:        " << (qqCrossSection / validEvents) 
             << " mb (" << fixed << setprecision(1) 
             << (100.0 * qqCrossSection / totalCrossSection) << "%)\n";
        cout << "  GG Cross Section:        " << scientific << (ggCrossSection / validEvents) 
             << " mb (" << fixed << (100.0 * ggCrossSection / totalCrossSection) << "%)\n";
        cout << "  GQ Cross Section:        " << scientific << (gqCrossSection / validEvents) 
             << " mb (" << fixed << (100.0 * gqCrossSection / totalCrossSection) << "%)\n";
        
        // Validate percentages
        double totalPercentage = (100.0 * qqCrossSection / totalCrossSection) + 
                                (100.0 * ggCrossSection / totalCrossSection) + 
                                (100.0 * gqCrossSection / totalCrossSection);
        cout << "  Total Percentage:        " << totalPercentage << "%";
        if (abs(totalPercentage - 100.0) < 0.1) {
            cout << " (OK)\n";
        } else {
            cout << " (Check for issues)\n";
        }

        // =====================================================================
        // WRITE OUTPUT
        // =====================================================================
        
        cout << "\nWriting output file...\n";
        outputFile->cd();
        outputFile->Write();
        
        cout << "\n=============================================================================\n";
        cout << "EVENT GENERATION COMPLETED SUCCESSFULLY\n";
        cout << "=============================================================================\n";
        cout << "Output File: " << outputFileName << "\n";
        cout << "\nTrees Created:\n";
        cout << "  QQ_Events       : " << setw(8) << treeQQ->GetEntries() << " events\n";
        cout << "  GG_Events       : " << setw(8) << treeGG->GetEntries() << " events\n";
        cout << "  GQ_Events       : " << setw(8) << treeGQ->GetEntries() << " events\n";
        cout << "  Combined_Events : " << setw(8) << treeCombined->GetEntries() << " events\n";
        cout << "  Resolved_Events : " << setw(8) << treeResolved->GetEntries() << " events\n";
        cout << "  Direct_Events   : " << setw(8) << treeDirect->GetEntries() << " events\n";
        
        cout << "\nData ready for analysis!\n";
        cout << "Each tree contains final state particles with complete 4-momentum info.\n";
        cout << "=============================================================================\n";
        
        // =====================================================================
        // CLEANUP
        // =====================================================================
        
        delete dataQQ;
        delete dataGG;
        delete dataGQ;
        delete dataCombined;
        delete dataResolved;
        delete dataDirect;
        
        cout << "\nEvent generation completed successfully!\n";
        
    } catch (const exception& e) {
        cerr << "\nERROR: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "\nUNKNOWN ERROR occurred during event generation!" << endl;
        return 1;
    }
    
    return 0;
}