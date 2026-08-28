

#include <string>

//TString sim_path = "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/EPCollisions.edm4eic.root";
TString sim_path = "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/EPCollisionsMuonGun.edm4eic.root.root";

void ResolutionAnalysis(TString infile=sim_path){
  // Set output file for the histograms
  TFile *ofile = TFile::Open("copper_eic_hybrid/ResolutionTesting.root","RECREATE");

  // Analysis code will go here
  // Set up input file chain
  TChain *mychain = new TChain("events");
  mychain->Add(infile);

  // Initialize reader
  TTreeReader tree_reader(mychain);

  // Get Particle Information
  TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
  TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
  TTreeReaderArray<double> partVertexX(tree_reader, "MCParticles.vertex.x");
  TTreeReaderArray<double> partVertexY(tree_reader, "MCParticles.vertex.y");
  TTreeReaderArray<double> partVertexZ(tree_reader, "MCParticles.vertex.z");
  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
  TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
  TTreeReaderArray<float> trackVertexX(tree_reader, "ReconstructedChargedParticles.referencePoint.x");
  TTreeReaderArray<float> trackVertexY(tree_reader, "ReconstructedChargedParticles.referencePoint.y");
  TTreeReaderArray<float> trackVertexZ(tree_reader, "ReconstructedChargedParticles.referencePoint.z");
  // Get Associations Between MCParticles and ReconstructedChargedParticles
  TTreeReaderArray<int> recoAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_rec.index");
  TTreeReaderArray<int> simuAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_sim.index");
    
  // Define Histograms
  TH1D *trackMomentumRes = new TH1D("trackMomentumRes","Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  trackMomentumRes->GetXaxis()->SetRangeUser(-1, 1);
  TH2D* trackMomResP = new TH2D("trackMomResP", "Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -1, 1, 150, 0, 150);
  trackMomResP->GetXaxis()->SetRangeUser(-1, 1);
  TH2D* trackMomResEta = new TH2D("trackMomResEta", "Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -1, 1, 120, -6, 6);
  trackMomResEta->GetXaxis()->SetRangeUser(-1, 1);
	TH2D* trackResVertexX = new TH2D("trackResVertexX", "Track Vertex X Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -1, 1, 120, -6, 6);
  TH1D *trackMomentumRes_e = new TH1D("trackMomentumRes_e","e^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_e = new TH2D("trackMomResP_e", "e^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 25);
  TH2D* trackMomResEta_e = new TH2D("trackMomResEta_e", "e^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);


  TH1D *trackMomentumRes_mu = new TH1D("trackMomentumRes_mu","#mu^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -1, 1);
  TH2D* trackMomResP_mu = new TH2D("trackMomResP_mu", "#mu^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -1, 1, 150, 0, 25);
  TH2D* trackMomResEta_mu = new TH2D("trackMomResEta_mu", "#mu^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -1, 1, 120, -6, 6);

  TH1D *trackMomentumRes_pi = new TH1D("trackMomentumRes_pi","#pi^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_pi = new TH2D("trackMomResP_pi", "#pi^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
  TH2D* trackMomResEta_pi = new TH2D("trackMomResEta_pi", "#pi^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

  TH1D *trackMomentumRes_K = new TH1D("trackMomentumRes_K","K^{#pm} Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_K = new TH2D("trackMomResP_K", "K^{#pm} Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
  TH2D* trackMomResEta_K = new TH2D("trackMomResEta_K", "K^{#pm} Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);

  TH1D *trackMomentumRes_p = new TH1D("trackMomentumRes_p","p Track Momentum Resolution; (P_{rec} - P_{MC})/P_{MC}", 400, -2, 2);
  TH2D* trackMomResP_p = new TH2D("trackMomResP_p", "p Track Momentum Resolution vs P; (P_{rec} - P_{MC})/P_{MC}; P_{MC}(GeV/c)", 400, -2, 2, 150, 0, 150);
  TH2D* trackMomResEta_p = new TH2D("trackMomResEta_p", "p Track Momentum Resolution vs #eta; (P_{rec} - P_{MC})/P_{MC}; #eta_{MC}", 400, -2, 2, 120, -6, 6);
  
  TH1D *matchedPartTrackDeltaEta = new TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructed Charged Particle; #Delta#eta", 100, -0.25, 0.25);
  TH1D *matchedPartTrackDeltaPhi = new TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3);
  TH1D *matchedPartTrackDeltaMom = new TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10);

  while(tree_reader.Next()) { // Loop over events

    for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over thrown particles
      {
	if(partGenStat[i] == 1) // Select stable thrown particles
	  {
	    int pdg = TMath::Abs(partPdg[i]);

	    if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212) // Look at charged particles (electrons, muons, pions, kaons, protons)
	      {
		TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);

		float trueEta = trueMom.PseudoRapidity();
		float truePhi = trueMom.Phi();

		// Loop over associations to find matching ReconstructedChargedParticle
		for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
		  {
		    if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
		      {
			TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

			// Check the distance between the thrown and reconstructed particle
			float deltaEta = trueEta - recMom.PseudoRapidity();
			float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recMom.Phi());
			float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
			float deltaMom = ((trueMom.Mag()) - (recMom.Mag()));

			float deltaPrimaryX = (partVertexX - trackVertexX);
			float deltaPrimaryY = (partVertexY- trackVertexY);
			float deltaPrimaryZ = (partVertexZ - trackVertexZ);

			double momRes = (recMom.Mag() - trueMom.Mag())/trueMom.Mag();
	
			trackMomentumRes->Fill(momRes); // Could also multiply by 100 and express as a percentage instead
			trackMomResP->Fill(momRes, trueMom.Mag());
			trackMomResEta->Fill(momRes, trueEta);

			if( pdg == 11){
			  trackMomentumRes_e->Fill(momRes);
			  trackMomResP_e->Fill(momRes, trueMom.Mag());
			  trackMomResEta_e->Fill(momRes, trueEta);
			}
			else if( pdg == 13){
			  trackMomentumRes_mu->Fill(momRes);
			  trackMomResP_mu->Fill(momRes, trueMom.Mag());
			  trackMomResEta_mu->Fill(momRes, trueEta);
			}
			else if( pdg == 211){
			  trackMomentumRes_pi->Fill(momRes);
			  trackMomResP_pi->Fill(momRes, trueMom.Mag());
			  trackMomResEta_pi->Fill(momRes, trueEta);
			}
			else if( pdg == 321){
			  trackMomentumRes_K->Fill(momRes);
			  trackMomResP_K->Fill(momRes, trueMom.Mag());
			  trackMomResEta_K->Fill(momRes, trueEta);
			}
			else if( pdg == 2212){
			  trackMomentumRes_p->Fill(momRes);
			  trackMomResP_p->Fill(momRes, trueMom.Mag());
			  trackMomResEta_p->Fill(momRes, trueEta);
			}
			  
			matchedPartTrackDeltaEta->Fill(deltaEta);
			matchedPartTrackDeltaPhi->Fill(deltaPhi);
			matchedPartTrackDeltaR->Fill(deltaR);
			matchedPartTrackDeltaMom->Fill(deltaMom);
			
		      }
		  }// End loop over associations
	      } // End PDG check
	  } // End stable particles condition
      } // End loop over thrown particles
  } // End loop over events

  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
}