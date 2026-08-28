
//Goal: track efficiency and resolution in C code. Need to change Tstring sim_path to be path of the eicrecon file
// Track Efficiency is defined as which particles generated have been truth-matched (Associated). Also takes only the charged partilces 
// within the eta region of the detector (|eta| < 4)
//https://eic.github.io/tutorial-analysis/03-analysis.html

#include <string>

//TString sim_path = "root://hpceph-xrootd.twgrid.org:1094//cephfs/epic//RECO/26.07.1/epic_craterlake/Bkg_Exact1S_2us/GoldCt/10um/DIS/NC/10x100/minQ2=1/pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.2638.eicrecon.edm4eic.root";

//TString sim_path = "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/EPCollisions.edm4eic.root";

//TString sim_path = "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/MuonsFullCraterlakeCutGDML.edm4eic.root";

TString sim_path = "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/MuonsFullCraterlakeCutGDMLMuon10GeV.edm4eic.root";

void TrackEfficiency(TString infile=sim_path){
  
  // Set output file for the histograms
  TFile *ofile = TFile::Open("Plots2808/TrackEfficiencyMuonsFullCraterlakeCutGDMLMuon10GeVNoweight.root","RECREATE");
  cout <<"Graphs of efficiencies only consider the particles within the eta range of the detector" << endl;
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
  // Get Reconstructed Track Information
  TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
  TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
  // Get Associations Between MCParticles and ReconstructedChargedParticles
  TTreeReaderArray<int> recoAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_rec.index");
  TTreeReaderArray<int> simuAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_sim.index");
  TTreeReaderArray<float> weightAssoc(tree_reader, "ReconstructedChargedParticleAssociations.weight");
  TTreeReaderArray<unsigned int> partParentsBegin(tree_reader, "MCParticles.parents_begin");

TTreeReaderArray<unsigned int> partParentsEnd(tree_reader, "MCParticles.parents_end");

TTreeReaderArray<int> parentIndices(tree_reader, "_MCParticles_parents.index");
  
  // Define Histograms
  TH1D *partEta = new TH1D("partEta","#eta of Thrown Charged Particles; #eta", 120, -10, 10);
  TH1D *matchedPartEta = new TH1D("matchedPartEta","#eta of Thrown Charged Particles That Have Matching Track; #eta", 120, -10, 10);
  TH1D* partMom = new TH1D("partMom", "Momentum of Thrown Charged Particles (truth); P(GeV/c)", 100, 0, 100);
  partMom->GetXaxis()->SetRangeUser(0, 40);
  TH1D* matchedPartMom = new TH1D("matchedPartMom", "Momentum of Thrown Charged Particles (truth), with matching track; P(GeV/c)", 100, 0, 100);
  matchedPartMom->GetXaxis()->SetRangeUser(0, 40);
  TH1D* partPhi = new TH1D("partPhi", "#phi of Thrown Charged Particles (truth); #phi(rad)", 320, -3.2, 3.2);
  TH1D* matchedPartPhi = new TH1D("matchedPartPhi", "#phi of Thrown Charged Particles (truth), with matching track; #phi(rad)", 320, -3.2, 3.2);
  TH1D* partMomT =  new TH1D("partMomT", "Transverse Momentum of Thrown Charged Particles (truth); P(GeV/c)", 100, 0, 5);
  TH1D* matchedPartMomT = new TH1D("matchedPartMomT", "Transverse Momentum of Thrown Charged Particles (truth), with matching track; P(GeV/c)", 100, 0, 5);

  TH1D* partElecMomAll = new TH1D("partElecMomAll", "P of Thrown Electron Particles (truth); P", 50, 0, 30);
  TH1D* partKaonMomAll = new TH1D("partKaonMomAll", "P of Thrown Kaon Particles (truth); P", 50, 0, 30);
  TH1D* partPionMomAll = new TH1D("partPionMomAll", "P of Thrown Pion Particles (truth); P", 50, 0, 30);
  TH1D* partMuonMomAll = new TH1D("partMuonMomAll", "P of Thrown Muon Particles (truth); P", 50, 0, 20);

  //Histograms for particles which are daughters
  TH1D* partElecMomDaughter = new TH1D("partElecMomDaughter", "Electrons which have a different Parent; P", 20, 0, 30);
  TH1D* partKaonMomDaughter = new TH1D("partKaonMomDaughter", "Kaons which have a different Parent; P", 20, 0, 30);
  TH1D* partPionMomDaughter = new TH1D("partPionMomDaughter", "Pions which have a different Parent; P", 20, 0, 30);
  TH1D* partMuonMomDaughter = new TH1D("partMuonMomDaughter", "Muons which have a different Parent; P", 20, 0, 30);

  TH1D* matchedPartElecMomDaughter = new TH1D("matchedPartElecMomDaughter", "Electrons which have a different Parent (matched); P", 20, 0, 30);
  TH1D* matchedPartKaonMomDaughter = new TH1D("matchedPartKaonMomDaughter", "Kaons which have a different Parent (matched); P", 20, 0, 30);
  TH1D* matchedPartPionMomDaughter = new TH1D("matchedPartPionMomDaughter", "Pions which have a different Parent (matched); P", 20, 0, 30);
  TH1D* matchedPartMuonMomDaughter = new TH1D("matchedPartMuonMomDaughter", "Muons which have a different Parent (matched); P", 20, 0, 30);

  //Histograms for particles which are not daughters
  TH1D* partElecMom = new TH1D("partElecMom", "Electrons which are not daughters ; P", 20, 0, 30);
  TH1D* partKaonMom= new TH1D("partKaonMom", "Kaons which are not daughters ; P", 20, 0, 30);
  TH1D* partPionMom = new TH1D("partPionMom", "Pions which are not daughters ; P", 20, 0, 30);
  TH1D* partMuonMom = new TH1D("partMuonMom", "Muons which are not daughters ; P", 20, 0, 30);

  TH1D* matchedPartElecMom = new TH1D("matchedPartElecMom", "Electrons which have a different Parent (matched); P", 20, 0, 30);
  TH1D* matchedPartKaonMom = new TH1D("matchedPartKaonMom", "Kaons which have a different Parent (matched); P", 20, 0, 30);
  TH1D* matchedPartPionMom = new TH1D("matchedPartPionMom", "Pions which have a different Parent (matched); P", 20, 0, 30);
  TH1D* matchedPartMuonMom = new TH1D("matchedPartMuonMom", "Muons which have a different Parent (matched); P", 20, 0, 30);


  TH2D* partPEta = new TH2D("partPEta", "P vs #eta of Thrown Charged Particles; P(GeV/c); #eta", 100, 0, 100, 120, -6, 6);
  TH2D* matchedPartPEta = new TH2D("matchedPartPEta", "P vs #eta of Thrown Charged Particles, with matching track; P(GeV/c); #eta", 100, 0, 100, 120, -6, 6);
  TH2D* partPhiEta = new TH2D("partPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);
  TH2D* matchedPartPhiEta = new TH2D("matchedPartPhiEta", "#phi vs #eta of Thrown Charged Particles; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);
  

  TH1D *matchedPartTrackDeltaEta = new TH1D("matchedPartTrackDeltaEta","#Delta#eta Between Matching Thrown and Reconstructed Charged Particle; #Delta#eta", 100, -0.25, 0.25);
  TH1D *matchedPartTrackDeltaPhi = new TH1D("matchedPartTrackDeltaPhi","#Detla #phi Between Matching Thrown and Reconstructed Charged Particle; #Delta#phi", 200, -0.2, 0.2);
  TH1D *matchedPartTrackDeltaR = new TH1D("matchedPartTrackDeltaR","#Delta R Between Matching Thrown and Reconstructed Charged Particle; #Delta R", 300, 0, 0.3);
  TH1D *matchedPartTrackDeltaMom = new TH1D("matchedPartTrackDeltaMom","#Delta P Between Matching Thrown and Reconstructed Charged Particle; #Delta P", 200, -10, 10);
  
  TH1D* matchedPartElecAll = new TH1D("matchedPartElec", "P of Thrown Electron Particles (matched); P", 50, 0, 30);
  TH1D* matchedPartKaonAll = new TH1D("matchedPartKaon", "P of Thrown Kaon Particles (matched); P", 50, 0, 30);
  TH1D* matchedPartPionAll = new TH1D("matchedPartPion", "P of Thrown Pion Particles (matched); P", 50, 0, 30);
  TH1D* matchedPartMuonAll = new TH1D("matchedPartMuon", "P of Thrown Muon Particles (matched); P", 50, 0, 20);
  
  // Define some histograms for our efficiencies
  TH1D *TrackEff_Eta = new TH1D("TrackEff_Eta", "Tracking efficiency as fn of #eta; #eta; Eff(%)", 120, -10, 10); 
  TrackEff_Eta->GetXaxis()->SetRangeUser(-3, 3);
  TH1D *TrackEff_Mom = new TH1D("TrackEff_Mom", "Tracking efficiency as fn of P; P(GeV/c); Eff(%)", 100, 0, 100); 
  TrackEff_Mom->GetXaxis()->SetRangeUser(0, 40);
  TH1D *TrackEff_Phi = new TH1D("TrackEff_Phi", "Tracking efficiency as fn of #phi; #phi(rad); Eff(%)", 320, -3.2, 3.2);
  TH1D *TrackEff_MomT = new TH1D("TrackEff_MomT", "Tracking efficiency as fn of PT; P(GeV/c); Eff(%)", 100, 0, 5); 
  TH1D *TrackEff_Elec_All = new TH1D("TrackEff_Elec_All", "Tracking efficiency for all Electrons as fn of P; P(GeV/c); Eff(%)", 50, 0, 30);
  TH1D *TrackEff_Kaon_All = new TH1D("TrackEff_Kaon_All", "Tracking efficiency for all Kaons as fn of P; P(GeV/c); Eff(%)", 50, 0, 30);
  TH1D *TrackEff_Pion_All = new TH1D("TrackEff_Pion_All", "Tracking efficiency for all Pions as fn of P; P(GeV/c); Eff(%)", 50, 0, 30);
  TH1D *TrackEff_Muon_All = new TH1D("TrackEff_Muon_All", "Tracking efficiency for all Muons as fn of P; P(GeV/c); Eff(%)", 50, 0, 20);


  TH1D *TrackEff_Elec_daughter = new TH1D("TrackEff_Elec_daughter", "Tracking efficiency for daughter Electrons as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  TH1D *TrackEff_Kaon_daughter = new TH1D("TrackEff_Kaon_daughter", "Tracking efficiency for daughter Kaons as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  TH1D *TrackEff_Pion_daughter = new TH1D("TrackEff_Pion_daughter", "Tracking efficiency for daughter Pions as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  TH1D *TrackEff_Muon_daughter = new TH1D("TrackEff_Muon_daughter", "Tracking efficiency for daughter Muons as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  
  
  TH1D *TrackEff_Elec = new TH1D("TrackEff_Elec", "Tracking efficiency for (not daughters) Electrons as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  TH1D *TrackEff_Kaon = new TH1D("TrackEff_Kaon", "Tracking efficiency for (not daughters) Kaons as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  TH1D *TrackEff_Pion = new TH1D("TrackEff_Pion", "Tracking efficiency for (not daughters) Pions as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  TH1D *TrackEff_Muon = new TH1D("TrackEff_Muon", "Tracking efficiency for (not daughters) Muons as fn of P; P(GeV/c); Eff(%)", 20, 0, 30);
  
  //Distribution of eta vs no. of Kaons (truthmatched)
  TH1D *partEtaElectron = new TH1D("EtaElectron", "Distribution of thrown Electron particles ; #eta", 50, -10, 10);
  TH1D *partEtaKaon = new TH1D("EtaKaon", "Distribution of thrown Kaon particles ; #eta", 50, -10, 10);
  TH1D *partEtaPion = new TH1D("EtaPion", "Distribution of thrown Pion particles ; #eta", 50, -10, 10);
  TH1D *partEtaMuon = new TH1D("EtaMuon", "Distribution of thrown Muon particles ; #eta", 50, -10, 10);

  TH1D *matchedEtaElectron = new TH1D("matchedEtaElectron", "Distribution of matched thrown Electron particles ; #eta", 50, -10, 10);
  TH1D *matchedEtaKaon = new TH1D("matchedEtaKaon", "Distribution of matched thrown Kaon particles ; #eta", 50, -10, 10);
  TH1D *matchedEtaPion = new TH1D("matchedEtaPion", "Distribution of matched thrown Pion particles ; #eta", 50, -10, 10);
  TH1D *matchedEtaMuon = new TH1D("matchedEtaMuon", "Distribution of matched thrown Muon particles ; #eta", 50, -10, 10);

  TH1D *TrackEff_Elec_eta = new TH1D("TrackEff_Elec_eta", "Tracking efficiency for Electrons as fn of #eta; P(GeV/c); Eff(%)", 50, -10, 10);
  TH1D *TrackEff_Kaon_eta = new TH1D("TrackEff_Kaon_eta", "Tracking efficiency for Kaons as fn of #eta; P(GeV/c); Eff(%)", 50, -10, 10);
  TH1D *TrackEff_Pion_eta = new TH1D("TrackEff_Pion_eta", "Tracking efficiency for Pions as fn of #eta; P(GeV/c); Eff(%)", 50, -10, 10);
  TH1D *TrackEff_Muon_eta = new TH1D("TrackEff_Muon_eta", "Tracking efficiency for Muons as fn of #eta; P(GeV/c); Eff(%)", 50, -10, 10);
  
  
  // 2D Efficiencies
  TH2D* TrackEff_PEta = new TH2D("TrackEff_PEta", "Tracking efficiency as fn of P and #eta; P(GeV/c); #eta", 100, 0, 100, 120, -6, 6);
  TH2D* TrackEff_PhiEta = new TH2D("TrackEff_PhiEta", "Tracking efficiency as fn of #phi and #eta; #phi(rad); #eta", 160, -3.2, 3.2, 120, -6, 6);

  // All charged particle histos
  TH1D *ChargedEta = new TH1D("ChargedEta", "#eta of all charged particles; #eta", 120, -6, 6);
  TH1D *ChargedPhi = new TH1D("ChargedPhi", "#phi of all charged particles; #phi (rad)", 120, -3.2, 3.2);
  TH1D *ChargedP = new TH1D("ChargedP", "P of all charged particles; P(GeV/c)", 100, 0, 100);
  
  int numberParticleOutsideEtaRange = 0;
  int totalNumberOfChargedParticles = 0;
  while(tree_reader.Next()) { // Loop over events
    int numberOfParticleForEvent = 0;
    for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over all MC particles
    //for (unsigned int i=0; i<10; i++)
      {
        
	if(partGenStat[i] == 1) // Select stable thrown particles
	  {
      numberOfParticleForEvent +=1;
	    int pdg = TMath::Abs(partPdg[i]);

	    if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212) // Look at charged particles (electrons, muons, pions, kaons, protons)
	      { 
		TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);
    
    ++totalNumberOfChargedParticles;
    float truePT =  sqrt(partMomX[i] * partMomX[i] + partMomY[i]*partMomY[i]);
		float trueEta = trueMom.PseudoRapidity();
		float truePhi = trueMom.Phi();
	    
	

    if (trueEta > 4 || trueEta < -4){
      ++numberParticleOutsideEtaRange;
      continue;
    }

    partEta->Fill(trueEta);
		partPhi->Fill(truePhi);
		partMom->Fill(trueMom.Mag());
    partMomT->Fill(truePT);
		partPEta->Fill(trueMom.Mag(), trueEta);
		partPhiEta->Fill(truePhi, trueEta);

    bool isDaughter = false;
    if (pdg == 11){
      partElecMomAll->Fill(trueMom.Mag());
      partEtaElectron->Fill(trueMom.PseudoRapidity());
      isDaughter = partParentsEnd[i] > partParentsBegin[i];
      
      if (isDaughter){  
        
        bool hasDifferentParent = false;
        // A particle can have more than one parent
        for (unsigned int relation = partParentsBegin[i];
             relation < partParentsEnd[i];
             ++relation) {

            int parentIndex = parentIndices[relation];

            if (TMath::Abs(partPdg[parentIndex]) != TMath::Abs(partPdg[i])) {
                hasDifferentParent = true;
                partElecMomDaughter->Fill(trueMom.Mag());

               // std::cout << "  Parent index: " << parentIndex
                 //     << ", parent PDG: " << partPdg[parentIndex]
                   //   << std::endl;
        }
      }
    }

        else if (!isDaughter){
            partElecMom->Fill(trueMom.Mag());
    }
  }

    if (pdg == 321){ //Positive or Neg Kaons
      partKaonMomAll->Fill(trueMom.Mag());
      partEtaKaon->Fill(trueMom.PseudoRapidity());
      isDaughter = partParentsEnd[i] > partParentsBegin[i];
      if (isDaughter){  
        
        bool hasDifferentParent = false;
        // A particle can have more than one parent
        for (unsigned int relation = partParentsBegin[i];
             relation < partParentsEnd[i];
             ++relation) {

            int parentIndex = parentIndices[relation];

            if (TMath::Abs(partPdg[parentIndex]) != TMath::Abs(partPdg[i])) {
                hasDifferentParent = true;
                partKaonMomDaughter->Fill(trueMom.Mag());

                //std::cout << "  Parent index: " << parentIndex
                  //    << ", parent PDG: " << partPdg[parentIndex]
                    //  << std::endl;
        }
      }
    }

        else if (!isDaughter){
            partKaonMom->Fill(trueMom.Mag());
    }
    }

    if (pdg == 211){ //Positive or Neg Pions
      partPionMomAll->Fill(trueMom.Mag());
      partEtaPion->Fill(trueMom.PseudoRapidity());
      isDaughter = partParentsEnd[i] > partParentsBegin[i];
      if (isDaughter){  
        
        bool hasDifferentParent = false;
        // A particle can have more than one parent
        for (unsigned int relation = partParentsBegin[i];
             relation < partParentsEnd[i];
             ++relation) {

            int parentIndex = parentIndices[relation];

            if (TMath::Abs(partPdg[parentIndex]) != TMath::Abs(partPdg[i])) {
                hasDifferentParent = true;
                partPionMomDaughter->Fill(trueMom.Mag());
              //  std::cout << "PION : 211" << endl;
                //std::cout << "  Parent index: " << parentIndex
                  //    << ", parent PDG: " << partPdg[parentIndex]
                    //  << std::endl;
        }
      }
    }
     else if (!isDaughter){
        partPionMom->Fill(trueMom.Mag());
    }
    }

    if (pdg == 13){
      partMuonMomAll->Fill(trueMom.Mag());
      partEtaMuon->Fill(trueMom.PseudoRapidity());
      isDaughter = partParentsEnd[i] > partParentsBegin[i];
      if (isDaughter){  
        
        bool hasDifferentParent = false;
        // A particle can have more than one parent
        for (unsigned int relation = partParentsBegin[i];
             relation < partParentsEnd[i];
             ++relation) {

            int parentIndex = parentIndices[relation];

            if (TMath::Abs(partPdg[parentIndex]) != TMath::Abs(partPdg[i])) {
                hasDifferentParent = true;
                partMuonMomDaughter->Fill(trueMom.Mag());

               // std::cout << "  Parent index: " << parentIndex
                 //     << ", parent PDG: " << partPdg[parentIndex]
                   //   << std::endl;
        }
      }
    }
    else if (!isDaughter){
        partMuonMom->Fill(trueMom.Mag());
    }
    }

		// Loop over associations to find matching ReconstructedChargedParticle
		for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
		  {
		    if(simuAssoc[j] == i and weightAssoc[j] > 0.01) // Find association index matching the index of the thrown particle we are looking at
		      { //Loop over the whole simu assoc to see if any of the entries is equal to the MC partilce
			TVector3 recMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); // recoAssoc[j] is the index of the matched ReconstructedChargedParticle

			// Check the distance between the thrown and reconstructed particle
			float deltaEta = trueEta - recMom.PseudoRapidity();
			float deltaPhi = TVector2::Phi_mpi_pi(truePhi - recMom.Phi());
			float deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
			float deltaMom = ((trueMom.Mag()) - (recMom.Mag()));

			matchedPartTrackDeltaEta->Fill(deltaEta);
			matchedPartTrackDeltaPhi->Fill(deltaPhi);
			matchedPartTrackDeltaR->Fill(deltaR);
			matchedPartTrackDeltaMom->Fill(deltaMom);

			matchedPartEta->Fill(trueEta); // Plot the thrown eta if a matched ReconstructedChargedParticle was found
			matchedPartPhi->Fill(truePhi);
			matchedPartMom->Fill(trueMom.Mag());
      matchedPartMomT->Fill(truePT);

			matchedPartPEta->Fill(trueMom.Mag(), trueEta);
			matchedPartPhiEta->Fill(truePhi, trueEta);

      if (pdg == 11){
      matchedPartElecAll->Fill(trueMom.Mag());
      matchedEtaElectron->Fill(trueMom.PseudoRapidity());
      if (isDaughter){
        matchedPartElecMomDaughter->Fill(trueMom.Mag());
        
        continue;
      }
      else if(!isDaughter){
        matchedPartElecMom->Fill(trueMom.Mag());
      }
    
    }

    if (pdg == 321){ //Positive or Neg Kaons
      matchedPartKaonAll->Fill(trueMom.Mag());
      matchedEtaKaon->Fill(trueMom.PseudoRapidity());

      if (isDaughter){
        matchedPartKaonMomDaughter->Fill(trueMom.Mag());
        continue;
      }
      else if(!isDaughter){
        matchedPartKaonMom->Fill(trueMom.Mag());
      }

    }

    if (pdg == 211){ //Positive or Neg Pions
      matchedPartPionAll->Fill(trueMom.Mag());
      matchedEtaPion->Fill(trueMom.PseudoRapidity());
      if (isDaughter){
        matchedPartPionMomDaughter->Fill(trueMom.Mag());
        continue;
      }
      else if(!isDaughter){
        matchedPartPionMom->Fill(trueMom.Mag());
      }
    }

    if (pdg == 13){
      matchedPartMuonAll->Fill(trueMom.Mag());
      matchedEtaMuon->Fill(trueMom.PseudoRapidity());
      if (isDaughter){
        matchedPartMuonMomDaughter->Fill(trueMom.Mag());
        continue;
      }
      else if(!isDaughter){
        matchedPartMuonMom->Fill(trueMom.Mag());
      }
    }
	
		      }
		  }// End loop over associations


	      } // End PDG check
	  } // End stable particles condition
      } // End loop over thrown particles
    // Loop over all charged particles and fill some histograms of kinematics quantities
    for(unsigned int k=0; k<trackMomX.GetSize(); k++){ // Loop over all charged particles, thrown or not
      
      TVector3 CPartMom(trackMomX[k], trackMomY[k], trackMomZ[k]);

      float CPartEta = CPartMom.PseudoRapidity();
      float CPartPhi = CPartMom.Phi();

      ChargedEta->Fill(CPartEta);
      ChargedPhi->Fill(CPartPhi);
      ChargedP->Fill(CPartMom.Mag());
      
    } // End loop over all charged particles


   // cout << "Number of Particle for Event is " << numberOfParticleForEvent << endl;
  
  
  }// End loop over events

  cout << "Number of charged Particles Outside the #eta range of the detector " << numberParticleOutsideEtaRange << endl;
  cout << "As a fraction of total charged particles: " << static_cast<double>(numberParticleOutsideEtaRange) / totalNumberOfChargedParticles << endl;
  // Take the ratio of the histograms above to get our efficiency plots
  TrackEff_Eta->Divide(matchedPartEta, partEta, 1, 1, "b");
  TrackEff_Mom->Divide(matchedPartMom, partMom, 1, 1, "b");
  TrackEff_MomT->Divide(matchedPartMomT, partMomT, 1, 1, "b");
  TrackEff_Phi->Divide(matchedPartPhi, partPhi, 1, 1, "b");
  TrackEff_PEta->Divide(matchedPartPEta, partPEta, 1, 1, "b");
  TrackEff_PhiEta->Divide(matchedPartPhiEta, partPhiEta, 1, 1, "b");

  TrackEff_Elec_All->Divide(matchedPartElecAll, partElecMomAll, 1, 1,"b");
  TrackEff_Kaon_All->Divide(matchedPartKaonAll, partKaonMomAll, 1, 1,"b");
  TrackEff_Pion_All->Divide(matchedPartPionAll, partPionMomAll, 1, 1,"b");
  TrackEff_Muon_All->Divide(matchedPartMuonAll, partMuonMomAll, 1, 1, "b");

  TrackEff_Elec_daughter->Divide(matchedPartElecMomDaughter, partElecMomDaughter, 1, 1, "b");
  TrackEff_Kaon_daughter->Divide(matchedPartKaonMomDaughter, partKaonMomDaughter, 1, 1, "b");
  TrackEff_Pion_daughter->Divide(matchedPartPionMomDaughter, partKaonMomDaughter, 1, 1, "b");
  TrackEff_Muon_daughter->Divide(matchedPartMuonMomDaughter, partMuonMomDaughter, 1, 1, "b");

  TrackEff_Elec->Divide(matchedPartElecMom, partElecMom, 1, 1, "b");
  TrackEff_Kaon->Divide(matchedPartKaonMom, partKaonMom, 1, 1, "b");
  TrackEff_Pion->Divide(matchedPartPionMom, partKaonMom, 1, 1, "b");
  TrackEff_Muon->Divide(matchedPartMuonMom, partMuonMom, 1, 1, "b");


  TrackEff_Elec_eta->Divide(matchedEtaElectron, partEtaElectron, 1, 1, "b");
  TrackEff_Kaon_eta->Divide(matchedEtaKaon, partEtaKaon, 1, 1, "b");
  TrackEff_Pion_eta->Divide(matchedEtaPion, partEtaPion, 1, 1, "b");
  TrackEff_Muon_eta->Divide(matchedEtaMuon, partEtaMuon, 1, 1, "b");


  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
}
