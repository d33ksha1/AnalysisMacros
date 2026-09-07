TString sim_path = "/data/cryoedm/users/pant/hybrid_svtob/epic/FinalBenchmarksNoCopper/ElectronProtonCollision.edm4eic.root";

void AcceptanceRatio(TString infile=sim_path){

    TChain *mychain = new TChain("events");
    mychain->Add(infile);

    TTreeReader tree_reader(mychain);

    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<int> recoAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_rec.index");
    TTreeReaderArray<int> simuAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_sim.index");


    TH2D* hAssociatedHitsVsMomentum = new TH2D("hAssociatedHitsVsMomentum","Associated RecHits vs Momentum;""MC particle momentum [GeV/c];"
    "Number of associated RecHits",
    100, 0.0, 20.0,       // momentum bins
    20, -0.5, 19.5        // hit-count bins
);
    

    int numberOfAssociated = 0;
    int numberOfSelected = 0;

    while (tree_reader.Next()) {
        for (unsigned int i =0; i<partGenStat.GetSize(); i++){ //loops over al MC particles generated
            if (partGenStat[i] == 1){ //only stable particles   
                int pdg = TMath::Abs(partPdg[i]);
                if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212){
                    ++numberOfSelected;
                    for(unsigned int j=0; j<simuAssoc.GetSize(); j++){
                        if (simuAssoc[j] == i){
                            ++numberOfAssociated;
                            break;
                        }
                    }
                }
            }
            }
    }
    double fraction = 0.0;
    fraction = static_cast<double>(numberOfAssociated) /
                static_cast<double>(numberOfSelected);
    cout << "Accpetance Ratio of the Generated, Stable particles " << fraction << endl;
}