#include <podio/ObjectID.h>
#include <podio/ROOTReader.h>
#include <podio/Frame.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4eic/MCRecoParticleAssociationCollection.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <THStack.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraph.h>
#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>


//std::string input_file = "/data/cryoedm/users/pant/hybrid_svtob/epic/FinalBenchmarksNoCopper/SmallEPCollisions.edm4eic.root";
//std::string input_file =  "/data/cryoedm/users/pant/hybrid_svtob/epic/FinalBenchmarksNoCopper/simMuFullTrackingRangeMomentumActual.edm4eic.root";
//std::string input_file = "/data/cryoedm/users/pant/copper_eic/epic/FinalBenchmarks/simMuFullCraterlakeEpicMain.edm4eic.root";
std::string input_file = "/data/cryoedm/users/pant/testing_hybrid/epic/FinalSimBenchmark/test.edm4eic.root";
//std::string input_file = "/data/cryoedm/users/pant/hybrid_svtob/epic/FinalBenchmarksNoCopper/EPCollisionsFullCraterlakeCorrectedActual.edm4eic.root
//std::string input_file = "//std::string input_file =  "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/EPCollisions.edm4eic.root";

//SimHit collections
std::vector<std::string> sim_coll_names{
	"VertexBarrelHits", "SiBarrelHits", "TrackerEndcapHits",
	"MPGDBarrelHits","BackwardMPGDEndcapHits","ForwardMPGDEndcapHits","OuterMPGDBarrelHits"
};

//TrackerHit RecHit (digitized) collections
std::vector<std::string> rec_coll_names{
	"SiBarrelVertexRecHits","SiBarrelTrackerRecHits", "SiEndcapTrackerRecHits",
	"MPGDBarrelRecHits","BackwardMPGDEndcapRecHits","ForwardMPGDEndcapRecHits","OuterMPGDBarrelRecHits"
};

//Collection IDs for RecHit collections.
std::vector<unsigned int> rec_coll_ids;

//------------------
//Template to access 'sign' of radius in (x,y) plane
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
//-----------------
//Function to set Collection IDs
vector<unsigned int> get_coll_ids(vector<string> coll_names){

	vector<unsigned int> coll_ids;

	podio::ROOTReader reader;
	reader.openFile(input_file.c_str());

	unsigned nEntries = reader.getEntries(podio::Category::Event);

	//Find the associated collection IDs
	for (unsigned i = 0; i < nEntries; ++i) {

		if(i > 0) break;

    		auto frameData = reader.readEntry(podio::Category::Event, i);
    		auto frame     = podio::Frame(std::move(frameData));

    		auto collectionNames = frame.getAvailableCollections();

		// Find associated collection ID for each track collection
		for(int iname = 0; iname < (int)coll_names.size(); iname++) {

			for (const auto& name : collectionNames) {
				const auto* coll = frame.get(name);
      				if (coll) {
        				unsigned int collectionID = coll->getID();
					if( name == coll_names[iname] )
					coll_ids.push_back(collectionID);
				}
			}

		} //Loop over track collections
	} //Loop over events

	return coll_ids;
}


//-----------------
// Declaration of "index_map", to be defined infra
std::unordered_map<std::string, int> index_map;

//-----------------

// Main function
void MissedMCParticles(){
    //Defining Histograms

    unsigned int numberTotalSimHits = 0;
    unsigned int numberSimHits = 0;
    unsigned int numberRawHits = 0;
    unsigned int numberRecoHits = 0;
    unsigned int numberMeasurements = 0;
    int numberMCChargedParticles = 0;
    int numberReconstructedChargedParticles = 0;
    int numberUsedByTracks = 0;
    int numberTotalStableParticles = 0;
    int numberGhostParticles = 0;

    TFile *ofile = TFile::Open("testing_hybrid/test.root","RECREATE");

    TH1D *missedElectron = new TH1D("missedElectrons","MC Particles Not Reconstruced; P (GeV/c)", 50, 0, 25);
    TH1D *missedKaon = new TH1D("missedKaons","MC Particles Not Reconstruced; P (GeV/c)", 50, 0, 25);
    TH1D *missedPion = new TH1D("missedPions","MC Particles Not Reconstruced; P (GeV/c)", 50, 0, 25);

    TH1D *MCElectron = new TH1D("MCElectron", "Momentum of MC Electrons; P(GeV/c)", 50, 0, 25);
    TH1D *MCKaon = new TH1D("MCKaon", "Momentum of MC Kaons; P(GeV/c)", 50, 0, 25);
    TH1D *MCPion = new TH1D("MCPion", "Momentum of MC Pions; P(GeV/c)", 50, 0, 25);


    TH2D *simHitPositions = new TH2D("simHitPositions", "SimHit positions;x [mm];y [mm]",200, -500, 500, 200, -500, 500);
    
    TH2D *missedRawHits = new TH2D("miseedRawHits","RawHits missed on the SiBarrel;x [mm];y [mm]",200, -500, 500, 200, -500, 500);
    //The SimHit was not translated successfully to RecoHit
    TH2D *missedRecoHits = new TH2D("missedRecoHits","RecoHits missed on the SiBarrel;x [mm];y [mm]",200, -500, 500, 200, -500, 500);
    //The SimHit was translated successfully to RecoHit but not to a Measurment
    TH2D *RecoHitsNoMeasurement = new TH2D("RecoHitsNoMeasurement","RecoHits missed on the SiBarrel;x [mm];y [mm]",200, -500, 500, 200, -500, 500);
     //The SimHit was translated successfully to RecoHit and Measurment but not used in the track
    TH2D *MeasuredHitsNotUsedForTracks = new TH2D("MeasuredHitsNotUsedForTracks","RecoHitsNotUsedForTracks; x [mm];y [mm]",200, -500, 500, 200, -500, 500);
    //The simHit was translated succesffuly to RecHit and then to a Measurement and then uses for tracks. 
    TH2D *MeasurementUsedForTracks = new TH2D("MeasurementUsedForTracks","MeasurementUsedForTracks; x [mm];y [mm]",200, -500, 500, 200, -500, 500);


    TH2D *simHitPositionsYZ = new TH2D("simHitPositionsYZ", "SimHit positions;y [mm];z [mm]",200, -500, 500, 200, -500, 500);
    TH2D *missedRawHitsYZ = new TH2D("miseedRawHitsYZ","RawHits missed on the SiBarrel;y [mm];z [mm]",200, -500, 500, 200, -500, 500);
    //The SimHit was not translated successfully to RecoHit
    TH2D *missedRecoHitsYZ = new TH2D("missedRecoHitsYZ","RecoHits missed on the SiBarrel;y [mm];z [mm]",200, -500, 500, 200, -500, 500);
    //The SimHit was translated successfully to RecoHit but not to a Measurment
    TH2D *RecoHitsNoMeasurementYZ = new TH2D("RecoHitsNoMeasuremenYZ","RecoHits missed on the SiBarrel;y [mm];z [mm]",200, -500, 500, 200, -500, 500);
     //The SimHit was translated successfully to RecoHit and Measurment but not used in the track
    TH2D *MeasuredHitsNotUsedForTracksYZ = new TH2D("MeasuredHitsNotUsedForTracksYZ","RecoHitsNotUsedForTracks; y [mm];z [mm]",200, -500, 500, 200, -500, 500);
    //The simHit was translated succesffuly to RecHit and then to a Measurement and then uses for tracks. 
    TH2D *MeasurementUsedForTracksYZ = new TH2D("MeasurementUsedForTracksYZ","MeasurementUsedForTracks; y [mm];z [mm]",200, -500, 500, 200, -500, 500);



    auto reconstructedXY = new TGraph();
    auto reconstructedZR = new TGraph();

    auto missedXY = new TGraph();
    auto missedZR = new TGraph();

    reconstructedXY->SetName("reconstructedHitsXY");
    reconstructedXY->SetTitle(
        "Particle hits;x [mm];y [mm]"
    );
    reconstructedXY->SetMarkerColor(kBlue);
    reconstructedXY->SetMarkerStyle(20);
    reconstructedXY->SetMarkerSize(0.7);
    reconstructedXY->SetLineWidth(0.0);

    reconstructedZR->SetName("reconstructedHitsZR");
    reconstructedZR->SetTitle(
        "Particle hits;z [mm];r [mm]"
    );
    reconstructedZR->SetMarkerColor(kBlue);
    reconstructedZR->SetMarkerStyle(20);
    reconstructedZR->SetMarkerSize(0.7);
    reconstructedZR->SetLineWidth(0.0);

    missedXY->SetName("missedHitsXY");
    missedXY->SetMarkerColor(kRed);
    missedXY->SetMarkerStyle(20);
    missedXY->SetMarkerSize(0.7);
    missedXY->SetLineWidth(0.0);

    missedZR->SetName("missedHitsZR");
    missedZR->SetMarkerColor(kRed);
    missedZR->SetMarkerStyle(20);
    missedZR->SetMarkerSize(0.7);
    missedZR->SetLineWidth(0.0);


    auto ghostXY = new TGraph();
    auto ghostZR = new TGraph();


    ghostXY->SetName("GhostHitsXY");
    ghostXY->SetTitle(
        "Particle hits used to reconstruct a Ghost;x [mm];y [mm]"
    );
    ghostXY->SetMarkerColor(kRed);
    ghostXY->SetMarkerStyle(20);
    ghostXY->SetMarkerSize(0.1);
    ghostXY->SetLineWidth(0.0);

    ghostZR->SetName("GhostHitsZR");
    ghostZR->SetTitle(
        "Particle hits used to reconstruct a Ghost;z[mm];r [mm]"
    );
    ghostZR->SetMarkerColor(kRed);
    ghostZR->SetMarkerStyle(20);
    ghostZR->SetMarkerSize(0.1);
    ghostZR->SetLineWidth(0.0);
    


    int pointNumber = 0;
    int pointNumberGhost = 0;

	//Track collection to use: real-seeded or truth-seeded
	std::string trk_coll = "CentralCKFTracks";

    podio::ROOTReader r;
	r.openFile(input_file);

	auto nevents = r.getEntries(podio::Category::Event);
	cout<<"---------------"<<endl;
	cout<<"Total number of events = "<<nevents<<"!"<<endl;
	cout<<"---------------"<<endl; 

	//Set collection Ids
	rec_coll_ids = get_coll_ids(rec_coll_names);

	//Create map between RecHit collection IDs and names
    	std::unordered_map<unsigned int, std::string> RecHitCollMap;
	for (int i = 0; i < (int)rec_coll_ids.size(); i++) {
		RecHitCollMap[rec_coll_ids[i]] = rec_coll_names[i];
    	}

	//Print out Map information
	cout << endl << "Created Map between Collection IDs and Names:"<<endl;
	for (const auto& pair : RecHitCollMap) {
        	cout << "Collection ID: "<< pair.first << " | Collection Name: " << pair.second << endl;
	}

	//Create map between RecHit collection names and vector index...
	//...making room for two indices for the 2DStrip MPGDs.
	size_t nRecColls = rec_coll_names.size();
	for (size_t i = 0; i < nRecColls; i++) {
		index_map[rec_coll_names[i]] = i;
	}


    for (unsigned int ievent =0; ievent< nevents; ievent++){ //Looping over all events

        auto f = podio::Frame(r.readNextEntry(podio::Category::Event));
        const auto& mc_particles = f.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto& reconstructed_charged_particles = f.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedParticles");
        const auto& recon_associations_charged = f.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedChargedParticleAssociations");
        const auto& reconstructed_particles = f.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
        const auto& reconstructed_particles_assoc = f.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
        numberReconstructedChargedParticles += reconstructed_charged_particles.size();
        //RawHit / SimHit association collection
        const auto& rawHitAssociations =f.get<edm4eic::MCRecoTrackerHitAssociationCollection>("CentralTrackingRawHitAssociations");
        const auto& measurements = f.get<edm4eic::Measurement2DCollection>("CentralTrackerMeasurements");
        const auto& track_coll = f.get<edm4eic::TrackCollection>("CentralCKFTracks");

        for(const auto particle : mc_particles){ //Go over all the MC particles
                if (particle.getGeneratorStatus() == 1){
                    numberTotalStableParticles ++;
                }

                if (std::abs(particle.getCharge()) < 1e-6 || particle.getGeneratorStatus() != 1) {
                    continue;
                }   
        
                numberMCChargedParticles ++;

                auto momentum = particle.getMomentum();
                auto pdg = TMath::Abs(particle.getPDG());
                double p = std::sqrt(momentum.x * momentum.x +momentum.y * momentum.y +momentum.z * momentum.z);
                
                if (pdg == 11){ //Electron
                        MCElectron->Fill(p);
                    }
                    if (pdg == 321){ //Kaon
                        MCKaon->Fill(p);
                    }
                    if (pdg == 211){ //Pion
                        MCPion->Fill(p);
                    }
                
                bool recordedParticle = false;
                
                edm4eic::ReconstructedParticle particle_reconstructed = edm4eic::ReconstructedParticle::makeEmpty();
                double bestWeight = 0;
                
                for (const auto recon_association : recon_associations_charged){
                    const auto sim = recon_association.getSim();
                    const auto rec = recon_association.getRec();

                     if (!sim.isAvailable() || !rec.isAvailable()) {
                        continue;
                    }
                    
                    if (recon_association.getSim() == particle && recon_association.getWeight() > 0.8){
                        recordedParticle = true;
                        particle_reconstructed = recon_association.getRec();
                        break;
                    }
                }

                if (recordedParticle){ //This plots the trajectory of particles 
                    
                    for (const auto& track : particle_reconstructed.getTracks()) {
                        if (!track.isAvailable()){
                            continue;
                        }
                        for (const auto& measurement : track.getMeasurements()) {
                            if (!measurement.isAvailable()){
                                continue;
                            }
                            // A measurement can be related to one or more TrackerHits
                            for (const auto& hit : measurement.getHits()) {
                                if (!hit.isAvailable()){
                                    continue;
                                }
                                const auto position = hit.getPosition();

                                const double x = position.x;
                                const double y = position.y;
                                const double z = position.z;
                                const double r = sgn<double>(x) * sqrt(x*x + y*y);

                                reconstructedXY->SetPoint(pointNumber, x, y);
                                reconstructedZR->SetPoint(pointNumber, z, r);

                                ++pointNumber;
                            }
                        }
                    }

                    continue;

                }  
                
                //The MC Particle was not reconstructed 

                if (pdg == 11){ //Electron
                    missedElectron->Fill(p);
                }
                if (pdg == 321){
                    missedKaon->Fill(p);
                }
                if (pdg == 211){
                    missedPion->Fill(p);
                }
            }


        for(int icoll = 0; icoll < (int)sim_coll_names.size(); icoll++){
			//SimHits
			auto coll_name_sim = sim_coll_names.at(icoll);
            auto coll_name_rec = rec_coll_names.at(icoll);
			auto& simHits = f.get<edm4hep::SimTrackerHitCollection>(coll_name_sim);
			auto num_simhits = simHits.size();
            numberTotalSimHits += simHits.size();
            const auto& recoHits = f.get<edm4eic::TrackerHitCollection>(coll_name_rec);
            
            for(const auto particle : mc_particles){ //Go over all the MC particles
                if (particle.getGeneratorStatus() == 1){
                    numberTotalStableParticles ++;
                }

                if (std::abs(particle.getCharge()) < 1e-6 || particle.getGeneratorStatus() != 1) {
                    continue;
                }   
            
                bool recordedParticle = false;
                
                for (const auto recon_association : recon_associations_charged){
                    if (recon_association.getSim() == particle && recon_association.getWeight() > 0.8){
                        recordedParticle = true;
                    }
                }

                if (recordedParticle){ //Dont care about reoconstructed particles 
                    continue;
                }  
                //Particle wasnt reconstructe

                
                for (const auto& simHit : simHits){ //Silicon Barrel Sim Hits for not Reconstructed 
                
                    if (!simHit.getParticle().isAvailable()){
                    continue;}

                    if (simHit.getParticle() != particle){
                            continue;
                        }
                    //The Sim Hit is with the particle
                    const auto position = simHit.getPosition();

                    const double x = position.x;
                    const double y = position.y;
                    const double z = position.z;
                    const double r =  sgn<double>(x) * std::sqrt(x*x + y*y);

                    missedXY->SetPoint(pointNumber, x, y);
                    missedZR->SetPoint(pointNumber, z, r);
                    simHitPositions->Fill(simHit.getPosition().x, simHit.getPosition().y);
                    simHitPositionsYZ->Fill(simHit.getPosition().y, simHit.getPosition().z);

                    numberSimHits ++;
                    bool hasRawHit = false;            
                    bool hasRecoHit = false;
                    bool hasMeasurement = false;
                    bool usedByTrack = false;

                    for (const auto& association : rawHitAssociations){
                        if (association.getSimHit() != simHit){
                            continue;
                        }


                        //The sim Hit has a raw hit found
                        hasRawHit = true;

                        for (const auto& recoHit : recoHits){
                            if (recoHit.getRawHit() != association.getRawHit()){
                                continue;
                            }
                        


                        //The Raw Hit has found a recHit
                        hasRecoHit = true;

                        //Number of matched SimHit - RecoHit paths 
                        ++numberRecoHits;

                        
                        
                        //Find any Measurment2D that has this TrackerHIt
                        for (const auto& measurement : measurements){ 
                            bool thisMeasurementContainsRecoHit = false; 
                            for (const auto& hit : measurement.getHits()){
                                if (hit == recoHit){
                                    thisMeasurementContainsRecoHit = true;
                                    break;
                                }
                            }
                            if (!thisMeasurementContainsRecoHit){
                                continue;
                            }

                            hasMeasurement = true;
                            numberMeasurements ++ ; //Increases the number of measurment Hits

                            for (const auto& track : track_coll){ //Does this measurement appear in any track?
                                //Each event has size of track_coll number of tracks 
                                bool usedByThisTrack = false;
                                for (const auto& trackMeasurement : track.getMeasurements()){ //looping 
                                    if (trackMeasurement == measurement){
                                        usedByThisTrack = true;
                                        break;
                                    } //if loop
                                } //looping over all track measurements

                                if (usedByThisTrack){
                                    usedByTrack = true;
                                    break;
                                }
                            } //looping over all the tracking 

                        } //looping over all the measurement

                    } //Looping over all recoHits
                    
                } //Looping all associations for the simHit to find rawHit
                    //This is for the SimHit
                    if (!hasRawHit){
                        missedRawHits->Fill(simHit.getPosition().x, simHit.getPosition().y);
                        missedRawHitsYZ->Fill(simHit.getPosition().y, simHit.getPosition().z);
                    }
                    
                    else if (!hasRecoHit) {
                        missedRecoHits -> Fill(simHit.getPosition().x, simHit.getPosition().y);
                        missedRecoHitsYZ->Fill(simHit.getPosition().y, simHit.getPosition().z);
                    }   
                    else if (!hasMeasurement) {
                        RecoHitsNoMeasurement-> Fill(simHit.getPosition().x, simHit.getPosition().y);
                        RecoHitsNoMeasurementYZ->Fill(simHit.getPosition().y, simHit.getPosition().z);
                    }
                    else if (!usedByTrack) {
                        MeasuredHitsNotUsedForTracks->Fill(simHit.getPosition().x, simHit.getPosition().y);
                        MeasuredHitsNotUsedForTracksYZ->Fill(simHit.getPosition().y, simHit.getPosition().z);
                    }
                    else if (usedByTrack){
                        ++numberUsedByTracks;
                        MeasurementUsedForTracks->Fill(simHit.getPosition().x, simHit.getPosition().y);
                        MeasurementUsedForTracksYZ->Fill(simHit.getPosition().y, simHit.getPosition().z);
                    }
                    
            }//Loop over all SimHits 
        } //Loop over all the MC particles
    }//Loop over all collections

        for(const auto recon_particle : reconstructed_charged_particles){ //Go over all the Recon partilces
            bool matchedParticle = false;
            for (const auto recon_association: recon_associations_charged){
                if (recon_association.getRec() == recon_particle && recon_association.getWeight() > 0.8){
                    matchedParticle = true;
                }
            }

            if (!matchedParticle){
                numberGhostParticles ++;
             //   cout << "The particle " << recon_particle.getPDG() << "is a ghost from Event " << ievent << "with charge " << recon_particle.getCharge() << endl;

                
                    for (const auto& track : recon_particle.getTracks()) {
                        if (!track.isAvailable()){
                            continue;
                        }
                        for (const auto& measurement : track.getMeasurements()) {
                            if (!measurement.isAvailable()){
                                continue;
                            }
                            // A measurement can be related to one or more TrackerHits
                            for (const auto& hit : measurement.getHits()) {
                                if (!hit.isAvailable()){
                                    continue;
                                }
                                const auto position = hit.getPosition();

                                const double x = position.x;
                                const double y = position.y;
                                const double z = position.z;
                                const double r = sgn<double>(x) * sqrt(x*x + y*y);

                                ghostXY->SetPoint(pointNumber, x, y);
                                ghostZR->SetPoint(pointNumber, z, r);

                                ++pointNumberGhost;
                            }
                        }
                    }



            }
        }

}  //Loop all over eents
    cout << "Total Number of Sim hits :" << numberTotalSimHits << endl;
    cout << "BELOW IS STATISTICS FOR MC PARTICLES FAILED TO BE RECONSTRUCTED" << endl;
    cout << "Number of Sim Hits in Si Barrel for MC particles failed to be reconstructed " << numberSimHits << endl;
    cout << "Number of Reco Hits in Si Barrels for MC particles failed to be reconstructed " << numberRecoHits << endl;
    cout << "Number of Measurement Hits in Si Barrels " << numberMeasurements << endl;
    cout << "Number of Measurement Hits used in SiBarrels for Tracks " << numberUsedByTracks <<endl;
    cout << "Number of Stable MC Particles (neutral and charged)" <<numberTotalStableParticles <<endl;
    cout << "Number of MC Charged Particles : " << numberMCChargedParticles << endl;
    cout << "Number of Recon Charged Particles: " << numberReconstructedChargedParticles << endl;
    cout << "Number of Ghost Particles (charge) " << numberGhostParticles << endl;
    
    reconstructedXY->Write();
    reconstructedZR->Write();
    missedXY->Write();
    missedZR->Write();
    ghostXY->Write();
    ghostZR->Write();
    ofile->Write();
    ofile->Close();


} 
