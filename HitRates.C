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
#include <string.h>

#include <vector>

std::string input_file = "/data/cryoedm/users/pant/hybrid_svtob/epic/FinalBenchmarksNoCopper/SmallEPCollisions.edm4eic.root";
//std::string input_file =  "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/MuonsFullCraterlakeCutGDML.edm4eic.root";

//std::string input_file = ""
//std::string input_file = "root://hpceph-xrootd.twgrid.org:1094//cephfs/epic//RECO/26.07.1/epic_craterlake/Bkg_Exact1S_2us/GoldCt/10um/DIS/NC/10x275/minQ2=1000/pythia8NCDIS_10x275_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_1.0532.eicrecon.edm4eic.root";
//SimHit collections
std::vector<std::string> sim_coll_names{
	"SiBarrelHits"
};

//TrackerHit RecHit (digitized) collections
std::vector<std::string> rec_coll_names{
	"SiBarrelTrackerRecHits"
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

vector<double> get_time_interval(vector<double> hits){
	std::sort(hits.begin(), hits.end());
    std::vector<double> TimeInterval;
    for (std::size_t i = 1; i < hits.size(); ++i) {
        TimeInterval.push_back(hits[i] - hits[i - 1]);
     }

	 return TimeInterval;
}

double get_mean(vector<double> array){
	double mean = 0;
	double sum = std::accumulate(array.begin(), array.end(), 0);
	mean = sum/array.size();

	return mean;
}
//-----------------
// Declaration of "index_map", to be defined infra
std::unordered_map<std::string, int> index_map;
// Main function
void HitRates(){
    //Defining Histograms
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



    for(int icoll = 0; icoll < (int)sim_coll_names.size(); icoll++){ //Go over each simHit collection
		//vector is destroyed at the end of a loop 
		std::vector<double> Hits;  //timing of each hit 
        std::vector<double> HitsForEvents; //number of hits per event
		std::vector<double> MeanTimeForEvents;

		auto coll_name_rec = rec_coll_names.at(icoll);
		
		for (unsigned int ievent =0; ievent< nevents; ievent++){ //Looping over all events

            auto f = podio::Frame(r.readNextEntry(podio::Category::Event));
            const auto& mc_particles = f.get<edm4hep::MCParticleCollection>("MCParticles");
            const auto& reconstructed_charged_particles = f.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedParticles");
            const auto& recon_associations_charged = f.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedChargedParticleAssociations");
            const auto& reconstructed_particles = f.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
            const auto& reconstructed_particles_assoc = f.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");       
            //RawHit / SimHit association collection
            const auto& rawHitAssociations =f.get<edm4eic::MCRecoTrackerHitAssociationCollection>("CentralTrackingRawHitAssociations");
            const auto& measurements = f.get<edm4eic::Measurement2DCollection>("CentralTrackerMeasurements");
            const auto& track_coll = f.get<edm4eic::TrackCollection>("CentralCKFTracks");

			std::vector<double> MeanTime;
            
            //Only care about RecHit
            const auto& recoHits = f.get<edm4eic::TrackerHitCollection>(coll_name_rec);  //rechits in that collection for that event 
            HitsForEvents.push_back(recoHits.size()); 
            
			for (const auto hit : recoHits){
				cout << "ID of the sensor" << hit.getCellID()<<endl;
                auto timeHit = hit.getTime();
				Hits.push_back(timeHit);
				MeanTime.push_back(timeHit);
            }

			double mean_timeinterval_for_event = get_mean(get_time_interval(MeanTime));
			MeanTimeForEvents.push_back(mean_timeinterval_for_event);

            } //Loop over all events

		
		
		std::vector<double> TimeIntervalHits = get_time_interval(Hits);
        double mean_time_between_hits = get_mean(TimeIntervalHits);


		double mean_hits_per_event = get_mean(HitsForEvents);
		double meantime_per_event = get_mean(MeanTimeForEvents);

		cout << "Mean Number of Hits per event for " << coll_name_rec << "is " << mean_hits_per_event << endl;
        cout << "Mean Time between each Hits for" << coll_name_rec << "is " << mean_time_between_hits << "ns" << endl;
		cout << "Mean time in between hits of each event for " << coll_name_rec << "is " << meantime_per_event << endl;

 
}  
}
   


