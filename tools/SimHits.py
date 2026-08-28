import argparse as ap
import ROOT
import numpy as np

from podio.reading import get_reader


#If running a single muon gun, then each event will just have one particle
#A more general case will have multiple particles being produced in a single event - proton electron collision
in_file_def = "/data/cryoedm/users/pant/hybrid_svtob/epic/FinalBenchmarksNoCopper/SmallEPCollisions.edm4eic.root";

reader = get_reader(in_file_def)

numberOfEvents = 0
numEventsNoSimHits = 0


for frame in reader.get("events"):
    particles    = frame.get("MCParticles");

    numberOfEvents +=1

    tracks       = frame.get("CentralCKFTracks");
    clusters_neg = frame.get("EcalEndcapNClusters");
    clusters_cen = frame.get("EcalBarrelClusters");
    clusters_pos = frame.get("EcalEndcapPClusters");
    associations = frame.get("CentralCKFTrackAssociations"); # can have getSim which associations.getSim() == daughter must be true to be valid
    rawhitassoc = frame.get("CentralTrackingRawHitAssociations");
    sim_collection_name = frame.get("SiBarrelHits");

    reco_associations = frame.get("ReconstructedChargedParticleAssociations");

    matched_simhits = 0;
    numberstableParticlesForEvent = 0;

    num_sim_hits = sim_collection_name.size()
    electronSimHitsMissed = 0;
    kaonSimHitsMissed = 0;
    pionSimHitsMissed = 0;
    electronSimHits = 0;
    kaonSimHits = 0;
    pionSimHits= 0;

    if num_sim_hits == 0:
        numEventsNoSimHits +=1
        continue

    for ihit in range(num_sim_hits): #Going through all the Sim Hits
        hit = sim_collection_name.at(ihit)
        object_id = hit.getObjectID()
        hit_collid = object_id.collectionID
        hit_index = object_id.index
        hit_cellid = hit.getCellID()
        quality = hit.getQuality()

        #Getting the MC Particle Asssocited with that sim hit

        hit_mcpart = hit.getParticle();
        mc_index = hit_mcpart.getObjectID().index;
        mc_pdg = abs(hit_mcpart.getPDG());
        mc_status = hit_mcpart.getGeneratorStatus();

        if (mc_status ==1 and quality == 0):
            if mc_pdg == 11:
                electronSimHits +=1
            if mc_pdg == 211:
                pionSimHits +=1
            if mc_pdg == 321:
                kaonSimHits +=1
            #If the sim hit is not matched then we want to store the information of what the particle is 
            #

   # print(f"Number of total Sim Hits : {num_sim_hits} | Number of matched Sim Hits : {matched_simhits}")



print(f"Total Number of Events: {numberOfEvents}")
print(f"Total Number of Events with No Sim Hits : {numEventsNoSimHits}")
print(f"Total Number of electro, pion, kaon sim hits : {electronSimHits}, {pionSimHits}, {kaonSimHits}")



reader = get_reader(in_file_def)
numberParticles = 0;
numberNeutralParticles = 0;

for frame in reader.get("events"):    
    particles = frame.get("MCParticles");
    recoAssoc = frame.get("ReconstructedChargedParticleAssociations_rec.index");
    simuAssoc = frame.get("ReconstructedChargedParticleAssoications_sim.index");
    weightAssoc = frame.get("ReconstructedChargedParticleAssociations")
    reconstructedParticles = frame.get("ReconstructedChargedParticleAssociations.weight")
    for particle in particles: #loop over all MC particles generated
        numberParticles +=1
        if particle.getCharge() == 0:
            numberNeutralParticles +=1
    numberGhost = 0
    for reco in reconstructedParticles: #Loop over all Reconstructed particles:
        ghost = True
        for i in range (recoAssoc.GetSize()):
            if reco == recoAssoc[i]:
                ghost = False
    
        if ghost == True:
            numberGhost +=1

print(f"Total Number of Particles: {numberParticles}");
print(f"Total Number of Neutral Particle : {numberNeutralParticles}");
print(f"Fraction of Neutral particles:  {numberNeutralParticles/numberParticles}");
print(f"Total Number of Ghost Particles :  {numberGhost}" )