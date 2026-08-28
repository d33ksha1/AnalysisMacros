#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>

#include <array>
#include <cmath>
#include <vector>
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

//std::string input_file = "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/MuonsFullCraterlake.edm4eic.root";
std::string input_file = "root://hpceph-xrootd.twgrid.org:1094//cephfs/epic//RECO/26.07.0/epic_craterlake/DIS/NC/10x100/minQ2=10/pythia8NCDIS_10x100_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.07*.eicrecon.edm4eic.root";
//std::string input_file = "/data/cryoedm/users/pant/copper_eic_hybrid/epic/FinalSimBenchmark/MuonsFullCraterlakeCutGDMLMuon25GeV.edm4eic.root";

void MoreTrackingPlots(){

        // Set output file for the histograms
    //TFile *ofile = TFile::Open("DeadCellCopper/TrackEfficiencyEPCollisions.root","RECREATE");

    // Analysis code will go here
    // Set up input file chain
    podio::ROOTReader r;
	r.openFile(input_file);
    auto nevents = r.getEntries(podio::Category::Event);

    // Edges defining the eta regions
    const std::vector<double> etaEdges = {-3.5, -1.0, 1.0, 3.5};

    const int nEtaBins = etaEdges.size() - 1;

    // One efficiency object for each eta region
    std::vector<TEfficiency*> efficiencies;

    for (int i = 0; i < nEtaBins; ++i) {
        auto* efficiency = new TEfficiency(
            Form("eff_eta_%d", i),
            Form("%.1f < #eta < %.1f;p_{T} [GeV];Tracking efficiency",
                 etaEdges[i], etaEdges[i + 1]),
            20, 0.0, 3.0   // 20 pT bins from 0 to 10 GeV
        );

        efficiencies.push_back(efficiency);
    }

    for (unsigned int ievent =0; ievent< nevents; ievent++){ // Loop over events
    auto f = podio::Frame(r.readNextEntry(podio::Category::Event));
    const auto& mc_particles = f.get<edm4hep::MCParticleCollection>("MCParticles");
    const auto& recon_associations_charged = f.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedChargedParticleAssociations");


    // Replace this section with your actual event/particle loop
    for (const auto& particle : mc_particles) {
        // Select particles that should enter the denominator
        if (particle.getCharge() == 0.0)
            continue;

        if (particle.getGeneratorStatus() != 1)
            continue;

        const auto momentum = particle.getMomentum();

        const double px = momentum.x;
        const double py = momentum.y;
        const double pz = momentum.z;

        const double pt = std::sqrt(px * px + py * py);
        const double p  = std::sqrt(px * px + py * py + pz * pz);

        // Protect against invalid eta values
        if (p <= std::abs(pz))
            continue;

        TVector3 trueMom(px,py,pz);
        float eta = trueMom.PseudoRapidity();

        bool reconstructed = false;

        for (const auto& association : recon_associations_charged) {
            if (association.getSim().getObjectID() == particle.getObjectID()) {
                reconstructed = true;
                break;
            }
        }

        // Find the particle's eta interval
        for (int i = 0; i < nEtaBins; ++i) {
            const bool insideEtaBin = eta >= etaEdges[i] && eta < etaEdges[i + 1];

            if (insideEtaBin) {
                // First argument: passed reconstruction selection?
                // Second argument: pT value
                efficiencies[i]->Fill(reconstructed, pt);
                break;
            }
        }
    }
}

    auto* canvas = new TCanvas(
        "efficiencyCanvas",
        "Efficiency versus transverse momentum",
        900,
        700
    );

    auto* legend = new TLegend(0.70, 0.75, 0.88, 0.88);

    const std::array<int, 6> colours = {
        kRed + 1,
        kBlue + 1,
        kGreen + 2,
        kMagenta + 1,
        kOrange + 7,
        kCyan + 2
    };

    for (int i = 0; i < nEtaBins; ++i) {
        efficiencies[i]->SetLineColor(colours[i % colours.size()]);
        efficiencies[i]->SetMarkerColor(colours[i % colours.size()]);
        efficiencies[i]->SetMarkerStyle(20 + i);
        efficiencies[i]->SetLineWidth(2);

        if (i == 0)
            efficiencies[i]->Draw("AP");
        else
            efficiencies[i]->Draw("P SAME");

        legend->AddEntry(
            efficiencies[i],
            Form("%.1f < #eta < %.1f", etaEdges[i], etaEdges[i + 1]),
            "lp"
        );
    }

    legend->Draw();

    canvas->Update();

    // Set common axis ranges after the first efficiency has been painted
    efficiencies[0]->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 1.05);
    efficiencies[0]->GetPaintedGraph()->GetXaxis()->SetTitle("p_{T} [GeV]");
    efficiencies[0]->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");

    canvas->Modified();
    canvas->Update();

    canvas->SaveAs("efficiency_vs_pt_eta_binsDIS.pdf");

    
}