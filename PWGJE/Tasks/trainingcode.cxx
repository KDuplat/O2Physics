
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"


#include "/home/kduplat/Documents/MotherfromSubprocess.h"
#include "/home/kduplat/Documents/Truejetfinder.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "TLorentzVector.h"
#include <TGraph.h>
#include "TCanvas.h"

    /**********Vertexing***********/
    #include <algorithm>

    #include "Common/Core/RecoDecay.h"
    #include "Common/Core/TrackSelectorPID.h"
    #include "Common/Core/trackUtilities.h"
    #include "Common/DataModel/PIDResponse.h"
    #include "Common/DataModel/TrackSelectionTables.h"
    #include "DCAFitter/DCAFitterN.h"
    #include "Framework/AnalysisDataModel.h"
    #include "Framework/AnalysisTask.h"
    #include "Framework/HistogramRegistry.h"
    #include "Framework/runDataProcessing.h"

    #include "PWGHF/Core/SelectorCuts.h"

    using namespace o2;
    using namespace o2::analysis;
    using namespace o2::aod;
    using namespace o2::framework;
    using namespace o2::framework::expressions;
    /***********Vertexing***********/

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h" //doit être ajouté pour le workflow en fin de tâche



struct MyTask
{

    HistogramRegistry histo;


    Configurable<int> motherbin{"motherbin", 20, "Number of bin"};
    Configurable<float> mothermin{"mothermin", 0, "minimum nb of mother in a jet"};
    Configurable<float> mothermax{"mothermax", 20, "maximum nb of mother in a jet"};
    Configurable<int> jetbin{"jetbin", 200, "Number of bin"};
    Configurable<float> jetmin{"jetmin", 0, "minimum track pT"};
    Configurable<float> jetmax{"jetmax", 50, "maximum track pT"};
    Configurable<int> jetptbin{"jetptbin", 200, "Number of bin"};
    Configurable<float> jetptmin{"jetptmin", 0, "minimum jet pT"};
    Configurable<float> jetptmax{"jetptmax", 50, "maximum jet pT"};
    Configurable<int> truejetbin{"truejetbin", 150, "Number of bin"};
    Configurable<float> truejetmin{"truejetmin", 0, "minimum nb of jets"};
    Configurable<float> truejetmax{"truejetmax", 150, "maximum nb of jets"};

    Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
    Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
    Configurable<float> trackEtaMin{"trackEtaMin", -0.8, "minimum track eta"};
    Configurable<float> trackEtaMax{"trackEtaMax", 0.8, "maximum track eta"};
    Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
    Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
    Configurable<float> jetptcut{"jetptcut", 5, ""};

    Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
    Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
    Configurable<double> bz{"bz", 5., "Magnetic field"};
    
    void init(InitContext const&)
    {   


        const AxisSpec specmother{motherbin,mothermin,mothermax};
        const AxisSpec specjet{jetbin,jetmin, jetmax};
        const AxisSpec spectruejet{truejetbin,truejetmin,truejetmax};
        const AxisSpec specjetpt{jetptbin,jetptmin, jetptmax};

        if (doprocessMother){

            histo.add("PDGcode", "First mother particle frequency in jets", HistType::kTH1F, {{8,0, 8}} );
            histo.get<TH1>(HIST("PDGcode"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PDGcode"))->GetXaxis()->SetTitle("Particle");
            auto hPDG = histo.get<TH1>(HIST("PDGcode"));
            hPDG->GetXaxis()->SetBinLabel(1, "u");
            hPDG->GetXaxis()->SetBinLabel(2, "d");
            hPDG->GetXaxis()->SetBinLabel(3, "s");
            hPDG->GetXaxis()->SetBinLabel(4, "c");
            hPDG->GetXaxis()->SetBinLabel(5, "b");
            hPDG->GetXaxis()->SetBinLabel(6, "t");
            hPDG->GetXaxis()->SetBinLabel(7, "g");
            hPDG->GetXaxis()->SetBinLabel(8, "other");

            histo.add("Purity", "Purity in jets", HistType::kTH1F, {{105, 0, 105}});
            histo.get<TH1>(HIST("Purity"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Purity"))->GetXaxis()->SetTitle("Purity (%)");

            histo.add("PMpercentage","Percentage of first mother in jets", HistType::kTH1F, {{200,0,101}});
            histo.get<TH1>(HIST("PMpercentage"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PMpercentage"))->GetXaxis()->SetTitle("%");

            histo.add("DCAhistoXYbbar", "Tracks with b XY DCA in jets", HistType::kTH1F, {{200,-60,60}});
            histo.get<TH1>(HIST("DCAhistoXYbbar"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DCAhistoXYbbar"))->GetXaxis()->SetTitle("DCA XY (cm)");

            histo.add("DCAhistoZbbar", "Tracks with b Z DCA in jets", HistType::kTH1F, {{200,-60,60}}); 
            histo.get<TH1>(HIST("DCAhistoZbbar"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DCAhistoZbbar"))->GetXaxis()->SetTitle("DCA Z (cm)");

            histo.add("DCAhistoXYother", "Tracks without b XY DCA in jets", HistType::kTH1F, {{200,-60,60}});
            histo.get<TH1>(HIST("DCAhistoXYother"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DCAhistoXYother"))->GetXaxis()->SetTitle("DCA XY (cm)");

            histo.add("DCAhistoZother", "Tracks without b Z DCA in jets", HistType::kTH1F, {{200,-60,60}}); 
            histo.get<TH1>(HIST("DCAhistoZother"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DCAhistoZother"))->GetXaxis()->SetTitle("DCA Z (cm)");

            histo.add("DCAhistoXY", "Tracks DCA XY in jets", HistType::kTH1F, {{200,-60,60}});
            histo.get<TH1>(HIST("DCAhistoXY"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DCAhistoXY"))->GetXaxis()->SetTitle("DCA XY (cm)");

            histo.add("DCAhistoZ", "Tracks DCA Z in jets", HistType::kTH1F, {{200,-60,60}});
            histo.get<TH1>(HIST("DCAhistoZ"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DCAhistoZ"))->GetXaxis()->SetTitle("DCA Z (cm)");

            histo.add("2DJetMother","Number of mother in tracks/Jet size", HistType::kTH2F, {specmother, specjet});
            histo.get<TH2>(HIST("2DJetMother"))->GetYaxis()->SetTitle("Jet size");
            histo.get<TH2>(HIST("2DJetMother"))->GetXaxis()->SetTitle("Number of mothers");

            


        }

        if (doprocessTruejet){

            histo.add("TrueJetsize", "Number of tracks in truejets", HistType::kTH1F, {spectruejet});
            histo.get<TH1>(HIST("TrueJetsize"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("TrueJetsize"))->GetXaxis()->SetTitle("Number of tracks");

            histo.add("trueTrackpt", "Transverse momentum of tracks in truejets", HistType::kTH1F, {{500, 0, 20}});
            histo.get<TH1>(HIST("trueTrackpt"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("trueTrackpt"))->GetXaxis()->SetTitle("Transverse momentum (GeV)");
            histo.add("trueTracketa", "Rapidity of tracks in truejets", HistType::kTH1F, {{50, -1, 1}});
            histo.get<TH1>(HIST("trueTracketa"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("trueTracketa"))->GetXaxis()->SetTitle("Rapidity");
            histo.add("trueTrackphi", "Phi of tracks in truejets", HistType::kTH1F, {{50, -1, 7}});
            histo.get<TH1>(HIST("trueTrackphi"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("trueTrackphi"))->GetXaxis()->SetTitle("rad");

            histo.add("2DMomentumTruejet","Truejet Pt/Truejet reconstructed Pt", HistType::kTH2F, {{100, 0, 20}, {100, 0, 20}});
            histo.get<TH2>(HIST("2DMomentumTruejet"))->GetYaxis()->SetTitle("First mother tranverse momentum");
            histo.get<TH2>(HIST("2DMomentumTruejet"))->GetXaxis()->SetTitle("Reconstructed tranverse momentum");

            histo.add("2DetaTruejet","Truejet rapidity/Truejet reconstructed rapidity", HistType::kTH2F, {{100, -2, 2}, {100, -8, 8}});
            histo.get<TH2>(HIST("2DetaTruejet"))->GetYaxis()->SetTitle("First mother rapidity");
            histo.get<TH2>(HIST("2DetaTruejet"))->GetXaxis()->SetTitle("Reconstructed rapidity");

            histo.add("2DphiTruejet","Truejet azimuth angle/Truejet reconstructed azimuth angle", HistType::kTH2F, {{100, -4, 4}, {100, -4, 4}});
            histo.get<TH2>(HIST("2DphiTruejet"))->GetYaxis()->SetTitle("First mother azimuth angle");
            histo.get<TH2>(HIST("2DphiTruejet"))->GetXaxis()->SetTitle("Reconstructed azimuth angle");

            histo.add("Jetsize", "Number of tracks in jets", HistType::kTH1F, {specjet});
            histo.get<TH1>(HIST("Jetsize"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Jetsize"))->GetXaxis()->SetTitle("Number of tracks");

            histo.add("jetTrackpt", "Transverse momentum of tracks in jets", HistType::kTH1F, {{500, 0, 20}});
            histo.get<TH1>(HIST("jetTrackpt"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("jetTrackpt"))->GetXaxis()->SetTitle("Momentum (GeV)");
            histo.add("jetTracketa", "Rapidity of tracks in jets", HistType::kTH1F, {{50, -1, 1}});
            histo.get<TH1>(HIST("jetTracketa"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("jetTracketa"))->GetXaxis()->SetTitle("Rapidity");
            histo.add("jetTrackphi", "Phi of tracks in jets", HistType::kTH1F, {{100, -1, 7}});
            histo.get<TH1>(HIST("jetTrackphi"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("jetTrackphi"))->GetXaxis()->SetTitle("Rad");

            histo.add("JetPt", "Transverse momentum of jets", HistType::kTH1F, {specjetpt});
            histo.get<TH1>(HIST("JetPt"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("JetPt"))->GetXaxis()->SetTitle("Momentum (GeV)");

            histo.add("Truejetmatching","Percentage of correspondance between true jet and jetfinder(Nb of match/jetsize)", HistType::kTH1F, {{200,0,101}});
            histo.get<TH1>(HIST("Truejetmatching"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Truejetmatching"))->GetXaxis()->SetTitle("%");

            histo.add("Nbtracksmatching","Truejet size when matching", HistType::kTH1F, {spectruejet});
            histo.get<TH1>(HIST("Nbtracksmatching"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Nbtracksmatching"))->GetXaxis()->SetTitle("Number of tracks");


            histo.add("PDGcode2", "Finalstate particle frequency in truejets", HistType::kTH1F, {{15,0, 15}} );
            histo.get<TH1>(HIST("PDGcode2"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PDGcode2"))->GetXaxis()->SetTitle("Particle");
            auto hPDG = histo.get<TH1>(HIST("PDGcode2"));
            hPDG->GetXaxis()->SetBinLabel(1, "p");
            hPDG->GetXaxis()->SetBinLabel(2, "pbar");
            hPDG->GetXaxis()->SetBinLabel(3, "K+");
            hPDG->GetXaxis()->SetBinLabel(4, "K-");
            hPDG->GetXaxis()->SetBinLabel(5, "e-");
            hPDG->GetXaxis()->SetBinLabel(6, "e+");
            hPDG->GetXaxis()->SetBinLabel(7, "pi+");
            hPDG->GetXaxis()->SetBinLabel(8, "pi-");
            hPDG->GetXaxis()->SetBinLabel(9, "mu+");
            hPDG->GetXaxis()->SetBinLabel(10, "mu-");
            hPDG->GetXaxis()->SetBinLabel(11, "sigma-");
            hPDG->GetXaxis()->SetBinLabel(12, "sigma+");
            hPDG->GetXaxis()->SetBinLabel(13, "ksi-");
            hPDG->GetXaxis()->SetBinLabel(14, "ksi+");
            hPDG->GetXaxis()->SetBinLabel(15, "Other");

            histo.add("PDGcode6", "First mother particle frequency in truejets", HistType::kTH1F, {{8,0, 8}} );
            histo.get<TH1>(HIST("PDGcode6"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PDGcode6"))->GetXaxis()->SetTitle("Particle");
            hPDG = histo.get<TH1>(HIST("PDGcode6"));
            hPDG->GetXaxis()->SetBinLabel(1, "u");
            hPDG->GetXaxis()->SetBinLabel(2, "d");
            hPDG->GetXaxis()->SetBinLabel(3, "s");
            hPDG->GetXaxis()->SetBinLabel(4, "c");
            hPDG->GetXaxis()->SetBinLabel(5, "b");
            hPDG->GetXaxis()->SetBinLabel(6, "t");
            hPDG->GetXaxis()->SetBinLabel(7, "g");
            hPDG->GetXaxis()->SetBinLabel(8, "other");

            histo.add("PDGcode3", "Finalstate particle frequency in jets", HistType::kTH1F, {{15,0, 15}} );
            histo.get<TH1>(HIST("PDGcode3"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PDGcode3"))->GetXaxis()->SetTitle("Particle");
            hPDG = histo.get<TH1>(HIST("PDGcode3"));
            hPDG->GetXaxis()->SetBinLabel(1, "p");
            hPDG->GetXaxis()->SetBinLabel(2, "pbar");
            hPDG->GetXaxis()->SetBinLabel(3, "K+");
            hPDG->GetXaxis()->SetBinLabel(4, "K-");
            hPDG->GetXaxis()->SetBinLabel(5, "e-");
            hPDG->GetXaxis()->SetBinLabel(6, "e+");
            hPDG->GetXaxis()->SetBinLabel(7, "pi+");
            hPDG->GetXaxis()->SetBinLabel(8, "pi-");
            hPDG->GetXaxis()->SetBinLabel(9, "mu+");
            hPDG->GetXaxis()->SetBinLabel(10, "mu-");
            hPDG->GetXaxis()->SetBinLabel(11, "sigma-");
            hPDG->GetXaxis()->SetBinLabel(12, "sigma+");
            hPDG->GetXaxis()->SetBinLabel(13, "ksi-");
            hPDG->GetXaxis()->SetBinLabel(14, "ksi+");
            hPDG->GetXaxis()->SetBinLabel(15, "Other");

            histo.add("Rdistribution","Truejet radius distribution", HistType::kTH1F, {{1000, 0,15}});
            histo.get<TH1>(HIST("Rdistribution"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Rdistribution"))->GetXaxis()->SetTitle("R");

             histo.add("Rdistributionjet","Jet radius distribution", HistType::kTH1F, {{200, 0,2}});
            histo.get<TH1>(HIST("Rdistributionjet"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Rdistributionjet"))->GetXaxis()->SetTitle("R");

            histo.add("2DMomentumjet","Jet Pt/Jet reconstructed Pt", HistType::kTH2F, {{100, 0, 20}, {100, 0, 20}});
            histo.get<TH2>(HIST("2DMomentumjet"))->GetYaxis()->SetTitle("Jet tranverse momentum");
            histo.get<TH2>(HIST("2DMomentumjet"))->GetXaxis()->SetTitle("Reconstructed tranverse momentum");

            histo.add("2Detajet","Jet Rapidity/Jet reconstructed rapidity", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});
            histo.get<TH2>(HIST("2Detajet"))->GetYaxis()->SetTitle("Jetrapidity");
            histo.get<TH2>(HIST("2Detajet"))->GetXaxis()->SetTitle("Reconstructed rapidity");

            histo.add("2Dphijet","Jet azimuth angle/Jet reconstructed azimuth angle", HistType::kTH2F, {{100, -4, 4}, {100, -4, 4}});
            histo.get<TH2>(HIST("2Dphijet"))->GetYaxis()->SetTitle("Jet azimuth angle");
            histo.get<TH2>(HIST("2Dphijet"))->GetXaxis()->SetTitle("Reconstructed azimuth angle");

        }

        if (doprocessParticleLevel){

            histo.add("Pjetsize", "Size of ParticleLevel jets",HistType::kTH1F, {specjet});
            histo.get<TH1>(HIST("Pjetsize"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Pjetsize"))->GetXaxis()->SetTitle("Number of tracks");

            histo.add("Pjetpt", "Momentum of Pjets", HistType::kTH1F, {specjetpt});
            histo.get<TH1>(HIST("Pjetpt"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Pjetpt"))->GetXaxis()->SetTitle("Momentum (GeV)");

            histo.add("PDGcode4", "Finalstate Particle frequency for ParticleLevel jets", HistType::kTH1F, {{15,0, 15}} );
            histo.get<TH1>(HIST("PDGcode4"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PDGcode4"))->GetXaxis()->SetTitle("Particle");
            auto hPDG = histo.get<TH1>(HIST("PDGcode4"));
            hPDG->GetXaxis()->SetBinLabel(1, "p");
            hPDG->GetXaxis()->SetBinLabel(2, "pbar");
            hPDG->GetXaxis()->SetBinLabel(3, "K+");
            hPDG->GetXaxis()->SetBinLabel(4, "K-");
            hPDG->GetXaxis()->SetBinLabel(5, "e-");
            hPDG->GetXaxis()->SetBinLabel(6, "e+");
            hPDG->GetXaxis()->SetBinLabel(7, "pi+");
            hPDG->GetXaxis()->SetBinLabel(8, "pi-");
            hPDG->GetXaxis()->SetBinLabel(9, "mu+");
            hPDG->GetXaxis()->SetBinLabel(10, "mu-");
            hPDG->GetXaxis()->SetBinLabel(11, "sigma-");
            hPDG->GetXaxis()->SetBinLabel(12, "sigma+");
            hPDG->GetXaxis()->SetBinLabel(13, "ksi-");
            hPDG->GetXaxis()->SetBinLabel(14, "ksi+");
            hPDG->GetXaxis()->SetBinLabel(15, "Other");

            histo.add("PDGcode5", "Primordial mother particle frequency in ParticleLevel jets", HistType::kTH1F, {{8,0, 8}} );
            histo.get<TH1>(HIST("PDGcode5"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PDGcode5"))->GetXaxis()->SetTitle("Particle");
            hPDG = histo.get<TH1>(HIST("PDGcode5"));
            hPDG->GetXaxis()->SetBinLabel(1, "u");
            hPDG->GetXaxis()->SetBinLabel(2, "d");
            hPDG->GetXaxis()->SetBinLabel(3, "s");
            hPDG->GetXaxis()->SetBinLabel(4, "c");
            hPDG->GetXaxis()->SetBinLabel(5, "b");
            hPDG->GetXaxis()->SetBinLabel(6, "t");
            hPDG->GetXaxis()->SetBinLabel(7, "g");
            hPDG->GetXaxis()->SetBinLabel(8, "other");

            histo.add("PtrueTrackpt", "Transverse momentum of tracks in truejets", HistType::kTH1F, {{500, 0, 20}});
            histo.get<TH1>(HIST("PtrueTrackpt"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PtrueTrackpt"))->GetXaxis()->SetTitle("Transverse momentum (GeV)");
            histo.add("PtrueTracketa", "Eta of tracks in truejets", HistType::kTH1F, {{50, -1, 1}});
            histo.get<TH1>(HIST("PtrueTracketa"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PtrueTracketa"))->GetXaxis()->SetTitle("Rapidity");
            histo.add("PtrueTrackphi", "Rapidity of tracks in truejets", HistType::kTH1F, {{50, -1, 7}});
            histo.get<TH1>(HIST("PtrueTrackphi"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("PtrueTrackphi"))->GetXaxis()->SetTitle("Rad");

            histo.add("P2DMomentum","Truejet Pt/Truejet Reconstructed Pt in ParticleLevel", HistType::kTH2F, {{100, 0, 10}, {100, 0, 10}});
            histo.get<TH2>(HIST("P2DMomentum"))->GetYaxis()->SetTitle("Primordial mother tranverse momentum");
            histo.get<TH2>(HIST("P2DMomentum"))->GetXaxis()->SetTitle("Reconstructed tranverse momentum");

        }

        if (doprocessMatching){
            histo.add("Matching", "Number of matching first mother between DetectorLevel jets and ParticleLevel jets",HistType::kTH1F, {{20, 0, 20}});
            histo.get<TH1>(HIST("Matching"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Matching"))->GetXaxis()->SetTitle("Number of matches in a Detector jet");

            /*histo.add("Matching2", "Number of matching of final state between Detector Level and Particle Level",HistType::kTH1F, {{20, 0, 20}});
            histo.get<TH1>(HIST("Matching2"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Matching2"))->GetXaxis()->SetTitle("Number of match in a Detector jet");*/

        }

        if (doprocessSV){
            histo.add("Declenght", "Decaylength for light partons", HistType::kTH1F, {{200,0,20}});
            histo.get<TH1>(HIST("Declenght"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Declenght"))->GetXaxis()->SetTitle("Distance (cm)");

            histo.add("DeclenghtB", "Decaylength for b quarks", HistType::kTH1F, {{200,0,20}});
            histo.get<TH1>(HIST("DeclenghtB"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DeclenghtB"))->GetXaxis()->SetTitle("Distance (cm)");

            histo.add("Declenghtxy", "Decay length xy for light partons", HistType::kTH1F, {{500,0,1}});
            histo.get<TH1>(HIST("Declenghtxy"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Declenghtxy"))->GetXaxis()->SetTitle("Distance (cm)");

            histo.add("DeclenghtBxy", "Decay length xy for b quarks", HistType::kTH1F, {{500,0,1}});
            histo.get<TH1>(HIST("DeclenghtBxy"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DeclenghtBxy"))->GetXaxis()->SetTitle("Distance (cm)");

            histo.add("Declenghtz", "Decay length z for light partons", HistType::kTH1F, {{100,-10,10}});
            histo.get<TH1>(HIST("Declenghtz"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Declenghtz"))->GetXaxis()->SetTitle("Distance (cm)");

            histo.add("DeclenghtBz", "Decay length z for b quarks", HistType::kTH1F, {{100,-10,10}});
            histo.get<TH1>(HIST("DeclenghtBz"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("DeclenghtBz"))->GetXaxis()->SetTitle("Distance (cm)");

            histo.add("Significance", "Significance b secondary vertex", HistType::kTH1F, {{200, 0, 30}});
            histo.get<TH1>(HIST("Significance"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Significance"))->GetXaxis()->SetTitle("significance");

             histo.add("SignificanceB", "Significance light and gluon secondary vertex", HistType::kTH1F, {{200, 0, 30}});
            histo.get<TH1>(HIST("SignificanceB"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("SignificanceB"))->GetXaxis()->SetTitle("significance");
        }
        
        if (doprocessMygraph){
            histo.add("Alltrackspt", "Transverse momentum of all tracks (in tracks cut)", HistType::kTH1F, {{500,0,20}});
            histo.get<TH1>(HIST("Alltrackspt"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Alltrackspt"))->GetXaxis()->SetTitle("Momentum(GeV)");

            histo.add("Nbjetincollision", "Number of jets in collisions", HistType::kTH1F, {{30,0,30}});
            histo.get<TH1>(HIST("Nbjetincollision"))->GetYaxis()->SetTitle("Count");
            histo.get<TH1>(HIST("Nbjetincollision"))->GetXaxis()->SetTitle("Number of jets");

        }
    }

    //Process that look for the DCA, PDGcode frequency and purity in jets by looking at the primordial mother
    int ijet2=0;
    int nbBjet=0;
    void processMother(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet, soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksDCA> const& tracks, aod::McParticles const& mcparticles){

        ijet2++;
        std::ofstream fichier;
        std::ofstream fichier2;
        if (ijet2==1)
        {
            fichier.open("Outputparticle.txt", std::ios::out | std::ios::trunc);
            fichier2.open("NoPM.txt", std::ios::out | std::ios::trunc);
        }else
        {
            fichier.open("Outputparticle.txt", std::ios::out | std::ios::app);
            fichier2.open("NoPM.txt", std::ios::out | std::ios::app);
        }

        //Array that contain all the mothers in a jet
        //First item: Id, Second item: Count
        std::vector<std::vector<int64_t>> mothers;
        int msize=0; 

        //First item:PDGcode, Second item: Count
        std::vector<int> dominantmother(2,0); 

        int totalcount=0;   //Total count of mothers with the statuscode needed
        fichier<<"Jet "<<ijet2<<":"<<std::endl;
        fichier2<<"Jet "<<ijet2<<":"<<std::endl;
        for (auto const& track : jet.tracks_as<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksDCA>>())
        {
            if(!track.has_mcParticle())
            {
                continue;
            }

            auto finalstateParticle=track.mcParticle();
            auto outputIparticle=finalstateParticle;

            histo.fill(HIST("DCAhistoXY"), track.dcaXY());
            histo.fill(HIST("DCAhistoZ"), track.dcaZ());

            if (getMotherfromSubprocess(mcparticles, finalstateParticle, outputIparticle ))
            {
                totalcount++;
                //******Outputparticle.txt*******//
                if (outputIparticle.pdgCode()==5)
                {
                    histo.fill(HIST("DCAhistoXYbbar"), track.dcaXY());
                    histo.fill(HIST("DCAhistoZbbar"), track.dcaZ());
                }else{
                    histo.fill(HIST("DCAhistoZother"), track.dcaZ());
                    histo.fill(HIST("DCAhistoXYother"), track.dcaXY());
                }
                
                if (abs(outputIparticle.getGenStatusCode())==51)
                {
                    auto tempo= outputIparticle;
                    tempo=tempo.mothers_first_as<aod::McParticles>();
                    if(tempo.pdgCode()!=21)
                    {
                        tempo=tempo.mothers_first_as<aod::McParticles>();
                        if(tempo.pdgCode()==21)
                        {
                            fichier<<"|Gluon Jet| ";
                        }
                        if(tempo.pdgCode()!=21)
                        {
                            TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(tempo.pdgCode());
                            if (particle != nullptr){
                                fichier<<"|"<<TDatabasePDG::Instance()->GetParticle(tempo.pdgCode())->GetName()<<" Jet| ";
                            }else{
                                fichier<<"|"<<tempo.pdgCode()<<" Jet| ";
                            }
                        }
                    }
                }

                TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(outputIparticle.pdgCode());
                if (particle != nullptr){
                    fichier<<particle->GetName()<<", StatusCode:"<<outputIparticle.getGenStatusCode();
                }else{
                    fichier<<outputIparticle.pdgCode()<<", StatusCode:"<<outputIparticle.getGenStatusCode();
                }
                //******Outputparticle.txt*******//
                
                int counter=0;
                for(int i=0; i<msize; i++)
                {
                    if (mothers[i][0]==outputIparticle.globalIndex())
                    {
                        mothers[i][2]++;
                        counter++;
                        break;
                    }
                }
                if (counter==0)
                {
                    mothers.push_back({outputIparticle.globalIndex(),outputIparticle.pdgCode(),1});
                    msize=mothers.size();
                }
            }else{

                auto mother=finalstateParticle;
                TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(mother.pdgCode());
                if (particle != nullptr){
                    fichier2<<particle->GetName()<<", StatusCode:"<<mother.getGenStatusCode()<<"-->";
                }else{
                    fichier2<<mother.pdgCode()<<", StatusCode:"<<mother.getGenStatusCode()<<"-->";
                }
                while(mother.has_mothers()){
                    mother=mother.mothers_first_as<aod::McParticles>();
                    particle = TDatabasePDG::Instance()->GetParticle(mother.pdgCode());
                    if (particle != nullptr){
                        fichier2<<particle->GetName()<<", StatusCode:"<<mother.getGenStatusCode()<<"-->";
                    }else{
                        fichier2<<mother.pdgCode()<<", StatusCode:"<<mother.getGenStatusCode()<<"-->";
                    }
                }
                
            }
            fichier<<std::endl;
            fichier2<<std::endl;
        }
        for(int i=0; i<msize; i++)
        {
            if (mothers[i][2]>dominantmother[1])
            {
                dominantmother[1]=mothers[i][2];
                dominantmother[0]=mothers[i][1];
            }
        }

        auto hPDG = histo.get<TH1>(HIST("PDGcode"));   //Fill the PDGcode histogram
        if(abs(dominantmother[0])==1){
            hPDG->Fill("u",1);
        }else if (abs(dominantmother[0])==2){
            hPDG->Fill("d",1);
        }else if (abs(dominantmother[0])==3){
            hPDG->Fill("s",1);
        }else if (abs(dominantmother[0])==4){
            hPDG->Fill("c",1);
        }else if (abs(dominantmother[0])==5){
            nbBjet++;
            hPDG->Fill("b",1);
        }else if (abs(dominantmother[0])==6){
            hPDG->Fill("t",1);
        }else if (abs(dominantmother[0])==21){
            hPDG->Fill("g",1);

        }else if (abs(dominantmother[0])!=0){ // In the case where no mother was found we don't increment
            hPDG->Fill("other",1);
        }

        float purity=dominantmother[1]*100./totalcount; //Calculate the purity of the jet
        if (totalcount>0)
        {
            histo.fill(HIST("Purity"), purity);
        }

        histo.fill(HIST("2DJetMother"), msize, jet.tracks().size());
        histo.get<TH2>(HIST("2DJetMother"))->SetOption("colz");
        float PMpercentage =totalcount*100./jet.tracks().size();
        histo.fill(HIST("PMpercentage"), PMpercentage);
        
        //float Bpurity=nbBjet*100./ijet2;
        //printf("Percentage of B jet:%f\n", Bpurity);


    }
    PROCESS_SWITCH(MyTask, processMother, "Process to find the mother with a Statuscode == 23 || 33", true);


    Filter Trackincut = (aod::track::pt>trackPtMin && aod::track::pt<trackPtMax && aod::track::eta>trackEtaMin && aod::track::eta<trackEtaMax && aod::track::phi>trackPhiMin && aod::track::phi<trackPhiMax);
    using Detectorjet = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;

    void processTruejet(Detectorjet const& jets,soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels>> const& tracks, aod::McParticles const& mcparticles){
        
        //Partition<soa::Join<aod::Tracks, aod::McTrackLabels>> Trackincut = (aod::track::pt>trackPtMin && aod::track::pt<trackPtMax && aod::track::eta>trackEtaMin && aod::track::eta<trackEtaMax && aod::track::phi>trackPhiMin && aod::track::phi<trackPhiMax);

        std::vector<std::vector<int64_t>>truejet=Truejetfinder(mcparticles,tracks);
        int tjetsize=truejet.size();

        for (int i=0; i<tjetsize; i++)
        {
            float jetpt=0;
            int subsize=truejet[i].size();
            
            auto track=tracks.rawIteratorAt(truejet[i][0]-tracks.offset());
            auto finalstate=track.mcParticle();
            auto primordialmother=finalstate;
            getMotherfromSubprocess(mcparticles, finalstate, primordialmother);

            float eta_jet = primordialmother.eta();
            float phi_jet = primordialmother.phi();
            float delta_eta=0;
            float delta_phi=0;
            float deltaR=0;
            float sumE=0;
            float sumpx=0;
            float sumpy=0;
            float sumpz=0;

            for(int j=0; j<subsize; j++)
            {
                auto track=tracks.rawIteratorAt(truejet[i][j]-tracks.offset());
                auto particle=track.mcParticle();

                histo.fill(HIST("trueTrackpt"), track.pt());
                histo.fill(HIST("trueTracketa"), track.eta());
                histo.fill(HIST("trueTrackphi"), track.phi());
                

                auto hPDG = histo.get<TH1>(HIST("PDGcode2"));
                int PDG=particle.pdgCode();
                if(PDG==321)       { hPDG->Fill("K+",1); }
                else if(PDG==-321)       { hPDG->Fill("K-",1); }
                else if(PDG==-211)  { hPDG->Fill("pi+",1); }
                else if(PDG==211)  { hPDG->Fill("pi-",1); }
                else if(PDG==11)   { hPDG->Fill("e-",1); }
                else if(PDG==-11)  { hPDG->Fill("e+",1); }
                else if(PDG==2212) { hPDG->Fill("p",1); }
                else if(PDG==-2212) { hPDG->Fill("pbar",1); }
                else if(PDG==-13)   { hPDG->Fill("mu+",1); }
                else if(PDG==13)   { hPDG->Fill("mu-",1); }
                else if(PDG==3112) { hPDG->Fill("sigma-",1); }
                else if(PDG==3222) { hPDG->Fill("sigma+",1); }
                else if(PDG==3312) { hPDG->Fill("ksi-",1); }
                else if(PDG==-3312) { hPDG->Fill("ksi+",1); }
                else{ hPDG->Fill("Other",1); }

                jetpt+=track.pt();

                float eta_track=track.eta();
                float phi_track=track.phi();

                delta_eta+=eta_track-eta_jet;
                delta_phi+=phi_track-phi_jet;

                sumE+=track.mcParticle().e();
                sumpx+=track.px();
                sumpy+=track.py();
                sumpz+=track.pz();
            }

            

            if (jetpt>jetptcut)
            {
                histo.fill(HIST("TrueJetsize"), truejet[i].size());
                int PMpdgcode=primordialmother.pdgCode();
                auto hPDG = histo.get<TH1>(HIST("PDGcode6"));   //Fill the PDGcode histogram
                if(abs(PMpdgcode)==1){
                    hPDG->Fill("u",1);
                }else if (abs(PMpdgcode)==2){
                    hPDG->Fill("d",1);
                }else if (abs(PMpdgcode)==3){
                    hPDG->Fill("s",1);
                }else if (abs(PMpdgcode)==4){
                    hPDG->Fill("c",1);
                }else if (abs(PMpdgcode)==5){
                    hPDG->Fill("b",1);
                }else if (abs(PMpdgcode)==6){
                    hPDG->Fill("t",1);
                }else if (abs(PMpdgcode)==21){
                    hPDG->Fill("g",1);

                }else if (abs(PMpdgcode)!=0){ // In the case where no mother was found we don't increment
                    hPDG->Fill("other",1);
                }


                deltaR= sqrt((delta_eta*delta_eta+delta_phi*delta_phi)/(subsize*subsize));
                histo.fill(HIST("Rdistribution"), deltaR);


                if (phi_jet> M_PI)
                {
                    phi_jet=phi_jet-2*M_PI;
                }

                TLorentzVector truejetvector(sumpx, sumpy, sumpz,sumE);
                float reconstructedphi=truejetvector.Phi();
                if(reconstructedphi > M_PI)
                {
                    reconstructedphi=reconstructedphi - 2*M_PI;
                }

                histo.fill(HIST("2DetaTruejet"), truejetvector.Eta(), eta_jet);
                histo.get<TH2>(HIST("2DetaTruejet"))->SetOption("colz");
                histo.fill(HIST("2DphiTruejet"), reconstructedphi, phi_jet);
                histo.get<TH2>(HIST("2DphiTruejet"))->SetOption("colz");
                
            }

            histo.fill(HIST("2DMomentumTruejet"), primordialmother.pt(),jetpt);
            histo.get<TH2>(HIST("2DMomentumTruejet"))->SetOption("colz");
        }

        
        for (auto const& jet : jets)
        {   
            histo.fill(HIST("Jetsize"),jet.tracks().size());
            histo.fill(HIST("JetPt"), jet.pt());
            float matchpercentage=0;    //Number of match / Number of 
            int subsize=0;
            int temposubsize;
            std::vector<int64_t> Idfound;  //Table that contain primordial mothers already matched
            int Idsize=0;
            int truecounter=0;      //Number of time a primordial mother is found
            float sumpt=0;

            float sumE=0;
            float sumpx=0;
            float sumpy=0;
            float sumpz=0;

            float eta_jet = jet.eta();
            float phi_jet = jet.phi();
            float delta_eta=0;
            float delta_phi=0;
            float deltaR=0;
            
            auto hPDG = histo.get<TH1>(HIST("PDGcode3"));
            for (auto const& track : jet.tracks_as<soa::Join<aod::Tracks, aod::McTrackLabels>>())
            {

                int PDG=track.mcParticle().pdgCode();
                if(PDG==321)       { hPDG->Fill("K+",1); }
                else if(PDG==-321)       { hPDG->Fill("K-",1); }
                else if(PDG==-211)  { hPDG->Fill("pi+",1); }
                else if(PDG==211)  { hPDG->Fill("pi-",1); }
                else if(PDG==11)   { hPDG->Fill("e-",1); }
                else if(PDG==-11)  { hPDG->Fill("e+",1); }
                else if(PDG==2212) { hPDG->Fill("p",1); }
                else if(PDG==-2212) { hPDG->Fill("pbar",1); }
                else if(PDG==-13)   { hPDG->Fill("mu+",1); }
                else if(PDG==13)   { hPDG->Fill("mu-",1); }
                else if(PDG==3112) { hPDG->Fill("sigma-",1); }
                else if(PDG==3222) { hPDG->Fill("sigma+",1); }
                else if(PDG==3312) { hPDG->Fill("ksi-",1); }
                else if(PDG==-3312) { hPDG->Fill("ksi+",1); }
                else{ hPDG->Fill("Other",1); }

                histo.fill(HIST("jetTrackpt"), track.pt());
                histo.fill(HIST("jetTracketa"), track.eta());
                histo.fill(HIST("jetTrackphi"), track.phi());

                if(!track.has_mcParticle())
                {
                    continue;
                }
                int match=0;
                float tempo;    //temporise the value for matchpercentage
                auto finalstateParticle=track.mcParticle();
                auto outputIparticle=finalstateParticle;
                if (getMotherfromSubprocess(mcparticles, finalstateParticle, outputIparticle) != false)
                {
                    truecounter++;
                    for (int i=0; i<tjetsize; i++)
                    {
                        auto trackinput=tracks.rawIteratorAt(truejet[i][0]-tracks.offset());
                        auto inputparticle=trackinput.mcParticle();
                        auto firstmother=inputparticle;
                        getMotherfromSubprocess(mcparticles, inputparticle, firstmother);

                        if (outputIparticle.globalIndex()==firstmother.globalIndex())
                        {
                            subsize=truejet[i].size();
                            int c=0;
                            for (int k=0; k<Idsize; k++)
                            {
                                if (outputIparticle.globalIndex()==Idfound[k])
                                {
                                    c++;
                                    break;
                                }
                            }
                            if (c==0)
                            {
                                Idfound.push_back(outputIparticle.globalIndex());
                                for (auto const & track2 :jet.tracks_as<soa::Join<aod::Tracks, aod::McTrackLabels>>())
                                {

                                    if(!track2.has_mcParticle())
                                    {
                                        continue;
                                    }
                                    for (int j=0; j<subsize; j++)
                                    {
                                        if (truejet[i][j]==track2.globalIndex())
                                        {
                                            match++;
                                            break;
                                        }
                                    }
                                }

                                break;
                            }
                        }
                    }
                    tempo= match*100./jet.tracks().size();
                    if (tempo>matchpercentage)
                    {
                        temposubsize=subsize;
                        matchpercentage=tempo;
                    }
                }
                sumE+=track.mcParticle().e();
                sumpx+=track.px();
                sumpy+=track.py();
                sumpz+=track.pz();
                sumpt+=track.pt();

                float eta_track=track.eta();
                float phi_track=track.phi();

                delta_eta+=eta_track-eta_jet;
                delta_phi+=phi_track-phi_jet;

            }
            histo.fill(HIST("Truejetmatching"), matchpercentage);
            if (truecounter>0) //Avoid to take into account tracks with no 23, 33, 51
            {
                histo.fill(HIST("Nbtracksmatching"), temposubsize);
            }

            deltaR= sqrt((delta_eta*delta_eta+delta_phi*delta_phi)/(jet.tracks().size()*jet.tracks().size()));
            histo.fill(HIST("Rdistributionjet"), deltaR);

            TLorentzVector jetvector(sumpx, sumpy, sumpz,sumE);
            histo.fill(HIST("2DMomentumjet"), sumpt, jet.pt());
            histo.get<TH2>(HIST("2DMomentumjet"))->SetOption("colz");
            
            float jetphi=jet.phi();
            if (jetphi > M_PI)
            {
                jetphi=jetphi-2*M_PI;
            }
            float reconstructedphi=jetvector.Phi();
            if(reconstructedphi > M_PI)
            {
                reconstructedphi=reconstructedphi - 2*M_PI;
            }

            histo.fill(HIST("2Detajet"), jetvector.Eta(), jet.eta());
            histo.get<TH2>(HIST("2Detajet"))->SetOption("colz");
            histo.fill(HIST("2Dphijet"), reconstructedphi, jetphi);
            histo.get<TH2>(HIST("2Dphijet"))->SetOption("colz");

        }

    }
    PROCESS_SWITCH(MyTask, processTruejet, "Process to compare True jet and Jetfinder", true);

    using ParticleJet = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

    
    void processParticleLevel(ParticleJet const& jets,soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks, aod::McParticles const& mcparticles){
        int nbjet=0;
        for (auto const& jet: jets)
        {
            nbjet++;
            std::ofstream fichier;
            //std::ofstream fichier2;
            if (nbjet==1)
            {
                fichier.open("MotherlistParticleLvl.txt", std::ios::out | std::ios::trunc);
                //fichier2.open("StatuscodeParticleLvl.txt", std::ios::out | std::ios::trunc);
            }else
            {
                fichier.open("MotherlistParticleLvl.txt", std::ios::out | std::ios::app);
                //fichier2.open("StatuscodeParticleLvl.txt", std::ios::out | std::ios::app);
            }
            
            //Array that contain all the mothers in a jet
            //First item: Id, Second item: Count
            std::vector<std::vector<int64_t>> mothers;
            int msize=0;

            std::vector<int> dominantmother(2,0); 

            histo.fill(HIST("Pjetsize"), jet.size());
            histo.fill(HIST("Pjetpt"), jet.pt());
            fichier<<"Jet "<<nbjet<<":"<<std::endl;
            auto hPDG = histo.get<TH1>(HIST("PDGcode4"));   //Fill the PDGcode histogram of final states
            auto hPDG2 = histo.get<TH1>(HIST("PDGcode5"));   //Fill the PDGcode histogram of primordial mother
            for (auto const& jetconstituent: jet.tracks_as<aod::McParticles>())
            {
                auto mother=jetconstituent;
                TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(mother.pdgCode());
                if (particle != nullptr) 
                {
                    fichier<<"|" <<particle->GetName()<<", Status code:"<<mother.getGenStatusCode()<<"|, ";
                }else{
                    fichier<<"|"<<static_cast<long int>(mother.pdgCode())<<", Status code:"<<mother.getGenStatusCode()<< "|, ";
                }

                while (mother.has_mothers()){
                    mother=mother.mothers_first_as<aod::McParticles>();

                    particle = TDatabasePDG::Instance()->GetParticle(mother.pdgCode());
                    if (particle != nullptr) 
                    {
                    fichier<<"|" <<particle->GetName()<<", Status code:"<<mother.getGenStatusCode()<<"|, ";
                    }else{
                        fichier<<"|"<<static_cast<long int>(mother.pdgCode())<<", Status code:"<<mother.getGenStatusCode()<< "|, ";
                    }
                }
                fichier<<std::endl;

                int PDG=jetconstituent.pdgCode();
                if(PDG==321)       { hPDG->Fill("K+",1); }
                else if(PDG==-321)       { hPDG->Fill("K-",1); }
                else if(PDG==-211)  { hPDG->Fill("pi+",1); }
                else if(PDG==211)  { hPDG->Fill("pi-",1); }
                else if(PDG==11)   { hPDG->Fill("e-",1); }
                else if(PDG==-11)  { hPDG->Fill("e+",1); }
                else if(PDG==2212) { hPDG->Fill("p",1); }
                else if(PDG==-2212) { hPDG->Fill("pbar",1); }
                else if(PDG==-13)   { hPDG->Fill("mu+",1); }
                else if(PDG==13)   { hPDG->Fill("mu-",1); }
                else if(PDG==3112) { hPDG->Fill("sigma-",1); }
                else if(PDG==3222) { hPDG->Fill("sigma+",1); }
                else if(PDG==3312) { hPDG->Fill("ksi-",1); }
                else if(PDG==-3312) { hPDG->Fill("ksi+",1); }
                else{ hPDG->Fill("Other",1); }

                auto finalstateParticle=jetconstituent;
                auto outputIparticle=finalstateParticle;

                if (getMotherfromSubprocess(mcparticles, finalstateParticle, outputIparticle ) != false)
                {
                    int counter=0;
                    for(int i=0; i<msize; i++)
                    {
                        if (mothers[i][0]==outputIparticle.globalIndex())
                        {
                            mothers[i][2]++;
                            counter++;
                            break;
                        }
                    }
                    if (counter==0)
                    {
                        mothers.push_back({outputIparticle.globalIndex(),outputIparticle.pdgCode(),1});
                        msize=mothers.size();
                    }
                }
            }
            fichier.close();

            for(int i=0; i<msize; i++)
            {
                if (mothers[i][2]>dominantmother[1])
                {
                    dominantmother[1]=mothers[i][2];
                    dominantmother[0]=mothers[i][1];
                }
            }

        
            if(abs(dominantmother[0])==1){
                hPDG2->Fill("u",1);
            }else if (abs(dominantmother[0])==2){
                hPDG2->Fill("d",1);
            }else if (abs(dominantmother[0])==3){
                hPDG2->Fill("s",1);
            }else if (abs(dominantmother[0])==4){
                hPDG2->Fill("c",1);
            }else if (abs(dominantmother[0])==5){
                hPDG2->Fill("b",1);
            }else if (abs(dominantmother[0])==6){
                hPDG2->Fill("t",1);
            }else if (abs(dominantmother[0])==21){
                hPDG2->Fill("g",1);
            }else if (abs(dominantmother[0])!=0){ // In the case where no mother was found we don't increment
                hPDG2->Fill("other",1);
            }
        }


        std::vector<std::vector<int64_t>> truejet=TruejetfinderParticle(mcparticles);
        int truejetsize=truejet.size();
        
        
        for (int i=0; i<truejetsize; i++)
        {
            float jetpt=0;
            int subsize=truejet[i].size();

            for(int j=0; j<subsize; j++)
            {
                auto particle=mcparticles.rawIteratorAt(truejet[i][j]-mcparticles.offset());

                histo.fill(HIST("PtrueTrackpt"), particle.pt());
                histo.fill(HIST("PtrueTracketa"), particle.eta());
                histo.fill(HIST("PtrueTrackphi"), particle.phi());
                jetpt+=particle.pt();
            }

            auto inputparticle=mcparticles.rawIteratorAt(truejet[i][0]-mcparticles.offset());
            auto PMparticle=inputparticle;
            getMotherfromSubprocess(mcparticles, inputparticle, PMparticle);

            histo.fill(HIST("P2DMomentum"), PMparticle.pt(),jetpt);
            histo.get<TH2>(HIST("P2DMomentum"))->SetOption("colz");

        } 

        

    }
    PROCESS_SWITCH(MyTask, processParticleLevel, "Process at ParticleLevel", false);

    void processMatching(ParticleJet const& Pjets, Detectorjet const& Djets, soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks, aod::McParticles const& mcparticles){
        
        for (auto const& Djet :Djets)
        {
            int nbtrackswithPM=0;
            
            std::vector<std::vector<int64_t>> mothers;
            int msize=0;
            float matching=0;

            std::vector<int64_t> dominantmother(2,0);
            
            for(auto const& dtrack: Djet.tracks_as<soa::Join<aod::Tracks, aod::McTrackLabels>>())
            {
                if(!dtrack.has_mcParticle())
                {
                    continue;
                }

                auto finalstateParticle=dtrack.mcParticle();
                auto outputIparticle=finalstateParticle;
                if (getMotherfromSubprocess(mcparticles, finalstateParticle, outputIparticle) != false)
                {
                    int counter=0;
                    for(int i=0; i<msize; i++)
                    {
                        if (mothers[i][0]==outputIparticle.globalIndex())
                        {
                            mothers[i][1]++;
                            counter++;
                            break;
                        }
                    }
                    if (counter==0)
                    {
                        mothers.push_back({outputIparticle.globalIndex(),1});
                        msize=mothers.size();
                    }
                }
            }

            for(int i=0; i<msize; i++)
            {
                if (mothers[i][1]>dominantmother[1])
                {
                    dominantmother[1]=mothers[i][1];
                    dominantmother[0]=mothers[i][0];
                }
            }

            for(auto const& Pjet: Pjets)
            {
                
                for(auto const& ptrack: Pjet.tracks_as<aod::McParticles>())
                {
                    auto finalstateParticle=ptrack;
                    auto outputIparticle=finalstateParticle;
                    if (getMotherfromSubprocess(mcparticles, finalstateParticle, outputIparticle) != false)
                    {
                        if (outputIparticle.globalIndex()==dominantmother[0])
                        {
                            nbtrackswithPM++;   //Number of tracks with the primordial mother
                        }
                    }
                }
                if(nbtrackswithPM!=0)
                {
                    matching=nbtrackswithPM*100./dominantmother[1];
                    //break;
                }
            }

            /*for(auto const& dtrack: Djet.tracks_as<soa::Join<aod::Tracks, aod::McTrackLabels>>())
            {
                if(!dtrack.has_mcParticle())
                {
                    continue;
                }

                auto finalstateParticle=dtrack.mcParticle();
                auto outputIparticle=finalstateParticle;
                if (getMotherfromSubprocess(mcparticles, finalstateParticle, outputIparticle) != false)
                {
                    for(auto const& Pjet: Pjets)
                    {
                        for(auto const& pparticle: Pjet.tracks_as<aod::McParticles>())
                        {
                            auto finalstateParticle2=pparticle;
                            auto outputIparticle2=pparticle;
                            if (getMotherfromSubprocess(mcparticles, finalstateParticle2, outputIparticle2) != false)
                            {
                                if(outputIparticle.globalIndex()==outputIparticle2.globalIndex())
                                {
                                    nbtrackswithPM++;
                                }
                            }
                        }

                    }
                }

                for(auto const& Pjet: Pjets)
                {
                    for(auto const& pparticle: Pjet.tracks_as<aod::McParticles>())
                    {
                        if(dtrack.mcParticle().globalIndex()==pparticle.globalIndex()){
                            nbtracks++;
                        }
                    }

                }
            }*/
            histo.fill(HIST("Matching"), matching);
        }

    }
    PROCESS_SWITCH(MyTask, processMatching, "Match detector level to particlelevel", false);


    using SelectedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksCov, aod::TracksDCA>>;

    template <typename T> int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }
    int nbSV=0;
    int nbSVb=0;
    
    int nbpoint=0;
    void processSV(Detectorjet::iterator const& jet, SelectedTracks const& tracks, aod::McParticles const& mcparticles, aod::Collisions const& collisions){

        o2::vertexing::DCAFitterN<2> df2;
        df2.setPropagateToPCA(propagateToPCA);
        df2.setMaxR(maxR);
        df2.setMaxDZIni(maxDZIni);
        df2.setMinParamChange(minParamChange);
        df2.setMinRelChi2Change(minRelChi2Change);
        df2.setUseAbsDCA(useAbsDCA);
        df2.setWeightedFinalPCA(useWeightedFinalPCA);
        df2.setBz(bz);
        
        int counterloop1=0;
        double finaldecaylenght=0;
        double finaldcaxy=0;
        double finaldcaz=0;

        int64_t Id1=0;
        int64_t Id2=0;

        double d1=0;
        double d2=0;
        double L=0;

        for (auto const& track1 : jet.tracks_as<SelectedTracks>()) {

            if(!track1.has_mcParticle() || track1.pt()<=1|| track1.dcaXY()>0.2|| track1.dcaZ()>17)
            {
                continue;
            }

            counterloop1++;
            int counterloop2=0;
            auto trackParVarPos1 = getTrackParCov(track1);
            auto collision= track1.collision();
            auto primaryVertex=getPrimaryVertex(collision);

            for (auto const& track2 : jet.tracks_as<SelectedTracks>()) {
                if(!track2.has_mcParticle() || track2.pt()<=1 || track2.dcaXY()>0.2|| track2.dcaZ()>17)
                {
                    continue;
                }
                counterloop2++;
                if (counterloop2>counterloop1){
                    auto trackParVarNeg1 = getTrackParCov(track2);
                    
                    // secondary vertex reconstruction and further 2-prong selections
                    if (df2.process(trackParVarPos1, trackParVarNeg1) == 0) {
                    continue;
                    }
                    //  get secondary vertex
                    const auto& secondaryVertex = df2.getPCACandidate();

                    double distX=primaryVertex.getX()-secondaryVertex[0];
                    double distY=primaryVertex.getY()-secondaryVertex[1];
                    double distZ=primaryVertex.getZ()-secondaryVertex[2];


                    //printf("PV: %f, SV: %f\n", primaryVertex.getX(), secondaryVertex[0]);
                    double decaylenght=sqrt(distX*distX+distY*distY+distZ*distZ);
                    double dcaxy=sqrt(distX*distX+distY*distY);
                    if (decaylenght>finaldecaylenght)
                    {
                        finaldecaylenght=decaylenght;
                        finaldcaxy=dcaxy;
                        /*if(track2.sign()==1)
                        {
                        finaldecaylenght=decaylenght;
                        finaldcaxy=dcaxy;
                        }else{
                        finaldecaylenght=-decaylenght;
                        finaldcaxy=-dcaxy;    
                        }*/
                        finaldcaz=distZ;
                        Id1=track1.globalIndex();
                        Id2=track2.globalIndex();

                        double dx=track1.x()-secondaryVertex[0];
                        double dy=track1.y()-secondaryVertex[1];
                        double dz=track1.z()-secondaryVertex[2];
                        d1=sqrt(dx*dx+dy*dy+dz*dz);

                        dx=track2.x()-secondaryVertex[0];
                        dy=track2.y()-secondaryVertex[1];
                        dz=track2.z()-secondaryVertex[2];
                        d2=sqrt(dx*dx+dy*dy+dz*dz);

                        double Lnorm= sqrt(secondaryVertex[0]*secondaryVertex[0]+secondaryVertex[1]*secondaryVertex[1]+secondaryVertex[2]*secondaryVertex[2]);
                        double Lp=secondaryVertex[0]*jet.px()+secondaryVertex[1]*jet.py()+secondaryVertex[2]*jet.pz();
                        L= Lnorm*sgn(Lp);
                    }
                }
            }
        }

        auto track1=tracks.rawIteratorAt(Id1-tracks.offset());
        auto inputpart1= track1.mcParticle();
        auto track2=tracks.rawIteratorAt(Id2-tracks.offset());
        auto inputpart2= track2.mcParticle();
        auto primordialmother1=inputpart1;
        auto primordialmother2=inputpart2;

        
        double error=sqrt(d1*d1+d2*d2);
        double significance=L/error;
        if (getMotherfromSubprocess(mcparticles,inputpart1,primordialmother1) || getMotherfromSubprocess(mcparticles,inputpart2,primordialmother2))
        {
            nbSV++;
            if(primordialmother1.pdgCode()==5 || primordialmother2.pdgCode()==5)
            {
                nbSVb++;
                histo.fill(HIST("DeclenghtB"), finaldecaylenght);
                histo.fill(HIST("DeclenghtBxy"), finaldcaxy);
                histo.fill(HIST("DeclenghtBz"), finaldcaz); 
                histo.fill(HIST("SignificanceB"), significance);

            }else{
                histo.fill(HIST("Declenght"), finaldecaylenght);
                histo.fill(HIST("Declenghtxy"), finaldcaxy);
                histo.fill(HIST("Declenghtz"), finaldcaz); 
                histo.fill(HIST("Significance"), significance);
            }
        }
        float purity=nbSVb*100./nbSV;
        float efficiency= nbSVb*100./nbBjet;
        printf("Purity: %f, Efficiency: %f\n", purity, efficiency);
        
        
    }
    PROCESS_SWITCH(MyTask, processSV, "Test SV", true);

    

    void processMygraph(aod::Collisions const& collisions, SelectedTracks const& tracks, Detectorjet const& jets,aod::McParticles const& mcparticles){
        for(auto const& track:tracks)
        {
        histo.fill(HIST("Alltrackspt"), track.pt());
        }

        //int jetcounter=0;
        int c=0;
        for(auto const& jet: jets)
        {
            /*if (jet.collisionId()==collision.globalIndex())
            {
            jetcounter++;
            }*/
            c++;
            printf("JET %d: \n", c);
            int trackcounter=0;
            for (auto const& track : jet.tracks_as<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksCov, aod::TracksDCA>>())
            {
                if(!track.has_mcParticle())
                {
                    continue;
                }
                trackcounter++;
                printf("track %d\n", trackcounter);
                auto finalstateParticle=track.mcParticle();
                auto mother=finalstateParticle;
                int countmother=0;
                while (mother.has_mothers()) {
                    countmother++;
                    printf("NewMother");
                    for (auto iMother = mother.mothersIds().front(); iMother <= mother.mothersIds().back(); ++iMother) 
                    {
                        auto mother2 = mcparticles.rawIteratorAt(iMother - mcparticles.offset());
                        if(c<10&&trackcounter==1)
                        {
                            printf("mother pdgcode: %d\n",mother2.pdgCode());
                        }
                        
                    }
                    if (countmother==10)
                    {
                        break;
                    }
                    printf("Test");
                    mother=mother.mothers_first_as<aod::McParticles>();
                }
            }
        }
        //histo.fill(HIST("Nbjetincollision"), jetcounter);
    }
    PROCESS_SWITCH(MyTask, processMygraph, "Process just to plot basics data", true);

};






WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MyTask>(cfgc)};//, TaskName{"trainingcode"})};
};

