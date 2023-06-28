#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/DataModel/MotherfromSubprocessTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h" //doit être ajouté pour le workflow en fin de tâche

struct FirstMotherExample
{
    HistogramRegistry histo;

    void init(InitContext const&)
    { 
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

    }

    void process (soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetFirstMother> const& jets){
        
        for(auto const& jet:jets){
            

            int firstmotherpdgcode=jet.firstmotherpdgCode();
            auto hPDG = histo.get<TH1>(HIST("PDGcode"));   //Fill the PDGcode histogram
            if(abs(firstmotherpdgcode)==1){
                hPDG->Fill("u",1);
            }else if (abs(firstmotherpdgcode)==2){
                hPDG->Fill("d",1);
            }else if (abs(firstmotherpdgcode)==3){
                hPDG->Fill("s",1);
            }else if (abs(firstmotherpdgcode)==4){
                hPDG->Fill("c",1);
            }else if (abs(firstmotherpdgcode)==5){
                hPDG->Fill("b",1);
            }else if (abs(firstmotherpdgcode)==6){
                hPDG->Fill("t",1);
            }else if (abs(firstmotherpdgcode)==21){
                hPDG->Fill("g",1);
            }else if (abs(firstmotherpdgcode)!=0){ // In the case where no mother was found we don't increment
                hPDG->Fill("other",1);
            }
        }

    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FirstMotherExample>(cfgc),
  };
}