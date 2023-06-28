// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
/// \brief Production of a table that contain the first dominant mother of a jet.
/// \author
/// \since

#include "MotherfromSubprocess.h"
#include "PWGJE/DataModel/MotherfromSubprocessTable.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/runDataProcessing.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"


using namespace o2;
using namespace o2::framework;


struct ProduceFirstMother {
    Produces<aod::ChargedMCDetectorLevelJetFirstMother> firstmother;
    using Detectorjet = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;

    void process(Detectorjet::iterator const& jet,soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks, aod::McParticles const& mcparticles){
        //Array that contain all the mothers in a jet
        //First item:ID, Second item: PDGcode, Third item: Counter
        std::vector<std::vector<int64_t>> mothers;
        int msize=0; 

        //First item:ID, Second item: PDGcode, Third item: Counter
        std::vector<int> dominantmother(3,0); 

        for(auto const& track: jet.tracks_as<soa::Join<aod::Tracks, aod::McTrackLabels>>())
        {
            if(!track.has_mcParticle())
            {
                continue;
            }

            auto finalstateParticle=track.mcParticle();
            auto outputIparticle=finalstateParticle;
            if (getMotherfromSubprocess(mcparticles, finalstateParticle, outputIparticle))
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

        for(int i=0; i<msize; i++)
        {
            if (mothers[i][2]>dominantmother[2])
            {
                dominantmother[2]=mothers[i][2];
                dominantmother[1]=mothers[i][1];
                dominantmother[0]=mothers[i][0];
            }
        }
        
        firstmother(dominantmother[0],dominantmother[1]);
        

    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ProduceFirstMother>(cfgc),
  };
}