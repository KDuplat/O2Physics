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

/// \file Truejetfinder.h
/// \brief Fonction to find jets based on the Track table based on the statuscode. It returns an array of track Id.
///
/// \author Kilian DUPLAT

#ifndef PWGJE_CORE_TRUEJETFINDER_H_
#define PWGJE_CORE_TRUEJETFINDER_H_

#include "/home/kduplat/Documents/MotherfromSubprocess.h"


/// \param particlesMC  table with MC particles
/// \param Tracks   table with Tracks

template <typename T,typename P>
    static std::vector<std::vector<int64_t>> Truejetfinder(const P& particlesMC,
                        const T& Tracks)
        {
            std::vector<std::vector<int64_t>>truejet;       //Jet found only base on the statuscode
            int truejetsize=0;


            for (auto const& track: Tracks)  // Loop that find all the truejets
            {   
                
                if(!track.has_mcParticle())
                {
                    continue;
                }

                auto finalstateParticle=track.mcParticle();
                auto outputIparticle=track.mcParticle();

                int counter=0;

                if (getMotherfromSubprocess(particlesMC, finalstateParticle, outputIparticle))
                {
                    for (int i=0; i<truejetsize; i++ )
                    {
                        auto trackinput=Tracks.rawIteratorAt(truejet[i][0]-Tracks.offset());
                        auto inputparticle=trackinput.mcParticle();
                        auto firstmother=inputparticle;
                        getMotherfromSubprocess(particlesMC, inputparticle, firstmother);

                        if (outputIparticle.globalIndex()==firstmother.globalIndex())
                        {
                            truejet[i].push_back(track.globalIndex());
                            counter++;
                            break;
                        }
                    }
                    if (counter==0){
                        truejet.push_back({track.globalIndex()});
                        truejetsize=truejet.size();
                    }

                }
            }
            return truejet;
        }

template <typename P>
    static std::vector<std::vector<int64_t>> TruejetfinderParticle(const P& particlesMC)
        {
            std::vector<std::vector<int64_t>>truejet;       //Jet found only base on the statuscode
            int truejetsize=0;

            for (auto const& particle: particlesMC)  // Loop that find all the truejets
            {   

                auto outputIparticle=particle;
                int counter=0;

                if (getMotherfromSubprocess(particlesMC, particle, outputIparticle))
                {
                    for (int i=0; i<truejetsize; i++ )
                    {
                        auto inputparticle=particlesMC.rawIteratorAt(truejet[i][0]-particlesMC.offset());
                        auto firstmother=inputparticle;
                        getMotherfromSubprocess(particlesMC, inputparticle, firstmother);

                        if (outputIparticle.globalIndex()==firstmother.globalIndex())
                        {
                            truejet[i].push_back(particle.globalIndex());
                            counter++;
                            break;
                        }
                    }
                    if (counter==0){
                        truejet.push_back({particle.globalIndex()});
                        truejetsize++;
                    }

                }
            }
            return truejet;
        }

#endif