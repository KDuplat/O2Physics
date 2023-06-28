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

/// \file getMotherfromSubprocess.h
/// \brief Fonction to find the first mother from a track based on the statuscode
///
/// \author Kilian DUPLAT

#ifndef PWGJE_CORE_MOTHERFROMSUBPROCESS_H_
#define PWGJE_CORE_MOTHERFROMSUBPROCESS_H_


/// \param particlesMC  table with MC particles
/// \param particle  MC particle
/// \param OutputParticle particle found with the expected status code (if there is one)
/// \return true if the mother particle is found, false otherwise
#include <cmath>

template <typename T,typename McParticleType>       
    static bool getMotherfromSubprocess(const T& particlesMC,
                        McParticleType particle,
                        McParticleType& OutputParticle)
    {
        auto& mother=particle;
        while (mother.has_mothers()){

            mother=mother.template mothers_first_as<T>();

            auto StatusCodeMother=mother.getGenStatusCode();


            if (abs(StatusCodeMother) == 23 || abs(StatusCodeMother) == 33) {

                OutputParticle=mother;
                return true;
            }

            if (abs(StatusCodeMother) == 51)//Ignore gluon jets
            {
                auto tempo= mother;
                tempo=tempo.template mothers_first_as<T>();
                if(tempo.pdgCode()!=21)
                {
                    tempo=tempo.template mothers_first_as<T>();
                    if(tempo.pdgCode()!=21){
                        OutputParticle=mother;
                        return true;
                    }
                }        
            }
        }
        return false;
    }

#endif

