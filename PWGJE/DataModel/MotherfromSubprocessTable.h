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

/// \file MotherfromSubprocessTable.h
/// \brief Produce the firstmother table.
///
/// \author Kilian DUPLAT

#ifndef PWGJE_CORE_MOTHERFROMSUBPROCESSTABLE_H_
#define PWGJE_CORE_MOTHERFROMSUBPROCESSTABLE_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace firstmother
{
    DECLARE_SOA_COLUMN(FirstMotherGlobalIndex, firstmotherglobalIndex, int64_t);
    DECLARE_SOA_COLUMN(FirstMotherPDGcode, firstmotherpdgCode, int);
} // namespace firstmother
    DECLARE_SOA_TABLE(ChargedMCDetectorLevelJetFirstMother, "AOD", "FIRSTMOTHER",
                  firstmother::FirstMotherGlobalIndex, firstmother::FirstMotherPDGcode);
} // namespace o2::aod


#endif //PWGJE_CORE_MOTHERFROMSUBPROCESSTABLE_H_