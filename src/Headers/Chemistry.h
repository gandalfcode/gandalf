//=================================================================================================
//  Chemistry.h
//  Contains main parent virtual class plus child classes for various hydrodynamics
//  algorithms that are implemented (e.g. SPH, Meshless Finite-Volume)
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_


#include <assert.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Particle.h"



enum ChemistrySpecies {                ///< Enumeration of all chemical species
  AB_H2,                               ///< Abundance of H2 (Molecular hydrogen)
  AB_HP,                               ///<       "  "   H?
  NSPECIES                             ///< Total number of chemical species
};


template<ndim>
struct ChemistryParticle {
  FLOAT abundances[NSPECIES];          ///< Abundances of all chemical species
  FLOAT av_mean;                       ///< Mean optical extinction
  FLOAT chi_mean;                      ///< Dust attentuation
  FLOAT dl;                            ///< 0.5*local cell size to compute self-shielding of cell
  FLOAT div_v;                         ///< 0.0
  FLOAT fshield_H2;                    ///< Self-shielding factor for H2
  FLOAT fshield_CO;                    ///<   ""           ""         CO
  FLOAT Tdust;                         ///< Dust temperature
  Particle<ndim> *hydropart;           ///< Pointer to hydro particle in main arrays

  ChemistryParticle() {
    for (k=0; k<NSPECIES; k++) abundances[k] = (FLOAT) 0.0;
    av_mean    = (FLOAT) 0.0;
    chi_mean   = (FLOAT) 1.0;
    div_v      = (FLOAT) 0.0;
    fshield_H2 = (FLOAT) 1.0;
    fshield_CO = (FLOAT) 1.0;
  }
};




//=================================================================================================
//  Class Chemistry
/// \brief   Main parent Chemistry class
/// \details
/// \author  D. A. Hubber, S. Walch, T. Balduin, P. Rohde
/// \date    19/01/2016
//=================================================================================================
template <int ndim>
class Chemistry
{
public:

  FLOAT initialAbundances[NSPECIES];
  ChemistryParticle* chemdata;

  Chemistry(Parameters *);
  ~Chemistry();

  void EvolveAbundances(ChemistryParticle<ndim> &, FLOAT, SimUnits &) = 0;

};




//=================================================================================================
//  Class SCOChemistry
/// \brief   ...
/// \details
/// \author  D. A. Hubber, S. Walch, T. Balduin, P. Rohde
/// \date    19/01/2016
//=================================================================================================
template <int ndim>
class SCOChemistry : public Chemistry
{
public:

  SCOChemistry(Parameters *params);
  ~SCOChemistry();

  void EvolveAbundances(ChemistryParticle<ndim> &, FLOAT, SimUnits &);

};


#endif
