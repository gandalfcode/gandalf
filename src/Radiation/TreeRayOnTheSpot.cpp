//=================================================================================================
//  TreeRayOnTheSpot.cpp
//  Contains routines for "On-the-spot" physics module in TreeRay.
//  Translated from original Fortran code written by R. Wunsch (2015).
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//            (C) 2015  R. Wunsch
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


#include "TreeRay.h"


//=================================================================================================
//  TreeRayOnTheSpot::TreeRayOnTheSpot
/// Constructor for TreeRay 'On-the-spot' physics module object.  Initialises variables that
/// were initialised in the 'tr_osInit.F90' in the original implementation by R. Wunsch.
//=================================================================================================
template <int ndim, int nfreq, template<int> class TreeCell>
TreeRayOnTheSpot<ndim,nfreq,TreeCell>::TreeRayOnTheSpot
 (int _bhNR) :
 bhNR(_bhNR)
{
  // conversion between eint and T
  TH2ToEint = boltz*TH2 / (AbarH2*tr_mH*(GammaH2 - 1.0));
  THpToEint = boltz*THp / (AbarHp*tr_mH*(GammaHp - 1.0));

  // recombination constant
  recombConst = AlphaStar * powf(Xhydro / tr_mH, 2)/ 3.0;

  // constant for converting photon flux to  energy density
  eflx2Erad = UVPhotonE / (lightSpeed);
  erad2Eflx = 1.0 / eflx2Erad;

  // DAVID : More things in tr_isInit such as reading in sources.  Ask Richard


  // Counts up the number of energy bands.  Check this is consistent with nfreq.
  //tr_nEb = tr_nEb + 1
  //  tr_iEbEUV = tr_nEb
  //  tr_mapEbSoln(tr_nEb) = EUVE_VAR
}



//=================================================================================================
//  TreeRayOnTheSpot::~TreeRayOnTheSpot
/// Destructor for TreeRay 'On-the-spot' physics module object.
//=================================================================================================
template <int ndim, int nfreq, template<int> class TreeCell>
TreeRayOnTheSpot<ndim,nfreq,TreeCell>::~TreeRayOnTheSpot()
{

}



//=================================================================================================
//  TreeRayOnTheSpot::~TreeRayOnTheSpot
/// Destructor for TreeRay 'On-the-spot' physics module object.
//=================================================================================================
/*template <int ndim, int nfreq, template<int> class TreeCell>
void TreeRayOnTheSpot<ndim,nfreq,TreeCell>::NodeContribution(void)
{
  // DAVID : What are these?
  FLOAT node_srcf = 0.0;
  FLOAT node_erad = 0.0;

  for (int i=0; i<ilNI; i++) {

    // Check if the node intersect with the ray
    if (intersectList[IIL(i,ins,iph,ith)] < 0) break;

    // Determine the ray index (ipix) and the weight of the intersection
    ipix = (int) (intersectList[IIL(i,ins,iph,ith)]);
    weight = (intersectList[IIL(i,ins,iph,ith)] - (FLOAT) ipix)/(FLOAT) 0.999;

    for (int j=0; j<bhNR; j++) {
      ir = radNodeMapIndex[IRNM(j,irf,iNodeSize)];
      if (ir < 0) break;
      rays[ipix][ir].srcF[0] += node_srcf*weight*radNodeMapValue[IRNM(j,irf,iNodeSize)];
      rays[ipix][ir].Erad[0] += node_erad*weight*radNodeMapValue[IRNM(j,irf,iNodeSize)];
    }
  }

  return;
}*/



//=================================================================================================
//  TreeRayOnTheSpot::~TreeRayOnTheSpot
/// Destructor for TreeRay 'On-the-spot' physics module object.
//=================================================================================================
template <int ndim, int nfreq, template<int> class TreeCell>
void TreeRayOnTheSpot<ndim,nfreq,TreeCell>::IntegrateRay
 (Rays *ray,                           ///< ..
  FLOAT eflux[nfreq])                  ///< ..
{
  int i;                               // ..
  int ir;                              // ..
  int is;                              // ..
  FLOAT Nphotons;                      // ..
  FLOAT Nrecomb;                       // ..
  FLOAT rho_2mean;                     // ..
  FLOAT dV;                            // ..
  FLOAT nH_av;                         // ..
  FLOAT recombContrib;                 // ..
  FLOAT recombContribSrc;              // ..
  FLOAT eflux_ray[bhNR+1];             // ..
  FLOAT eRadFluxSrc[bhNR+1];           // ..
  FLOAT eFluxSrc[bhNR+1];              // ..
  FLOAT recombSrc[bhNR+1];             // ..
  //real,dimension(0:tr_bhNR) :: eflux_ray, eRadFluxSrc, eFluxSrc, recombSrc

  FLOAT tr_osEflx2Erad = 1.0;  // DAVID : What are these?
  FLOAT tr_osAlphaStar = 1.0;

  // Radiant flux of energy of a source: energy (# of photons) per unit space angle,
  // in the direction of the point-of-calulation, for each source at a given point on the ray
  // initially, fill with fluxes of sources
  for (i=0; i<bhNR+1; i++) eRadFluxSrc[i] = ray[i].srcF[0]/(4.0*pi);

  // Energy (# of photons) per unit area, in the direction of the point-of-calulation,
  // from all sources added together, at each point along the ray
  for (i=0; i<bhNR+1; i++) eflux_ray[i] = 0.0;


  //---------------------------------------------------------------------------------------------
  for (ir=bhNR-1; ir>=0; ir--) {

    // Average particle density between point ir and ir+1
    // DAVID : Add physical constants later (osXhydro/mH)
    nH_av = 0.5*(ray[ir+1].rho + ray[ir].rho);
    //nH_av = (0.5*tr_osXhydro/tr_mH)*(rho_ray(ir+1) + rho_ray(ir))

    // Contribution to recombination for this ray
    recombContrib = min(2.0, eflux_ray[ir+1]*tr_osEflx2Erad / (ray[ir+1].erad[0] + 1.0E-99));

    // Calculate recombination rates belonging to sources
    //-------------------------------------------------------------------------------------------
    for (is=ir+1; is<bhNR+1; is++) {

      // DAVID : What is tr_bhRayR??
      // volume of the element between point ir and ir+1, on the cone coming from source "is"
      //dV = ((tr_bhRayR(is)-tr_bhRayR(ir))**3 - (tr_bhRayR(is)-tr_bhRayR(ir+1))**3)/3.0

      if (eRadFluxSrc[is] > 0.0) {
        // eFluxSrc includes values for ir+1 point
        recombContribSrc = recombContrib*min(1.0, eFluxSrc[is]/(eflux_ray[ir+1] + 1.0e-99));
        recombSrc[is] = tr_osAlphaStar*nH_av*nH_av*dV*recombContribSrc;
        eRadFluxSrc[is] = max((FLOAT) 0.0, eRadFluxSrc[is] - recombSrc[is]);
        //eFluxSrc[is] = eRadFluxSrc[is]/(tr_bhRayR(is)-tr_bhRayR(ir))**2    DAVID : tr_bhRayR
        eflux_ray[ir] = eflux_ray[ir] + eFluxSrc[is];
      }

    }
    //-------------------------------------------------------------------------------------------

  }
  //---------------------------------------------------------------------------------------------

  eflux[0] = 0.5*(eflux_ray[0] + eflux_ray[1]);

  return;
}


//=================================================================================================
//  TreeRayOnTheSpot::~TreeRayOnTheSpot
/// Destructor for TreeRay 'On-the-spot' physics module object.
//=================================================================================================
template <int ndim, int nfreq, template<int> class TreeCell>
void TreeRayOnTheSpot<ndim,nfreq,TreeCell>::FinaliseCell
 (TreeCell<ndim> &cell,                ///< ..
  FLOAT phFluxInt[nfreq],              ///< ..
  FLOAT **eflux,                       ///< ..
  FLOAT **cdMaps)                      ///< ..
{
  FLOAT EradErr = 0.0;
  cell.erdEUV[0] = eflx2Erad*phFluxInt[0];

  // DAVID : Forgot what these are
  /*tr_bhLocEradTot = tr_bhLocEradTot + solnPoint(EUVE_VAR)
  tr_bhLocMionTot = tr_bhLocMionTot + solnPoint(DENS_VAR)*solnPoint(IHP_SPEC)*vol_poc
  EradErr = 2.0*abs(solnPoint(EUVE_VAR) - solnPoint(EUVO_VAR)) /
    (solnPoint(EUVE_VAR) + solnPoint(EUVO_VAR) + 1d-99)*/

  if (EradErr > bhLocRelErr) bhLocRelErr = EradErr;

  return;
}


template class TreeRayOnTheSpot<1, 1, TreeRayCell>;
template class TreeRayOnTheSpot<2, 1, TreeRayCell>;
template class TreeRayOnTheSpot<3, 1, TreeRayCell>;
