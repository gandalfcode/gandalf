//=============================================================================
//  KDRadiationTree.h
//  ..
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
//=============================================================================




#ifndef _KD_RADIATION_TREE_H_
#define _KD_RADIATION_TREE_H_


#include <map>
#include <string>
#include <list>
#include "CodeTiming.h"
#include "DomainBox.h"
#include "Precision.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "SphParticle.h"
using namespace std;



//=============================================================================
//  Struct KDRadTreeCell
/// KD radiation-tree cell data structure.
//=============================================================================
template <int ndim>
struct KDRadTreeCell {
  bool uniform;                     ///< Is cell sufficiently uniform?
  int id;
  int c1;                           ///< First child cell
  int c2;                           ///< Second child cell
  int cnext;                        ///< i.d. of next cell if not opened
  int k_divide;                     ///< Dimension along which cell is split
  int level;                        ///< Level of cell on tree
  int ifirst;                       ///< i.d. of first particle in cell
  int ilast;                        ///< i.d. of last particle in cell
  int N;                            ///< No. of particles contained in cell
  int cexit[2][ndim];               ///< Left and right exit cells (per dim)
  FLOAT bbmin[ndim];                ///< Minimum extent of bounding box
  FLOAT bbmax[ndim];                ///< Maximum extent of bounding box
  FLOAT r[ndim];                    ///< COM position of cell
  FLOAT rcell[ndim];                ///< Geometric centre of cell
  FLOAT v[ndim];                    ///< Velocity of COM of cell
  FLOAT m;                          ///< Mass contained in cell
  FLOAT rho;                        ///< Average density in cell
  FLOAT sigma_rho;                  ///< Variance of density
  FLOAT temp;                       ///< Average temperature in cell
  FLOAT uphoton;                    ///< Photon energy density
  FLOAT volume;                     ///< Cell volume
  FLOAT opacity;                    ///< Opacity of cell
};



//=============================================================================
//  Class KDRadiationTree
/// \brief   Controls propagation of radiation through KD-tree
/// \details Controls propagation of radiation through KD-tree
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    25/04/2014
//=============================================================================
template <int ndim, template<int> class ParticleType>
class KDRadiationTree
{
 public:

  KDRadiationTree(int);
  ~KDRadiationTree();


  //---------------------------------------------------------------------------
  void BuildTree(int, int, ParticleType<ndim> *);
  void AllocateMemory(void);
  void DeallocateMemory(void);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void DivideTreeCell(int, int, ParticleType<ndim> *, KDRadTreeCell<ndim> &);
  int FindCell(int, FLOAT *);
  void OptimiseTree(void);
  FLOAT QuickSelect(int, int, int, int, ParticleType<ndim> *);
  void StockTree(KDRadTreeCell<ndim> &, ParticleType<ndim> *);
  void StockCellProperties(KDRadTreeCell<ndim> &, ParticleType<ndim> *);


  //---------------------------------------------------------------------------
  bool allocated_tree;
  int gmax;
  int gtot;
  int lmax;
  int ltot;
  int ltot_old;
  int Ncell;
  int Ncellmax;
  int Nleafmax;
  int Ntot;
  int Ntotold;
  int Ntotmax;
  int Ntotmaxold;
  int Nthreads;                     ///< No. of OpenMP threads
  int *ids;                         ///< Particle ids
  int *inext;                       ///< Linked list for grid search
  KDRadTreeCell<ndim> *radcell;

};
#endif
