//=================================================================================================
//  KDRadiationTree.h
//  Class for controlling and propagating radiation on a KD-tree.
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
#include "Particle.h"
using namespace std;



//=================================================================================================
//  Struct KDRadTreeCell
/// KD radiation-tree cell data structure.
//=================================================================================================
template <int ndim, int nfreq>
struct KDRadTreeCell {
  bool uniform;                     ///< Is cell sufficiently uniform?
  int id;                           ///< i.d. number of cell
  int c1;                           ///< First child cell
  int c2;                           ///< Second child cell
  int cnext;                        ///< i.d. of next cell if not opened
  int k_divide;                     ///< Dimension along which cell is split
  int level;                        ///< Level of cell on tree
  int ifirst;                       ///< i.d. of first particle in cell
  int ilast;                        ///< i.d. of last particle in cell
  int N;                            ///< No. of particles contained in cell
  int Nphoton;                      ///< No. of photon packets intercepted
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
  FLOAT volume;                     ///< Cell volume
  FLOAT rmax;                       ///< ..
  FLOAT lsum[nfreq];                ///< Summation of photon path lengths
  FLOAT opacity[nfreq];             ///< ...
  FLOAT tau;
};



//=================================================================================================
//  Struct MonoIonTreeCell
/// KD radiation-tree cell data structure.
//=================================================================================================
template <int ndim, int nfreq>
struct MonoIonTreeCell : public KDRadTreeCell<ndim,nfreq>
{
  FLOAT Xion;
  FLOAT Xold;
  //FLOAT tau;

  MonoIonTreeCell() {
    Xion = 0.99999;
    Xold = 0.99999;
    //tau = 0.0;
  }

};



//=================================================================================================
//  Class KDRadiationTree
/// \brief   Controls propagation of radiation through KD-tree
/// \details Controls propagation of radiation through KD-tree
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    25/04/2014
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
class KDRadiationTree
{
 public:


  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  KDRadiationTree(int);
  ~KDRadiationTree();


  // Function prototypes
  //-----------------------------------------------------------------------------------------------
  void BuildTree(int, int, ParticleType<ndim> *);
  void AllocateMemory(void);
  void DeallocateMemory(void);
  int ComputeGatherCellList(const FLOAT *, const FLOAT, const int, int *);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void DivideTreeCell(int, int, ParticleType<ndim> *, CellType<ndim,nfreq> &);
  int FindAdjacentCell(const int, const int, const FLOAT *);
  int FindCell(const int, const int, const FLOAT *);
  int FindRayExitFace(CellType<ndim,nfreq> &, const FLOAT *, const FLOAT *, const FLOAT *, FLOAT &);
  void OptimiseTree(void);
  FLOAT QuickSelect(int, int, int, int, ParticleType<ndim> *);
  void StockTree(CellType<ndim,nfreq> &, ParticleType<ndim> *, bool);
  void StockCellProperties(CellType<ndim,nfreq> &, ParticleType<ndim> *, bool);
  void SumRadiationField(const int, CellType<ndim,nfreq> &);


  // Variables
  //-----------------------------------------------------------------------------------------------
  bool allocated_tree;                 ///< Is tree memory allocated?
  int gmax;                            ///< Max. no. of leaf cells
  int gtot;                            ///< Total no. of leaf cells
  int lmax;                            ///< Max. no. of tree levels
  int ltot;                            ///< Total. no. of tree levels
  int ltot_old;                        ///< Previous no. of tree levels
  int Ncell;                           ///< No. of tree cells
  int Ncellmax;                        ///< Max. no. of tree cells
  int Nleafmax;                        ///< Max. no. of particles per leaf cell
  int Ntot;                            ///< Total no. of particles in tree
  int Ntotold;                         ///< Previous no. of particles in tree
  int Ntotmax;                         ///< Max. no. of particles allowed in tree
  int Ntotmaxold;                      ///< Prev. value of Ntotmax
  int Nthreads;                        ///< No. of OpenMP threads
  int *ids;                            ///< Particle ids
  int *inext;                          ///< Linked list for grid search
  CellType<ndim,nfreq> *radcell;       ///< Array of tree cells

};
#endif
