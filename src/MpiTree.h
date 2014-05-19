//=============================================================================
//  MpiTree.h
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




#ifndef _MPI_TREE_H_
#define _MPI_TREE_H_


#include <map>
#include <string>
#include <list>
#include "CodeTiming.h"
#include "DomainBox.h"
#include "Precision.h"
#include "Parameters.h"
using namespace std;



//=============================================================================
//  Struct MpiTreeCell
/// ..
//=============================================================================
template <int ndim>
struct MpiTreeCell {
  bool uniform;                     ///< Is cell sufficiently uniform?
  int id;                           ///< ..
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
  FLOAT m;                          ///< Mass contained in cell
  FLOAT volume;                     ///< Cell volume
};



//=============================================================================
//  Class MpiTree
/// \brief   ..
/// \details ..
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    25/04/2014
//=============================================================================
template <int ndim, template<int> class ParticleType>
class MpiTree
{
 public:

  MpiTree(int);
  ~MpiTree();


  //---------------------------------------------------------------------------
  void BuildTree(int, int, ParticleType<ndim> *);
  void AllocateMemory(void);
  void DeallocateMemory(void);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void DivideTreeCell(int, int, ParticleType<ndim> *, KDRadTreeCell<ndim> &);
  int FindCell(int, FLOAT *);
  FLOAT QuickSelect(int, int, int, int, ParticleType<ndim> *);


  //---------------------------------------------------------------------------
  bool allocated_tree;
  int gtot;
  int lmax;
  int ltot;
  int Ncell;
  int Ncellmax;
  int Ntot;
  int Ntotmax;
  int Nthreads;                     ///< No. of OpenMP threads
  int *ids;                         ///< Particle ids
  int *inext;                       ///< Linked list for grid search
  MpiTreeCell<ndim> *mpicell;       ///< ..

};
#endif
