//=================================================================================================
//  NbodySystemTree.h
//  Contains definitions for Nearest neighbour tree class and data structures.
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


#ifndef _NBODY_SYSTEM_TREE_H_
#define _NBODY_SYSTEM_TREE_H_


#include <iostream>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Nbody.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "BinaryOrbit.h"
#include "Parameters.h"
using namespace std;



//=================================================================================================
//  Structure NNTreeCell
/// \brief   Class definition for N-body nearest neighbour tree data structure.
/// \details Class definition for N-body nearest neighbour tree data structure.
/// \author  D. A. Hubber, G. Rosotti
/// \date    10/06/2013
//=================================================================================================
template <int ndim>
struct NNTreeCell {
  int ichild1;                                ///< i.d. of child 1
  int ichild2;                                ///< i.d. of child 2
  int inearest;                               ///< i.d. of nearest node (for building)
  int iparent;                                ///< i.d. of parent node
  int Ncomp;                                  ///< No. of components in node
  int Nstar;                                  ///< No. of stars in node
  int Nchildlist;                             ///< No. of systems in list
  FLOAT radius;                              ///< Radius of bounding sphere of node
  FLOAT rsqdnearest;                         ///< Distance squared to nearest node
  FLOAT r[ndim];                             ///< Position of centre of mass
  FLOAT rpos[ndim];                          ///< Centre of position of node
  FLOAT v[ndim];                             ///< Velocity of centre of mass
  FLOAT a[ndim];                             ///< Acceleration of centre of mass
  FLOAT adot[ndim];                          ///< Jerk of node
  FLOAT a2dot[ndim];                         ///< 2nd derivative of node
  FLOAT a3dot[ndim];                         ///< 3rd derivative of node
  FLOAT m;                                   ///< Mass contained in node
  FLOAT h;                                   ///< Smoothing length of node
  FLOAT gpot;                                ///< Gravitational potential at node centre
  FLOAT gpe;                                 ///< Total gpe of node components
  FLOAT gpe_internal;                        ///< Total internal gpe of node components
  FLOAT tcross;                              ///< Crossing time of node components
  NbodyParticle<ndim>* childlist[Ncompmax];   ///< List of child components
};



//=================================================================================================
//  Class NbodySystemTree
/// \brief   Class definition for N-body nearest neighbour tree.
/// \details Class definition for N-body nearest neighbour tree.  Used to identify and construct
///          sub-systems for efficient integration of bound mutiple systems and high-accuracy
///          integration of close encounters.
/// \author  D. A. Hubber, G. Rosotti
/// \date    10/06/2013
//=================================================================================================
template <int ndim>
class NbodySystemTree
{
 public:

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  NbodySystemTree();
  ~NbodySystemTree();


  // Other function definitions
  //-----------------------------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void ComputeNewBinaryOrbit(const int, const int, const int, NNTreeCell<ndim> &,
                             NNTreeCell<ndim> &, NNTreeCell<ndim> &);
  void CreateNbodySystemTree(Nbody<ndim> *);
  void BuildSubSystems(int, int, FLOAT, Nbody<ndim> *);
  void FindBinarySystems(Nbody<ndim> *);
  void FindPerturberLists(Nbody<ndim> *);
  void OutputBinaryProperties(Nbody<ndim> *);
  void RestockTreeNodes(Nbody<ndim> *);


  // Class variables and main arrays for nearest neighbour tree and binaries
  //-----------------------------------------------------------------------------------------------
  bool allocated_tree;              ///< Is NN-tree memory allocated?
  int Nbinary;                      ///< No. of binary stars
  int Nnode;                        ///< No. of nodes of NN-tree
  int Nnodemax;                     ///< Max. no. of nodes on NN-tree.
  int Norbit;                       ///< No. of binary orbits
  int Norbitmax;                    ///< Max. no. of binary orbits
  int Nquadruple;                   ///< No. of quadruple systems
  int Ntriple;                      ///< No. of triple systems
  FLOAT gpehard;                   ///< Grav. energy limit hard sub-systems
  FLOAT gpesoft;                   ///< Grav. energy limit soft sub-systems
  struct NNTreeCell<ndim> *NNtree;  ///< Main NN-tree array
  struct BinaryOrbit *orbit;        ///< Main binary star array

};
#endif
