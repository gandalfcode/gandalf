//=============================================================================
//  NbodySystemTree.h
//  Contains definitions for Nearest neighbour tree class and data structures.
//=============================================================================


#ifndef _NBODY_SYSTEM_TREE_H_
#define _NBODY_SYSTEM_TREE_H_


#include <iostream>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "MergeList.h"
#include "Nbody.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "Parameters.h"
using namespace std;



//=============================================================================
//  Structure NNTreeCell
/// \brief   Class definition for N-body nearest neighbour tree data structure.
/// \details Class definition for N-body nearest neighbour tree data structure.
/// \author  D. A. Hubber, G. Rosotti
/// \date    10/06/2013
//=============================================================================
template <int ndim>
struct NNTreeCell {
  int ichild1;                     ///< i.d. of child 1
  int ichild2;                     ///< i.d. of child 2
  int inearest;                    ///< i.d. of nearest node (for building)
  int iparent;                     ///< i.d. of parent node
  int Ncomp;                       ///< No. of components in node
  int Nstar;                       ///< No. of stars in node
  int Nchildlist;                  ///< No. of systems in list
  DOUBLE rsqdnearest;              ///< Distance squared to nearest node
  DOUBLE r[ndim];                  ///< Position of node
  DOUBLE v[ndim];                  ///< Velocity of node
  DOUBLE a[ndim];                  ///< Acceleration of node
  DOUBLE adot[ndim];               ///< Jerk of node
  DOUBLE a2dot[ndim];              ///< 2nd derivative of node
  DOUBLE a3dot[ndim];              ///< 3rd derivative of node
  DOUBLE m;                        ///< Mass contained in node
  DOUBLE h;                        ///< Smoothing length of node
  DOUBLE gpot;                     ///< Gravitational potential at node centre
  DOUBLE gpe;                      ///< Total gpe of node components
  DOUBLE gpe_internal;             ///< Total internal gpe of node components
  DOUBLE tcross;                   ///< Crossing time of node components
  NbodyParticle<ndim>* childlist[Ncompmax];  ///< List of child components
  MergeList<NbodyParticle<ndim> *> clist;    ///< (To be deleted??)
};



//=============================================================================
//  Class NbodySystemTree
/// \brief   Class definition for N-body nearest neighbour tree.
/// \details Class definition for N-body nearest neighbour tree.  Used to 
///          identify and construct sub-systems for efficient integration of 
///          bound mutiple systems and high-accuracy integration of close 
///          encounters.
/// \author  D. A. Hubber, G. Rosotti
/// \date    10/06/2013
//=============================================================================
template <int ndim>
class NbodySystemTree
{
protected:
  typedef typename MergeList<NbodyParticle<ndim>* >::iterator NbodyListIterator;
 public:

  NbodySystemTree();
  ~NbodySystemTree();

  void AllocateMemory(int);
  void DeallocateMemory(void);
  void CreateNbodySystemTree(Nbody<ndim> *);
  void BuildSubSystems(Nbody<ndim> *);
  void FindPerturberLists(Nbody<ndim> *);

  bool allocated_tree;              ///< Is NN-tree memory allocated?
  int Nnode;                        ///< No. of nodes of NN-tree
  int Nnodemax;                     ///< Max. no. of nodes on NN-tree.
  DOUBLE gpefrac;                   ///< Grav. energy limit for sub-system
  struct NNTreeCell<ndim> *NNtree;  ///< Main NN-tree array

};
#endif
