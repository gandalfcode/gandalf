//=============================================================================
//  NbodySystemTree.h
//  Header file containing class definition for ..
//=============================================================================


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
#include "Parameters.h"
using namespace std;



//=============================================================================
//  Structure NNTreeCell
/// \brief   ..
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    10/06/2013
//=============================================================================
template <int ndim>
struct NNTreeCell {
  int ichild1;                     ///< ..
  int ichild2;                     ///< ..
  int inearest;                    ///< ..
  int iparent;                     ///< ..
  int Ncomp;                       ///< ..
  DOUBLE rsqdnearest;              ///< ..
  DOUBLE r[ndim];                  ///< ..
  DOUBLE v[ndim];                  ///< ..
  DOUBLE a[ndim];                  ///< ..
  DOUBLE adot[ndim];               ///< ..
  DOUBLE m;                        ///< ..
  DOUBLE h;                        ///< ..
  DOUBLE gpot;                     ///< ..
  DOUBLE gpe;                      ///< ..
  DOUBLE gpe_internal;             ///< ..
  DOUBLE tcross;
};



//=============================================================================
//  Class NbodySystemTree
/// \brief   ..
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    10/06/2013
//=============================================================================
template <int ndim>
class NbodySystemTree
{
 public:

  NbodySystemTree();
  ~NbodySystemTree();

  void AllocateMemory(int);
  void DeallocateMemory(void);
  void CreateNbodySystemTree(Nbody<ndim> *);
  void StockNbodySystemTree(Nbody<ndim> *);
  void BuildSubSystems(Nbody<ndim> *);
  void FindPerturberLists(Nbody<ndim> *);

  bool allocated_tree;
  int Nnode;
  int Nnodemax;
  struct NNTreeCell<ndim> *NNtree;

};
#endif
