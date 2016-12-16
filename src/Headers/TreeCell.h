/*
 * TreeCell.h
 *
 *  Created on: 5 Dec 2016
 *      Author: rosotti
 */

#ifndef SRC_HEADERS_TREECELL_H_
#define SRC_HEADERS_TREECELL_H_


//=================================================================================================
//  Struct TreeCellBase
/// Base tree cell data structure which contains all data elements common to all trees.
//=================================================================================================
template <int ndim>
struct TreeCellBase {
  int cnext;                           ///< i.d. of next cell if not opened
  int copen;                           ///< i.d. of first child cell
  int id;                              ///< Cell id
  int level;                           ///< Level of cell on tree
  int ifirst;                          ///< i.d. of first particle in cell
  int ilast;                           ///< i.d. of last particle in cell
  int N;                               ///< No. of particles in cell
  int Nactive;                         ///< No. of active particles in cell
  int cexit[2][ndim];                  ///< Left and right exit cells (per dim)
  FLOAT cdistsqd;                      ///< Minimum distance to use COM values
  FLOAT mac;                           ///< Multipole-opening criterion value
  Box<ndim> bb ;                       ///< Bounding box
  Box<ndim> hbox;                      ///< Bounding box for smoothing volume
  Box<ndim> vbox ;                     ///< Velocity space bounding box
  FLOAT rcell[ndim];                   ///< Geometric centre of cell bounding box
  FLOAT r[ndim];                       ///< Position of cell COM
  FLOAT v[ndim];                       ///< Velocity of cell COM
  FLOAT m;                             ///< Mass contained in cell
  FLOAT rmax;                          ///< Radius of bounding sphere
  FLOAT hmax;                          ///< Maximum smoothing length inside cell
  FLOAT maxsound;                      ///< Maximum sound speed inside the cell
  FLOAT drmaxdt;                       ///< Rate of change of bounding sphere
  FLOAT dhmaxdt;                       ///< Rate of change of maximum h
  FLOAT q[5];                          ///< Quadrupole moment tensor
  FLOAT amin;                          ///< Minimum grav accel of particles in the cell
  #ifdef MPI_PARALLEL
  double worktot;                      ///< Total work in cell
#endif
};

//=================================================================================================
//  Struct MultipoleMoment
/// Structure to hold the multipole moment data.
//=================================================================================================
template<int ndim>
struct MultipoleMoment {
  MultipoleMoment()
  {
    for (int k=0; k<ndim; k++) r[k] = 0 ;
    for (int k=0; k<5; k++) q[k] = 0 ;
    m = 0;
    id = 0;
  }

  explicit MultipoleMoment(const TreeCellBase<ndim>& cell)
  {
    for (int k=0; k<ndim; k++) r[k] = cell.r[k] ;
    for (int k=0; k<5; k++) q[k] = cell.q[k] ;
    m = cell.m;
    id = cell.id;
  }

  FLOAT r[ndim];                       ///< Position of cell COM
  FLOAT m;                             ///< Mass contained in cell
  FLOAT q[5];                          ///< Quadrupole moment tensor
  int id ;
};




#endif /* SRC_HEADERS_TREECELL_H_ */
