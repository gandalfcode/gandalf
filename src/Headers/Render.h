//=================================================================================================
//  Render.h
//  Contains class and function definitions for generating rendered plots in
//  the python front-end.
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


#ifndef _RENDER_H_
#define _RENDER_H_


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include "Hydrodynamics.h"
#include "Particle.h"
#include "SphSnapshot.h"
#include "Exception.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  Class RenderBase
/// \brief   Parent class for generating rendered images in python.
/// \details Parent class for generating rendered images in python.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
class RenderBase
{
public:

  RenderBase() {};
  virtual ~RenderBase() {};

  static RenderBase* RenderFactory(int ndim, SimulationBase* sim);

#ifdef GANDALF_SNAPSHOT_SINGLE_PRECISION
  virtual int CreateColumnRenderingGrid(const int, const int, const string, const string, const string,
                                const string, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                                float* values, const int Ngrid, SphSnapshotBase &,
                                const string, float& scaling_factor)=0;
  virtual int CreateSliceRenderingGrid(const int, const int, const string, const string,
                                       const string, const string, const string, const SNAPFLOAT,
                                       const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                                       float* values, const int Ngrid, SphSnapshotBase &,
                                       const string , float& scaling_factor)=0;
#else
  virtual int CreateColumnRenderingGrid(const int, const int, const string, const string, const string,
                                const string, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                                double* values, const int Ngrid, SphSnapshotBase &,
                                const string, double& scaling_factor)=0;
  virtual int CreateSliceRenderingGrid(const int, const int, const string, const string,
                                       const string, const string, const string, const SNAPFLOAT,
                                       const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                                       double* values, const int Ngrid, SphSnapshotBase &,
                                       const string , double& scaling_factor)=0;
#endif

#ifdef GANDALF_SNAPSHOT_SINGLE_PRECISION
  static const bool single=true;
#else
  static const bool single=false;
#endif

};



//=================================================================================================
//  Class Render
/// \brief   Class for generating rendered images in python.
/// \details Class for generating rendered images in python.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Render : public RenderBase
{
 public:

  // Constructor and Destructor
  //-----------------------------------------------------------------------------------------------
  Render(SimulationBase* sim);
  ~Render();

  // Subroutine prototypes
  //-----------------------------------------------------------------------------------------------
#ifdef GANDALF_SNAPSHOT_SINGLE_PRECISION
  virtual int CreateColumnRenderingGrid(const int, const int, const string, const string, const string,
                                const string, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                                float* values, const int Ngrid, SphSnapshotBase &,
                                const string, float& scaling_factor);
  virtual int CreateSliceRenderingGrid(const int, const int, const string, const string, const string,
                               const string, const string, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                               const SNAPFLOAT, const SNAPFLOAT,
                               float* values, const int Ngrid,
                               SphSnapshotBase &, const string, float& scaling_factor);
#else
  virtual int CreateColumnRenderingGrid(const int, const int, const string, const string, const string,
                                const string, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                                double* values, const int Ngrid, SphSnapshotBase &,
                                const string, double& scaling_factor);
  virtual int CreateSliceRenderingGrid(const int, const int, const string, const string, const string,
                               const string, const string, const SNAPFLOAT, const SNAPFLOAT, const SNAPFLOAT,
                               const SNAPFLOAT, const SNAPFLOAT,
                               double* values, const int Ngrid,
                               SphSnapshotBase &, const string, double& scaling_factor);
#endif


  Hydrodynamics<ndim>* hydro;          ///< Pointer to Hydrodynamics object to be rendered


};
#endif
