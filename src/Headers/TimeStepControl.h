//=================================================================================================
//  TimeStepControl.h
//  Contains class & functions related to controlling the time-step heirachy
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

#ifndef _TIME_STEP_CONTROL_H_
#define _TIME_STEP_CONTROL_H_

//=================================================================================================
//  Class TimeStepController
/// \brief Handles individual particle time-steps
/// \details Handles the time-step information for particles on time step levels [0, nlevelmax].
///          The controller handles a maximum of 32 time-step levels.
/// \author  R. A. Booth
/// \date    26/09/2016
//=================================================================================================
class TimeStepController {
public:
  TimeStepController(int num_levels=32)
    : _tCurrent(0),
      _nlevelmax(num_levels)
  {
    assert(_nlevelmax <= 32) ;

    _levelfactor[0] = 1 ;
    for (int i=1; i < 32; ++i)
      _levelfactor[i] = 2*_levelfactor[i-1] ;
  } ;

  void reset_timestep_hierarchy(double dt, int levelmax) {
    assert(_nlevelmax <= 32) ;
    _nlevelmax = levelmax ;
    reset_timestep_hierarchy(dt) ;
  }
  void reset_timestep_hierarchy(double dt_min) {

    for (unsigned int i=0; i <= _nlevelmax; ++i)
      _dtLevel[i] = dt_min * _levelfactor[_nlevelmax - i] ;

    _dtCurrent = dt_min ;
    _nstepCurrent = 0 ;
  }
  void update_step() {
    _tCurrent += _dtCurrent ;
    _nstepCurrent += 1 ;
  }


  int compute_timestep_level(double dt) const {
    return max((int) (log(_dtLevel[0]/dt)/log(2)) + 1, 0) ;
  }
  double tlast(int level) const ;
  double nlast(int nlevel) const ;

  double timestep(int level) const {
    return _dtLevel[level] ;
  }
  double time() const {
    return _tCurrent ;
  }
  unsigned int step_number() const {
    return _nstepCurrent ;
  }

  bool is_resync_point() const {
    return _nstepCurrent == _levelfactor[_nlevelmax] ;
  }

  bool is_active(int level) const {
    return (_nstepCurrent % _levelfactor[_nlevelmax-level]) == 0;
  }


private:
  double _tCurrent, _dtCurrent ;
  double _dtLevel[32] ;
  unsigned int _levelfactor[32] ;
  unsigned int _nlevelmax, _nstepCurrent ;
};



#endif//_TIME_STEP_CONTROL_H_
