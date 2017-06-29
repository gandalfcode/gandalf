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
///          The controller handles a maximum of 32 time-step levels. The timestep levels differ by
///          a factor of two in length and are arranged so that particles remain synchronised.
///          Particles on level 0 have the largest timestep, and those on nlevelmax have the
///          smallest.
/// \author  R. A. Booth
/// \date    26/09/2016
//=================================================================================================
class TimeStepController {
public:
  TimeStepController(unsigned int num_levels=__maxlevels)
    : _tCurrent(0),
      _nstepCurrent(0),
      _nlevelmax(num_levels),
      _tot_num_steps(0)
  {
    assert( _nlevelmax <= __maxlevels) ;

    _levelfactor[0] = 1 ;
    for (unsigned int i=1; i < __maxlevels; ++i)
      _levelfactor[i] = 2*_levelfactor[i-1] ;
  } ;

  // rest_time_step_hierachy
  /// \brief Reset the time-step hierarchy with a new step size and maximum
  ///        time-step level. Must only be called when synchronised.
  void reset_timestep_hierarchy(double dt_max, int levelmax) {
    assert(_nlevelmax <= __maxlevels) ;
    _nlevelmax = levelmax ;
    reset_timestep_hierarchy(dt_max) ;
  }
  void reset_timestep_hierarchy(double dt_max) {
    assert(is_resync_point()) ;

    for (unsigned int i=0; i < __maxlevels; ++i)
      _dtLevel[i] = dt_max / _levelfactor[i] ;

    _nstepCurrent = 0 ;
  }

  // update_step
  /// \brief Increment the time and counters
  void update_step() {
    _tCurrent += _dtLevel[_nlevelmax] ;
    _nstepCurrent++ ;
    _tot_num_steps++ ;
  }

  // change_max_level
  /// \brief    Change the number of active time-step levels.
  /// \details  Updates the number of active levels without
  ///           breaking the level synchronisation. If the requested
  ///           level is not synchronised then we will only reduce
  ///           the number of levels as far as possible.
  /// \returns  The number of active levels after the reduction
  ///           has been applied. Guaranteed to be at least the
  ///           requested number of levels.
  int change_num_active_levels(unsigned int level) {

    // We can always increase the number of levels
    if (level >= _nlevelmax) {
      unsigned int ratio =
        _levelfactor[level] / _levelfactor[_nlevelmax] ;

      _nstepCurrent *= ratio ;
      _nlevelmax = level ;

      return level ;
    }
    else {
      unsigned int ratio =
              _levelfactor[_nlevelmax] / _levelfactor[level] ;

      // Find the best synchronised level
      while (_nstepCurrent % ratio) {
        level++ ;
        ratio /= 2 ;
      }
      assert(ratio >= 1) ;

      _nstepCurrent /= ratio ;
      _nlevelmax = level ;
      return level ;
    }
  }

  // compute_time_step_level
  /// \brief Compute the level necessary to ensure the time step is smaller
  ///        than the requested value.
  int compute_timestep_level(double dt) const {
    return max((int) (log(_dtLevel[0]/dt)/log(2)) + 1, 0) ;
  }


  // tlast
  // \brief The last time this level was active
  double tlast(int level) const {
    unsigned int dn = _nstepCurrent - nlast(level) ;
    return current_time() - dn * _dtLevel[_nlevelmax] ;
  }
  // nlast
  // \brief The last step number when this level was active
  double nlast(int level) const {
    unsigned int nstepLevel = _levelfactor[level] ;
    return (_nstepCurrent / nstepLevel) * nstepLevel ;
  }

  // timestep
  /// \brief The time step of the given level
  double timestep(int level) const {
    return _dtLevel[level] ;
  }
  // current_time
  /// \brief The current time
  double current_time() const {
    return _tCurrent ;
  }
  // step_number
  /// \brief The number of steps since the last resync point
  unsigned int step_number() const {
    return _nstepCurrent ;
  }
  // total_step_number
  /// \brief The total number of steps in taken in the simulation
  unsigned long long total_step_number() const {
    return _tot_num_steps ;
  }

  // is_resync_point
  /// \brief Determine whether all time-step levels are synchronised
  bool is_resync_point() const {
    return is_synchronised(0);
  }
  // is_synchronised
  /// \brief Determine whether a timestep level is synchronised (active).
  bool is_synchronised(int level) const {
    return is_active(level) ;
  }
  // is_active
  /// \brief Determine whether particles on a timestep level are active
  bool is_active(int level) const {
    return (_nstepCurrent % _levelfactor[_nlevelmax-level]) == 0;
  }

  // min_synchronised_level
  /// \brief Determine the lowest timestep level that is currently
  ///        synchronised.
  int min_synchronised_level() const {
    int level = _nlevelmax;
    while (level > 0 && is_active(level))
      level-- ;

    return ++level ;
  }


private:
  static const unsigned int __maxlevels = 32 ;

  double _tCurrent;
  double _dtLevel[__maxlevels] ;
  unsigned int _levelfactor[__maxlevels] ;
  unsigned int _nlevelmax, _nstepCurrent ;
  unsigned long long int _tot_num_steps ;
};



#endif//_TIME_STEP_CONTROL_H_
