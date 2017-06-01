//=================================================================================================
//  Flags.h
//  Class for handling bit fl data structures
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

#ifndef _FLAGS_H_
#define _FLAGS_H_




enum flags {

    none = 0,
    dead = 1 << 0,
    active = 1 << 1,
    end_timestep = 1 << 2, // Particle has reached the end of this timestep
    potmin = 1 << 3,

    update_density = 1 << 4, // For meshless
    bad_gradients  = 1 << 5,  // For meshless


    x_periodic_lhs = 1 << 7,
    y_periodic_lhs = 1 << 8,
    z_periodic_lhs = 1 << 9,

    x_periodic_rhs = 1 << 10,
    y_periodic_rhs = 1 << 11,
    z_periodic_rhs = 1 << 12,

  x_periodic = x_periodic_lhs | x_periodic_rhs,
  y_periodic = y_periodic_lhs | y_periodic_rhs,
  z_periodic = z_periodic_lhs | z_periodic_rhs,


    x_mirror_lhs = 1 << 13,
    y_mirror_lhs = 1 << 14,
    z_mirror_lhs = 1 << 15,

    x_mirror_rhs = 1 << 16,
    y_mirror_rhs = 1 << 17,
    z_mirror_rhs = 1 << 18,

  x_mirror = x_mirror_lhs | x_mirror_rhs,
  y_mirror = y_mirror_lhs | y_mirror_rhs,
  z_mirror = z_mirror_lhs | z_mirror_rhs,

  periodic_boundary = x_periodic | y_periodic | z_periodic,
  mirror_boundary   = x_mirror   | y_mirror   | z_mirror,

  any_boundary = periodic_boundary | mirror_boundary,
};


const int periodic_bound_flags[3][2] = {
   {x_periodic_lhs, x_periodic_rhs},
   {y_periodic_lhs, y_periodic_rhs},
   {z_periodic_lhs, z_periodic_rhs},
};


const int mirror_bound_flags[3][2] = {
   {x_mirror_lhs, x_mirror_rhs},
   {y_mirror_lhs, y_mirror_rhs},
   {z_mirror_lhs, z_mirror_rhs},
};


//=================================================================================================
//  Class type_flag
/// \brief  ...
/// \author R. A. Booth
/// \date   31/3/2016
//=================================================================================================
class type_flag {
private:
  unsigned int _flag;


public:
  type_flag(unsigned int flag = none) : _flag(flag) {}

  unsigned int& set(unsigned int flag, bool value) {
    if (value)
      return set(flag) ;
    else
      return unset(flag) ;
  }
  unsigned int& set(unsigned int flag) {
    return _flag |= flag;
  }
  unsigned int& unset(unsigned int flag) {
    return _flag &= ~flag;
  }
  unsigned int get() const {
    return _flag;
  }
  void reset() {
    _flag = none;
  }

  bool check(unsigned int flag) const {
    return (_flag & flag) ;
  }
  bool is_dead() const {
    return _flag & dead;
  }
  bool is_boundary() const {
    return _flag & any_boundary;
  }
  bool is_periodic() const {
    return _flag & periodic_boundary;
  }
  bool is_mirror() const {
    return _flag & mirror_boundary;
  }

};



#endif /* _FLAGS_H_ */
