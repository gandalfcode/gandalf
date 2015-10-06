//=============================================================================
//  formatted_output.h
//  Helper class to output in custom way floats, integers and booleans
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

#include <ostream>
#include <iomanip>
using namespace std;

class formatted_output
{
    private:
        ostream& stream_obj;
        int width_float;
        int width_integer;
        int precision;
    public:
        formatted_output(ostream& obj, int w_f, int w_i, int p):
          stream_obj(obj), width_float(w_f), width_integer(w_i), precision (p) {}

        formatted_output& operator<<(const int& output)
        {
          stream_obj << setw(width_integer) << output;
          return *this;
        }

        formatted_output& operator<<(const float& output)
        {
          stream_obj << setw(width_float) << setprecision(precision) << output;
          return *this;
        }

        formatted_output& operator<<(const double& output)
        {
           stream_obj << setw(width_float) << setprecision(precision) << output;
           return *this;
        }

        formatted_output& operator<<(const bool& output)
        {
          char boolean;
          if (output) boolean='T';
          else boolean='F';
          stream_obj << boolean;
          return *this;
        }

        formatted_output& operator<<(ostream& (*func)(ostream&))
        {
            func(stream_obj);
            return *this;
        }

        template<typename T>
        formatted_output& operator<<(const T& output)
        {
            stream_obj << output;

            return *this;
        }

        void set_width_integer(int width) {width_integer=width;};
};
