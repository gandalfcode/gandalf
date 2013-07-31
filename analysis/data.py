#==============================================================================
#  data.py
#  ..
#
#  This file is part of GANDALF :
#  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
#  https://github.com/gandalfcode/gandalf
#  Contact : gandalfcode@gmail.com
#
#  Copyright (C) 2013  D. A. Hubber, G Rosotti
#
#  GANDALF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  GANDALF is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License (http://www.gnu.org/licenses) for more details.
#==============================================================================
class Data:
    '''Generic class to communicate the data with the plotting part.
    Simple container with x_data, y_data and render_data members,
    which are supposed to be numpy arrays.
    '''
    def __init__(self, x_data, y_data, render_data = None):
        self.x_data = x_data
        self.y_data = y_data
        self.render_data = render_data