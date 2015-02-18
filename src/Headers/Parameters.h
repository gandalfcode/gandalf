//=============================================================================
//  Parameters.h
//  Class definition for Parameters class.  Contains all data structures and
//  functions for reading and storing simulation parameters before they
//  are processed by the Simulation class.
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


#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_


#include <map>
#include <string>


//=============================================================================
//  Class Parameters
/// \brief   Class for continaing all simulation parameter information.
/// \details Class for continaing all simulation parameter information.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
class Parameters
{
 public:

  // Constructors and destructor
  //---------------------------------------------------------------------------
  Parameters();
  ~Parameters();
  Parameters(const Parameters&);

  // Other function prototypes
  //---------------------------------------------------------------------------
  void CheckInvalidParameters(void);
  void ReadParamsFile(std::string);
  void ParseLine(std::string);
  void SetDefaultValues(void);
  std::string GetParameter (std::string);
  void SetParameter(std::string , std::string);
  void PrintParameters(void);
  void RecordParametersToFile(void);
  std::string TrimWhiteSpace(std::string);

  // Maps containing parameters for all available data types
  //---------------------------------------------------------------------------
  std::map <std::string,int> intparams;             ///< Integer parameters
  std::map <std::string,double> floatparams;        ///< Float parameters
  std::map <std::string,std::string> stringparams;  ///< String parameters

};

#endif
