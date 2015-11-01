//=================================================================================================
//  Parameters.cpp
//  Contains all functions for calculating default values and reading in
//  new values from the simulation parameter file.
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


#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Exception.h"
#include "Parameters.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  Parameters::Parameters
/// Constructor for Parameters class
//=================================================================================================
Parameters::Parameters()
{
  SetDefaultValues();
}



//=================================================================================================
//  Parameters::~Parameters
/// Parameters destructor
//=================================================================================================
Parameters::~Parameters()
{
}


//=================================================================================================
//  Parameters::Parameters
/// Alternative Parameters constructor
//=================================================================================================
Parameters::Parameters(const Parameters& other)
{
  this->intparams = other.intparams;
  this->stringparams = other.stringparams;
  this->floatparams = other.floatparams;
}



//=================================================================================================
//  Parameters::ReadParamsFile
/// Read and parse parameter file 'filename'.  If file doesn't exist, or
/// file does not contain a simulation run id, then quit program here.
//=================================================================================================
void Parameters::ReadParamsFile
 (string filename)                     ///< [in] Parameters file to be read
{
  ifstream inputfile;                  // Input file stream
  std::string line;                    // Parameter file line

  debug1("[Parameters::ReadParamsFile]");

  // Set-up all parameters and assign default values
  SetDefaultValues();

  // If parameter file can be opened, parse each line in turn.
  // Else, quit program with exception
  inputfile.open(filename.c_str(), ios::in);
  if (inputfile.is_open()) {
    while ( inputfile.good() ) {
      getline(inputfile, line);
      ParseLine(line);
    }
  }
  else {
    string message = "The specified parameter file: " + filename +
      " does not exist, aborting";
    ExceptionHandler::getIstance().raise(message);
  }
  inputfile.close();

  // Now verify that parameters file contains a run id in order to generate
  // labelled output files.  If not defined, then quit program with exception.
  if (stringparams["run_id"] == "") {
    string message = "The parameter file: " + filename +
      " does not contain a run id string, aborting";
    ExceptionHandler::getIstance().raise(message);
  }

  // Check if any deactivated options have been selected or not, and also
  // verify certain key parameters are valid.
  CheckInvalidParameters();

  return;
}



//=================================================================================================
//  Parameters::ParseLine
/// Parse a single line read from the parameters file.
/// Identifies if the line is in the form 'Comments : variable = value', or
/// 'variable = value', and if so, stores value in memory.
/// If line begins with a hash charatcer '#', then ignore line as a comment.
///================================================================================================
void Parameters::ParseLine
 (string paramline)                      ///< [in] Line from parameters file to be parsed.
{
  // First, trim all white space from line
  paramline = TrimWhiteSpace(paramline);

  int colon_pos = paramline.find(':');   // Position of colon in string
  int equal_pos = paramline.find('=');   // Position of equals in string
  int hash_pos  = paramline.find('#');   // Position of hash in string
  int length    = paramline.length();    // Length of string

  // Ignore line if it is a comment (i.e. begins with a hash character)
  if (hash_pos == 0 || length == 0) return;

  // If line is not in the correct format (either equals is not present,
  // or equals is before the colon) then skip line and return.
  if (equal_pos >= length || (colon_pos < length && colon_pos > equal_pos)) return;

  // Extract variable name and value from line
  std::string var_name = paramline.substr(colon_pos+1, equal_pos-colon_pos-1);
  std::string var_value = paramline.substr(equal_pos+1, length-equal_pos-1);

  // Finally, set parameter value in memory
  SetParameter(var_name, var_value);

  return;
}



//=================================================================================================
//  Parameters::SetDefaultValues
/// Record all parameter variable names in memory also setting default values.
//=================================================================================================
void Parameters::SetDefaultValues(void)
{
  debug1("[Parameters::SetDefaultValues]");

  // Main simulation algorithm parameters
  //-----------------------------------------------------------------------------------------------
  intparams["ndim"] = 3;
  stringparams["sim"] = "sph";
  stringparams["nbody"] = "hermite4";

  // Simulation id, filename and output time parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["ic"] = "box";
  stringparams["run_id"] = "";
  stringparams["in_file"] = "";
  stringparams["in_file_form"] = "su";
  stringparams["out_file_form"] = "su";
  floatparams["tend"] = 1.0;
  floatparams["tmax_wallclock"] = 9.99e20;
  floatparams["dt_snap"] = 0.2;
  floatparams["tsnapfirst"] = 0.2;
  intparams["Nstepsmax"] = 99999999;
  intparams["noutputstep"] = 128;
  intparams["ndiagstep"] = 1024;
  intparams["nrestartstep"] = 512;
  intparams["litesnap"] = 0;
  floatparams["dt_litesnap"] = 0.2;
  floatparams["tlitesnapfirst"] = 0.0;

  // Unit and scaling parameters
  //-----------------------------------------------------------------------------------------------
  intparams["dimensionless"] = 0;
  stringparams["rinunit"] = "";
  stringparams["minunit"] = "";
  stringparams["tinunit"] = "";
  stringparams["vinunit"] = "";
  stringparams["ainunit"] = "";
  stringparams["rhoinunit"] = "";
  stringparams["sigmainunit"] = "";
  stringparams["pressinunit"] = "";
  stringparams["finunit"] = "";
  stringparams["Einunit"] = "";
  stringparams["mominunit"] = "";
  stringparams["angmominunit"] = "";
  stringparams["angvelinunit"] = "";
  stringparams["dmdtinunit"] = "";
  stringparams["Linunit"] = "";
  stringparams["kappainunit"] = "";
  stringparams["Binunit"] = "";
  stringparams["Qinunit"] = "";
  stringparams["Jcurinunit"] = "";
  stringparams["uinunit"] = "";
  stringparams["dudtinunit"] = "";
  stringparams["tempinunit"] = "";
  stringparams["routunit"] = "pc";
  stringparams["moutunit"] = "m_sun";
  stringparams["toutunit"] = "myr";
  stringparams["voutunit"] = "km_s";
  stringparams["aoutunit"] = "km_s2";
  stringparams["rhooutunit"] = "g_cm3";
  stringparams["sigmaoutunit"] = "m_sun_pc2";
  stringparams["pressoutunit"] = "Pa";
  stringparams["foutunit"] = "N";
  stringparams["Eoutunit"] = "J";
  stringparams["momoutunit"] = "m_sunkm_s";
  stringparams["angmomoutunit"] = "m_sunkm2_s";
  stringparams["angveloutunit"] = "rad_s";
  stringparams["dmdtoutunit"] = "m_sun_yr";
  stringparams["Loutunit"] = "L_sun";
  stringparams["kappaoutunit"] = "m2_kg";
  stringparams["Boutunit"] = "tesla";
  stringparams["Qoutunit"] = "C";
  stringparams["Jcuroutunit"] = "C_s_m2";
  stringparams["uoutunit"] = "J_kg";
  stringparams["dudtoutunit"] = "J_kg_s";
  stringparams["tempoutunit"] = "K";

  // Integration scheme and timestep parameters
  //-----------------------------------------------------------------------------------------------
  floatparams["accel_mult"] = 0.3;
  floatparams["courant_mult"] = 0.15;
  floatparams["nbody_mult"] = 0.1;
  floatparams["subsys_mult"] = 0.05;
  intparams["Nlevels"] = 1;
  intparams["level_diff_max"] = 1;
  intparams["sph_single_timestep"] = 0;
  intparams["nbody_single_timestep"] = 0;

  // SPH parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["sph_integration"] = "lfkdk";
  stringparams["kernel"] = "m4";
  intparams["tabulated_kernel"] = 1;
  floatparams["h_fac"] = 1.2;
  floatparams["h_converge"] = 0.01;

  // Thermal physics parameters
  //-----------------------------------------------------------------------------------------------
  intparams["hydro_forces"] = 1;
  stringparams["gas_eos"] = "energy_eqn";
  stringparams["energy_integration"] = "null";
  floatparams["energy_mult"] = 0.4;
  floatparams["gamma_eos"] = 1.66666666666666;
  floatparams["temp0"] = 1.0;
  floatparams["mu_bar"] = 1.0;
  floatparams["rho_bary"] = 1.0e-14;
  floatparams["eta_eos"] = 1.4;
  stringparams["radws_table"] = "eos.bell.cc.dat";
  floatparams["temp_ambient"] = 5.0;

  // Artificial viscosity parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["avisc"] = "mon97";
  stringparams["acond"] = "none";
  stringparams["time_dependent_avisc"] = "none";
  floatparams["alpha_visc"] = 1.0;
  floatparams["alpha_visc_min"] = 0.1;
  floatparams["beta_visc"] = 2.0;

  // Meshless Finite-Volume parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["riemann_solver"] = "exact";
  stringparams["slope_limiter"] = "springel2009";
  intparams["zero_mass_flux"] = 0;
  intparams["static_particles"] = 0;

  // Gravity parameters
  //-----------------------------------------------------------------------------------------------
  intparams["self_gravity"] = 0;
  intparams["kgrav"] = 1;
  stringparams["grav_kernel"] = "mean_h";
  stringparams["external_potential"] = "none";
  floatparams["avert"] = -0.5;
  floatparams["rplummer_extpot"] = 1.0;
  floatparams["mplummer_extpot"] = 1.0;

  // Neighbour searching and tree-gravity parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["neib_search"] = "kdtree";
  stringparams["gravity_mac"] = "geometric";
  stringparams["multipole"] = "quadrupole";
  intparams["Nleafmax"] = 6;
  intparams["ntreebuildstep"] = 1;
  intparams["ntreestockstep"] = 1;
  floatparams["thetamaxsqd"] = 0.1;
  floatparams["macerror"] = 0.0001;

  // N-body parameters
  //-----------------------------------------------------------------------------------------------
  intparams["sub_systems"] = 0;
  stringparams["sub_system_integration"] = "hermite4";
  intparams["Npec"] = 1;
  intparams["nbody_softening"] = 1;
  intparams["perturbers"] = 0;
  intparams["binary_stats"] = 0;
  intparams["nsystembuildstep"] = 1;
  floatparams["gpefrac"] = 5.0e-2;
  floatparams["gpesoft"] = 2.0e-2;
  floatparams["gpehard"] = 1.0e-3;

  // Sink particle parameters
  //-----------------------------------------------------------------------------------------------
  intparams["sink_particles"] = 0;
  intparams["create_sinks"] = 0;
  intparams["smooth_accretion"] = 0;
  intparams["fixed_sink_mass"] = 0;
  intparams["extra_sink_output"] = 0;
  floatparams["rho_sink"] = 1.e-12;
  floatparams["alpha_ss"] = 0.01;
  floatparams["sink_radius"] = 2.0;
  floatparams["smooth_accrete_frac"] = 0.01;
  floatparams["smooth_accrete_dt"] = 0.01;
  stringparams["sink_radius_mode"] = "hmult";

  // Radiation algortihm parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["radiation"] = "none";
  intparams["Nraditerations"] = 2;
  intparams["Nradlevels"] = 1;
  intparams["nradstep"] = 1;
  floatparams["Nphotonratio"] = 8;
  floatparams["mu_ion"] = 0.678;
  floatparams["temp_ion"] = 1e4;
  floatparams["arecomb"] = 2.7e-13;
  floatparams["Ndotmin"] = 1e47;
  floatparams["NLyC"] = 1e47;

  // TreeRay algorithm parameters
  //-----------------------------------------------------------------------------------------------
  intparams["on_the_spot"] = 0;
  intparams["nside"] = 4;
  intparams["ilNR"] = 50;
  intparams["ilNTheta"] = 25;
  intparams["ilNPhi"] = 50;
  intparams["ilNNS"] = 20;
  intparams["ilFinePix"] = 4;
  floatparams["maxDist"] = 1.0e99;
  floatparams["rayRadRes"] = 1.0;
  floatparams["relErr"] = 0.01;
  stringparams["errControl"] = "erad_tot";

  // Boundary conditions parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["boundary_lhs[0]"] = "open";
  stringparams["boundary_rhs[0]"] = "open";
  stringparams["boundary_lhs[1]"] = "open";
  stringparams["boundary_rhs[1]"] = "open";
  stringparams["boundary_lhs[2]"] = "open";
  stringparams["boundary_rhs[2]"] = "open";
  floatparams["boxmin[0]"] = -9.9e30;
  floatparams["boxmin[1]"] = 9.9e30;
  floatparams["boxmin[2]"] = -9.9e30;
  floatparams["boxmax[0]"] = 9.9e30;
  floatparams["boxmax[1]"] = -9.9e30;
  floatparams["boxmax[2]"] = 9.9e30;

  // Ewald periodic gravity parameters
  //-----------------------------------------------------------------------------------------------
  intparams["ewald"] = 1;
  intparams["gr_bhewaldseriesn"] = 10;
  intparams["in"] = 500;
  intparams["nEwaldGrid"] = 16;
  floatparams["ewald_mult"] = 1.0;
  floatparams["ixmin"] = 1.0e-8;
  floatparams["ixmax"] = 5.0;
  floatparams["EFratio"] = 1.0;

  // Initial conditions parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["particle_distribution"] = "cubic_lattice";
  intparams["smooth_ic"] = 0;
  intparams["com_frame"] = 0;
  intparams["regularise_particle_ics"] = 0;
  intparams["Nreg"] = 1;
  intparams["field_type"] = 1;
  intparams["gridsize"] = 64;
  intparams["Nhydro"] = 0;
  intparams["Nhydromax"] = -1;
  intparams["Nstar"] = 0;
  intparams["Nstarmax"] = -1;
  intparams["Nlattice1[0]"] = 1;
  intparams["Nlattice1[1]"] = 1;
  intparams["Nlattice1[2]"] = 1;
  intparams["Nlattice2[0]"] = 1;
  intparams["Nlattice2[1]"] = 1;
  intparams["Nlattice2[2]"] = 1;
  floatparams["vfluid1[0]"] = 0.0;
  floatparams["vfluid1[1]"] = 0.0;
  floatparams["vfluid1[2]"] = 0.0;
  floatparams["vfluid2[0]"] = 0.0;
  floatparams["vfluid2[1]"] = 0.0;
  floatparams["vfluid2[2]"] = 0.0;
  floatparams["rhofluid1"] = 1.0;
  floatparams["rhofluid2"] = 1.0;
  floatparams["press1"] = 1.0;
  floatparams["press2"] = 1.0;
  floatparams["rexplosion"] = 0.2;
  floatparams["amp"] = 0.1;
  floatparams["lambda"] = 0.5;
  floatparams["kefrac"] = 0.0;
  floatparams["radius"] = 1.0;
  floatparams["angvel"] = 0.0;
  floatparams["mcloud"] = 1.0;
  floatparams["mplummer"] = 1.0;
  floatparams["rplummer"] = 1.0;
  floatparams["rstar"] = 0.1;
  floatparams["cdmfrac"] = 0.0;
  floatparams["gasfrac"] = 0.0;
  floatparams["starfrac"] = 1.0;
  floatparams["m1"] = 0.5;
  floatparams["m2"] = 0.5;
  floatparams["m3"] = 0.5;
  floatparams["m4"] = 0.5;
  floatparams["abin"] = 1.0;
  floatparams["abin2"] = 0.1;
  floatparams["ebin"] = 0.0;
  floatparams["ebin2"] = 0.0;
  floatparams["phirot"] = 0.0;
  floatparams["thetarot"] = 0.0;
  floatparams["psirot"] = 0.0;
  floatparams["vmachbin"] = 1.0;
  floatparams["alpha_turb"] = 0.1;
  floatparams["power_turb"] = -4.0;
  floatparams["asound"] = 1.0;
  floatparams["zmax"] = 1.0;
  floatparams["thermal_energy"] = 1.0;

  // Random number generator parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["rand_algorithm"] = "none";
  intparams["randseed"] = 0;

  // MPI parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["mpi_decomposition"] = "kdtree";
  intparams["pruning_level_min"] = 6;
  intparams["pruning_level_max"] = 6;

  // Python parameters
  //-----------------------------------------------------------------------------------------------
  floatparams["dt_python"] = 8.0;

  // Dust Parameters
  //-----------------------------------------------------------------------------------------------
  stringparams["dust_forces"] = "none" ;
  stringparams["drag_law"] = "none" ;
  floatparams["drag_coeff"] = 0 ;
  floatparams["dust_mass_factor"] = 1 ;


  return;
}



//=================================================================================================
//  Parameters::GetParameter
/// Gets (string) value for given parameter
//=================================================================================================
string Parameters::GetParameter
 (string key)                          ///< [in] Parameter key to be searched
{
  if (intparams.count(key) == 1) {
    std::stringstream value;
    value << intparams[key];
    return value.str();
  }
  else if (floatparams.count(key) == 1) {
    std::stringstream value;
    value << floatparams[key];
    return value.str();
  }
  else if (stringparams.count(key) == 1) {
    return stringparams[key];
  }
  else {
    string msg = "Error: the parameter " + key + " is not recognized ";
    ExceptionHandler::getIstance().raise(msg);
  }

  return NULL;
}



//=================================================================================================
//  Parameters::SetParameter
/// Set parameter value in memory.  Checks in turn if parameter is a
/// string, float or integer before recording value.
//=================================================================================================
void Parameters::SetParameter
 (string key,                          ///< [in] Parameter key to be searched
  string value)                        ///< [in] Parameter value if key found
{
  if (intparams.count(key) == 1) {
    std::stringstream(value) >> intparams[key];
  }
  else if (floatparams.count(key) == 1) {
    std::stringstream(value) >> floatparams[key];
  }
  else if (stringparams.count(key) == 1) {
    stringparams[key] = value;
  }
  else {
    cout << "Warning: parameter " << key << " was not recognized" << endl;
  }

  return;
}



//=================================================================================================
//  Parameters::PrintParameters
/// Prints all parameters stored in memory to screen.
//=================================================================================================
void Parameters::PrintParameters(void)
{
  debug1("[Parameters::PrintParameters]");

  // Print all integer parameters
  std::map <std::string, int>::iterator it;
  for (it=intparams.begin(); it != intparams.end(); ++it) {
    std::cout << it->first << " = " << it->second << std::endl;
  }

  // Print all float parameters
  std::map <std::string, double>::iterator it2;
  for (it2=floatparams.begin(); it2 != floatparams.end(); ++it2) {
    std::cout << it2->first << " = " << it2->second << std::endl;
  }

  // Print all string parameters
  std::map <std::string, std::string>::iterator it3;
  for (it3=stringparams.begin(); it3 != stringparams.end(); ++it3) {
    std::cout << it3->first << " = " << it3->second << std::endl;
  }

  return;
}



//=================================================================================================
//  Parameters::RecordParametersToFile
/// Writes all recorded parameters to file named 'run_id.param'
//=================================================================================================
void Parameters::RecordParametersToFile(void)
{
  string filename = stringparams["run_id"] + ".param";  // Output filename
  ofstream outfile;                                     // Output file stream

  debug1("[Parameters::RecordParametersToFile]");

  outfile.open(filename.c_str());

  // Write all integer parameters
  std::map <std::string, int>::iterator it;
  for (it=intparams.begin(); it != intparams.end(); ++it) {
    outfile << it->first << " = " << it->second << endl;
  }

  // Write all float parameters
  std::map <std::string, double>::iterator it2;
  for (it2=floatparams.begin(); it2 != floatparams.end(); ++it2) {
    outfile << it2->first << " = " << it2->second << endl;
  }

  // Write all string parameters
  std::map <std::string, std::string>::iterator it3;
  for (it3=stringparams.begin(); it3 != stringparams.end(); ++it3) {
    outfile << it3->first << " = " << it3->second << endl;
  }

  outfile.close();

  return;
}



//=================================================================================================
//  Parameters::TrimWhiteSpace
/// Trims string of all white space
//=================================================================================================
std::string Parameters::TrimWhiteSpace
 (std::string instr)                   ///< [in] Input string to be trimmed
{
  string outstr;                       // Final string without any whitespace

  // Loop over all characters and ignore any white-space characters
  for (unsigned int i=0; i<instr.length(); i++) {
    if (instr[i] != ' ') outstr += instr[i];
  }

  return outstr;
}



//=================================================================================================
//  Parameters::CheckInvalidParameters
/// Check if any features deemed unstable and/or buggy are switched on in the
/// parameters file.  If so, then throw an exception with an error message.
//=================================================================================================
void Parameters::CheckInvalidParameters(void)
{
  bool errorflag = false;              ///< Flag for any invlaid parameters

  debug2("[Parameters::CheckInvalidParameters]");

  // SM2012 SPH simulation specific errors
  //-----------------------------------------------------------------------------------------------
  if (stringparams["sim"] == "sm2012sph") {

    // Saitoh & Makino (2012) currently deactivated while development of MPI
    cout << "Saitoh & Makino (2012) SPH algorithm currently disabled" << endl;
    errorflag = true;

  }
  //-----------------------------------------------------------------------------------------------


  // If any single error has been registered, throw exception and terminate.
  if (errorflag) {
    string message = "Aborting simulation due to invalid parameters";
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}
