//=============================================================================
// Parameters.cpp
// Contains all functions for calculating default values and reading in 
// new values from the simulation parameter file.
//=============================================================================


#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Exception.h"
#include "Parameters.h"
#include "Debug.h"
using namespace std;


//=============================================================================
//  Parameters::Parameters
/// Constructor for Parameters class
//=============================================================================
Parameters::Parameters()
{
  SetDefaultValues();
}



//=============================================================================
//  Parameters::~Parameters
/// Parameters destructor
//=============================================================================
Parameters::~Parameters()
{
}



Parameters::Parameters(const Parameters& other) {
  this->intparams = other.intparams;
  this->stringparams = other.stringparams;
  this->floatparams = other.floatparams;
}



//=============================================================================
//  Parameters::ReadParamsFile
/// Read and parse parameter file 'filename'.  If file doesn't exist, or 
/// file does not contain a simulation run id, then quit program here.
//=============================================================================
void Parameters::ReadParamsFile
(string filename                    ///< [in] Parameters file to be read
 )
{
  ifstream inputfile;               // Input file stream
  std::string line;                 // Parameter file line

  debug1("[Parameters::ReadParamsFile]");

  // Set-up all parameters and assign default values
  SetDefaultValues();

  // If parameter file can be opened, parse each line in turn.
  // Else, quit program with exception
  inputfile.open(filename.c_str(), ios::in);
  if (inputfile.is_open()) {
    while ( inputfile.good() ) {
      getline(inputfile, line);
      ParseLine (line);
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

  // Record parameters to file 'run_id.param'
  RecordParametersToFile();

  return;
}



//=============================================================================
//  Parameters::ParseLine
/// Parse a single line read from the parameters file.
/// Identifies if the line is in the form 'Comments : variable = value', or 
/// 'variable = value', and if so, stores value in memory.
/// If line begins with a hash charatcer '#', then ignore line as a comment.
///============================================================================
void Parameters::ParseLine
(string paramline                       ///< [in] Line from parameters file 
 )                                      ///< to be parsed.
{
  // First, trim all white space from line
  paramline = TrimWhiteSpace(paramline);

  int colon_pos = paramline.find(':');  // Position of colon in string
  int equal_pos = paramline.find('=');  // Position of equals in string
  int hash_pos = paramline.find('#');   // Position of hash in string
  int length = paramline.length();      // Length of string

  // Ignore line if it is a comment (i.e. begins with a has character)
  if (hash_pos == 0) return;

  // If line is not in the correct format (either equals is not present, or 
  // equals is before the colon) then skip line and return
  if (equal_pos == std::string::npos || colon_pos >= equal_pos) return;

  // Extract variable name and value from line
  std::string var_name = paramline.substr(colon_pos+1,equal_pos-colon_pos-1);
  std::string var_value = paramline.substr(equal_pos+1,length-equal_pos-1);

  // Finally, set parameter value in memory
  SetParameter(var_name,var_value);

  return;
}



//=============================================================================
//  Parameters::SetDefaultValues
/// Record all parameter variable names in memory also setting default values.
//=============================================================================
void Parameters::SetDefaultValues(void)
{
  debug1("[Parameters::SetDefaultValues]");

  // Simulation id, filename and output time parameters
  // --------------------------------------------------------------------------
  stringparams["run_id"] = "";
  stringparams["in_file_form"] = "ascii";
  stringparams["out_file_form"] = "ascii";
  floatparams["tend"] = 1.0;
  floatparams["dt_snap"] = 0.1;
  intparams["Nstepsmax"] = 9999999;
  intparams["noutputstep"] = 32;

  // Initial conditions parameters
  // --------------------------------------------------------------------------
  stringparams["ic"] = "random_cube";
  stringparams["particle_distribution"] = "random";
  intparams["smooth_ic"] = 0;
  intparams["Npart"] = 256;
  intparams["ndim"] = 3;
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
  floatparams["amp"] = 0.1;
  floatparams["lambda"] = 0.5;
  floatparams["kefrac"] = 0.0;
  floatparams["radius"] = 1.0;
  floatparams["angvel"] = 0.0;
  floatparams["mcloud"] = 1.0;

  // Integration scheme and timestep parameters
  // --------------------------------------------------------------------------
  stringparams["sph_integration"] = "lfkdk";
  floatparams["accel_mult"] = 0.3;
  floatparams["courant_mult"] = 0.15;
  floatparams["nbody_mult"] = 0.075;
  intparams["Nlevels"] = 1;
  intparams["sph_single_timestep"] = 0;
  intparams["nbody_single_timestep"] = 0;

  // SPH parameters
  // --------------------------------------------------------------------------
  stringparams["sph"] = "gradh";
  stringparams["kernel"] = "m4";
  intparams["tabulated_kernel"] = 0;
  stringparams["neib_search"] = "grid";
  floatparams["h_fac"] = 1.2;
  floatparams["h_converge"] = 0.01;

  // Thermal physics parameters
  // --------------------------------------------------------------------------
  intparams["hydro_forces"] = 1;
  stringparams["gas_eos"] = "energy_eqn";
  stringparams["energy_integration"] = "PEC";
  floatparams["energy_mult"] = 0.4;
  floatparams["gamma_eos"] = 1.66666666666666;
  floatparams["temp0"] = 1.0;
  floatparams["mu_bar"] = 1.0;
  floatparams["rho_bary"] = 1.0e-14;
  floatparams["eta_eos"] = 1.4;

  // Artificial viscosity parameters
  // --------------------------------------------------------------------------
  stringparams["avisc"] = "mon97";
  stringparams["acond"] = "none";
  floatparams["alpha_visc"] = 1.0;
  floatparams["beta_visc"] = 2.0;

  // Riemann solver parameters
  // --------------------------------------------------------------------------
  stringparams["riemann_solver"] = "exact";
  stringparams["slope_limiter"] = "mine";
  intparams["riemann_order"] = 1;

  // Gravity parameters
  // --------------------------------------------------------------------------
  intparams["self_gravity"] = 0;
  stringparams["grav_kernel"] = "mean_h";

  // N-body parameters
  // --------------------------------------------------------------------------
  stringparams["nbody"] = "lfkdk";
  intparams["nbody_softening"] = 0;
  intparams["sub_systems"] = 0;

  // Boundary conditions parameters
  // --------------------------------------------------------------------------
  stringparams["x_boundary_lhs"] = "open";
  stringparams["x_boundary_rhs"] = "open";
  stringparams["y_boundary_lhs"] = "open";
  stringparams["y_boundary_rhs"] = "open";
  stringparams["z_boundary_lhs"] = "open";
  stringparams["z_boundary_rhs"] = "open";
  floatparams["boxmin[0]"] = 0.0;
  floatparams["boxmin[1]"] = 0.0;
  floatparams["boxmin[2]"] = 0.0;
  floatparams["boxmax[0]"] = 1.0;
  floatparams["boxmax[1]"] = 1.0;
  floatparams["boxmax[2]"] = 1.0;

  // Unit and scaling parameters
  // --------------------------------------------------------------------------
  intparams["dimensionless"] = 0;
  stringparams["rinunit"] = "";
  stringparams["minunit"] = "";
  stringparams["tinunit"] = "";
  stringparams["vinunit"] = "";
  stringparams["ainunit"] = "";
  stringparams["rhoinunit"] = "";
  stringparams["Einunit"] = "";
  stringparams["mominunit"] = "";
  stringparams["angmominunit"] = "";
  stringparams["angvelinunit"] = "";
  stringparams["uinunit"] = "";
  stringparams["dudtinunit"] = "";
  stringparams["tempinunit"] = "";
  stringparams["routunit"] = "pc";
  stringparams["moutunit"] = "m_sun";
  stringparams["toutunit"] = "myr";
  stringparams["voutunit"] = "km_s";
  stringparams["aoutunit"] = "km_s2";
  stringparams["rhooutunit"] = "g_cm3";
  stringparams["Eoutunit"] = "J";
  stringparams["momoutunit"] = "m_sunkm_s";
  stringparams["angmomoutunit"] = "m_sunkm2_s";
  stringparams["angveloutunit"] = "rad_s";
  stringparams["uoutunit"] = "J_kg";
  stringparams["dudtoutunit"] = "J_kg_s";
  stringparams["tempoutunit"] = "K";

  return;
}



//=============================================================================
//  Parameters::SetParameter
/// Set parameter value in memory.  Checks in turn if parameter is a 
/// string, float or integer before recording value.
//=============================================================================
string Parameters::GetParameter (string key)
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
  else if (stringparams.count(key) == 1)
    return stringparams[key];
  else {
    string msg = "Error: the parameter " + key + " is not recognized ";
    ExceptionHandler::getIstance().raise(msg);
  }

  return NULL;
}



//=============================================================================
//  Parameters::SetParameter
/// Set parameter value in memory.  Checks in turn if parameter is a
/// string, float or integer before recording value.
//=============================================================================
void Parameters::SetParameter
(string key,                        ///< [in] Parameter key to be searched
 string value)                      ///< [in] Parameter value if key found
{
  if (intparams.count(key) == 1)
    std::stringstream(value) >> intparams[key];
  else if (floatparams.count(key) == 1)
    std::stringstream(value) >> floatparams[key];
  else if (stringparams.count(key) == 1)
    stringparams[key] = value;
  else cout << "Warning: parameter " << key << "was not recognized" << endl;

  return;
}



//=============================================================================
//  Parameters::PrintParameters
/// Prints all parameters stored in memory to screen.
//=============================================================================
void Parameters::PrintParameters(void)
{
  debug1("[Parameters::PrintParameters]");

  // Print all integer parameters
  std::map <std::string, int>::iterator it;
  for (it=intparams.begin(); it != intparams.end(); ++it) {
    std::cout << it->first << " = " << it->second << std::endl;
  }

  // Print all float parameters
  std::map <std::string, float>::iterator it2;
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



//=============================================================================
//  Parameters::RecordParametersToFile
/// Writes all recorded parameters to file named 'run_id.param'
//=============================================================================
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
  std::map <std::string, float>::iterator it2;
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



//=============================================================================
//  Parameters::TrimWhiteSpace
/// Trims string of all white space
//=============================================================================
std::string Parameters::TrimWhiteSpace(std::string instr)
{
  int i;                             // Character counter
  string outstr;                     // Resultant string without whitespace

  // Loop over all characters and ignore any white-space characters
  for (i=0; i < instr.length(); i++) {
    if (instr[i] != ' ') outstr += instr[i];
  }

  return outstr;
}
