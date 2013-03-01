// ============================================================================
// Parameters.cpp
// ============================================================================


#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Exception.h"
#include "Parameters.h"
#include "Debug.h"
using namespace std;


// ============================================================================
// Parameters::Parameters
// ============================================================================
Parameters::Parameters()
{
  SetDefaultValues();
}



// ============================================================================
// Parameters::~Parameters
// ============================================================================
Parameters::~Parameters()
{
}



// ============================================================================
// Parameters::ReadParamsFile
// ============================================================================
void Parameters::ReadParamsFile(std::string filename)
{
  ifstream inputfile;
  std::string line;

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

  // Now verify that parameters file contains a (unique) run id.
  // If not defined, then quit program with exception
  if (stringparams["run_id"] == "") {
	string message = "The parameter file: " + filename +
	  " does not contain a run id string, aborting";
	ExceptionHandler::getIstance().raise(message);
  }

  // Record parameters to file
  RecordParametersToFile();

  return;
}



// ============================================================================
// Parameters::ParseLine
// ============================================================================
void Parameters::ParseLine(std::string paramline)
{
  int colon_pos;
  int equal_pos;
  int length;

  colon_pos = paramline.find(':');
  equal_pos = paramline.find('=');
  length = paramline.length();

  if (colon_pos == std::string::npos || equal_pos == std::string::npos || 
      colon_pos >= equal_pos) return;

  std::string var_name = paramline.substr(colon_pos+1,equal_pos-colon_pos-1);
  std::string var_value = paramline.substr(equal_pos+1,length-equal_pos-1);

  //std::remove(var_name.begin(), var_name.end(), ' ');
  //std::remove(var_value.begin(), var_value.end(), ' ');
  trim2(var_name);
  trim2(var_value);

  SetParameter(var_name,var_value);

  return;
}



// ============================================================================
// Parameters::SetDefaultValues
// ============================================================================
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
  intparams["Npart"] = 100;
#if defined(FIXED_DIMENSIONS)
  intparams["ndim"] = NDIM;
#else
  intparams["ndim"] = 3;
#endif
  intparams["Nlattice1[0]"] = 16;
  intparams["Nlattice1[1]"] = 16;
  intparams["Nlattice1[2]"] = 16;
  intparams["Nlattice2[0]"] = 16;
  intparams["Nlattice2[1]"] = 16;
  intparams["Nlattice2[2]"] = 16;
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

  // Integration scheme and timestep parameters
  // --------------------------------------------------------------------------
  stringparams["sph_integration"] = "lfkdk";
  floatparams["accel_mult"] = 0.3;
  floatparams["courant_mult"] = 0.15;
  intparams["Nlevels"] = 1;
  intparams["sph_single_timestep"] = 0;
  intparams["nbody_single_timestep"] = 0;

  // SPH parameters
  // --------------------------------------------------------------------------
  stringparams["sph"] = "gradh";
  stringparams["kernel"] = "m4";
  stringparams["neib_search"] = "bruteforce";
  floatparams["h_fac"] = 1.2;
  floatparams["h_converge"] = 0.005;

  // Artificial viscosity parameters
  // --------------------------------------------------------------------------
  stringparams["avisc"] = "mon97";
  stringparams["acond"] = "none";
  floatparams["alpha_visc"] = 1.0;
  floatparams["beta_visc"] = 2.0;

  // Thermal physics parameters
  // --------------------------------------------------------------------------
  intparams["hydro_forces"] = 1;
  stringparams["gas_eos"] = "isothermal";
  stringparams["energy_integration"] = "PEC";
  floatparams["energy_mult"] = 0.2;
  floatparams["gamma_eos"] = 1.6666666666666;
  floatparams["temp0"] = 1.0;
  floatparams["mu_bar"] = 1.0;

  // Gravity parameters
  // --------------------------------------------------------------------------
  intparams["self_gravity"] = 0;
  stringparams["grav_kernel"] = "mean_h";

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
  floatparams["boxmax[0]"] = 0.0;
  floatparams["boxmax[1]"] = 0.0;
  floatparams["boxmax[2]"] = 0.0;

  // Unit and scaling parameters
  // --------------------------------------------------------------------------
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

  PrintParameters();

  return;
}



// ============================================================================
// Parameters::SetParameter
// ============================================================================
void Parameters::SetParameter(std::string key, std::string value)
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



// ============================================================================
// Parameters::PrintParameters
// ============================================================================
void Parameters::PrintParameters(void)
{
  debug1("[Parameters::PrintParameters]");

  std::map <std::string, int>::iterator it;
  for (it=intparams.begin(); it != intparams.end(); ++it) {
    std::cout << it->first << " = " << it->second << std::endl;
  }

  std::map <std::string, float>::iterator it2;
  for (it2=floatparams.begin(); it2 != floatparams.end(); ++it2) {
    std::cout << it2->first << " = " << it2->second << std::endl;
  }

  std::map <std::string, std::string>::iterator it3;
  for (it3=stringparams.begin(); it3 != stringparams.end(); ++it3) {
    std::cout << it3->first << " = " << it3->second << std::endl;
  }

}



// ============================================================================
// Parameters::RecordParametersToFile
// ============================================================================
void Parameters::RecordParametersToFile(void)
{
  string filename = stringparams["run_id"] + ".param";
  ofstream outfile;

  debug1("[Parameters::RecordParametersToFile]");

  outfile.open(filename.c_str());

  std::map <std::string, int>::iterator it;
  for (it=intparams.begin(); it != intparams.end(); ++it) {
    outfile << it->first << " = " << it->second << endl;
  }

  std::map <std::string, float>::iterator it2;
  for (it2=floatparams.begin(); it2 != floatparams.end(); ++it2) {
    outfile << it2->first << " = " << it2->second << endl;
  }

  std::map <std::string, std::string>::iterator it3;
  for (it3=stringparams.begin(); it3 != stringparams.end(); ++it3) {
    outfile << it3->first << " = " << it3->second << endl;
  }

  outfile.close();

}



// ============================================================================
// Parameters::trim2
// ============================================================================
void Parameters::trim2(std::string& str)
{
  std::string::size_type pos = str.find_last_not_of(' ');
  if (pos != std::string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if (pos != std::string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}
