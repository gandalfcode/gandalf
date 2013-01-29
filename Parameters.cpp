// ============================================================================
// PARAMETERS.CPP
// ============================================================================


#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Parameters.h"
using namespace std;


// ============================================================================
// PARAMETERS::PARAMETERS
// ============================================================================
Parameters::Parameters()
{
  SetDefaultValues();
}


// ============================================================================
// PARAMETERS::~PARAMETERS
// ============================================================================
Parameters::~Parameters()
{
}



// ============================================================================
// PARAMETERS::PARSELINE
// ============================================================================
void Parameters::ReadParamsFile(std::string filename)
{
  ifstream inputfile;
  std::string line;

  inputfile.open(filename.c_str(), ios::in);
  if (inputfile.is_open()) {
    while ( inputfile.good() ) {
      getline(inputfile, line);
      ParseLine (line);
    }
  }
  else {
    cout << "The specified parameter file: " << filename 
	 << " does not exist, aborting" << endl;
    exit(-1);
  }
  inputfile.close();
  return;
}


// ============================================================================
// PARAMETERS::PARSELINE
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

  /*
  std::cout << "colon_pos : " << colon_pos << std::endl;
  std::cout << "equal_pos : " << equal_pos << std::endl;
  std::cout << "length    : " << length << std::endl;
  std::cout << "var_name  : " << var_name << std::endl;
  std::cout << "var_value : " << var_value << std::endl;
  */

  SetParameter(var_name,var_value);

  return;
}



// ============================================================================
// PARAMETERS::
// ============================================================================
void Parameters::SetDefaultValues(void)
{

  intparams["ndim"]=3;

  intparams["Npart"]=100;

  floatparams["h_fac"]=1.2;
  floatparams["h_converge"]=0.005;
  floatparams["gamma_eos"]=1.6666666666666;
  floatparams["temp0"]=1.0;
  floatparams["mu_bar"]=1.0;
  floatparams["accel_mult"]=0.3;
  floatparams["courant_mult"]=0.15;
  floatparams["alpha_visc"]=0.1;
  floatparams["beta_visc"]=0.2;

  floatparams["boxmin[0]"]=0.0;
  floatparams["boxmin[1]"]=0.0;
  floatparams["boxmin[2]"]=0.0;
  floatparams["boxmax[0]"]=0.0;
  floatparams["boxmax[1]"]=0.0;
  floatparams["boxmax[2]"]=0.0;

  stringparams["x_boundary_lhs"]="open";
  stringparams["x_boundary_rhs"]="open";
  stringparams["y_boundary_lhs"]="open";
  stringparams["y_boundary_rhs"]="open";
  stringparams["z_boundary_lhs"]="open";
  stringparams["z_boundary_rhs"]="open";

  stringparams["ic"]="random_cube";

  stringparams["sph"]="gradh";

  stringparams["neib_search"]="bruteforce";

  stringparams["kernel"]="m4";

  stringparams["sph_integration"]="lfkdk";

  stringparams["gas_eos"]="isothermal";

  stringparams["out_file_form"]="sf";

  stringparams["rinunit"]="";
  stringparams["minunit"]="";
  stringparams["tinunit"]="";

  stringparams["routunit"]="pc";
  stringparams["moutunit"]="m_sun";
  stringparams["toutunit"]="myr";

  return;
}



// ============================================================================
// PARAMETERS::
// ============================================================================
void Parameters::SetParameter(std::string key, std::string value)
{
  std::cout << "@" << key << "@" << std::endl;
  std::cout << "@" << value << "@" << std::endl;
  std::cout << intparams.count(key) << "  " << stringparams.count(key) << std::endl;
  if (intparams.count(key) == 1) std::stringstream(value) >> intparams[key];
  else if (floatparams.count(key) == 1) std::stringstream(value) >> floatparams[key];
  else if (stringparams.count(key) == 1) stringparams[key] = value;
  else cout << "Warning: parameter " << key << "was not recognized" << endl;

  return;
}



// ============================================================================
// PARAMETERS::
// ============================================================================
void Parameters::PrintParameters(void)
{

  std::map <std::string, int>::iterator it;
  for (it=intparams.begin(); it != intparams.end(); ++it){
    std::cout << it->first << " " << it->second << std::endl;
  }

  std::map <std::string, std::string>::iterator it2;
  for (it2=stringparams.begin(); it2 != stringparams.end(); ++it2){
    std::cout << it2->first << " " << it2->second << std::endl;
  }

}




void Parameters::trim2(std::string& str)
{
  std::string::size_type pos = str.find_last_not_of(' ');
  if(pos != std::string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != std::string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}
