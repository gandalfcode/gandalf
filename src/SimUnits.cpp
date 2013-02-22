// ============================================================================
// SimUnits.cpp
// ============================================================================


#include <math.h>
#include <map>
#include <string>
#include "Constants.h"
#include "Exception.h"
#include "SimUnits.h"
using namespace std;



// ============================================================================
// SimUnit::SimUnit
// ============================================================================
SimUnit::SimUnit()
{
  inscale = 1.0;
  inSI = 1.0;
  outcgs = 1.0;
  outscale = 1.0;
  outSI = 1.0;
}



// ============================================================================
// LengthUnit::OutputScale
// ============================================================================
double SimUnit::OutputScale(string unit_string)
{
  return inscale*inSI/SIUnit(unit_string);
}



// ============================================================================
// LengthUnit::SIUnit
// ============================================================================
double LengthUnit::SIUnit(string unit)
{
  if (unit == "mpc") return 1.0E6*r_pc;
  else if (unit == "kpc") return 1.0E3*r_pc;
  else if (unit == "pc") return r_pc;
  else if (unit == "au") return r_au;
  else if (unit == "r_sun") return r_sun;
  else if (unit == "r_earth") return r_earth;
  else if (unit == "km") return 1000.0;
  else if (unit == "m") return 1.0;
  else if (unit == "cm") return 0.01;
  else if (unit == "") return 1.0;
  else {
    string message = "Parameter error : Unrecognised unit = " + unit;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// MassUnit::SIUnit
// ============================================================================
double MassUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sun") return m_sun;
  else if (unit_string == "m_jup") return m_jup;
  else if (unit_string == "m_earth") return m_earth;
  else if (unit_string == "kg") return 1.0;
  else if (unit_string == "g") return 1.0e-3;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// TimeUnit::SIUnit
// ============================================================================
double TimeUnit::SIUnit(string unit_string)
{
  if (unit_string == "gyr") return 1000.0*myr;
  else if (unit_string == "myr") return myr;
  else if (unit_string == "yr") return yr;
  else if (unit_string == "day") return day;
  else if (unit_string == "s") return 1.0;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// VelocityUnit::SIUnit
// ============================================================================
double VelocityUnit::SIUnit(string unit_string)
{
  if (unit_string == "km_s") return 1000.0;
  else if (unit_string == "au_yr") return r_au/yr;
  else if (unit_string == "m_s") return 1.0;
  else if (unit_string == "cm_s") return 0.01;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}


// ============================================================================
// AccelerationUnit::SIUnit
// ============================================================================
double AccelerationUnit::SIUnit(string unit_string)
{
  if (unit_string == "km_s2") return 1000.0;
  else if (unit_string == "au_yr2") return r_au/yr/yr;
  else if (unit_string == "m_s2") return 1.0;
  else if (unit_string == "cm_s2") return 0.01;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// DensityUnit::SIUnit
// ============================================================================
double DensityUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sun_pc3") return r_pc*r_pc*r_pc;
  else if (unit_string == "kg_m3") return 1.0;
  else if (unit_string == "g_cm3") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// EnergyUnit::SIUnit
// ============================================================================
double EnergyUnit::SIUnit(string unit_string)
{
  if (unit_string == "J") return 1.0;
  else if (unit_string == "erg") return 1.0e-7;
  else if (unit_string == "GJ") return 1.0e12;
  else if (unit_string == "10^40erg") return 1.0e33;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// MomentumUnit::SIUnit
// ============================================================================
double MomentumUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sunkm_s") return m_sun*1000.0;
  else if (unit_string == "m_sunau_yr") return m_sun*r_au/yr;
  else if (unit_string == "kgm_s") return 1.0;
  else if (unit_string == "gcm_s") return 1.0e-5;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// AngularMomentumUnit::SIUnit
// ============================================================================
double AngularMomentumUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sunkm2_s") return m_sun*1000.0*1000.0;
  else if (unit_string == "m_sunau2_yr") return m_sun*r_au*r_au/yr;
  else if (unit_string == "kgm2_s") return 1.0;
  else if (unit_string == "gcm2_s") return 1.0e-7;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// AngularVelocityUnit::SIUnit
// ============================================================================
double AngularVelocityUnit::SIUnit(string unit_string)
{
  if (unit_string == "rad_s") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// SpecificEnergyUnit::SIUnit
// ============================================================================
double SpecificEnergyUnit::SIUnit(string unit_string)
{
  if (unit_string == "J_kg") return 1.0;
  else if (unit_string == "erg_g") return 1.0e-4;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// SpecificEnergyRateUnit::SIUnit
// ============================================================================
double SpecificEnergyRateUnit::SIUnit(string unit_string)
{
  if (unit_string == "J_kg_s") return 1.0;
  else if (unit_string == "erg_g_s") return 1.0e-4;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



// ============================================================================
// TemperatureUnit::SIUnit
// ============================================================================
double TemperatureUnit::SIUnit(string unit_string)
{
  if (unit_string == "K") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}




// ============================================================================
// SimUnits::SimUnits
// ============================================================================
SimUnits::SimUnits()
{
  ReadInputUnits = false;
}



// ============================================================================
// SimUnits::~SimUnits
// ============================================================================
SimUnits::~SimUnits()
{
}



// ============================================================================
// Units::SetupUnits
// ============================================================================
void SimUnits::SetupUnits(Parameters &params)
{

  if (!ReadInputUnits) {
    params.stringparams["rinunit"] = params.stringparams["routunit"];
    params.stringparams["minunit"] = params.stringparams["moutunit"];
    params.stringparams["tinunit"] = params.stringparams["toutunit"];
    params.stringparams["vinunit"] = params.stringparams["voutunit"];
    params.stringparams["ainunit"] = params.stringparams["aoutunit"];
    params.stringparams["rhoinunit"] = params.stringparams["rhooutunit"];
    params.stringparams["sigmainunit"] = params.stringparams["sigmaoutunit"];
    params.stringparams["pressinunit"] = params.stringparams["pressoutunit"];
    params.stringparams["finunit"] = params.stringparams["foutunit"];
    params.stringparams["Einunit"] = params.stringparams["Eoutunit"];
    params.stringparams["mominunit"] = params.stringparams["momoutunit"];
    params.stringparams["angmominunit"] = params.stringparams["angmomoutunit"];
    params.stringparams["angvelinunit"] = params.stringparams["angveloutunit"];
    params.stringparams["dmdtinunit"] = params.stringparams["dmdtoutunit"];
    params.stringparams["Linunit"] = params.stringparams["Loutunit"];
    params.stringparams["kappainunit"] = params.stringparams["kappaoutunit"];
    params.stringparams["Binunit"] = params.stringparams["Boutunit"];
    params.stringparams["Qinunit"] = params.stringparams["Qoutunit"];
    params.stringparams["Jcurinunit"] = params.stringparams["Jcuroutunit"];
    params.stringparams["uinunit"] = params.stringparams["uoutunit"];
    params.stringparams["dudtinunit"] = params.stringparams["dudtoutunit"];
    params.stringparams["tempinunit"] = params.stringparams["tempoutunit"];
  }

  // Length units
  // --------------------------------------------------------------------------
  r.inSI = r.SIUnit(params.stringparams["rinunit"]);
  r.outSI = r.SIUnit(params.stringparams["routunit"]);
  r.inscale = 1.0;
  r.outscale = r.inscale*r.inSI/r.outSI;
  r.outcgs = 100.0*r.outSI;

  // Mass units
  // --------------------------------------------------------------------------
  m.inSI = m.SIUnit(params.stringparams["minunit"]);
  m.outSI = m.SIUnit(params.stringparams["moutunit"]);
  m.inscale = 1.0;
  m.outscale = m.inscale*m.inSI/m.outSI;
  m.outcgs = 1000.0*m.outSI;

  // Time units
  // --------------------------------------------------------------------------
  t.inSI = t.SIUnit(params.stringparams["tinunit"]);
  t.outSI = t.SIUnit(params.stringparams["toutunit"]);
  t.inscale = pow(r.inscale*r.inSI,1.5)/sqrt(m.inscale*m.inSI*G_const);
  t.inscale /= t.inSI;
  t.outscale = pow(r.outscale*r.outSI,1.5)/sqrt(m.outscale*m.outSI*G_const);
  t.outscale /= t.outSI;
  t.outcgs = t.outSI;

  // Velocity units
  // --------------------------------------------------------------------------
  v.inSI = v.SIUnit(params.stringparams["vinunit"]);
  v.outSI = v.SIUnit(params.stringparams["voutunit"]);
  v.inscale = r.inscale*r.inSI/(t.inscale*t.inSI);
  v.inscale /= v.inSI;
  v.outscale = r.outscale*r.outSI/(t.outscale*t.outSI);
  v.outscale /= v.outSI;
  v.outcgs = 100.0*v.outSI;

  // Acceleration units
  // --------------------------------------------------------------------------
  a.inSI = a.SIUnit(params.stringparams["ainunit"]);
  a.outSI = a.SIUnit(params.stringparams["aoutunit"]);
  a.inscale = (r.inscale*r.inSI)/(t.inscale*t.inSI*t.inscale*t.inSI);
  a.inscale = a.inscale / a.inSI;
  a.outscale = (r.outscale*r.outSI)/(t.outscale*t.outSI*t.outscale*t.outSI);
  a.outscale = a.outscale / a.outSI;
  a.outcgs = 100.0*a.outSI;

  // Density units
  // --------------------------------------------------------------------------
  rho.inSI = rho.SIUnit(params.stringparams["rhoinunit"]);
  rho.outSI = rho.SIUnit(params.stringparams["rhooutunit"]);
  rho.inscale = (m.inscale*m.inSI) / pow(r.inscale*r.inSI,3);
  rho.inscale = rho.inscale / rho.inSI;
  rho.outscale = (m.outscale*m.outSI) / pow(r.outscale*r.outSI,3);
  rho.outscale = rho.outscale / rho.outSI;
  rho.outcgs = 1.0e-3*rho.outSI;

  // Energy units
  // --------------------------------------------------------------------------
  E.inSI = E.SIUnit(params.stringparams["Einunit"]);
  E.outSI = E.SIUnit(params.stringparams["Eoutunit"]);
  E.inscale = m.inscale*m.inSI*pow(r.inscale*r.inSI,2)/pow(t.inscale*t.inSI,2);
  E.inscale = E.inscale / E.inSI;
  E.outscale = m.outscale*m.outSI*pow(r.outscale*r.outSI,2)/
    pow(t.outscale*t.outSI,2);
  E.outscale = E.outscale / E.outSI;
  E.outcgs = 1.0e7*E.outSI;

  // Momentum units
  // --------------------------------------------------------------------------
  mom.inSI = mom.SIUnit(params.stringparams["mominunit"]);
  mom.outSI = mom.SIUnit(params.stringparams["momoutunit"]);
  mom.inscale = m.inscale*m.inSI*r.inscale*r.inSI/(t.inscale*t.inSI);
  mom.inscale = mom.inscale / mom.inSI;
  mom.outscale = m.outscale*m.outSI*r.outscale*r.outSI/(t.outscale*t.outSI);
  mom.outscale = mom.outscale / mom.outSI;
  mom.outcgs = 1.0e5*mom.outSI;

  // Angular momentum units
  // --------------------------------------------------------------------------
  angmom.inSI = angmom.SIUnit(params.stringparams["angmominunit"]);
  angmom.outSI = angmom.SIUnit(params.stringparams["angmomoutunit"]);
  angmom.inscale = m.inscale*m.inSI*pow(r.inscale*r.inSI,2)/(t.inscale*t.inSI);
  angmom.inscale = angmom.inscale / angmom.inSI;
  angmom.outscale = m.outscale*m.outSI*pow(r.outscale*r.outSI,2)/
    (t.outscale*t.outSI);
  angmom.outscale = angmom.outscale / angmom.outSI;
  angmom.outcgs = 1.0e7*angmom.outSI;

  // Angular velocity units
  // --------------------------------------------------------------------------
  angvel.inSI = angvel.SIUnit(params.stringparams["angvelinunit"]);
  angvel.outSI = angvel.SIUnit(params.stringparams["angveloutunit"]);
  angvel.inscale = 1.0/(t.inscale*t.inSI);
  angvel.inscale = angvel.inscale / angvel.inSI;
  angvel.outscale = 1.0/(t.outscale*t.outSI);
  angvel.outscale = angvel.outscale / angvel.outSI;
  angvel.outcgs = angvel.outSI;

  // Specific internal energy units
  // --------------------------------------------------------------------------
  u.inSI = u.SIUnit(params.stringparams["uinunit"]);
  u.outSI = u.SIUnit(params.stringparams["uoutunit"]);
  u.inscale = pow(r.inscale*r.inSI,2)/pow(t.inscale*t.inSI,2);
  u.inscale = u.inscale / u.inSI;
  u.outscale = pow(r.outscale*r.outSI,2)/pow(t.outscale*t.outSI,2);
  u.outscale = u.outscale / u.outSI;
  u.outcgs = 1.0e4*u.outSI;

  // Rate of change of specific internal energy units
  // --------------------------------------------------------------------------
  dudt.inSI = dudt.SIUnit(params.stringparams["dudtinunit"]);
  dudt.outSI = dudt.SIUnit(params.stringparams["dudtoutunit"]);
  dudt.inscale = pow(r.inscale*r.inSI,2)/pow(t.inscale*t.inSI,3);
  dudt.inscale = dudt.inscale / dudt.inSI;
  dudt.outscale = pow(r.outscale*r.outSI,2)/pow(t.outscale*t.outSI,3);
  dudt.outscale = dudt.outscale / dudt.outSI;
  dudt.outcgs = 1.0e4*dudt.outSI;

  // Temperature units
  // --------------------------------------------------------------------------
  temp.inSI = temp.SIUnit(params.stringparams["tempinunit"]);
  temp.outSI = temp.SIUnit(params.stringparams["tempoutunit"]);
  temp.inscale = (m_hydrogen*u.inscale*u.inSI)/k_boltzmann;
  temp.inscale = temp.inscale / temp.inSI;
  temp.outscale = (m_hydrogen*u.outscale*u.outSI)/k_boltzmann;
  temp.outscale = temp.outscale / temp.outSI;
  temp.outcgs = temp.outSI;

  cout << "r.inscale : " << r.inscale 
       << "    r.outscale : " << r.outscale << endl;
  cout << "m.inscale : " << m.inscale 
       << "    m.outscale : " << m.outscale << endl;
  cout << "t.inscale : " << t.inscale 
       << "    t.outscale : " << t.outscale << endl;

  return;
}
