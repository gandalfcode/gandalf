//=============================================================================
//  SimUnits.cpp
//  Contains functions for computing all scaling factors for converting 
//  between dimensionless and physical units.
//=============================================================================


#include <math.h>
#include <map>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Debug.h"
#include "Exception.h"
#include "SimUnits.h"
using namespace std;



//=============================================================================
//  SimUnit::SimUnit
/// Constructor for SimUnit class.  Defaults all scaling variables to unity 
/// in case of a dimensionless simulation.
//=============================================================================
SimUnit::SimUnit()
{
  inscale = 1.0;
  inSI = 1.0;
  outcgs = 1.0;
  outscale = 1.0;
  outSI = 1.0;
}



//=============================================================================
//  SimUnit::OutputScale
/// Compute scaling factor to convert between dimensionless units and the 
/// requested unit (unit_string).
//=============================================================================
double SimUnit::OutputScale(string unit_string)
{
  return inscale*inSI/SIUnit(unit_string);
}



//=============================================================================
//  LengthUnit::SIUnit
/// Return numerical value requested length unit in SI units
//=============================================================================
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



//=============================================================================
//  LengthUnit::LatexLabel
/// Return latex string of requested length unit for external plotting
//=============================================================================
string LengthUnit::LatexLabel(string unit)
{
  if (unit == "mpc") return "Mpc";
  else if (unit == "kpc") return "kpc";
  else if (unit == "pc") return "pc";
  else if (unit == "au") return "AU";
  else if (unit == "r_sun") return "R_{\\\\odot}";
  else if (unit == "r_earth") return "R_{\\\\oplus}";
  else if (unit == "km") return "km";
  else if (unit == "m") return "m";
  else if (unit == "cm") return "cm";
  else return "";
}



//=============================================================================
//  MassUnit::SIUnit
/// Return numerical value requested mass unit in SI units
//=============================================================================
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



//=============================================================================
//  MassUnit::LatexLabel
/// Return latex string of requested mass unit for external plotting
//=============================================================================
string MassUnit::LatexLabel(string unit_string)
{
  if (unit_string == "m_sun") return "M_{\\odot}";
  else if (unit_string == "m_jup") return "M_{J}";
  else if (unit_string == "m_earth") return "M_{\\oplus}";
  else if (unit_string == "kg") return "kg";
  else if (unit_string == "g") return "g";
  else return "";
}



//=============================================================================
//  TimeUnit::SIUnit
/// Return numerical value requested time unit in SI units
//=============================================================================
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



//=============================================================================
//  TimeUnit::LatexLabel
/// Return latex string of requested time unit for external plotting
//=============================================================================
string TimeUnit::LatexLabel(string unit_string)
{
  if (unit_string == "gyr") return "Gyr";
  else if (unit_string == "myr") return "Myr";
  else if (unit_string == "yr") return "yr";
  else if (unit_string == "day") return "days";
  else if (unit_string == "s") return "s";
  else return "";
}



//=============================================================================
//  VelocityUnit::SIUnit
/// Return numerical value requested velocity unit in SI units
//=============================================================================
double VelocityUnit::SIUnit(string unit_string)
{
  if (unit_string == "km_s") return 1000.0;
  else if (unit_string == "au_yr") return r_au/yr;
  else if (unit_string == "m_s") return 1.0;
  else if (unit_string == "cm_s") return 0.01;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  VelocityUnit::LatexLabel
/// Return latex string of requested velocity unit for external plotting
//=============================================================================
string VelocityUnit::LatexLabel(string unit_string)
{
  if (unit_string == "km_s") return "km\\,s^{-1}";
  else if (unit_string == "au_yr") return "AU\\,yr^{-1}";
  else if (unit_string == "m_s") return "m/,s^{-1}";
  else if (unit_string == "cm_s") return "cm\\,s^{-1}";
  else return "";
}



//=============================================================================
//  AccelerationUnit::SIUnit
/// Return numerical value requested acceleration unit in SI units
//=============================================================================
double AccelerationUnit::SIUnit(string unit_string)
{
  if (unit_string == "km_s2") return 1000.0;
  else if (unit_string == "au_yr2") return r_au/yr/yr;
  else if (unit_string == "m_s2") return 1.0;
  else if (unit_string == "cm_s2") return 0.01;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  AccelerationUnit::LatexLabel
/// Return latex string of requested acceleration unit for external plotting
//=============================================================================
string AccelerationUnit::LatexLabel(string unit_string)
{
  if (unit_string == "km_s2") return "km\\,s^{-2}";
  else if (unit_string == "au_yr2") return "AU\\,yr^{-2}";
  else if (unit_string == "m_s2") return "m\\,s^{-2}";
  else if (unit_string == "cm_s2") return "cm\\,s^{-2}";
  else return "";
}



//=============================================================================
//  DensityUnit::SIUnit
/// Return numerical value requested density unit in SI units
//=============================================================================
double DensityUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sun_pc3") return m_sun/(r_pc*r_pc*r_pc);
  else if (unit_string == "kg_m3") return 1.0;
  else if (unit_string == "g_cm3") return 1000.0;
  else if (unit_string == "") return 1.0;
  else {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  DensityUnit::LatexLabel
/// Return latex string of requested density unit for external plotting
//=============================================================================
string DensityUnit::LatexLabel(string unit_string)
{
  if (unit_string == "m_sun_pc3") return "M_{\\odot}\\,pc^{-3}";
  else if (unit_string == "kg_m3") return "kg\\,m^{-3}";
  else if (unit_string == "g_cm3") return "g\\,cm^{-3}";
  else return "";
}



//=============================================================================
//  PressureUnit::SIUnit
/// Return numerical value requested pressure unit in SI units
//=============================================================================
double PressureUnit::SIUnit(string unit_string)
{
  if (unit_string == "Pa") return 1.0;
  else if (unit_string == "bar") return 1.0e5;
  else if (unit_string == "") return 1.0;
  else {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  PressureUnit::LatexLabel
/// Return latex string of requested pressure unit for external plotting
//=============================================================================
string PressureUnit::LatexLabel(string unit_string)
{
  if (unit_string == "Pa") return "Pa";
  else if (unit_string == "bar") return "bar";
  else return "";
}



//=============================================================================
//  ForceUnit::SIUnit
/// Return numerical value of requested force unit in SI units
//=============================================================================
double ForceUnit::SIUnit(string unit_string)
{
  if (unit_string == "N") return 1.0;
  else if (unit_string == "dyn") return 1.0e-5;
  else if (unit_string == "") return 1.0;
  else {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  ForceUnit::LatexLabel
/// Return latex string of requested force unit for external plotting
//=============================================================================
string ForceUnit::LatexLabel(string unit_string)
{
  if (unit_string == "N") return "N";
  else if (unit_string == "dyn") return "dyn";
  else return "";
}



//=============================================================================
//  EnergyUnit::SIUnit
/// Return numerical value requested energy unit in SI units
//=============================================================================
double EnergyUnit::SIUnit(string unit_string)
{
  if (unit_string == "J") return 1.0;
  else if (unit_string == "erg") return 1.0e-7;
  else if (unit_string == "GJ") return 1.0e12;
  else if (unit_string == "10^40erg") return 1.0e33;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  EnergyUnit::LatexLabel
/// Return latex string of requested energy unit for external plotting
//=============================================================================
string EnergyUnit::LatexLabel(string unit_string)
{
  if (unit_string == "J") return "J";
  else if (unit_string == "erg") return "erg";
  else if (unit_string == "GJ") return "GJ";
  else if (unit_string == "10^40erg") return "10^{40}\\,erg";
  else return "";
}



//=============================================================================
//  MomentumUnit::SIUnit
/// Return numerical value requested momentum unit in SI units
//=============================================================================
double MomentumUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sunkm_s") return m_sun*1000.0;
  else if (unit_string == "m_sunau_yr") return m_sun*r_au/yr;
  else if (unit_string == "kgm_s") return 1.0;
  else if (unit_string == "gcm_s") return 1.0e-5;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  MomentumUnit::LatexLabel
/// Return latex string of requested momentum unit for external plotting
//=============================================================================
string MomentumUnit::LatexLabel(string unit_string)
{
  if (unit_string == "m_sunkm_s") return "M_{\\odot}\\,km\\,s^{-1}";
  else if (unit_string == "m_sunau_yr") return "M_{\\odot}\\,AU\\,yr^{-1}";
  else if (unit_string == "kgm_s") return "kg\\,m\\,s^{-1}";
  else if (unit_string == "gcm_s") return "g\\,cm\\,s^{-1}";
  else return "";
}



//=============================================================================
//  AngularMomentumUnit::SIUnit
/// Return numerical value requested angular momentum unit in SI units
//=============================================================================
double AngularMomentumUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sunkm2_s") return m_sun*1000.0*1000.0;
  else if (unit_string == "m_sunau2_yr") return m_sun*r_au*r_au/yr;
  else if (unit_string == "kgm2_s") return 1.0;
  else if (unit_string == "gcm2_s") return 1.0e-7;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  AngularMomentumUnit::LatexLabel
/// Return latex string of requested angular mom. unit for external plotting
//=============================================================================
string AngularMomentumUnit::LatexLabel(string unit_string)
{
  if (unit_string == "m_sunkm2_s") return "M_{\\odot}\\,km^2\\,s^{-1}";
  else if (unit_string == "m_sunau2_yr") return "M_{\\odot}\\,AU^{2}\\,yr^{-1}";
  else if (unit_string == "kgm2_s") return "kg\\,m^{2}\\,s^{-1}";
  else if (unit_string == "gcm2_s") return "g\\,cm^{2}\\,s^{-1}";
  else return "";
}



//=============================================================================
//  AngularVelocityUnit::SIUnit
/// Return numerical value requested angular velocity unit in SI units
//=============================================================================
double AngularVelocityUnit::SIUnit(string unit_string)
{
  if (unit_string == "rad_s") return 1.0;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  AngularVelocityUnit::LatexLabel
/// Return latex string of requested angular vel. unit for external plotting
//=============================================================================
string AngularVelocityUnit::LatexLabel(string unit_string)
{
  if (unit_string == "rad_s") return "rad\\,s^{-1}";
  else return "";
}



//=============================================================================
//  MassRateUnit::SIUnit
/// Return numerical value requested Mass rate unit in SI units
//=============================================================================
double MassRateUnit::SIUnit(string unit_string)
{
  if (unit_string == "m_sun_myr") return m_sun/myr;
  else if (unit_string == "m_sun_yr") return m_sun/yr;
  else if (unit_string == "kg_s") return 1.0;
  else if (unit_string == "g_s") return 1.0e-3;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
}



//=============================================================================
//  MassRateUnit::LatexLabel
/// Return latex string of requested mass rate unit for external plotting
//=============================================================================
string MassRateUnit::LatexLabel(string unit_string)
{
  if (unit_string == "m_sun_myr") return "M\\_{\\odot}\\,Myr^{-1}";
  else if (unit_string == "m_sun_yr") return "M\\_{\\odot}\\,yr^{-1}";
  else if (unit_string == "kg_s") return "kg\\,s^{-1}";
  else if (unit_string == "g_s") return "g\\,s^{-1}";
  else return "";
}



//=============================================================================
//  SpecificEnergyUnit::SIUnit
/// Return numerical value requested specific energy unit in SI units
//=============================================================================
double SpecificEnergyUnit::SIUnit(string unit_string)
{
  if (unit_string == "J_kg") return 1.0;
  else if (unit_string == "erg_g") return 1.0e-4;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
  return -1;
}



//=============================================================================
//  SpecificEnergyUnit::LatexLabel
/// Return latex string of requested specific energy unit for external plotting
//=============================================================================
string SpecificEnergyUnit::LatexLabel(string unit_string)
{
  if (unit_string == "J_kg") return "J\\,kg^{-1}";
  else if (unit_string == "erg_g") return "erg\\,g^{-1}";
  else return "";
}



//=============================================================================
//  SpecificEnergyRateUnit::SIUnit
/// Return numerical value requested specific energy rate unit in SI units
//=============================================================================
double SpecificEnergyRateUnit::SIUnit(string unit_string)
{
  if (unit_string == "J_kg_s") return 1.0;
  else if (unit_string == "erg_g_s") return 1.0e-4;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
  return -1;
}



//=============================================================================
//  SpecificEnergyRateUnit::LatexLabel
/// Return latex string of requested dudt unit for external plotting
//=============================================================================
string SpecificEnergyRateUnit::LatexLabel(string unit_string)
{
  if (unit_string == "J_kg_s") return "J\\,kg^{-1}\\,s^{-1}";
  else if (unit_string == "erg_g_s") return "erg\\,g^{-1}\\,s^{-1}";
  else return "";
}



//=============================================================================
//  TemperatureUnit::SIUnit
/// Return numerical value requested temperature unit in SI units
//=============================================================================
double TemperatureUnit::SIUnit(string unit_string)
{
  if (unit_string == "K") return 1.0;
  else if (unit_string == "") return 1.0;
  else if (unit_string != "") {
    string message = "Parameter error : Unrecognised unit = " + unit_string;
    ExceptionHandler::getIstance().raise(message);
  }
  return -1;
}



//=============================================================================
//  TemperatureUnit::LatexLabel
/// Return latex string of requested temperature unit for external plotting
//=============================================================================
string TemperatureUnit::LatexLabel(string unit_string)
{
  if (unit_string == "K") return "K";
  else return "";
}



//=============================================================================
//  SimUnits::SimUnits
/// Constructor for SimUnits class
//=============================================================================
SimUnits::SimUnits()
{
  ReadInputUnits = false;
}



//=============================================================================
//  SimUnits::~SimUnits
/// Destructor for SimUnits class
//=============================================================================
SimUnits::~SimUnits()
{
}



//=============================================================================
//  Units::SetupUnits
/// Calculate all scaling variables based on the chosen units from the 
/// parameters file (and possibly input snapshot).
//=============================================================================
void SimUnits::SetupUnits(Parameters &params)
{
  debug1("[SimUnits::SetupUnits]");

  // If we are using dimensionless units, then return immediately
  dimensionless = params.intparams["dimensionless"];
  if (dimensionless) return;

  // If not input units have been read from the snapshot file, then assume
  // units are the same as the output units in parameters file.
  // --------------------------------------------------------------------------
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

  // Pressure units
  // --------------------------------------------------------------------------
  press.inSI = press.SIUnit(params.stringparams["pressinunit"]);
  press.outSI = press.SIUnit(params.stringparams["pressoutunit"]);
  press.inscale = (m.inscale*m.inSI) / 
    (r.inscale*r.inSI*t.inscale*t.inSI*t.inscale*t.inSI);
  press.inscale = press.inscale / press.inSI;
  press.outscale = (m.outscale*m.outSI) / 
    (r.outscale*r.outSI*t.outscale*t.outSI*t.outscale*t.outSI);
  press.outscale = press.outscale / press.outSI;
  press.outcgs = 0.1*press.outSI;

  // Force units
  // --------------------------------------------------------------------------
  f.inSI = f.SIUnit(params.stringparams["forceinunit"]);
  f.outSI = f.SIUnit(params.stringparams["forceoutunit"]);
  f.inscale = (m.inscale*m.inSI*r.inscale*r.inSI)/
    (t.inscale*t.inscale*t.inSI*t.inSI);
  f.inscale = f.inscale / f.inSI;
  f.outscale = (m.outscale*m.outSI*r.outscale*r.outSI)/
    (t.outscale*t.outscale*t.outSI*t.outSI);
  f.outscale = f.outscale / f.outSI;
  f.outcgs = 1.0e5*f.outSI;

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

  // Mass rate units
  // --------------------------------------------------------------------------
  dmdt.inSI = dmdt.SIUnit(params.stringparams["dmdtinunit"]);
  dmdt.outSI = dmdt.SIUnit(params.stringparams["dmdtoutunit"]);
  dmdt.inscale = m.inscale*m.inSI/(t.inscale*t.inSI);
  dmdt.inscale = dmdt.inscale / dmdt.inSI;
  dmdt.outscale = m.outscale*m.outSI/(t.outscale*t.outSI);
  dmdt.outscale = dmdt.outscale / dmdt.outSI;
  dmdt.outcgs = 1.0e3*dmdt.outSI;

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

  return;
}
