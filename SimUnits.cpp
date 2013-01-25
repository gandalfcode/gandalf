// ============================================================================
// SimUnits.cpp
// ============================================================================


#include <math.h>
#include <map>
#include <string>
#include "Constants.h"
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
    cout << "Parameter error : Unrecognised unit = " << unit << endl;
    exit(0);
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
    cout << "Parameter error : Unrecognised unit = " << unit_string << endl;
    exit(0);
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
    cout << "Parameter error : Unrecognised unit = " << unit_string << endl;
    exit(0);
  }
}


/*
// ============================================================================
// TimeUnit::SIUnit
// ============================================================================
double TimeUnit::SIUnit(string unit_string)
{
  if (unit_string == "km_s") return 1000.0;
  else if (unit_string != "") {
    cout << "Parameter error : Unrecognised unit = " << unit_string << endl;
    exit(0);
  }
}
*/


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
  }

  
  r.inSI = r.SIUnit(params.stringparams["rinunit"]);
  r.outSI = r.SIUnit(params.stringparams["routunit"]);
  r.inscale = 1.0;
  r.outscale = r.inscale*r.inSI/r.outSI;
  r.outcgs = 100.0*r.outSI;

  m.inSI = m.SIUnit(params.stringparams["minunit"]);
  m.outSI = m.SIUnit(params.stringparams["moutunit"]);
  m.inscale = 1.0;
  m.outscale = m.inscale*m.inSI/m.outSI;
  m.outcgs = 1000.0*m.outSI;

  t.inSI = t.SIUnit(params.stringparams["tinunit"]);
  t.outSI = t.SIUnit(params.stringparams["toutunit"]);
  t.inscale = pow(r.inscale*r.inSI,1.5)/sqrt(m.inscale*m.inSI*G_const);
  t.inscale /= t.inSI;
  t.outscale = pow(r.outscale*r.outSI,1.5)/sqrt(m.outscale*m.outSI*G_const);
  t.outscale /= t.outSI;
  t.outcgs = t.outSI;

  cout << "r.inscale : " << r.inscale << "    r.outscale : " << r.outscale << endl;

  cout << "m.inscale : " << m.inscale << "    m.outscale : " << m.outscale << endl;


  cout << "t.inscale : " << t.inscale << "    t.outscale : " << t.outscale << endl;

  return;
}
