// ============================================================================
// Render.h
// ============================================================================


#ifndef _RENDER_H_
#define _RENDER_H_


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include "SphParticle.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include "Exception.h"
#include "Parameters.h"
#include "SphKernel.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// Clas Render
// ============================================================================
class Render
{
 public:

  // Constructor and Destructor
  // --------------------------------------------------------------------------
  Render();
  ~Render();

  // Subroutine prototypes
  // --------------------------------------------------------------------------
  int CreateRenderingGrid(int, int, string, string, string, string, float &,
		       float, float, float, float, float, float *, float *,
			   SphSnapshot &, SphKernel *);


  // ..
  // --------------------------------------------------------------------------



};


#endif
