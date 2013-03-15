#include "SphKernel.h"


SphKernel* KernelFactory (int ndimaux, string KernelName) {
  if (KernelName == "m4")
    return new M4Kernel(ndimaux, KernelName);
  else if (KernelName == "quintic")
    return new QuinticKernel(ndimaux,KernelName);
  else {
    string message = "Unrecognised kernel: " + KernelName;
    ExceptionHandler::getIstance().raise(message);
  }
  return NULL;
}


TabulatedKernel::TabulatedKernel(int ndimaux, string KernelName, int resaux)
    {
    res = resaux;

    kernel = KernelFactory (ndimaux, KernelName);

    kernrange = kernel->kernrange;
    kernrangesqd = kernel->kernrangesqd;
    invkernrange = kernel->invkernrange;
    resinvkernrange = res/kernrange;
    resinvkernrangesqd = res/kernrangesqd;
    kernnorm = kernel->kernnorm;
#if !defined (FIXED_DIMENSIONS)
    ndim = kernel->ndim;
    ndimpr = kernel->ndimpr;
#endif

    //allocate memory
    tableW0 = new FLOAT[res];
    tableW1 = new FLOAT[res];
    tableWomega = new FLOAT[res];
    tableWzeta = new FLOAT[res];
    tableWgrav = new FLOAT[res];
    tableWpot = new FLOAT[res];
    tableW0_s2 = new FLOAT[res];
    tableWomega_s2 = new FLOAT[res];
    tableWzeta_s2 = new FLOAT[res];
    tableLOS = new FLOAT[res];


    //initialize the tables
    initializeTable(tableW0,&SphKernel::w0);
    initializeTable(tableW1,&SphKernel::w1);
    initializeTable(tableWomega,&SphKernel::womega);
    initializeTable(tableWzeta,&SphKernel::wzeta);
    initializeTable(tableWgrav,&SphKernel::wgrav);
    initializeTable(tableWpot,&SphKernel::wpot);
    initializeTableSqd(tableW0_s2,&SphKernel::w0);
    initializeTableSqd(tableWomega_s2,&SphKernel::womega);
    initializeTableSqd(tableWzeta_s2,&SphKernel::wzeta);
    initializeTableLOS();

    //deallocates the kernel now that we don't need it anymore
    delete kernel;
  }

void TabulatedKernel::initializeTableLOS() {
    const FLOAT step = kernel->kernrange/res; //step in the tabulated variable
    const int intsteps = 4000; //how many steps for each integration
    FLOAT sum;
    FLOAT s; //distance from the center
    for (int i=0; i<res; i++) {
      FLOAT impactparameter = i*step;
      FLOAT impactparametersqd = impactparameter*impactparameter;
      sum=0;
      FLOAT dist = sqrt(kernrangesqd-impactparametersqd); //half-length of the integration path
      FLOAT intstep = dist/intsteps;
      for (int j=0; j<intsteps; j++) {
        FLOAT position = intstep*j;
        s = sqrt(position*position+impactparametersqd); //compute distance from the center
        sum += kernel->w0(s)*intstep;
      }
      tableLOS[i] = 2*sum; //multiply by 2 because we integrated only along half of the path
    }
  }
