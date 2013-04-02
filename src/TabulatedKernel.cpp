#include "SphKernel.h"

template <int ndim>
SphKernel<ndim>* SphKernel<ndim>::KernelFactory (string KernelName) {
  if (KernelName == "m4")
    return new M4Kernel<ndim>(KernelName);
  else if (KernelName == "quintic")
    return new QuinticKernel<ndim>(KernelName);
  else if (KernelName == "gaussian")
    return new GaussianKernel<ndim>(KernelName);
  else {
    string message = "Unrecognised kernel: " + KernelName;
    ExceptionHandler::getIstance().raise(message);
  }
  return NULL;
}

template <int ndim>
TabulatedKernel<ndim>::TabulatedKernel(string KernelName, int resaux):
  SphKernel<ndim>()
    {
    res = resaux;

    kernel = SphKernel<ndim>::KernelFactory (KernelName);

    this->kernrange = kernel->kernrange;
    this->kernrangesqd = kernel->kernrangesqd;
    this->invkernrange = kernel->invkernrange;
    resinvkernrange = res/this->kernrange;
    resinvkernrangesqd = res/this->kernrangesqd;
    this->kernnorm = kernel->kernnorm;
//#if !defined (FIXED_DIMENSIONS)
//    ndim = kernel->ndim;
//    ndimpr = kernel->ndimpr;
//#endif

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
    initializeTable(tableW0,&SphKernel<ndim>::w0);
    initializeTable(tableW1,&SphKernel<ndim>::w1);
    initializeTable(tableWomega,&SphKernel<ndim>::womega);
    initializeTable(tableWzeta,&SphKernel<ndim>::wzeta);
    initializeTable(tableWgrav,&SphKernel<ndim>::wgrav);
    initializeTable(tableWpot,&SphKernel<ndim>::wpot);
    initializeTableSqd(tableW0_s2,&SphKernel<ndim>::w0);
    initializeTableSqd(tableWomega_s2,&SphKernel<ndim>::womega);
    initializeTableSqd(tableWzeta_s2,&SphKernel<ndim>::wzeta);
    initializeTableLOS();

    //deallocates the kernel now that we don't need it anymore
    delete kernel;
  }

template <int ndim>
void TabulatedKernel<ndim>::initializeTableLOS() {
    const FLOAT step = kernel->kernrange/res; //step in the tabulated variable
    const int intsteps = 4000; //how many steps for each integration
    FLOAT sum;
    FLOAT s; //distance from the center
    for (int i=0; i<res; i++) {
      FLOAT impactparameter = i*step;
      FLOAT impactparametersqd = impactparameter*impactparameter;
      sum=0;
      FLOAT dist = sqrt(this->kernrangesqd-impactparametersqd); //half-length of the integration path
      FLOAT intstep = dist/intsteps;
      for (int j=0; j<intsteps; j++) {
        FLOAT position = intstep*j;
        s = sqrt(position*position+impactparametersqd); //compute distance from the center
        sum += kernel->w0(s)*intstep;
      }
      tableLOS[i] = 2*sum; //multiply by 2 because we integrated only along half of the path
    }
  }

template class TabulatedKernel<1>;
template class TabulatedKernel<2>;
template class TabulatedKernel<3>;
