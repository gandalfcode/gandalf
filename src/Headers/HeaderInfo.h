/*
 * HeaderInfo.h
 *
 *  Created on: Jul 3, 2013
 *      Author: giovanni
 */

#ifndef HEADERINFO_H_
#define HEADERINFO_H_

#include "Precision.h"

struct HeaderInfo {
  int Nhydro;
  int Nstar;
  int Ndust;
  int ndim;
  DOUBLE t;

  HeaderInfo(): Nhydro(0), Nstar(0), Ndust(0), ndim(0), t(0) {};

};


#endif /* HEADERINFO_H_ */
