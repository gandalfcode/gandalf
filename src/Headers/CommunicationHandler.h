//=================================================================================================
//  CommunicationHandler.h
//  Classes that handle the communication of the particle data structures
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#ifndef _COMMUNICATION_HANDLER_H_
#define _COMMUNICATION_HANDLER_H_

template <int ndim>
class GradhSphCommunicationHandler {

  struct GradhSphExportParticle {
    GradhSphExportParticle () {};

    GradhSphExportParticle (const GradhSphParticle<ndim>& p) {
      iorig = p.iorig;
      itype = p.itype;
      ptype = p.ptype;
      levelneib = p.levelneib;
      level = p.level;
      for (int k=0; k<ndim; k++) {
        r[k] = p.r[k];
        v[k] = p.v[k];
      }
      m = p.m;
      rho = p.rho;
      h_dust = p.h_dust;
      u = p.u;

      pfactor = p.pfactor;
      alpha = p.alpha;

      invomega = p.invomega;
      zeta = p.zeta;
    }

    int iorig;
    int itype;
    int ptype;
    int levelneib;
    int level;
    FLOAT r[ndim];
    FLOAT v[ndim];
    FLOAT m;
    FLOAT rho;
    FLOAT h_dust;
    FLOAT u;

    FLOAT pfactor;
    FLOAT alpha;

    FLOAT invomega;
    FLOAT zeta;

  };

public:
  typedef GradhSphExportParticle DataType;

  void ReceiveParticle (void* pointer, GradhSphParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    DataType& p = *reinterpret_cast<DataType*>(pointer);
    //GradhSph<ndim>* sph = static_cast<GradhSph<ndim>*>(hydro);

    p2.iorig = p.iorig;
    p2.itype = p.itype;
    p2.ptype = p.ptype;
    p2.levelneib = p.levelneib;
    p2.level = p.level;
    for (int k=0; k<ndim; k++) {
      p2.r[k]=p.r[k];
      p2.v[k]=p.v[k];
    }
    p2.m = p.m;
    p2.rho = p.rho;
    p2.h_dust = p.h_dust;
    p2.u = p.u;

    p2.pfactor = p.pfactor;
    p2.alpha = p.alpha;

    p2.invomega = p.invomega;
    p2.zeta = p.zeta;

    //Recompute things we store only for optimization
    p2.invrho = 1./p2.rho;
    p2.h = hydro->h_fac*powf(p2.m*p2.invrho, Sph<ndim>::invndim);
    p2.invh = 1./p2.h;
    p2.hrangesqd = hydro->kernfacsqd*hydro->kernrange*hydro->kernrange*p2.h*p2.h;
    p2.hfactor = pow(p2.invh,ndim+1);
    p2.sound = hydro->eos->SoundSpeed(p2);
    p2.active = true;



  }

};

template <int ndim>
class MeshlessCommunicationHandler {

  struct MeshlessExportParticle {
    MeshlessExportParticle () {};

    MeshlessExportParticle (const MeshlessFVParticle<ndim>& p) {
      ExceptionHandler::getIstance().raise("not implemented");
    }

    FLOAT rho;
    int iorig;

  };

public:
  typedef MeshlessExportParticle DataType;

  void ReceiveParticle (void* pointer, MeshlessFVParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    DataType& p = *reinterpret_cast<DataType*>(pointer);
    ExceptionHandler::getIstance().raise("not implemented");


  }
};

template <int ndim>
class SM2012CommunicationHandler {

  struct SM2012ExportParticle {
    SM2012ExportParticle () {};

    SM2012ExportParticle (const SM2012SphParticle<ndim>& p) {
      ExceptionHandler::getIstance().raise("not implemented");
    }

    int iorig;
    int itype;
    int ptype;
    int levelneib;
    int level;
    FLOAT r[ndim];
    FLOAT v[ndim];
    FLOAT m;
    FLOAT h;
    FLOAT h_dust;
    FLOAT u;
    FLOAT rho;


  };

public:
  typedef SM2012ExportParticle DataType;
  void ReceiveParticle (void* pointer, SM2012SphParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    DataType& p = *reinterpret_cast<DataType*>(pointer);
    ExceptionHandler::getIstance().raise("not implemented");


  }
};



#endif
