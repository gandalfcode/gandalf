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

  struct GradSphForcesParticle {

    GradSphForcesParticle (const GradhSphParticle<ndim>& p) {
      for (int k=0; k<ndim; k++) {
        a[k] = p.a[k];
        agrav[k] = p.agrav[k];
      }
      gpot = p.gpot;
      dudt = p.dudt;
      div_v = p.div_v;
      levelneib = p.levelneib;
      iorig = p.iorig;
    }

    FLOAT a[ndim];
    FLOAT agrav[ndim];
    FLOAT gpot;
    FLOAT dudt;
    FLOAT div_v;
    int levelneib;
    int iorig;

  };

  struct GradhSphExportParticle {

    GradhSphExportParticle (const GradhSphParticle<ndim>& p) {
      iorig = p.iorig;
      flags = p.flags.get();
      ptype = p.ptype;
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
    unsigned int flags;
    int ptype;
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
  typedef GradSphForcesParticle ReturnDataType;

  void ReceiveParticleAccelerations (ReturnDataType* pointer, GradhSphParticle<ndim>& p2) {

    const ReturnDataType& p = *pointer;

    for (int k=0; k<ndim; k++) {
      p2.a[k] += p.a[k];
      p2.agrav[k] += p.agrav[k];
    }

    p2.gpot += p.gpot;
    p2.dudt += p.dudt;
    p2.div_v += p.div_v;
    p2.levelneib = max(p.levelneib,p2.levelneib);

  }

  void ReceiveParticle (void* pointer, GradhSphParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    DataType& p = *reinterpret_cast<DataType*>(pointer);
    //GradhSph<ndim>* sph = static_cast<GradhSph<ndim>*>(hydro);

    p2.iorig = p.iorig;
    p2.flags = type_flag(p.flags);
    p2.ptype = p.ptype;
    p2.level = p.level;
    for (int k=0; k<ndim; k++) {
      p2.r[k]=p.r[k];
      p2.v[k]=p.v[k];
      p2.a[k]=0;
      p2.agrav[k]=0;
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
    p2.flags.set_flag(active);

    p2.gpot=0;
    p2.div_v=0;
    p2.dudt = 0;
    p2.levelneib=0;


  }

};

template <int ndim>
class MeshlessCommunicationHandler {

  struct MeshlessExportParticle {

    MeshlessExportParticle (const MeshlessFVParticle<ndim>& p) {
      ExceptionHandler::getIstance().raise("not implemented");
    }

    FLOAT rho;
    int iorig;

  };

  struct MeshlessForcesParticle {
    MeshlessForcesParticle (const MeshlessFVParticle<ndim>& p) {
      ExceptionHandler::getIstance().raise("not implemented");
    }

    int iorig;
  };

public:
  typedef MeshlessExportParticle DataType;
  typedef MeshlessForcesParticle ReturnDataType;

  void ReceiveParticleAccelerations (ReturnDataType* pointer, MeshlessFVParticle<ndim>& p2) {
    ExceptionHandler::getIstance().raise("not implemented");
  }

  void ReceiveParticle (void* pointer, MeshlessFVParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    DataType& p = *reinterpret_cast<DataType*>(pointer);
    ExceptionHandler::getIstance().raise("not implemented");


  }
};

template <int ndim>
class SM2012CommunicationHandler {

  struct SM2012ExportParticle {

    SM2012ExportParticle (const SM2012SphParticle<ndim>& p) {
      ExceptionHandler::getIstance().raise("not implemented");
    }

    int iorig;
    unsigned int flags;
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

  struct SM2012ForcesParticle {
    SM2012ForcesParticle(const SM2012SphParticle<ndim>& p) {
      ExceptionHandler::getIstance().raise("not implemented");
    }

    int iorig;

  };

public:
  typedef SM2012ExportParticle DataType;
  typedef SM2012ForcesParticle ReturnDataType;

  void ReceiveParticleAccelerations (ReturnDataType* pointer, SM2012SphParticle<ndim>& p2) {
    ExceptionHandler::getIstance().raise("not implemented");
  }

  void ReceiveParticle (void* pointer, SM2012SphParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    DataType& p = *reinterpret_cast<DataType*>(pointer);
    ExceptionHandler::getIstance().raise("not implemented");


  }
};

template <int ndim>
class TreeCommunicationHandler {

  struct TreeCellStreamlined {

    TreeCellStreamlined (int Nactive, int exported_particles) {

        ifirst = exported_particles;
        ilast = exported_particles+Nactive-1;
        N = Nactive;

      }


    int ifirst;
    int ilast;
    int N;



    };

public:

  typedef TreeCellStreamlined DataType;

  void ReceiveCell (void* pointer, TreeCellBase<ndim>& c2,int Ntot) {
    const DataType& c = *reinterpret_cast<DataType*>(pointer);

    c2.ifirst = c.ifirst + Ntot;
    c2.ilast = c.ilast + Ntot;

    c2.Nactive = c.N;
    c2.N = c.N;


  }

  template <template <int> class ParticleType>
  void ReconstructProperties (TreeCellBase<ndim>& c, ParticleType<ndim>* partdata, FLOAT kernrange) {

    c.m = 0;
    c.hmax = 0;
    for (int k=0; k<ndim; k++) {
      c.bbmin[k] = big_number;
      c.bbmax[k] = -big_number;
      c.hboxmin[k] = big_number;
      c.hboxmax[k] = -big_number;
      c.r[k] = 0;
    }


    for (int i=c.ifirst; i<=c.ilast; i++) {
      const ParticleType<ndim>& p = partdata[i];
      for (int k=0; k<ndim; k++) {
        c.r[k] += p.m*p.r[k];
        if (c.bbmin[k] > p.r[k]) c.bbmin[k] = p.r[k];
        if (c.bbmax[k] < p.r[k]) c.bbmax[k] = p.r[k];
        if (p.r[k] - kernrange*p.h < c.hboxmin[k])
          c.hboxmin[k] = p.r[k] - kernrange*p.h;
        if (p.r[k] + kernrange*p.h > c.hboxmax[k])
          c.hboxmax[k] = p.r[k] + kernrange*p.h;
      }
      c.m += p.m;
      c.hmax = max(c.hmax,p.h);
    }

    FLOAT dr[ndim];
    for (int k=0; k<ndim; k++) {
      c.r[k] /= c.m;
      c.rcell[k] = (FLOAT) 0.5*(c.bbmin[k] + c.bbmax[k]);
      dr[k] = (FLOAT) 0.5*(c.bbmax[k] - c.bbmin[k]);
    }
    c.rmax = sqrt(DotProduct(dr,dr,ndim));

  }

};


#endif
