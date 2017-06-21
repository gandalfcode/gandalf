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

#include "Particle.h"
#include "Precision.h"
#include "Exception.h"

template <int ndim>
class GradhSphCommunicationHandler {

  struct GradSphForcesParticle {

    GradSphForcesParticle (const GradhSphParticle<ndim>& p) {
      for (int k=0; k<ndim; k++) {
        a[k] = p.a[k];
        atree[k] = p.atree[k];
      }
      gpot = p.gpot;
      dudt = p.dudt;
      div_v = p.div_v;
      levelneib = p.levelneib;
      iorig = p.iorig;
    }

    FLOAT a[ndim];
    FLOAT atree[ndim];
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

  void ReceiveParticleAccelerations (const ReturnDataType* pointer, GradhSphParticle<ndim>& p2) {

    const ReturnDataType& p = *pointer;

    for (int k=0; k<ndim; k++) {
      p2.a[k] += p.a[k];
      p2.atree[k] += p.atree[k];
    }

    p2.gpot += p.gpot;
    p2.dudt += p.dudt;
    p2.div_v += p.div_v;
    p2.levelneib = max(p.levelneib,p2.levelneib);

  }

  void ReceiveParticle (const void* pointer, GradhSphParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    const DataType& p = *reinterpret_cast<const DataType*>(pointer);
    //GradhSph<ndim>* sph = static_cast<GradhSph<ndim>*>(hydro);

    p2.iorig = p.iorig;
    p2.flags = type_flag(p.flags);
    p2.ptype = p.ptype;
    p2.level = p.level;
    for (int k=0; k<ndim; k++) {
      p2.r[k]=p.r[k];
      p2.v[k]=p.v[k];
      p2.a[k]=0;
      p2.atree[k]=0;

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
    p2.h = hydro->h_fac*pow(p2.m/p2.rho, Sph<ndim>::invndim);
    p2.hrangesqd = hydro->kernrange*hydro->kernrange*p2.h*p2.h;
    p2.hfactor = pow(1/p2.h,ndim+1);
    p2.sound = hydro->eos->SoundSpeed(p2);
    p2.flags.set(active);

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
        iorig = p.iorig;
        flags = p.flags.get();
        ptype = p.ptype;
        level = p.level;
        nstep = p.nstep;
    	for (int k=0; k<ndim; k++) {
    		r[k] = p.r[k];
    		v[k] = p.v[k];
    		for (int kk=0; kk<ndim; kk++) {
    			B[k][kk]=p.B[k][kk];
    		}
    	}
    	for (int ivar=0; ivar<ndim+2; ivar++) {
    		Wprim[ivar]=p.Wprim[ivar];
    		alpha_slope[ivar]=p.alpha_slope[ivar];
    		for (int k=0; k<ndim; k++) {
        		grad[ivar][k] = p.grad[ivar][k];
    		}
    	}
    	m = p.m;
    	ndens = p.ndens;
    	vsig_max = p.vsig_max;
    	sound = p.sound;
    }

    int iorig;
    unsigned int flags;
    int ptype;
    int level;
    int nstep;

    FLOAT r[ndim];
    FLOAT v[ndim];
    FLOAT m;
    FLOAT ndens;
    FLOAT Wprim[ndim+2];
    FLOAT grad[ndim+2][ndim];
    FLOAT B[ndim][ndim];
    FLOAT alpha_slope[ndim+2];
    FLOAT vsig_max;
    FLOAT sound;
  };

  struct MeshlessForcesParticle {
    MeshlessForcesParticle (const MeshlessFVParticle<ndim>& p) {

      for (int k=0; k<ndim; k++) {
    	  rdmdt[k]=p.rdmdt[k];
    	  a[k]=p.a[k];
    	  atree[k]=p.atree[k];
      }
      for (int ivar=0; ivar<ndim+2; ivar++) {
    	  dQ[ivar]=p.dQ[ivar];
    	  dQdt[ivar]=p.dQdt[ivar];
      }
      iorig = p.iorig;
      gpot = p.gpot;
      vsig_max = p.vsig_max;

    }

    int iorig;
    FLOAT dQ[ndim+2];
    FLOAT dQdt[ndim+2];
    FLOAT rdmdt[ndim];
    FLOAT a[ndim];
    FLOAT atree[ndim];
    FLOAT gpot;
    FLOAT vsig_max;
  };

public:
  typedef MeshlessExportParticle DataType;
  typedef MeshlessForcesParticle ReturnDataType;

  void ReceiveParticleAccelerations (const ReturnDataType* pointer, MeshlessFVParticle<ndim>& p2) {
	    const ReturnDataType& p = *pointer;

	    for (int k=0; k<ndim; k++) {
	      p2.rdmdt[k] += p.rdmdt[k];
	      p2.a[k] += p.a[k];
	      p2.atree[k] += p.atree[k];
	    }

	    for (int ivar=0; ivar<ndim+2; ivar++) {
	    	  p2.dQ[ivar]+=p.dQ[ivar];
	    	  p2.dQdt[ivar]+=p.dQdt[ivar];
	    }

	    p2.gpot += p.gpot;
	    p2.vsig_max = max((FLOAT) p2.vsig_max,p.vsig_max);
  }

  void ReceiveParticle (const void* pointer, MeshlessFVParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    const DataType& p = *reinterpret_cast<const DataType*>(pointer);
    p2.iorig = p.iorig;
    p2.flags = type_flag(p.flags);
    p2.ptype = p.ptype;
    p2.level = p.level;
    p2.nstep = p.nstep;

    for (int k=0; k<ndim; k++) {
      p2.r[k]=p.r[k];
      p2.v[k]=p.v[k];
      p2.a[k]=0.0;
      p2.atree[k]=0;
	  for (int kk=0; kk<ndim; kk++) {
			p2.B[k][kk]=p.B[k][kk];
	  }
	  p2.rdmdt[k]=0.0;
  }

	for (int ivar=0; ivar<ndim+2; ivar++) {
		p2.Wprim[ivar]=p.Wprim[ivar];
		p2.alpha_slope[ivar]=p.alpha_slope[ivar];
		for (int k=0; k<ndim; k++) {
			p2.grad[ivar][k] = p.grad[ivar][k];
		}
		p2.dQ[ivar]=0.0;
		p2.dQdt[ivar]=0.0;
}

  p2.m = p.m;
  p2.ndens = p.ndens;
  p2.gpot = 0.0;
  p2.vsig_max=p.vsig_max;
  p2.sound=p.sound;

  //Recompute h dependent stuff
  p2.rho = p2.ndens*p2.m;
  p2.h = hydro->h_fac*pow(1/p2.ndens, (FLOAT)(1.)/ndim);
  p2.hfactor = pow(1/p2.h, ndim+1);
  p2.hrangesqd = hydro->kernrange*hydro->kernrange*p2.h*p2.h;
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

  void ReceiveParticleAccelerations (const ReturnDataType* pointer, SM2012SphParticle<ndim>& p2) {
    ExceptionHandler::getIstance().raise("not implemented");
  }

  void ReceiveParticle (const void* pointer, SM2012SphParticle<ndim>& p2, Hydrodynamics<ndim>* hydro) {
    //const DataType& p = *reinterpret_cast<const DataType*>(pointer);
    ExceptionHandler::getIstance().raise("not implemented");


  }
};

template <int ndim>
class TreeCommunicationHandler {

  struct TreeCellStreamlined {

    TreeCellStreamlined (const TreeCellBase<ndim>& c, int Nactive, int exported_particles) {

        ifirst = exported_particles;
        ilast = exported_particles+Nactive-1;
        N = Nactive;
        amin = c.amin;
      }


    FLOAT amin ;
    int ifirst;
    int ilast;
    int N;

    };

public:

  typedef TreeCellStreamlined DataType;

  void ReceiveCell (const void* pointer, TreeCellBase<ndim>& c2,int Ntot) {
    const DataType& c = *reinterpret_cast<const DataType*>(pointer);

    c2.ifirst = c.ifirst + Ntot;
    c2.ilast = c.ilast + Ntot;

    c2.Nactive = c.N;
    c2.N = c.N;
    c2.worktot=0.;

    // Only need amin as amin / macfactor are in a union
    c2.amin = c.amin ;
  }

  template <template <int> class ParticleType>
  void ReconstructProperties (TreeCellBase<ndim>& c, ParticleType<ndim>* partdata, FLOAT kernrange) {

    c.m = 0;
    c.hmax = 0;
    c.maxsound = 0.0f;
    c.amin = big_number;
    for (int k=0; k<ndim; k++) {
      c.bb.min[k] = big_number;
      c.bb.max[k] = -big_number;
      c.hbox.min[k] = big_number;
      c.hbox.max[k] = -big_number;
      c.vbox.min[k] = big_number;
      c.vbox.max[k] = -big_number;
      c.r[k] = 0;
    }


    for (int i=c.ifirst; i<=c.ilast; i++) {
      const ParticleType<ndim>& p = partdata[i];
      for (int k=0; k<ndim; k++) {
        c.r[k] += p.m*p.r[k];
        if (c.bb.min[k] > p.r[k]) c.bb.min[k] = p.r[k];
        if (c.bb.max[k] < p.r[k]) c.bb.max[k] = p.r[k];
        if (p.r[k] - kernrange*p.h < c.hbox.min[k])
          c.hbox.min[k] = p.r[k] - kernrange*p.h;
        if (p.r[k] + kernrange*p.h > c.hbox.max[k])
          c.hbox.max[k] = p.r[k] + kernrange*p.h;
        if (p.v[k] > c.vbox.max[k]) c.vbox.max[k] = p.v[k];
        if (p.v[k] < c.vbox.min[k]) c.vbox.min[k] = p.v[k];
      }
      c.m += p.m;
      c.hmax = max(c.hmax,p.h);
      c.maxsound = max(c.maxsound, p.sound);
      c.amin = min(c.amin,
                   sqrt(DotProduct(partdata[i].atree,partdata[i].atree,ndim)));
    }

    FLOAT dr[ndim];
    for (int k=0; k<ndim; k++) {
      c.r[k] /= c.m;
      //c.rcell[k] = (FLOAT) 0.5*(c.bb.min[k] + c.bb.max[k]);
      dr[k] = (FLOAT) 0.5*(c.bb.max[k] - c.bb.min[k]);
    }
    c.rmax = sqrt(DotProduct(dr,dr,ndim));

  }

};


#endif
