//=================================================================================================
//  TreeRay.cpp
//  ...
//
//=================================================================================================


#include "TreeRay.h"
#include "chealpix.h"



//=================================================================================================
//  TreeRay::TreeRay
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
TreeRay<ndim,nfreq,ParticleType,TreeCell>::TreeRay(bool _ots, int _nSide, int _ilNR, int _ilNTheta,
                                                   int _ilNPhi, int _ilNNS, int _ilFinePix,
                                                   FLOAT _maxDist, FLOAT _rayRadRes, FLOAT _relErr,
                                                   string _errControl, DomainBox<ndim> &simbox,
                                                   SimUnits *units, Parameters *params):
  Radiation<ndim>(),
  onTheSpot(_ots),
  nSide(_nSide),
  ilNR(_ilNR),
  ilNTheta(_ilNTheta),
  ilNPhi(_ilNPhi),
  ilNNS(_ilNNS),
  ilFinePix(_ilFinePix),
  maxDist(_maxDist),
  rayRadRes(_rayRadRes),
  relErr(_relErr),
  errControl(_errControl)
{
  debug2("[TreeRay::TreeRay]");

  // Create tree object
  tree = new OctTree<ndim,ParticleType,TreeCell>(params->intparams["Nleafmax"],
                                                 params->floatparams["thetamaxsqd"],
                                                 (FLOAT) 2.0,
                                                 params->floatparams["macerror"],
                                                 params->stringparams["gravity_mac"],
                                                 params->stringparams["multipole"]);

  // Set important variables
  NTBLevels = tree->lmax + 1;
  nPix      = 12*nSide*nSide;
  ilNI      = nPix;
  ilNSSampFac = 1.6;

  // Generate intersection list
  GenerateIntersectList();


  // Initialise On-the-spot variables
  //-----------------------------------------------------------------------------------------------
  if (onTheSpot) {

    // conversion between eint and T
    /*tr_osTH2ToEint = tr_boltz*tr_osTH2 / (tr_osAbarH2*tr_mH*(tr_osGammaH2-1.0))
    tr_osTHpToEint = tr_boltz*tr_osTHp / (tr_osAbarHp*tr_mH*(tr_osGammaHp-1.0))

    // recombination constant
    tr_osRecombConst = tr_osAlphaStar * (tr_osXhydro / tr_mH)**2 / 3.0

    // constant for converting photon flux to  energy density
    tr_osEflx2Erad = tr_osUVPhotonE / (tr_lightSpeed)
    tr_osErad2Eflx = 1.0 / tr_osEflx2Erad*/

    // DAVID : More things in tr_isInit such as reading in sources.  Ask Richard
  }


  // Allocate arrays for rays
  //call Grid_getMinCellSize(tr_bhMinCellSize)
  minCellSize = 0.0f;
  for (int k=0; k<ndim; k++) minCellSize = max(minCellSize, simbox.boxsize[k]);
  minCellSize = minCellSize/powf(2.0, tree->lmax);
  max_ray_length = sqrtf(DotProduct(simbox.boxsize, simbox.boxsize, ndim));
  bhNR = rayRadRes*((int) (sqrtf(2.0*max_ray_length / minCellSize)) + 1);

  cout << "boxsize : " << simbox.boxsize[0] << "   " << simbox.boxsize[1] << "   " << simbox.boxsize[2] << endl;
  cout << "minCellSize : " << minCellSize << "    max_ray_length : " << max_ray_length
       << "     bhNR : " << bhNR << endl;

  rays = new Rays*[nPix];
  for (int i=0; i<nPix; i++) rays[i] = new Rays[bhNR+1];


  rayR   = new FLOAT[bhNR+1];
  rayR2  = new FLOAT[bhNR+1];
  rayRi  = new FLOAT[bhNR+1];
  rayR2i = new FLOAT[bhNR+1];

  for (int ir=0; ir<bhNR+1; ir++) {
    rayR[ir]   = 0.5*minCellSize*(FLOAT) ir*(FLOAT) ir/(rayRadRes*rayRadRes);
    rayR2[ir]  = rayR[ir]*rayR[ir];
    rayRi[ir]  = 1.0/(rayR[ir] + 1.0e-99);
    rayR2i[ir] = 1.0/(rayR2[ir] + 1.0e-99);
  }
  // DAVID Add extra arrays here later perhaps!!??

  // ..
  GenerateRadNodeMapping();

}



//=================================================================================================
//  TreeRay::~TreeRay
/// ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
TreeRay<ndim,nfreq,ParticleType,TreeCell>::~TreeRay()
{
  delete[] rayR2i;
  delete[] rayRi;
  delete[] rayR2;
  delete[] rayR;

  for (int i=nPix-1; i>= 0; i--) delete[] rays[i];
  delete[] rays;
}



//=================================================================================================
//  TreeRay::UpdateRadiationField
/// ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::UpdateRadiationField
 (int Nhydro,                          ///< [in] No. of SPH particle
  int Nnbody,                          ///< [in] No. of N-body particles
  int Nsink,                           ///< [in] No. of sink particles
  SphParticle<ndim> *sph_gen,          ///< [in] Generic SPH particle data array
  NbodyParticle<ndim> **nbodydata,     ///< [in] N-body data array
  SinkParticle<ndim> *sinkdata)        ///< [in] Sink data array
{
  int it;
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* >(sph_gen);


  debug2("[TreeRay::UpdateRadiationField]");

  // Build, re-build or re-stock the tree
  //-----------------------------------------------------------------------------------------------
  tree->Ntot       = Nhydro;
  tree->Ntotmax    = max(tree->Ntotmax,tree->Ntot);
  tree->BuildTree(0, Nhydro - 1, Nhydro, tree->Ntotmax, partdata, 0.0);

  // Add sink particles to the tree
  AddRadiationSources(Nsink, sinkdata);

  // Now compute all tree cell properties, including point sources (i.e. sinks)
  tree->StockTree(tree->celldata[0], partdata);



  //-----------------------------------------------------------------------------------------------
  for (it=0; it<Nitmax; it++) {

    // reset variables to store error - for iterations (From TreeRay_bhInitFieldVar)
    bhLocRelErr = 0.0;
    bhMaxRelEradErr = 0.0;
    bhLocEradTot = 0.0;
    bhLocMionTot = 0.0;


    for (c=0; c<tree->Ncell; c++) {
      OsTreeRayCell<ndim> &cell = tree->celldata[c];
      cell.erdEUVold[0] = cell.erdEUV[0];
      cell.erdEUV[0] = 0.0;

      TreeWalk(cell);
    }

    // Need to add exit conditions (total energy or individual cells; TreeRay_bhTreeWalkEnd.F90)


  }
  //-----------------------------------------------------------------------------------------------

  for (c=0; c<tree->Ncell; c++) {
    OsTreeRayCell<ndim> &cell = tree->celldata[c];

  }




  return;
}



//=================================================================================================
//  TreeRay::GenerateIntersectList
//  ...
/// For now, ignored any MPI coding (i.e. to check which processors things should be on)
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::GenerateIntersectList(void)
{
  int i1d;                                     // ..
  int iclosest;                                // ..
  int il;                                      // ..
  int ir;                                      // ..
  int max_il = 0;                              // ..
  int node_count;                              // ..
  int ray_count[nPix];                         // ..
  long int ifpix;                              // ..
  long int ipix;                               // ..
  long int ispix;                              // ..
  FLOAT dotprod;                               // ..
  FLOAT ilNSSampFacI = 1.0/ilNSSampFac;        // ..
  FLOAT max_dotprod;                           // ..
  FLOAT nshalf;                                // ..
  FLOAT phfpix;                                // ..
  FLOAT phi;                                   // ..
  FLOAT r[ndim];                               // ..
  FLOAT rad;                                   // ..
  FLOAT rnode[ndim];                           // ..
  FLOAT theta;                                 // ..
  FLOAT thfpix;                                // ..

  debug2("[TreeRay::GenerateIntersectList]");

  intersectList = new FLOAT[ilNI*ilNNS*(ilNPhi+1)*(ilNTheta+1)];
  for (int i=0; i<ilNI*ilNNS*(ilNPhi+1)*(ilNTheta+1); i++) intersectList[i] = -1.0;


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(i1d,ins,iph,ith,rnode,theta)
  for (int ith=0; ith<ilNTheta + 1; ith++) {
    theta = pi*(FLOAT) ith/(FLOAT) ilNTheta;

    //---------------------------------------------------------------------------------------------
    for (int iph=0; iph<ilNPhi + 1; iph++) {
      phi = twopi*(FLOAT) iph/(FLOAT) ilNPhi;
      i1d = ith*(ilNPhi + 1) + iph;

      rnode[0] = sin(theta)*cos(phi);
      rnode[1] = sin(theta)*sin(phi);
      rnode[2] = cos(theta);

      //-------------------------------------------------------------------------------------------
      for (int ins=0; ins<ilNNS; ins++) {
        nshalf = 0.5*ilNSSampFac*(FLOAT) (ins + 1)/(FLOAT) ilNNS;
        node_count = 0;
        for (ipix=0; ipix<nPix; ipix++) ray_count[ipix] = 0;

        //-----------------------------------------------------------------------------------------
        for (ipix=0; ipix<nPix; ipix++) {
          for (ispix=0; ispix<ilFinePix*ilFinePix; ispix++) {
            ifpix = ilFinePix*ilFinePix*ipix + ispix;
            pix2ang_nest(nSide*ilFinePix, ifpix, &thfpix, &phfpix);

            for (ir=0; ir<ilNR+1; ir++) {
              rad = 3.0*(FLOAT) ir/(FLOAT) ilNR;
              r[0] = rad*sin(thfpix)*cos(phfpix);
              r[1] = rad*sin(thfpix)*sin(phfpix);
              r[2] = rad*cos(thfpix);

              if ((fabs(r[0] - rnode[0]) < nshalf) && (fabs(r[1] - rnode[1]) < nshalf) &&
                  (fabs(r[2] - rnode[2]) < nshalf)) {
                node_count++;
                ray_count[ipix]++;
              }
            }
          }
        }
        //-----------------------------------------------------------------------------------------

        // Store list of intersections between the node and rays into the array
        il = 0;
        for (ipix=0; ipix<nPix; ipix++) {
          if (ray_count[ipix] > 0) {
            if (il > max_il) max_il = il;
            if (il >= ilNI) {
              cout << "GenIntersectList: Too many intersection. Increase ilNI" << endl;
              exit(0);
            }
            intersectList[IIL(il,ins,iph,ith)] =
              (FLOAT) 0.999*(FLOAT) ray_count[ipix]/(FLOAT) node_count + (FLOAT) ipix;
            il++;
          }
        }

        // In case no intersection is found, just add the closest ray
        if (node_count == 0) {
          max_dotprod = (FLOAT) 0.0;
          for (ipix=0; ipix<nPix; ipix++) {
            pix2ang_nest(nSide, ipix, &thfpix, &phfpix);
            r[0] = sin(thfpix)*cos(phfpix);
            r[1] = sin(thfpix)*sin(phfpix);
            r[2] = cos(thfpix);
            dotprod = DotProduct(r, rnode, ndim);
            if (dotprod > max_dotprod) {
              max_dotprod = dotprod;
              iclosest = ipix;
            }
          }
          intersectList[IIL(0,ins,iph,ith)] = (FLOAT) iclosest + (FLOAT) 0.999;
        }


      }
      //-------------------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  // Write table to file (for debugging)
  FILE * pFile;
  pFile = fopen ("gandalfTreeRay.dat","w");

  fprintf(pFile, "# %d %d %d\n", ilNNS, ilNTheta, ilNPhi);
  for (int ins=0; ins<ilNNS; ins++) {
    FLOAT ns = ilNSSampFac*(FLOAT) (ins + 1)/(FLOAT) ilNNS;

    for (int ith=0; ith<ilNTheta+1; ith++) {
      theta = pi*(FLOAT) ith/(FLOAT) ilNTheta;

      for (int iph=0; iph<ilNPhi+1; iph++) {
        phi = twopi*(FLOAT) iph/(FLOAT) ilNPhi;

        rnode[0] = sin(theta)*cos(phi);
        rnode[1] = sin(theta)*sin(phi);
        rnode[2] = cos(theta);
        fprintf(pFile, "%12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  ",
                ns, theta, phi, rnode[0], rnode[1], rnode[2]);
        for (int il=0; il<ilNI; il++) {
          fprintf(pFile, "%15.8e  ",intersectList[IIL(il, ins, iph, ith)]);
        }
        fprintf(pFile, "\n");
      }
      fprintf(pFile, "\n");

    }
    fprintf(pFile, "\n");
  }
  fclose(pFile);

  return;
}



//=================================================================================================
//  TreeRay::GenerateRadNodeMapping
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::GenerateRadNodeMapping(void)
{
  int irStore;
  FLOAT ndSize = minCellSize;
  FLOAT r;
  FLOAT rMin;
  FLOAT rNodeCen;
  FLOAT totWeight;
  FLOAT weight;

  debug2("[TreeRay::GenerateRadNodeMapping]");


  radNodeMapIndex = new int[bhNR*nFineR*bhNR*NL];
  radNodeMapValue = new FLOAT[bhNR*nFineR*bhNR*NL];


  for (int i=0; i<bhNR*nFineR*bhNR*NL; i++) radNodeMapIndex[i] = -1;

  // DAVID : What is NTBLevels???  nint??
  //-----------------------------------------------------------------------------------------------
  for (int il=NTBLevels-1; il >=0; il--) {

    //---------------------------------------------------------------------------------------------
    for (int irf=0; irf<nFineR*bhNR; irf++) {

      rNodeCen = (FLOAT) 0.5*minCellSize*(FLOAT) (irf + 1)*(FLOAT) (irf + 1) /
        (nFineR*nFineR*rayRadRes*rayRadRes);
      totWeight = 0.0;
      irStore = 0;

      for (int ir=0; ir<bhNR; ir++) {
        r = (FLOAT) 0.5*minCellSize*(FLOAT) (ir + 1)*(FLOAT) (ir + 1)/(rayRadRes*rayRadRes);

        weight = NodeKernel(r, rNodeCen, ndSize);
        if (weight > (FLOAT) 0.0) {
          totWeight += weight;
          radNodeMapIndex[IRNM(irStore,irf,il)] = ir;
          radNodeMapValue[IRNM(irStore,irf,il)] = weight;
          irStore++;
        }
      }

      // Takes the nearest point in case no one is found
      if (totWeight == 0.0) {
        rMin = 1.0e99;

        for (int ir=0; ir<ilNR; ir++) {
          r = (FLOAT) 0.5*minCellSize*(FLOAT) (ir + 1)*(FLOAT) (ir + 1)/(rayRadRes*rayRadRes);
          if (fabs(r - rNodeCen) < rMin) {
            rMin = fabs(r - rNodeCen);
            irStore = ir;
          }
        }
        radNodeMapIndex[IRNM(0,irf,il)] = irStore;
        radNodeMapValue[IRNM(0,irf,il)] = 1.0;
        totWeight = 1.0;
      }

      // Renormalize
      for (int ir=0; ir<ilNR; ir++) {
        if (radNodeMapIndex[IRNM(ir,irf,il)] > 0) {
          radNodeMapValue[IRNM(ir,irf,il)] = radNodeMapValue[IRNM(ir,irf,il)] / totWeight;
        }
      }


      for (int ir=0; ir<bhNR; ir++) {
        FLOAT r = 0.5*rayRadRes*minCellSize*powf(radNodeMapIndex[IRNM(ir,irf,il)],2);
        cout << "RNM: " << il << "   " << ndSize << "   " << irf << "   " << rNodeCen
             << "    " << ir << "    " << r << "    " << radNodeMapIndex[IRNM(ir,irf,il)]
             << "    " << radNodeMapValue[IRNM(ir,irf,il)] << "   " << totWeight << endl;
        if (radNodeMapIndex[IRNM(ir,irf,il)] < 0) break;
      }

    }
    //---------------------------------------------------------------------------------------------

    ndSize *= 2.0;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  TreeRay::NodeKernel
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
FLOAT TreeRay<ndim,nfreq,ParticleType,TreeCell>::NodeKernel
 (const FLOAT r,
  const FLOAT rNodeCen,
  const FLOAT ndSize)
{
  static const FLOAT a1 = 1.21;
  static const FLOAT a2 = -6.35;
  static const FLOAT M = 1.184;
  static const FLOAT A = 0.40802854426;

  const FLOAT rnorm = fabs((r - rNodeCen)/ndSize);

  if (rnorm < (FLOAT) 0.5) {
    return a1*powf(rnorm,3) - ((FLOAT) 0.5*a1 + (FLOAT) 4.0*(M-A))*rnorm*rnorm + M;
  }
  else if (rnorm < (FLOAT) 0.5*sqrt((FLOAT) 3.0)) {
    return a2*powf(rnorm,3) - (0.5*a2*(1.0 + 2*sqrt(3.)) - 2.0*A*(sqrt(3.) + 2.0))*rnorm*rnorm -
      (9.0*a2/4.0 - 0.5*a2*(sqrt(3.0) + 6.0) + 2.0*A*(3.0 + 2.0*sqrt(3.0)))*rnorm +
      3.0*sqrt(3.0)*a2/4.0 - 3.0*a2*(1.0 + 2.0*sqrt(3.0))/8.0 + 1.5*A*(sqrt(3.0) + 2.0);
  }
  else {
    return -1.0;
  }
}



//=================================================================================================
//  TreeRay::TreeWalk
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::TreeWalk
 (TreeCell<ndim> &cell)                ///< [in] Pointer to cell
  //const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
{
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  long int ipix;                       // ..
  long int iR;                         // ..
  int iF;                              // ..
  FLOAT cellEdgeSize;                  // ..
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell


  FLOAT **eflux;


  eflux = new FLOAT*[nPix];
  for (ipix=0; ipix<nPix; ipix++) eflux[ipix] = new FLOAT[nfreq];


  for (ipix=0; ipix<nPix; ipix++) {
    for (iR=0; iR<ilNR; iR++) {
      rays[ipix][iR].rho    = 0.0;
      rays[ipix][iR].mass   = 0.0;
      rays[ipix][iR].volume = 0.0;
      for (iF=0; iF<nfreq; iF++) {
        rays[ipix][iR].srcF[iF] = 0.0;
        rays[ipix][iR].erad[iF] = 0.0;
        rays[ipix][iR].Erad[iF] = 0.0;
      }
    }
  }

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < tree->Ncell) {

    for (k=0; k<ndim; k++) dr[k] = tree->celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    if (drsqd > tree->celldata[cc].cdistsqd && tree->celldata[cc].N > 0) {
      cellEdgeSize = tree->rootCellSize/powf(2.0,tree->celldata[cc].level);
      NodeContribution(cell, tree->celldata[cc], dr, drsqd, cellEdgeSize);
      cc = tree->celldata[cc].cnext;
    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (!(drsqd > tree->celldata[cc].cdistsqd) && tree->celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree->celldata[cc].copen != -1)
        cc = tree->celldata[cc].copen;

      else {
        cellEdgeSize = tree->rootCellSize/powf(2.0,tree->celldata[cc].level);
        NodeContribution(cell, tree->celldata[cc], dr, drsqd, cellEdgeSize);
        cc = tree->celldata[cc].cnext;
      }
    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = tree->celldata[cc].cnext;
    }

  };
  //===============================================================================================


  for (ipix=0; ipix<nPix; ipix++) {
    rays[ipix][0].rho    = cell.m/cell.volume;
    rays[ipix][0].mass   = cell.m;
    rays[ipix][0].volume = cell.volume;
    for (iF=0; iF<nfreq; iF++) {
      rays[ipix][0].Erad[iF] = rays[ipix][0].erad[iF]*cell.volume;
      rays[ipix][0].erad[iF] = cell.srcEUV[iF];
    }
  }


  for (ipix=0; ipix<nPix; ipix++) {
    FixRay(rays[ipix]);
    IntegrateRay(rays[ipix],eflux[ipix]);
  }

  FLOAT **cdMaps;  // DEFINE THIS SAME AS RAYS
  FinaliseCell(cell, eflux, cdMaps);


  for (ipix=0; ipix<nPix; ipix++) delete[] eflux[ipix];
  delete[] eflux;


  return;
}



//=================================================================================================
//  TreeRay::NodeContribution
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::NodeContribution
 (TreeCell<ndim> &targetNode,              ///< [inout] ..
  TreeCell<ndim> &contributingNode,        ///< [in] ..
  FLOAT dr[ndim],
  FLOAT drsqd,
  FLOAT cellEdgeSize)
{
  int i;                                   // ..
  int iNodeSize;                           // ..
  int ins;                                 // ..
  int iph;                                 // ..
  int ipix;                                // ..
  int ir;                                  // ..
  int irf;                                 // ..
  int ith;                                 // ..
  int j;                                   // ..
  FLOAT drmag;                             // ..
  FLOAT invdrmag;                          // ..
  FLOAT node_vol;                          // ..
  FLOAT ns;                                // Node angular size
  FLOAT phi;                               // ..
  FLOAT theta;                             // ..
  FLOAT tot_weight = 0.0;                  // ..
  FLOAT tot_w2;                            // ..
  FLOAT weight;                            // ..
  FLOAT zor;                               // ..

  // Check if the distance does not exceed maxDist
  if (drsqd > maxDist*maxDist) return;
  drmag    = sqrt(drsqd) + small_number;
  invdrmag = (FLOAT) 1.0/drmag;

  // Node angular size - assumes that the node is cubic
  ns = cellEdgeSize*invdrmag;
  ins = (int) (ns*ilNNS*ilNSSampFacI - (FLOAT) 0.5);
  if ((ins < 0) || (ins >= ilNNS)) {
    cout << "Wrong ins: " << ins << "   " << ilNNS << "   " << ns
         << "    " << (int) (ns*ilNNS - (FLOAT) 0.5) << endl;
    exit(0);
  }

  // Node angular position
  zor   = max(-1.0, min(1.0, dr[2]*invdrmag));
  theta = acos(zor);
  phi   = fmod(atan2(dr[1], dr[0]) + twopi, twopi);
  ith   = (int) (theta*ilNTheta/pi + 0.5);
  iph   = (int) (0.5*phi*ilNPhi/pi + 0.5);
  if ((ith < 0) || (ith > ilNTheta)) {
    cout << "NC Wrong ith: " << ith << "    " << ilNTheta << "    " <<  theta
         << "    " << (int) (theta*ilNTheta/pi + 0.5) << endl;
    //call Driver_abortFlash("TreeRay_bhBotNodeContrib: incorrect theta index")
    exit(0);
  }
  if ((iph < 0) || (iph > ilNPhi)) {
    cout << "Wrong iph: " << iph << "    " << ilNPhi << "     "
         << phi << "    " << (int) (0.5*phi*ilNPhi/pi + 0.5) << endl;
    exit(0);
  }


  // Node mass, volume and column density
  //node_mass = node(tr_bhIM)
  node_vol = powf(cellEdgeSize,3);

  // Node index in the rays-array
  irf = (int) (nFineR*rayRadRes*sqrt(2.0*drmag/minCellSize) - 0.5);
  if ((irf < 0) || (irf >= bhNR*nFineR)) {
    cout << "Wrong irf: " << irf << "   " << bhNR << "    " << drmag
         << minCellSize << "    " << sqrt(2.0*drmag/minCellSize) << endl;
    exit(0);
  }

  // Node size index for RadNodeMap
  iNodeSize = NTBLevels - NINT(log(cellEdgeSize/minCellSize)/log(2.0)) - 1;

  // Map node to the rays
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<ilNI; i++) {
  //do i = 1, tr_ilNI

    // Check if the node intersect with the ray
    if (intersectList[IIL(i,ins,iph,ith)] < 0) break;

    // Determine the ray index (ipix) and the weight of the intersection
    ipix        = (int) (intersectList[IIL(i,ins,iph,ith)]);
    weight      = (intersectList[IIL(i,ins,iph,ith)] - ipix)/0.999;
    tot_weight += weight;

    tot_w2 = 0.0;
    for (j=0; j<bhNR; j++) {
    //do j = 1,tr_bhNR
      ir = radNodeMapIndex[IRNM(j,irf,iNodeSize)];
      if (ir < 0) break;
      tot_w2 += radNodeMapValue[IRNM(j,irf,iNodeSize)];

      rays[ipix][ir].volume += node_vol*weight*radNodeMapValue[IRNM(j,irf,iNodeSize)];
      rays[ipix][ir].mass   += contributingNode.m*weight*radNodeMapValue[IRNM(j,irf,iNodeSize)];

    }
    //if (abs(tot_w2 - 1.0) > 1.e-3) print *, "BNC tw2 = ", tot_w2, tr_bhRadNodeMapVal(:,irf,iNodeSize)
  }
  //if (abs(tot_weight - 1.0) > 1.e-3) print *, "NC tw = ", tot_weight, tr_intersectList(:,ins,iph,ith)


  // Sum all contributions for on-the-spot approximation
  //-----------------------------------------------------------------------------------------------
  if (onTheSpot) {

    // DAVID : What are these?
    FLOAT node_srcf = 0.0;
    FLOAT node_erad = 0.0;

    for (int i=0; i<ilNI; i++) {

      // Check if the node intersect with the ray
      if (intersectList[IIL(i,ins,iph,ith)] < 0) break;

      // Determine the ray index (ipix) and the weight of the intersection
      ipix = (int) (intersectList[IIL(i,ins,iph,ith)]);
      weight = (intersectList[IIL(i,ins,iph,ith)] - (FLOAT) ipix)/0.999;

      for (int j=0; j<bhNR; j++) {
        ir = radNodeMapIndex[IRNM(j,irf,iNodeSize)];
        if (ir < 0) break;
        rays[ipix][ir].srcF[0] += node_srcf*weight*radNodeMapValue[IRNM(j,irf,iNodeSize)];
        rays[ipix][ir].Erad[0] += node_erad*weight*radNodeMapValue[IRNM(j,irf,iNodeSize)];
      }
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  TreeRay::FixRay
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::FixRay
 (Rays *ray)
{
  int ir;
  int ir2;
  int irh;
  int kfreq;
  FLOAT q;

  // If there are still empty points, interpolate along individual rays
  //-----------------------------------------------------------------------------------------------
  for (ir=1; ir<bhNR+1; ir++) {

    // Test for the empty point - zero volume
    //---------------------------------------------------------------------------------------------
    if (ray[ir].volume <= 0.0) {

      // Go along the ray and search for a valid point
      irh = -1;
      for (ir2=ir+1; ir2<bhNR+1; ir2++) {

        // valid point found, set right point for interpolation to it
        if (ray[ir2].volume > 0.0) {
          irh = ir2;
          break;
        }
      }

      // Worst case scenarion: no valid point found
      // Just copy into ir the last valid point
      if (irh == -1) {
        ray[ir].volume = ray[ir-1].volume;
        ray[ir].mass   = ray[ir-1].mass;
        ray[ir].rho    = ray[ir-1].rho;
        for (kfreq=0; kfreq<nfreq; kfreq++) ray[ir].srcF[kfreq] = 0.0;
        for (kfreq=0; kfreq<nfreq; kfreq++) ray[ir].erad[kfreq] = 0.0;
        for (kfreq=0; kfreq<nfreq; kfreq++) ray[ir].Erad[kfreq] = 0.0;
      }

      // Linear interpolation between the last valid point and the
      // first valid point at larger r
      else {

        // DAVID : What is tr_bhRayR??  Radial position?
        //q = (tr_bhRayR(irh) - tr_bhRayR(ir)) / (tr_bhRayR(irh) - tr_bhRayR(ir-1))
        ray[ir].volume = (1.0 - q)*ray[ir-1].volume + q*ray[irh].volume;
        ray[ir].mass   = (1.0 - q)*ray[ir-1].volume + q*ray[irh].mass;
        ray[ir].rho    = ray[ir].mass/ray[ir].volume;
        for (kfreq=0; kfreq<nfreq; kfreq++) {
          ray[ir].srcF[kfreq] = (1.0 - q)*ray[ir-1].srcF[kfreq] + q*ray[irh].srcF[kfreq];
          ray[ir].erad[kfreq] = (1.0 - q)*ray[ir-1].erad[kfreq] + q*ray[irh].erad[kfreq];
          ray[ir].Erad[kfreq] = (1.0 - q)*ray[ir-1].Erad[kfreq] + q*ray[irh].Erad[kfreq];
        }
      }

    }
    // everything fine, just set rho_ray
    //---------------------------------------------------------------------------------------------
    else {
      ray[ir].rho = ray[ir].mass/ray[ir].volume;
      for (kfreq=0; kfreq<nfreq; kfreq++) ray[ir].erad[kfreq] = ray[ir].Erad[kfreq]/ray[ir].volume;
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  TreeRay::IntegrateRay
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::IntegrateRay
 (Rays *ray,
  FLOAT eflux[nfreq])
{

  //===============================================================================================
  if (onTheSpot) {

    int i;
    int ir;
    int is;
    FLOAT Nphotons;
    FLOAT Nrecomb;
    FLOAT rho_2mean;
    FLOAT dV;
    FLOAT nH_av;
    FLOAT recombContrib;
    FLOAT recombContribSrc;
    FLOAT eflux_ray[bhNR+1];
    FLOAT eRadFluxSrc[bhNR+1];
    FLOAT eFluxSrc[bhNR+1];
    FLOAT recombSrc[bhNR+1];
    //real,dimension(0:tr_bhNR) :: eflux_ray, eRadFluxSrc, eFluxSrc, recombSrc

    FLOAT tr_osEflx2Erad = 1.0;  // DAVID : What are these?
    FLOAT tr_osAlphaStar = 1.0;

    // Radiant flux of energy of a source: energy (# of photons) per unit space angle,
    // in the direction of the point-of-calulation, for each source at a given point on the ray
    // initially, fill with fluxes of sources
    for (i=0; i<bhNR+1; i++) eRadFluxSrc[i] = ray[i].srcF[0]/(4.0*pi);

    // Energy (# of photons) per unit area, in the direction of the point-of-calulation,
    // from all sources added together, at each point along the ray
    for (i=0; i<bhNR+1; i++) eflux_ray[i] = 0.0;


    //---------------------------------------------------------------------------------------------
    for (ir=bhNR-1; ir>=0; ir--) {

      // Average particle density between point ir and ir+1
      // DAVID : Add physical constants later (osXhydro/mH)
      nH_av = 0.5*(ray[ir+1].rho + ray[ir].rho);
      //nH_av = (0.5*tr_osXhydro/tr_mH)*(rho_ray(ir+1) + rho_ray(ir))

      // Contribution to recombination for this ray
      recombContrib = min(2.0, eflux_ray[ir+1]*tr_osEflx2Erad / (ray[ir+1].erad[0] + 1.0E-99));

      // Calculate recombination rates belonging to sources
      //-------------------------------------------------------------------------------------------
      for (is=ir+1; is<bhNR+1; is++) {

        // DAVID : What is tr_bhRayR??
        // volume of the element between point ir and ir+1, on the cone coming from source "is"
        //dV = ((tr_bhRayR(is)-tr_bhRayR(ir))**3 - (tr_bhRayR(is)-tr_bhRayR(ir+1))**3)/3.0

        if (eRadFluxSrc[is] > 0.0) {
          // eFluxSrc includes values for ir+1 point
          recombContribSrc = recombContrib*min(1.0, eFluxSrc[is]/(eflux_ray[ir+1] + 1.0e-99));
          recombSrc[is] = tr_osAlphaStar*nH_av*nH_av*dV*recombContribSrc;
          eRadFluxSrc[is] = max(0.0, eRadFluxSrc[is] - recombSrc[is]);
          //eFluxSrc[is] = eRadFluxSrc[is]/(tr_bhRayR(is)-tr_bhRayR(ir))**2    DAVID : tr_bhRayR
          eflux_ray[ir] = eflux_ray[ir] + eFluxSrc[is];
        }

      }
      //-------------------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------

    eflux[0] = 0.5*(eflux_ray[0] + eflux_ray[1]);

  }
  //===============================================================================================

  return;
}



//=================================================================================================
//  TreeRay::FinaliseCell
//  ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
void TreeRay<ndim,nfreq,ParticleType,TreeCell>::FinaliseCell
 (TreeCell<ndim> &cell,
  FLOAT **eflux,
  FLOAT **cdMaps)
{
  FLOAT phFluxInt[nfreq];

  for (int k=0; k<nfreq; k++) {
    phFluxInt[k] = 0.0;
    for (int i=0; i<nPix; i++) {
      phFluxInt[k] += eflux[ipix][k];
    }
  }


  if (onTheSpot) {
    cell.erdEUV[0] = Eflx2Erad*phFluxInt[0];

    tr_bhLocEradTot = tr_bhLocEradTot + solnPoint(EUVE_VAR)
    tr_bhLocMionTot = tr_bhLocMionTot + solnPoint(DENS_VAR)*solnPoint(IHP_SPEC)*vol_poc
    EradErr = 2.0*abs(solnPoint(EUVE_VAR) - solnPoint(EUVO_VAR)) &
    &       / (solnPoint(EUVE_VAR) + solnPoint(EUVO_VAR) + 1d-99)

    if (EradErr > tr_bhLocRelErr) tr_bhLocRelErr = EradErr
  }


  return;
}



//template class OctTree<3,GradhSphParticle,OsTreeRayCell>;
template class TreeRay<3,1,GradhSphParticle,OsTreeRayCell>;
