#include "Precision.h"
#include "Particle.h"
#include "Debug.h"
#include "SphNeighbourSearch.h"
#include "Radiation.h"
#include "Sinks.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
using namespace std;

//===================================================================
//	MultipleSourceIonisation::MultipleSourceIonisation
//	MultipleSourceIonisation class constructor
//====================================================================

template <int ndim, template<int> class ParticleType>
MultipleSourceIonisation<ndim,ParticleType>::MultipleSourceIonisation(
	SphNeighbourSearch<ndim> *sphneibaux ,
	float mu_baraux,
	float mu_ionaux,
	float temp0aux,
	float temp_ionaux,
        double Ndotminaux,
	float gamma_eosaux,
	float scaleaux,
	float tempscaleaux,
	double rad_contaux)
{
  sphneib = sphneibaux;
  mu_bar=mu_baraux;
  mu_ion=mu_ionaux;
  temp0=temp0aux;
  temp_ion=temp_ionaux;
  Ndotmin=Ndotminaux;
  gamma_eos=gamma_eosaux;
  scale=scaleaux;
  tempscale=tempscaleaux;
  rad_cont=rad_contaux;
}

//====================================================================
//	MultipleSourceIonisation::~MultipleSourceIonisation
//	MultipleSourceIonisation class destructor
//====================================================================

template <int ndim, template<int> class ParticleType>
MultipleSourceIonisation<ndim,ParticleType>::~MultipleSourceIonisation()
{
}



//====================================================================
//	MultipleSourceIonisation::UpdateRadiationFieldMMS
//	Calculates the internal energy of particles due to ionising
//	radiation.
//====================================================================

template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::UpdateRadiationField
(int N,
int nos,
int aux,
SphParticle<ndim>  * sphgen,
NbodyParticle<ndim> ** ndata,
SinkParticle<ndim> * sphaux)
{
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sphgen);

  ionisation_intergration(nos,N,ndata,sphdata,scale,tempscale,sphneib,
			  temp0,mu_bar,mu_ion,temp_ion,Ndotmin,1./gamma_eos);

  return;
}



/////////////////////////////////////////////////////////////////////
//Works out the fraction of ionisation each star is responsible for//
/////////////////////////////////////////////////////////////////////
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::probs(
	int &nos,
	ionpar *ionisedsph,
	int *sinkid,
	int &testpart,
	double *ndot)

	{
	int pp;
	double fluxcontrole[nos]; 	//To store temporary radii
	double sum=0;			//Controle to stop devision by zero

	//Loop over sources and add photon flux at current location
	for (pp=0;pp<nos;pp++) 						//Loop over all sources
		{
		if (ionisedsph[ionisedsph[testpart].neigh[pp]].ionised[pp]==1)	//If the particle is ionised by the source we are testing
			fluxcontrole[pp]=ndot[pp]-ionisedsph[ionisedsph[testpart].neigh[pp]].photons[pp];
		else
			fluxcontrole[pp]=0;				//If the particle is not ionised then contribution will be zero
		sum=sum+fluxcontrole[pp];				//Add to controle parameter to ensure we carry out the scaling
		}

	//Scale so total fraction of used photons is one
	for (pp=0;pp<nos;pp++)
		{
		if (sum>0)						//Do we have to scale
			ionisedsph[testpart].prob[pp]=fluxcontrole[pp]/sum;	//Scale to one
		else
			ionisedsph[testpart].prob[pp]=fluxcontrole[pp];	//Pass through zeros as no photons are recived
		}
	}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Recursive function that works out number of photons used in ionisation along the path from the source to the particle//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <int ndim, template<int> class ParticleType>
double MultipleSourceIonisation<ndim, ParticleType>::lost(
	ionpar *ionisedsph,
	int *sinkid,
	double *ndot,
	int &N,
	int &pp,
	int &testpart,
	int &nos,
	int &change)

	{
	double absorbed=0,d1,d2;   //Absorbed is the total amount absorbed, d1 is the distance between the particle and the source, d2 is the distance between the neighbour and the source

	//If the particle has not been checked
	if(ionisedsph[testpart].checked[pp]==0)
		{
		//Works out amount of photons lost along the path
		if (ionisedsph[testpart].sink==0)			//Is the particle a sink (this stops recursion when we get to the source)
			{
			probs(nos,ionisedsph,sinkid,testpart,ndot);	//Call the probs function to work out ionisation fraction of wach source
			if (ionisedsph[testpart].neigh[pp]!=N)		//Does the particle have a neighbour for this sourcetF
				{
				d1=sqrt(pow(ionisedsph[testpart].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[testpart].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[testpart].z-ionisedsph[sinkid[pp]].z,2.));
				d2=sqrt(pow(ionisedsph[ionisedsph[testpart].neigh[pp]].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[ionisedsph[testpart].neigh[pp]].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[ionisedsph[testpart].neigh[pp]].z-ionisedsph[sinkid[pp]].z,2.));

				if (ionisedsph[ionisedsph[testpart].neigh[pp]].sink==0)
					{
					absorbed=((pow((ionisedsph[testpart].rho+ionisedsph[ionisedsph[testpart].neigh[pp]].rho)/2.,2.))/3.)*(pow(d1,3.)-pow(d2,3.))*(ionisedsph[testpart].prob[pp])+lost(ionisedsph,sinkid,ndot,N,pp,ionisedsph[testpart].neigh[pp],nos,change);
					}
				else
					{
					absorbed=((pow((ionisedsph[testpart].rho),2.))/3.)*(pow(d1,3.)-pow(d2,3.))*(ionisedsph[testpart].prob[pp])+lost(ionisedsph,sinkid,ndot,N,pp,ionisedsph[testpart].neigh[pp],nos,change);
					}
					//absorbed=rho**2.*[d1**3-d2**3]/3*transmitonfrac+absorbed(previous particle in chain)
				}
			else
				{
				absorbed=ndot[pp];		//All photons used up as we cant have links to the dummy partcile
				}
			}
		ionisedsph[testpart].photons[pp]=absorbed;

		if((ndot[pp]-absorbed)>0)
			{
			if (ionisedsph[testpart].ionised[pp]==0) 			//Record if the particle is changing state (for convergence)
				{
				change=change+1;
				}
			ionisedsph[testpart].ionised[pp]=1; 				//Set particle as source ionised
			}
		else
			{
			if (ionisedsph[testpart].ionised[pp]==1)			//Record if the particle is changing state (for convergence)
				{
				change=change+1;
				}
			ionisedsph[testpart].ionised[pp]=0; 				//Set particle as not source ionised
			}
		ionisedsph[testpart].checked[pp]=1;
		return absorbed;
		}
	else
		{
		return ionisedsph[testpart].photons[pp];
		}


	}

//////////////////////////////////////////////////////////////
//Works out if a particle is ionised using partner functions//
//////////////////////////////////////////////////////////////
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::photoncount(
	ionpar *ionisedsph,
	int *sinkid,
	double *ndot,
	int &N,
	int &nos,
	int &testpart,
	int &change)

	{
	int pp;

	//Deturmines if a particle is ionised
	ionisedsph[testpart].fionised=0;						//Set as not ionised initially
	for (pp=0;pp<nos;pp++)							//Looping over all sources
		{
		lost(ionisedsph,sinkid,ndot,N,pp,testpart,nos,change);	//Is the ammount lost smaller than the ammount availbile (Calls lost function)
		}

	ionisedsph[testpart].fionised=0;	//Set particle as unionised
	//Check to see if the particle should be ionised
	for (pp=0;pp<nos;pp++)
		{
		if (ionisedsph[testpart].ionised[pp]==1){ionisedsph[testpart].fionised=1;} //Set particle to ionised
		}
	}

//////////////////////////
//Main contole roughtine//
//////////////////////////
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::ionisation_intergration(
	int newnos,				//Number of ionising sources
	int N,					//Number of SPH particles
	NbodyParticle<ndim> ** ndata,           //Source Data
	SphParticle<ndim> * sphgen,		//SPH particle data
	double scale,				//Scaling
	double tempscale,			//Temperature scaling
	SphNeighbourSearch<ndim> * sphneib,	//Neighbour Search Roughtine
	double tn,				//Neutral gas temperature
	double mu_bar,				//Average neutral gas mass
	double mu_ion,				//Average ionised gas mass
	double ti,				//Ionised gas temperature
	double Ndotmin,			//Minimum Ionising output
	double gammam1)				//1/gamma
{

//Casts particle arrays
ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sphgen);

struct timeval start, end;
gettimeofday(&start, NULL);

int ii,jj,kk,pp,tt; //Integer allocation for loops
int debug=2,maxneigh=N; //Debug mode controler and maximum number of neighbours allowed
float delta=0;

//Check that the stellar.dat file is present
bool check;
struct stat buf;
stat("stellar.dat", &buf);
check=S_ISREG(buf.st_mode);
if(check==0)
	{
	cout<<"Stellar.dat is not present in run directory, ionisation will not be included"<<endl;
	return;
	}

//Checks if there are currently any sinks in gandalf and if not exits
if (newnos==0)
	{
	cout<<"No stars"<<endl;
	return;
	}

int nos=0; //Alocation of nos variable
int *newnosid=new int[newnos]; //Alocation to store which sinks are active sources

//Deturmines which sinks are active sources based on user choices
for(ii=0;ii<newnos;ii++)
	{
	if(ndata[ii]->NLyC>=Ndotmin)
		{
		nos=nos+1;
		newnosid[nos-1]=ii;
		}
	}

//Checks if the sinks are of large enough size
if (nos==0)
{
cout<<"No stars of suitable mass"<<endl;
return;
}

cout<<"# of sources followed is "<<nos<<". ";

N=N+nos;	 //Inreases N to accomidate sinks

if (debug==1){
gettimeofday(&end, NULL);

delta = (((end.tv_sec  - start.tv_sec) * 1000000 +
         end.tv_usec - start.tv_usec) / 1.e6)-delta;
cout<<delta<<"s to ";
cout<<"Starting"<<endl; //Debug message
}

//Fill ionisedsph particle array
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create sink id table and ndot table
int *sinkid=new int[nos];
double *ndot=new double[nos];

//Create the ionisedsph array and resize to the number of particles (Number of simulation particiles, inluding sinks. Also the dummy particle.)
ionpar *ionisedsph=new ionpar[N+1];

if(ionisation_fraction.size()==0){
	ionisation_fraction.resize(N+1);
	for (ii;ii<N+1;ii++){
		ionisation_fraction[ii].resize(newnos);
		for (pp=0;pp<nos;pp++){
			ionisation_fraction[ii][pp]=0;
			}
		}
	}

ionisation_fraction.resize(N+1);

//Add ionisedsph particle data
#pragma omp parallel for private(ii,jj)
for (ii=0;ii<N-nos;ii++)
	{
	ionisedsph[ii].x=sphdata[ii].r[0];	//Particle x from gandalf
	ionisedsph[ii].y=sphdata[ii].r[1];	//Particle y from gandalf
	ionisedsph[ii].z=sphdata[ii].r[2];	//Particle z from gandalf
	ionisedsph[ii].rho=sphdata[ii].rho;	//Particle density from gandalf
	ionisedsph[ii].h=sphdata[ii].h;	//Particle h from gandalf
	ionisedsph[ii].sink=0;			//Is the particle a sink
	ionisedsph[ii].t=tn;			//Neutral gas temp
	ionisedsph[ii].fionised=0;		//Is the particle ionised at all
	ionisedsph[ii].neighstorcont=0;	//Controle varible for building of neighstore
	ionisedsph[ii].neighstor=new int[200];	//Array containing references to all particles that consider this one a neighbour
	ionisedsph[ii].angle=new double[nos];
	ionisedsph[ii].neigh=new int[nos];
	ionisedsph[ii].photons=new double[nos];
	ionisedsph[ii].prob=new double[nos];
	ionisedsph[ii].checked=new int[nos];
	ionisedsph[ii].ionised=new int[nos];
	ionisedsph[ii].rad_pre_acc=new double[ndim];
	for(jj=0;jj<ndim;jj++)
		{
		ionisedsph[ii].rad_pre_acc[jj]=0;
		}
	ionisation_fraction[ii].resize(newnos); //Resize operation
	//Correctly filling new spaces with 0
	for(jj=0;jj<(newnos-ionisation_fraction[ii].size());jj++)
		{
		ionisedsph[ii].ionised[ionisation_fraction[ii].size()+jj-1]=0;
		}
	//Copying in current values for active sources
	for(jj=0;jj<nos;jj++)
		{
		ionisedsph[ii].ionised[jj]=ionisation_fraction[ii][newnosid[jj]];
		}
	for(jj=0;jj<nos;jj++)
		{
		ionisedsph[ii].angle[jj]=2.*pi;
		ionisedsph[ii].neigh[jj]=N;
		ionisedsph[ii].checked[jj]=0;
		ionisedsph[ii].prob[jj]=0;
		}
	for(jj=0;jj<200;jj++)
		{
		ionisedsph[ii].neighstor[jj]=N;
		}

	}
//Add sink propertys to sink particles
#pragma omp parallel for private(ii,jj)
for (ii=0;ii<nos;ii++)
	{
	ionisedsph[N-nos+ii].t=ti;						//Set stars to ionised gas temp for smoothing
	ionisedsph[N-nos+ii].sink=1;						//Is the particle a sink
	ionisedsph[N-nos+ii].x=ndata[newnosid[ii]]->r[0];				//Source x from gandalf
	ionisedsph[N-nos+ii].y=ndata[newnosid[ii]]->r[1];				//Source y from gandalf
	ionisedsph[N-nos+ii].z=ndata[newnosid[ii]]->r[2];				//Source z from gandalf
	ionisedsph[N-nos+ii].fionised=1;					//As this is source it is marked as ionised
	ionisedsph[N-nos+ii].neighstorcont=0;					//Controle varible for building of neighstore
	ionisedsph[N-nos+ii].neighstor=new int[200];				//Array containing references to all particles that consider this one a neighbour
	sinkid[ii]=N-nos+ii;
	ndot[ii]=pow(2.4e-24,2.)*ndata[newnosid[ii]]->NLyC/(4.*pi*2.6e-13)*scale; 	//Ndot table s^-1
	ionisedsph[N-nos+ii].angle=new double[nos];
	ionisedsph[N-nos+ii].neigh=new int[nos];
	ionisedsph[N-nos+ii].photons=new double[nos];
	ionisedsph[N-nos+ii].checked=new int[nos];
	ionisedsph[N-nos+ii].prob=new double[nos];
	ionisedsph[N-nos+ii].ionised=new int[nos];
	ionisedsph[N-nos+ii].rad_pre_acc=new double[3];
	for(jj=0;jj<3;jj++)
		{
		ionisedsph[N-nos+ii].rad_pre_acc[jj]=0;
		}
	for(jj=0;jj<nos;jj++)
		{
		ionisedsph[N-nos+ii].angle[jj]=2.*pi;
		ionisedsph[N-nos+ii].neigh[jj]=N;
		ionisedsph[N-nos+ii].checked[jj]=0;
		ionisedsph[N-nos+ii].prob[jj]=0;
		ionisedsph[N-nos+ii].ionised[jj]=0;
		}
	ionisedsph[N-nos+ii].ionised[ii]=1;					//Set the sink location in the ionised array to 1
	for(jj=0;jj<200;jj++)
		{
		ionisedsph[N-nos+ii].neighstor[jj]=N;
		}
	}

//Add controle ionisedsph to which all particles are initally linked
//Large values are used so it is never linked in the furture
ionisedsph[N].x=100000.;
ionisedsph[N].y=100000.;
ionisedsph[N].z=100000.;
ionisedsph[N].rho=1e100;
ionisedsph[N].sink=0;
ionisedsph[N].angle=new double[nos];
ionisedsph[N].neigh=new int[nos];
ionisedsph[N].photons=new double[nos];
ionisedsph[N].checked=new int[nos];
ionisedsph[N].prob=new double[nos];
ionisedsph[N].ionised=new int[nos];
ionisedsph[N].rad_pre_acc=new double[ndim];
for(jj=0;jj<ndim;jj++)
	{
	ionisedsph[N].rad_pre_acc[jj]=0;
	}
for(jj=0;jj<nos;jj++)
	{
	ionisedsph[N].angle[jj]=2.*pi;
	ionisedsph[N].neigh[jj]=N;
	ionisedsph[N].checked[jj]=0;
	ionisedsph[N].prob[jj]=0;
	ionisedsph[N].ionised[jj]=0;
	}


//Debug message
if (debug==1){
gettimeofday(&end, NULL);

delta = (((end.tv_sec  - start.tv_sec) * 1000000 +
         end.tv_usec - start.tv_usec) / 1.e6)-delta;
cout<<delta<<"s to ";
cout<<"Particle arrays created"<<endl;
}

//Find the closest source neighbour in chain for each particle
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int current_paricle_nn[maxneigh],Nneigb;		//Working particle neighbour list,Nneigb is the number of neighbours found by neighbour roughtine
double dot,mag,angletest; 				//holding variables
double distanceii,distancejj,temp_radius,temp_radius2;	//Distance between the source of the test particle ii and the neighbour particle jj

//Begin neighbour find
#pragma omp parallel for private(current_paricle_nn,ii,jj,tt,Nneigb,temp_radius,temp_radius2,dot,mag,angletest,pp,distanceii,distancejj) //Initiate openmp
for (ii=0;ii<N;ii++) 				//Loop over all particles
	{

	Nneigb = sphneib->GetGatherNeighbourList(sphdata[ii].r,sphdata[ii].h*2.,sphgen,N,maxneigh,current_paricle_nn); //Find NNnumber nearest neighbours

	//Checks which sources lie within the smoothing length
	for (tt=0;tt<nos;tt++){
		if(sqrt(pow((ionisedsph[sinkid[tt]].x-ionisedsph[ii].x),2)+pow((ionisedsph[sinkid[tt]].y-ionisedsph[ii].y),2)+pow((ionisedsph[sinkid[tt]].z-ionisedsph[ii].z),2))<=ionisedsph[ii].h*2.)
			{
			current_paricle_nn[Nneigb]=sinkid[tt];
			Nneigb=Nneigb+1;
			}
		}

	//For each of these neighbours write my id to there neighstorcont array
	for(jj=0;jj<Nneigb;jj++)
		{
		//If there is still room in the neighstorecont array (This should always be true but this ensures no seg fault)
		if (ionisedsph[current_paricle_nn[jj]].neighstorcont<200)
			{
			ionisedsph[current_paricle_nn[jj]].neighstor[ionisedsph[current_paricle_nn[jj]].neighstorcont]=ii; //Write id
			ionisedsph[current_paricle_nn[jj]].neighstorcont+=1;					     //Advance id counter
			}
		else
			{
			temp_radius=sqrt(pow(ionisedsph[ii].x-ionisedsph[current_paricle_nn[jj]].x,2.)+pow(ionisedsph[ii].y-ionisedsph[current_paricle_nn[jj]].y,2.)+pow(ionisedsph[ii].z-ionisedsph[current_paricle_nn[jj]].z,2.));
			for(tt=0;tt<200;tt++)
				{
				temp_radius2=sqrt(pow(ionisedsph[ii].x-ionisedsph[ionisedsph[current_paricle_nn[jj]].neighstor[tt]].x,2.)+pow(ionisedsph[ii].y-ionisedsph[ionisedsph[current_paricle_nn[jj]].neighstor[tt]].y,2.)+pow(ionisedsph[ii].z-ionisedsph[ionisedsph[current_paricle_nn[jj]].neighstor[tt]].z,2.));
				if(temp_radius>temp_radius2)
					{
						ionisedsph[current_paricle_nn[jj]].neighstor[tt]=ii; //Write id
						break;
					}
				}
			}

		for(pp=0;pp<nos;pp++)
			{
			//Work out the distances for both test and candidate particle
			distanceii=sqrt(pow(ionisedsph[ii].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[ii].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[ii].z-ionisedsph[sinkid[pp]].z,2.));
			distancejj=sqrt(pow(ionisedsph[current_paricle_nn[jj]].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[current_paricle_nn[jj]].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[current_paricle_nn[jj]].z-ionisedsph[sinkid[pp]].z,2.));
			//If the candidate particle is closer than the test particle it is a candidate (Also has controle so a particle cant be its own neighbour)
			if (distancejj<distanceii and ii!=current_paricle_nn[jj])
				{

				//Use the dot product to work out the angle between the conneting line and the neighbour particle
				dot=((ionisedsph[ii].x-ionisedsph[sinkid[pp]].x)*(ionisedsph[current_paricle_nn[jj]].x-ionisedsph[sinkid[pp]].x))+((ionisedsph[ii].y-ionisedsph[sinkid[pp]].y)*(ionisedsph[current_paricle_nn[jj]].y-ionisedsph[sinkid[pp]].y))+((ionisedsph[ii].z-ionisedsph[sinkid[pp]].z)*(ionisedsph[current_paricle_nn[jj]].z-ionisedsph[sinkid[pp]].z));
				mag=(sqrt(pow(ionisedsph[ii].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[ii].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[ii].z-ionisedsph[sinkid[pp]].z,2.)))*(sqrt(pow(ionisedsph[current_paricle_nn[jj]].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[current_paricle_nn[jj]].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[current_paricle_nn[jj]].z-ionisedsph[sinkid[pp]].z,2.)));
				angletest=acos(dot/mag);

				//If the partcle is a sink set angle to be max (Stops non-relavant sinks becoming a neighbour)
				if (ionisedsph[current_paricle_nn[jj]].sink==1) angletest=2.*pi;

				//If the neighbour is the relavant source set angletest to be neg hence particle will be neighbour
				if (current_paricle_nn[jj]==sinkid[pp]) angletest=-1e50;

				//If Neighbour is closest so far to the connecting line then set current neighbour to be the neighbour
				if (angletest<ionisedsph[ii].angle[pp])
					{
					ionisedsph[ii].angle[pp]=angletest;				//Set new comparison angle to be that of the neighbour
					ionisedsph[ii].neigh[pp]=current_paricle_nn[jj];	//Write particle id to neigh array
					}
				}
			}
		}
	}

//Debug message
if (debug==1){
gettimeofday(&end, NULL);

delta = (((end.tv_sec  - start.tv_sec) * 1000000 +
         end.tv_usec - start.tv_usec) / 1.e6)-delta;
cout<<delta<<"s to ";
cout<<"neigbour step one compleate"<<endl;
}

//Initiate openmp
#pragma omp parallel for private(current_paricle_nn,dot,mag,angletest,jj,pp,distanceii,distancejj,ii,tt,Nneigb,temp_radius,temp_radius2)
for (ii=0;ii<N;ii++) 				//Loop over all particles
	{

	for (pp=0;pp<nos;pp++)		//For each source
		{
		//For each particle that considers me a neighbour
		for (jj=0;jj<ionisedsph[ii].neighstorcont;jj++)
			{
			//Work out the distances for both the test and candidate particle
			distanceii=sqrt(pow(ionisedsph[ii].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[ii].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[ii].z-ionisedsph[sinkid[pp]].z,2.));
			distancejj=sqrt(pow(ionisedsph[ionisedsph[ii].neighstor[jj]].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[ionisedsph[ii].neighstor[jj]].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[ionisedsph[ii].neighstor[jj]].z-ionisedsph[sinkid[pp]].z,2.));
			//If the candidate particle is closer than the test particle it is a candidate (Also has controle so a particle cant be its own neighbour)
			if (distancejj<distanceii and ii!=ionisedsph[ii].neighstor[jj])
				{

				//Use the dot product to work out the angle between the conneting line and the candidate particle
				dot=((ionisedsph[ii].x-ionisedsph[sinkid[pp]].x)*(ionisedsph[ionisedsph[ii].neighstor[jj]].x-ionisedsph[sinkid[pp]].x))+((ionisedsph[ii].y-ionisedsph[sinkid[pp]].y)*(ionisedsph[ionisedsph[ii].neighstor[jj]].y-ionisedsph[sinkid[pp]].y))+((ionisedsph[ii].z-ionisedsph[sinkid[pp]].z)*(ionisedsph[ionisedsph[ii].neighstor[jj]].z-ionisedsph[sinkid[pp]].z));
				mag=(sqrt(pow(ionisedsph[ii].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[ii].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[ii].z-ionisedsph[sinkid[pp]].z,2.)))*(sqrt(pow(ionisedsph[ionisedsph[ii].neighstor[jj]].x-ionisedsph[sinkid[pp]].x,2.)+pow(ionisedsph[ionisedsph[ii].neighstor[jj]].y-ionisedsph[sinkid[pp]].y,2.)+pow(ionisedsph[ionisedsph[ii].neighstor[jj]].z-ionisedsph[sinkid[pp]].z,2.)));
				angletest=acos(dot/mag);

				//If the partcle is a sink set angle to be max (Stops non-relavant sinks becoming a neighbour)
				if (ionisedsph[ionisedsph[ii].neighstor[jj]].sink==1) angletest=2.*pi;

				//If the neighbour is the relavant source set angletest to be neg hence particle will be neighbour (Ensures link to sources)
				if (ionisedsph[ii].neighstor[jj]==sinkid[pp]) angletest=-1e50;

				//If candidate is closest so far to the connecting line then set current candidate to be the neighbour
				if (angletest<ionisedsph[ii].angle[pp])
					{
					ionisedsph[ii].angle[pp]=angletest;				//Set new comparison angle to be that of the neighbour
					ionisedsph[ii].neigh[pp]=ionisedsph[ii].neighstor[jj];	//Write particle id to neigh array
					}
				}
			}

		}
	}

if (debug==1){
gettimeofday(&end, NULL);

delta = (((end.tv_sec  - start.tv_sec) * 1000000 +
         end.tv_usec - start.tv_usec) / 1.e6)-delta;
cout<<delta<<"s to ";
cout<<"All neighbours found"<<endl; //Debug message
}

//Begin working out if particle are ionised
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int change=N;              //Variable for checking convergence
int finalcheck=0;	//Finalcheck is a controle that initials the final check of the solution.
while (change!=0 or finalcheck==0)	//loop until no changes are made (We have converged)
	{

	if (debug==1){
		cout<<"The number of adjustments made is "<<change<<endl;	//Debug message
		}

	//Set state of finalcheck (If no changes made in last iteration then carry out final check)
	if (change==0){finalcheck=1;}

	change=0;	//Set change to zero (This will increase as changes are made in the calculation)

	//Open mp initalised
	#pragma omp parallel for schedule(dynamic)	//Begin open mp in dynamic mode i.e one argument is passed at a time rather than a block
	//For each particle
	for (ii=1;ii<N;ii++)
		{
		//Call fucntion to deturmine if test particle is ionised
		photoncount(ionisedsph,sinkid,ndot,N,nos,ii,change);
		}
	#pragma omp parallel for private(ii,pp)
	for (ii=1;ii<N;ii++)
		{
		for(pp=0;pp<nos;pp++)
			{
			ionisedsph[ii].checked[pp]=0;
			}
		}
	}
if (debug==1){
gettimeofday(&end, NULL);

delta = (((end.tv_sec  - start.tv_sec) * 1000000 +
         end.tv_usec - start.tv_usec) / 1.e6)-delta;
cout<<delta<<"s to ";
cout<<"Iterations compleate"<<endl; //Debug message
}

//Smooth the temperature
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double t,rad,s,w,invmu;

#pragma omp parallel for private(jj,ii,t,rad,s,w,invmu,current_paricle_nn,Nneigb) //Initalise openmp
for (ii=0;ii<N;ii++)								//For each particle
	{
	if (ionisedsph[ii].fionised==1)
		{
		Nneigb = sphneib->GetGatherNeighbourList(sphdata[ii].r,sphdata[ii].h*3.,sphgen,N,maxneigh,current_paricle_nn); //Find NNnumber nearest neighbours

		for (jj=0;jj<Nneigb;jj++)			//For each of the neighbours
			{
			if (ionisedsph[current_paricle_nn[jj]].fionised==0)
				{
				rad=(sqrt(pow(ionisedsph[current_paricle_nn[jj]].x-ionisedsph[ii].x,2.)+pow(ionisedsph[current_paricle_nn[jj]].y-ionisedsph[ii].y,2.)+pow(ionisedsph[current_paricle_nn[jj]].z-ionisedsph[ii].z,2.)));  //Ditance between particles
				s=rad/(ionisedsph[ii].h*1.5);						//Work out s for smoothing kernal

				//Work out w for the kernal
				if (s<1) w=1-(3./2.)*pow(s,2.)+(3./4.)*pow(s,3.);
				else if (s<2) w=(1./4.)*pow(2-s,3.);
				else w=0;

				if(ionisedsph[current_paricle_nn[jj]].t<ti*w)
					{
					ionisedsph[current_paricle_nn[jj]].t=ti*w;
					}
				}
				else
				{
					ionisedsph[current_paricle_nn[jj]].t=ti;
				}
			}
		}
	}

if (debug==1){
gettimeofday(&end, NULL);

delta = (((end.tv_sec  - start.tv_sec) * 1000000 +
         end.tv_usec - start.tv_usec) / 1.e6)-delta;

cout<<delta<<"s to ";
cout<<"Temperatures smoothed"<<endl;
}

double theta,thi,photon_acceleration;
#pragma omp parallel for private(jj,ii,t,rad,s,w,invmu,current_paricle_nn,Nneigb,theta,thi,photon_acceleration) //Initalise openmp
for (ii=0;ii<N;ii++)
	{
	//If the particle is ionised then its temperature must be the ionised temperature
	if (ionisedsph[ii].fionised==1)
		ionisedsph[ii].t=ti;

	//If the particles temp is less than the neutral temp because of smoothing set it back to the neutral temp
	if (ionisedsph[ii].t<tn)
		ionisedsph[ii].t=tn;

	invmu=(((ionisedsph[ii].t-tn)/mu_ion)+((ti-ionisedsph[ii].t)/mu_bar))/(ti-tn);//Work out corrected inverted mean gas particle mass
	ionisedsph[ii].u=ionisedsph[ii].t/tempscale/gammam1*invmu;			//Work out the internal energy

	//Set particle ionisation state
	if (ionisedsph[ii].t==tn)
		{
			sphdata[ii].flags.unset(ionised);
		//sphdata[ii].ionstate=0;
		}
	else if (ionisedsph[ii].fionised==1)
		{
			sphdata[ii].flags.set(ionised);
		//sphdata[ii].ionstate=2;
		}
	else
		{
			sphdata[ii].flags.set(ionised);
		//sphdata[ii].ionstate=1;
		}

	sphdata[ii].u=ionisedsph[ii].u;				//Write new internal energy to gandalf

	//Working out radiation pressure (NOT COMPLEATE)
	photon_acceleration=3.4455561764e-34*rad_cont*ionisedsph[ii].rho/(mu_bar*mu_bar);

	for(jj=0;jj<nos;jj++)
		{
		//Copying ionised state over to holding array
		ionisation_fraction[ii][newnosid[jj]]=ionisedsph[ii].ionised[jj];
		if(ionisedsph[ii].fionised==1)
			{
			theta=atan((ionisedsph[ii].y-ionisedsph[sinkid[jj]].y)/(ionisedsph[ii].z-ionisedsph[sinkid[jj]].z));
			thi=atan((ionisedsph[ii].y-ionisedsph[sinkid[jj]].y)/(ionisedsph[ii].x-ionisedsph[sinkid[jj]].x));
			ionisedsph[ii].rad_pre_acc[0]=ionisedsph[ii].rad_pre_acc[0]+ionisedsph[ii].prob[jj]*photon_acceleration*sin(theta)*cos(thi);
			ionisedsph[ii].rad_pre_acc[1]=ionisedsph[ii].rad_pre_acc[1]+ionisedsph[ii].prob[jj]*photon_acceleration*sin(theta)*sin(thi);
			ionisedsph[ii].rad_pre_acc[2]=ionisedsph[ii].rad_pre_acc[2]+ionisedsph[ii].prob[jj]*photon_acceleration*cos(theta);
			}
		}

	//sphdata[ii].rad_pres[0]=0;//ionisedsph[ii].rad_pre_acc[0];
	//sphdata[ii].rad_pres[1]=0;//ionisedsph[ii].rad_pre_acc[1];
	//sphdata[ii].rad_pres[2]=0;//ionisedsph[ii].rad_pre_acc[2];
	}


//Memory De-allocation
delete [] sinkid;
delete [] ndot;

#pragma omp parallel for private(ii)
for(ii=0;ii<N;ii++)
	{
	delete [] ionisedsph[ii].angle;
	delete [] ionisedsph[ii].neigh;
	delete [] ionisedsph[ii].photons;
	delete [] ionisedsph[ii].neighstor;
	delete [] ionisedsph[ii].checked;
	delete [] ionisedsph[ii].prob;
	delete [] ionisedsph[ii].ionised;
	delete [] ionisedsph[ii].rad_pre_acc;
	}

delete [] ionisedsph;

gettimeofday(&end, NULL);

delta = ((end.tv_sec  - start.tv_sec) * 1000000 +
         end.tv_usec - start.tv_usec) / 1.e6;
cout<<"The time taken to calculate ionisation temperatures = "<<delta<<" s"<<endl;

if (debug==1 or debug==2)
	{
	ofstream myfile;
	myfile.open ("timing.dat",ios::app);
	myfile <<delta<<"\n";
	}

}


template class MultipleSourceIonisation<1,GradhSphParticle>;
template class MultipleSourceIonisation<2,GradhSphParticle>;
template class MultipleSourceIonisation<3,GradhSphParticle>;
template class MultipleSourceIonisation<1,SM2012SphParticle>;
template class MultipleSourceIonisation<2,SM2012SphParticle>;
template class MultipleSourceIonisation<3,SM2012SphParticle>;
