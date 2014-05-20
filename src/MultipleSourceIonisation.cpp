#include "Precision.h"
#include "SphParticle.h"
#include "Debug.h"
//#include "ionisingmodule_MS.cpp"
#include "SphNeighbourSearch.h"
#include "Radiation.h"
#include "Sinks.h"
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
        float Ndotminaux,
	float gamma_eosaux,
	float scaleaux,
	float tempscaleaux)
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
void MultipleSourceIonisation<ndim, ParticleType>::probs(int &nos,vector<particle> &sph,vector<int> &sinkid,int &testpart,vector<double> &ndot)
	{
	int pp;
	vector<double> fluxcontrole; 	//To store temporary radii
	double sum=0;			//Controle to stop devision by zero
	fluxcontrole.resize(nos);	//Resizing array to number of sources
	sph[testpart].prob.resize(nos);	//^^^
	
	//Loop over sources and add photon flux at current location
	for (pp=0;pp<nos;pp++) 						//Loop over all sources
		{
		if (sph[sph[testpart].neigh[pp]].ionised[pp]==1)	//If the particle is ionised by the source we are testing
			fluxcontrole[pp]=ndot[pp]-sph[sph[testpart].neigh[pp]].photons[pp];
		else
			fluxcontrole[pp]=0;				//If the particle is not ionised then contribution will be zero
		sum=sum+fluxcontrole[pp];				//Add to controle parameter to ensure we carry out the scaling
		}

	//Scale so total fraction of used photons is one			
	for (pp=0;pp<nos;pp++)
		{		
		if (sum>0)						//Do we have to scale
			sph[testpart].prob[pp]=fluxcontrole[pp]/sum;	//Scale to one
		else
			sph[testpart].prob[pp]=fluxcontrole[pp];	//Pass through zeros as no photons are recived
		}
	}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Recursive function that works out number of photons used in ionisation along the path from the source to the particle//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <int ndim, template<int> class ParticleType>
double MultipleSourceIonisation<ndim, ParticleType>::lost(vector<particle> &sph,vector<int> &sinkid,vector<double> &ndot,int &N,int &pp,int &testpart,int &nos,int &change)
	{
	double absorbed=0,d1,d2;   //Absorbed is the total amount absorbed, d1 is the distance between the particle and the source, d2 is the distance between the neighbour and the source

	//Works out amount of photons lost along the path
	if (sph[testpart].sink==0)			//Is the particle a sink (this stops recursion when we get to the source)
		{
		probs(nos,sph,sinkid,testpart,ndot);	//Call the probs function to work out ionisation fraction of wach source
		if (sph[testpart].neigh[pp]!=N)		//Does the particle have a neighbour for this sourcetF
			{			
			d1=sqrt(pow(sph[testpart].x-sph[sinkid[pp]].x,2.)+pow(sph[testpart].y-sph[sinkid[pp]].y,2.)+pow(sph[testpart].z-sph[sinkid[pp]].z,2.));
			d2=sqrt(pow(sph[sph[testpart].neigh[pp]].x-sph[sinkid[pp]].x,2.)+pow(sph[sph[testpart].neigh[pp]].y-sph[sinkid[pp]].y,2.)+pow(sph[sph[testpart].neigh[pp]].z-sph[sinkid[pp]].z,2.));				
			
			if (sph[sph[testpart].neigh[pp]].sink==0)
				{
				absorbed=((pow((sph[testpart].rho+sph[sph[testpart].neigh[pp]].rho)/2.,2.))/3.)*(pow(d1,3.)-pow(d2,3.))*(sph[testpart].prob[pp])+lost(sph,sinkid,ndot,N,pp,sph[testpart].neigh[pp],nos,change);		
				}
			else
				{
				absorbed=((pow((sph[testpart].rho),2.))/3.)*(pow(d1,3.)-pow(d2,3.))*(sph[testpart].prob[pp])+lost(sph,sinkid,ndot,N,pp,sph[testpart].neigh[pp],nos,change);		
				}			
				//absorbed=rho**2.*[d1**3-d2**3]/3*transmitonfrac+absorbed(previous particle in chain)			
			}
		else
			{
			absorbed=ndot[pp];		//All photons used up as we cant have links to the dummy partcile
			}
		}
	sph[testpart].photons[pp]=absorbed;
	

	if((ndot[pp]-absorbed)>0)
		{
		if (sph[testpart].ionised[pp]==0) 			//Record if the particle is changing state (for convergence)
			{
			change=change+1;
			}
		sph[testpart].ionised[pp]=1; 				//Set particle as source ionised
		}
	else
		{
		if (sph[testpart].ionised[pp]==1)			//Record if the particle is changing state (for convergence)
			{
			change=change+1;
			}
		sph[testpart].ionised[pp]=0; 				//Set particle as not source ionised
		}
	return absorbed;
	}

//////////////////////////////////////////////////////////////
//Works out if a particle is ionised using partner functions//
//////////////////////////////////////////////////////////////
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::photoncount(vector<particle> &sph,vector<int> &sinkid,vector<double> &ndot,int &N,int &nos,int &testpart,int &change)
	{
	int pp;
	
	//Deturmines if a particle is ionised
	sph[testpart].fionised=0;						//Set as not ionised initially
	for (pp=0;pp<nos;pp++)							//Looping over all sources
		{
		lost(sph,sinkid,ndot,N,pp,testpart,nos,change);	//Is the ammount lost smaller than the ammount availbile (Calls lost function)
		}

	sph[testpart].fionised=0;	//Set particle as unionised	
	//Check to see if the particle should be ionised
	for (pp=0;pp<nos;pp++)
		{
		if (sph[testpart].ionised[pp]==1){sph[testpart].fionised=1;} //Set particle to ionised
		}
	}

//////////////////////////
//Main contole roughtine//
//////////////////////////
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::ionisation_intergration(
	int nos,				//Number of ionising sources
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


//////////TEMP/////////
Ndotmin=1e48;

//Casts particle arrays
ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sphgen);

struct timeval start, end;
gettimeofday(&start, NULL);

//Checks if there are currently any sinks in gandalf and if not exits
if (nos==0)
{
cout<<"No stars"<<endl;
return;
} 

int ii,jj,kk,pp,tt; //Integer allocation for loops
int debug=0,maxneigh=500; //Debug mode controler and maximum number of neighbours allowed
N=N+nos;	 //Inreases N to accomidate sinks

if (debug==1){
cout<<"Starting"<<endl; //Debug message
}

//Fill sph particle array
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create sink id table				
vector<int> sinkid;
sinkid.resize(nos);

//Allocate and resize all varibles relating to number of sources
vector<int> ionised;			//Array to store which sources particle is ionised from
vector<int> neigh;	                //Array to hold the closest source in the neighbour train
vector<double> prob;			//Array to hold the absorption probabilitys
vector<double> ndot;			//Array to hold the photon output of each source
vector<double> photons;			//Array to hold lost photons up to this evaluation point
ionised.resize(nos);
neigh.resize(nos);
prob.resize(nos);
ndot.resize(nos);
photons.resize(nos);

//Fill the variables with starting values
#pragma omp parallel for private(ii)
for (ii=0;ii<nos;ii++)
	{
	ionised[ii]=0;  //Each will start as zero and be corrected on particle fill
	neigh[ii]=N;	//Each neighbour will start as the dummy particle
	prob[ii]=0;	//Each will start as zero and be corrected on particle fill
	}

//Create the sph array and resize to the number of particles (Number of simulation particiles, inluding sinks. Also the dummy particle.)
vector<particle> sph;
sph.resize(N+1);

//Add sph particle data
#pragma omp parallel for private(ii)
for (ii=0;ii<N-nos;ii++)
	{
	sph[ii].x=sphdata[ii].r[0];	//Particle x from gandalf
	sph[ii].y=sphdata[ii].r[1];	//Particle y from gandalf
	sph[ii].z=sphdata[ii].r[2];	//Particle z from gandalf
	sph[ii].rho=sphdata[ii].rho;	//Particle density from gandalf
	sph[ii].h=sphdata[ii].h;	//Particle h from gandalf
	sph[ii].sink=0;			//Is the particle a sink 
	sph[ii].ionised=ionised;	
	sph[ii].neigh=neigh;
	sph[ii].prob=prob;
	sph[ii].photons=photons;
	sph[ii].t=tn;			//Neutral gas temp
	sph[ii].fionised=0;		//Is the particle ionised at all
	sph[ii].neighstorcont=0;	//Controle varible for building of neighstore
	sph[ii].neighstor.resize(200);	//Array containing references to all particles that consider this one a neighbour
	}

//Add sink propertys to sink particles
#pragma omp parallel for private(ii)
for (ii=0;ii<nos;ii++)
	{
	sph[N-nos+ii].ionised=ionised;	
	sph[N-nos+ii].neigh=neigh;
	sph[N-nos+ii].prob=prob;
	sph[N-nos+ii].photons=photons;
	sph[N-nos+ii].t=ti;						//Set stars to ionised gas temp for smoothing
	sph[N-nos+ii].ionised[ii]=1;					//Set the sink location in the ionised array to 1
	sph[N-nos+ii].sink=1;						//Is the particle a sink
	sph[N-nos+ii].x=ndata[ii]->r[0];				//Source x from gandalf
	sph[N-nos+ii].y=ndata[ii]->r[1];				//Source y from gandalf
	sph[N-nos+ii].z=ndata[ii]->r[2];				//Source z from gandalf
	sph[N-nos+ii].fionised=1;					//As this is source it is marked as ionised
	sph[N-nos+ii].neighstorcont=0;					//Controle varible for building of neighstore
	sph[N-nos+ii].neighstor.resize(200);				//Array containing references to all particles that consider this one a neighbour
	sinkid[ii]=N-nos+ii;						//Building sink id table
	ndot[ii]=pow(2.4e-24,2.)*Ndotmin/(4.*pi*2.6e-13)*scale; 	//Ndot table s^-1
	}

//Add controle sph to which all particles are initally linked
//Large values are used so it is never linked in the furture
sph[N].x=100000.;
sph[N].y=100000.;
sph[N].z=100000.;
sph[N].rho=1e100;
sph[N].sink=0;
sph[N].ionised=ionised;
sph[N].neigh=neigh;
sph[N].prob=prob;

//Debug message
if (debug==1){
cout<<"Particle arrays created"<<endl; 
}

//Find the closest source neighbour in chain for each particle
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int current_paricle_nn[maxneigh],Nneigb;		//Working particle neighbour list,Nneigb is the number of neighbours found by neighbour roughtine
int a;							//a is a controle so we dont check a particle twice
vector<double> angle; 					//Controle angle
double dot,mag,angletest; 				//holding variables
double distanceii,distancejj;				//Distance between the source of the test particle ii and the neighbour particle jj

//Begin neighbour find
#pragma omp parallel for private(current_paricle_nn,ii,jj,tt,Nneigb) //Initiate openmp
for (ii=0;ii<N;ii++) 				//Loop over all particles
	{
	
	//Is the particle a sink, if it is not we find the nearest NNnumber neighbours, if it is we only find those within h.
	if (sph[ii].sink==0)
		{
		Nneigb = sphneib->GetGatherNeighbourList(sphdata[ii].r,sphdata[ii].h*2.,sphgen,N,maxneigh,current_paricle_nn); //Find NNnumber nearest neighbours		
		if (Nneigb==-1){cout<<"Too many neighbours found 1 :("<<endl;}//Degug message	
		//Checks which sources lie within the smoothing length		
		for (tt=0;tt<nos;tt++){
			if(sqrt(pow((sph[sinkid[tt]].x-sph[ii].x),2)+pow((sph[sinkid[tt]].y-sph[ii].y),2)+pow((sph[sinkid[tt]].z-sph[ii].z),2))<=sph[ii].h*2.)
				{
				current_paricle_nn[Nneigb]=sinkid[tt];
				Nneigb=Nneigb+1;
				}
			}

		//For each of these neighbours write my id to there neighstorcont array		
		for(jj=0;jj<Nneigb;jj++)
			{
				//If there is still room in the neighstorecont array (This should always be true but this ensures no seg fault)
				if (sph[current_paricle_nn[jj]].neighstorcont<200)
					{
					sph[current_paricle_nn[jj]].neighstor[sph[current_paricle_nn[jj]].neighstorcont]=ii; //Write id
					sph[current_paricle_nn[jj]].neighstorcont+=1;					     //Advance id counter
					}
			}
		}
	else
		{
		Nneigb = sphneib->GetGatherNeighbourList(sphdata[ii].r,sphdata[ii].h,sphgen,N,maxneigh,current_paricle_nn); //Find NNnumber nearest neighbours		
		if (Nneigb==-1){cout<<"Too many neighbours found 2 :("<<endl;}//Degug message		
		
		//Checks which sources lie within the smoothing length		
		for (tt=0;tt<nos;tt++){
			if(sqrt(pow((sph[sinkid[tt]].x-sph[ii].x),2)+pow((sph[sinkid[tt]].y-sph[ii].y),2)+pow((sph[sinkid[tt]].z-sph[ii].z),2))<=sph[ii].h*2.)
				{
				current_paricle_nn[Nneigb]=sinkid[tt];  //Write id of source
				Nneigb=Nneigb+1;			//Advance number of neighbours
				}
			}

		//For each of these neighbours write my id to there neighstorcont array
		for(jj=0;jj<Nneigb;jj++)
			{
				//If there is still room in the neighstorecont array (This should always be true but this ensures no seg fault)
				if (sph[current_paricle_nn[jj]].neighstorcont<200)
					{
					sph[current_paricle_nn[jj]].neighstor[sph[current_paricle_nn[jj]].neighstorcont]=ii;	//Write id
					sph[current_paricle_nn[jj]].neighstorcont+=1;						//Advance id counter
					}
				else
					{
					//Debug error message to inform user that neighstore is full (This should not happen)
					cout<<"Neighstore full seg fault avoided, check code for mistakes"<<endl;
					}
			}
		}
	}

//Debug message
if (debug==1){
cout<<"neigbour step one compleate"<<endl;     
}

//Initiate openmp
#pragma omp parallel for private(current_paricle_nn,a,angle,dot,mag,angletest,jj,pp,distanceii,distancejj,ii,tt,Nneigb)
for (ii=0;ii<N;ii++) 				//Loop over all particles
	{
	a=0;	

	//If this particle has dummy neighbours it has not been checked so set to 1 to initiate check (This stops multiple checks)					
	for (pp=0;pp<nos;pp++)
		{
		if (sph[ii].neigh[pp]==N)
			{
			a=1;
			}
		}

	//If the neighbours are to be found (a=1)
	if (a==1)
		{

		angle.resize(nos); 					//Resize controle angle
		
		//Set controle angle to maximum for each source
		for (pp=0;pp<nos;pp++)
			{
			angle[pp]=2.*pi;
			}
		
		for (pp=0;pp<nos;pp++)		//For each source
			{
			//For each particle that considers me a neighbour
			for (jj=0;jj<sph[ii].neighstorcont;jj++)
				{
				//Work out the distances for both the test and candidate particle
				distanceii=sqrt(pow(sph[ii].x-sph[sinkid[pp]].x,2.)+pow(sph[ii].y-sph[sinkid[pp]].y,2.)+pow(sph[ii].z-sph[sinkid[pp]].z,2.));
				distancejj=sqrt(pow(sph[sph[ii].neighstor[jj]].x-sph[sinkid[pp]].x,2.)+pow(sph[sph[ii].neighstor[jj]].y-sph[sinkid[pp]].y,2.)+pow(sph[sph[ii].neighstor[jj]].z-sph[sinkid[pp]].z,2.));		
				//If the candidate particle is closer than the test particle it is a candidate (Also has controle so a particle cant be its own neighbour)
				if (distancejj<distanceii and ii!=sph[ii].neighstor[jj])
					{
					
					//Use the dot product to work out the angle between the conneting line and the candidate particle
					dot=((sph[ii].x-sph[sinkid[pp]].x)*(sph[sph[ii].neighstor[jj]].x-sph[sinkid[pp]].x))+((sph[ii].y-sph[sinkid[pp]].y)*(sph[sph[ii].neighstor[jj]].y-sph[sinkid[pp]].y))+((sph[ii].z-sph[sinkid[pp]].z)*(sph[sph[ii].neighstor[jj]].z-sph[sinkid[pp]].z));
					mag=(sqrt(pow(sph[ii].x-sph[sinkid[pp]].x,2.)+pow(sph[ii].y-sph[sinkid[pp]].y,2.)+pow(sph[ii].z-sph[sinkid[pp]].z,2.)))*(sqrt(pow(sph[sph[ii].neighstor[jj]].x-sph[sinkid[pp]].x,2.)+pow(sph[sph[ii].neighstor[jj]].y-sph[sinkid[pp]].y,2.)+pow(sph[sph[ii].neighstor[jj]].z-sph[sinkid[pp]].z,2.)));
					angletest=acos(dot/mag);

					//If the partcle is a sink set angle to be max (Stops non-relavant sinks becoming a neighbour)
					if (sph[sph[ii].neighstor[jj]].sink==1) angletest=2.*pi;

					//If the neighbour is the relavant source set angletest to be neg hence particle will be neighbour (Ensures link to sources)
					if (sph[ii].neighstor[jj]==sinkid[pp]) angletest=-1e50;
				
					//If candidate is closest so far to the connecting line then set current candidate to be the neighbour
					if (angletest<angle[pp])
						{
						angle[pp]=angletest;				//Set new comparison angle to be that of the neighbour
						sph[ii].neigh[pp]=sph[ii].neighstor[jj];	//Write particle id to neigh array
						}
					}
				}

			Nneigb = sphneib->GetGatherNeighbourList(sphdata[ii].r,sphdata[ii].h*2.,sphgen,N,maxneigh,current_paricle_nn); //Find NNnumber nearest neighbours		
			if (Nneigb==-1){cout<<"Too many neighbours found 3 :("<<endl;}//Degug message	

			//Checks which sources lie within the smoothing length		
			for (tt=0;tt<nos;tt++){
				if(sqrt(pow((sph[sinkid[tt]].x-sph[ii].x),2)+pow((sph[sinkid[tt]].y-sph[ii].y),2)+pow((sph[sinkid[tt]].z-sph[ii].z),2))<=sph[ii].h*2.)
					{
					current_paricle_nn[Nneigb]=sinkid[tt];
					Nneigb=Nneigb+1;
					}
				}

			//For each of my neighbours
			for (jj=0;jj<Nneigb;jj++)
				{
				//Work out the distances for both test and candidate particle
				distanceii=sqrt(pow(sph[ii].x-sph[sinkid[pp]].x,2.)+pow(sph[ii].y-sph[sinkid[pp]].y,2.)+pow(sph[ii].z-sph[sinkid[pp]].z,2.));
				distancejj=sqrt(pow(sph[current_paricle_nn[jj]].x-sph[sinkid[pp]].x,2.)+pow(sph[current_paricle_nn[jj]].y-sph[sinkid[pp]].y,2.)+pow(sph[current_paricle_nn[jj]].z-sph[sinkid[pp]].z,2.));		
				//If the candidate particle is closer than the test particle it is a candidate (Also has controle so a particle cant be its own neighbour)
				if (distancejj<distanceii and ii!=current_paricle_nn[jj])
					{
				
					//Use the dot product to work out the angle between the conneting line and the neighbour particle
					dot=((sph[ii].x-sph[sinkid[pp]].x)*(sph[current_paricle_nn[jj]].x-sph[sinkid[pp]].x))+((sph[ii].y-sph[sinkid[pp]].y)*(sph[current_paricle_nn[jj]].y-sph[sinkid[pp]].y))+((sph[ii].z-sph[sinkid[pp]].z)*(sph[current_paricle_nn[jj]].z-sph[sinkid[pp]].z));
					mag=(sqrt(pow(sph[ii].x-sph[sinkid[pp]].x,2.)+pow(sph[ii].y-sph[sinkid[pp]].y,2.)+pow(sph[ii].z-sph[sinkid[pp]].z,2.)))*(sqrt(pow(sph[current_paricle_nn[jj]].x-sph[sinkid[pp]].x,2.)+pow(sph[current_paricle_nn[jj]].y-sph[sinkid[pp]].y,2.)+pow(sph[current_paricle_nn[jj]].z-sph[sinkid[pp]].z,2.)));
					angletest=acos(dot/mag);

					//If the partcle is a sink set angle to be max (Stops non-relavant sinks becoming a neighbour)
					if (sph[current_paricle_nn[jj]].sink==1) angletest=2.*pi;

					//If the neighbour is the relavant source set angletest to be neg hence particle will be neighbour
					if (current_paricle_nn[jj]==sinkid[pp]) angletest=-1e50;
				
					//If Neighbour is closest so far to the connecting line then set current neighbour to be the neighbour
					if (angletest<angle[pp])
						{
						angle[pp]=angletest;				//Set new comparison angle to be that of the neighbour
						sph[ii].neigh[pp]=current_paricle_nn[jj];	//Write particle id to neigh array
						}
					}
				}
			}
		}
	} 

if (debug==1){
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

	//Open mp initalised for main iterations but disabled for final check to eliminate any race conditions (THIS CHECK WILL BE REMOVED IN FIANL VERSION)
	if (finalcheck==0)
		{
		#pragma omp parallel for schedule(dynamic)	//Begin open mp in dynamic mode i.e one argument is passed at a time rather than a block
		//For each particle
		for (ii=1;ii<N;ii++)
			{
			//Call fucntion to deturmine if test particle is ionised
			photoncount(sph,sinkid,ndot,N,nos,ii,change);
			}
		}	
        else
		{
		//cout<<"Final check run"<<endl;
		//For each particle
		//for (ii=1;ii<N;ii++)
		//	{	
		//	//Call fucntion to deturmine if test particle is ionised
		//	photoncount(sph,sinkid,ndot,N,nos,ii,change);
		//	}
		}	
	}
if (debug==1){
cout<<"Iterations compleate"<<endl; //Debug message
}

//Smooth the temperature
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double top,bottom,t,rad,s,w,invmu;
float radcont;
#pragma omp parallel for private(jj,ii,top,bottom,t,rad,s,w,radcont,invmu) //Initalise openmp
for (ii=0;ii<N;ii++)								//For each particle
	{
	top=0.;								//Intialise top (smoothing hold variable)
	radcont=1e5;
	for (jj=0;jj<sph[ii].neighstorcont;jj++)			//For each of the neighbours
		{
		if (sph[sph[ii].neighstor[jj]].fionised==1) t=ti;		//If the particle is ionised set temp to ionised temp
		if (sph[sph[ii].neighstor[jj]].fionised==0) t=tn;		//If the particle is not ionised set temp to neutral temp
		rad=(sqrt(pow(sph[sph[ii].neighstor[jj]].x-sph[ii].x,2.)+pow(sph[sph[ii].neighstor[jj]].y-sph[ii].y,2.)+pow(sph[sph[ii].neighstor[jj]].z-sph[ii].z,2.)));  //Ditance between particles
		s=rad/(sph[sph[ii].neighstor[jj]].h*2.);						//Work out s for smoothing kernal
		
		//Work out w for the kernal		
		if (s<1) w=1-(3./2.)*pow(s,2.)+(3./4.)*pow(s,3.);
		else if (s<2) w=(1./4.)*pow(2-s,3.);
		else w=0;
		
		//If the particle is the current closest ionised particle within the smoothing length then use to smooth
		if (t==ti and rad<radcont and s<2)
			{
			top=t*w;	//Work out the smoothed temperature
			radcont=rad;	//Set the radcont to the currect distance between particles
			}
		}
	
	//If there are no close ionised particles that satisfy the conditions then the particle is neutral else set the new particle temperature
	if (radcont==1e5)
		{
		sph[ii].t=tn;
		}
	else 
		{
		sph[ii].t=top;
		}
	
	//If the particle is ionised then its temperature must be the ionised temperature
	if (sph[ii].fionised==1)
		sph[ii].t=ti;
	
	//If the particles temp is less than the neutral temp because of smoothing set it back to the neutral temp
	if (sph[ii].t<tn)
		sph[ii].t=tn;
	
	invmu=(((sph[ii].t-tn)/mu_ion)+((ti-sph[ii].t)/mu_bar))/(ti-tn);//Work out corrected inverted mean gas particle mass						
	sph[ii].u=sph[ii].t/tempscale/gammam1*invmu;			//Work out the internal energy
	
	//Change the master sph paticle u if the temp should not be neutral	
	if(sph[ii].t!=tn)
		{
		sphdata[ii].u=sph[ii].u;				//Write new internal energy to gandalf
		}
	}
if (debug==1){
cout<<"Temperatures smoothed"<<endl;
}

float delta;
gettimeofday(&end, NULL);

delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;
cout<<"The time taken to calculate ionisation temperatures = "<<delta<<" s"<<endl;
}




template class MultipleSourceIonisation<1,GradhSphParticle>;
template class MultipleSourceIonisation<2,GradhSphParticle>;
template class MultipleSourceIonisation<3,GradhSphParticle>;
template class MultipleSourceIonisation<1,SM2012SphParticle>;
template class MultipleSourceIonisation<2,SM2012SphParticle>;
template class MultipleSourceIonisation<3,SM2012SphParticle>;
template class MultipleSourceIonisation<1,GodunovSphParticle>;
template class MultipleSourceIonisation<2,GodunovSphParticle>;
template class MultipleSourceIonisation<3,GodunovSphParticle>;

