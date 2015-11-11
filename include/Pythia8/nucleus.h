///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef NUCLEUS_H
#define NUCLEUS_H
#define protonMass   0.938272046212
#define neutronMass  0.939565378212
#define muonMass     0.105658371535
#define electronMass 0.000510998928
#define chargedPionMass 0.1395701835
#define unchargedPionMass 0.13497666
#define chargedKaonMass   0.49367716
#define unchargedKaonMass 0.49761424
#define nucleonMass (protonMass + neutronMass)/2.
//#include <cmath>
#include "Pythia.h"
#include "timer.h"
#include <vector>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
extern ofstream pathInNucleiOutput;
extern ofstream softCollisionsNumberOutput;
extern ofstream deltaPtOutput;
extern ofstream x1File;
extern  double getRandomFromFile();
//extern  double getRandom();
//extern  double ZMWLGaussIntegration5();
//extern  const gsl_rng_type *gslRandomGeneratorType;
//extern  gsl_rng *gslRandomGenerator;
extern ifstream coordinateFile;
extern int verbose;
extern double sNN;
extern double kEnergyLoss;
 //class Pythia8::Pythia;
//This class holds the information for a target nucleus
namespace Pythia8{
struct index {
  int i;
  int j;
} ;

class nucleus{
public:
	nucleus(const int Z, const int A);
	~nucleus(){};
	int    Z              () const { return _Z;                     }  ///< returns atomic number of nucleus
	int    A              () const { return _A;                     }  ///< returns nucleon number of nucleus

	double getNuclearDensity(const double r) const
	{ return (_A < 10) ? nuclearGaussDensity(r) : nuclearWoodSaxonDensity(r);}
	double renormalizedNuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;
	double renormalizedNuclearThicknessGauss5(double impactParameter, double zMin, double zMax)const;
	double ZMWLGaussIntegration5(double zMin, double zMax);
	void setPrecisionOfNuclearThicknessCalculation(double precision){_precisionOfNuclearThicknessCalculation = precision;}
	double getNuclearRadius  () const { return _r0 ; }
	void setInitialProjectileLabMomentum(double pz){_initialProjectileLabMomentum = pz;};
	double getInitialMomentum(void){return _initialProjectileLabMomentum;};
	int getId(void){return _id;};
	void setId(int id){_id = id;};
	//double getRandom(void){ return _pythia->rndm.flat();}
	std::string getNucleusName(){
		std::string	ElementName;

		switch (_Z) {
				case 4:
					ElementName = "Be";
				  break;
				case 26:
					ElementName = "Fe";
				  break;
				case 74:
					ElementName = "W";
				  break;
				case 1:
					if(_A == 1)	ElementName = "p";
					if(_A == 2)ElementName = "D";
				  break;
				case 7:
					ElementName = "N";
					break;
				case 36:
				  ElementName = "Kr";
				  break;
				default:
				  cout<<"Error no found such element. "<<endl;
				  ElementName = "";
				  break;
				}
				return ElementName;

	}

private:
	double woodSaxonSkinDepth() const { return 0.53;  }  ///< returns surface (0.53 fm for Au)
	double woodSaxonRadius() const { return 1.2 * pow(_A, 1. / 3.); }  ///< returns Wood-Saxon nuclear radius [fm] (Fermi model)

	double nuclearThicknessGauss5 (const double impactParameter, double zMin, double zMax) const;  ///< calculates nuclear thickness function for given distance b in impact parameter space (Eq. 4 in KN, PRC 60)
	double nuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;
	double getZMWL(double r);
	//double ZMWLGaussIntegration5(double zMin, double zMax);
	//	Pythia8::Pythia* _pythia;

	int    _Z;                      ///< atomic number of nucleus
	int    _A;                      ///< nucleon number of nucleus
	double _r0;				// nuclear radius
	double _rho0;
	double _precisionOfNuclearThicknessCalculation;
	double _initialProjectileLabMomentum;
	int _id;
	double nuclearWoodSaxonDensity(const double r) const
	{ return 3./4.*_A/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/(1+M_PIl*M_PIl*woodSaxonSkinDepth()*woodSaxonSkinDepth()/getNuclearRadius()/getNuclearRadius())/(1. + exp((r - getNuclearRadius()) / woodSaxonSkinDepth())); } ///< Wood-Saxon nuclear density
	double nuclearGaussDensity(const double r) const
	{ return exp(-r*r/getNuclearRadius()/getNuclearRadius())*_A/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/2*M_2_SQRTPIl;}
};

//Pythia* pythia;
class hardpingParticle : public Particle{
public:
	hardpingParticle():
		Particle(0),
		//
		//  Particle(particleId, 0, 0, 0, 0, 0, 0, 0, 0., 0., 0., 0., 0., 0., 9.),
		//
		_transferred4Momentum(0),
		_transferredCM4Momentum(0),
		_virtualPhotonEnergy(0),
		_x1(0),
		_x2(0),
		_hadronEnergyFraction(0),
		_motherParticleHistoryIndex(0),
		_hadronFormationLength(0),
		_history(0),
		_hardInteraction(false),
		_softInteraction(false),
		_outOfNucleus(false),
		_scatteringOnparticle(0),
		_numberOfCurrentGeneration(0),
		_indexInGeneration(0),
		_thetaHardping(0),
		_phiHardping(0),
		_thetaHardping2(0),
		_phiHardping2(0),
		_lastHard(false),
		_numberOfSoftCollisions(0),
		_numberOfHardCollisions(0),
		_preHadronFormationLength(0),
		_leftHadronFormationLength(0),
		_leftPreHadronFormationLength(0),
		_residualHadronFormationLength(0),
		_totalPathInNucleus(0),
		_energyLoss(0),
		_verbose(verbose),
		  _sinY1(0),
		  _cosY1(0),
		  _sinX2(0),
		  _cosX2(0)


	{
		_pythiaParticle = new Particle(0);
		_history = new vector <unsigned int>;
		_hadronNucleonCrossSection = 0;
		_preHadronNucleonCrossSection = 0;

	}
	hardpingParticle(int particleId):
	//	Particle(0),
		//
		  Particle(particleId, 211, 211, 211, 211, 211, 211, 211, 121., 111., 321., 222., 1., 3., 9.),
		  //Particle(particleId, 211, 211, 211, 211, 211, 211, 0, 0., 0., 0., 0., 0., 0., 9.),
		//
		_transferred4Momentum(0),
		_transferredCM4Momentum(0),
		_virtualPhotonEnergy(0),
		_x1(0),
		_x2(0),
		_hadronEnergyFraction(0),
		_motherParticleHistoryIndex(0),
		_hadronFormationLength(0),
		_history(0),
		_hardInteraction(false),
		_softInteraction(false),
		_outOfNucleus(false),
		_scatteringOnparticle(0),
		_numberOfCurrentGeneration(0),
		_indexInGeneration(0),
		_thetaHardping(0),
		_phiHardping(0),
		_thetaHardping2(0),
		_phiHardping2(0),
		_lastHard(false),
		_numberOfSoftCollisions(0),
		_numberOfHardCollisions(0),
		_preHadronFormationLength(0),
		_leftHadronFormationLength(0),
		_leftPreHadronFormationLength(0),
		_residualHadronFormationLength(0),
		_totalPathInNucleus(0),
		_energyLoss(0),
		_verbose(verbose),
		  _sinY1(0),
		  _cosY1(0),
		  _sinX2(0),
		  _cosX2(0)


	{
		_pythiaParticle = new Particle(0);
		_history = new vector <unsigned int>;
		_hadronNucleonCrossSection = 0;
		_preHadronNucleonCrossSection = 0;

	}
	//todo Написать конструктор, копирующий hardpingParticle;
	hardpingParticle(Particle &particle):
		Particle(particle),
		_transferred4Momentum(0),
		_transferredCM4Momentum(0),
		_virtualPhotonEnergy(0),
		_x1(0),
		_x2(0),
		_hadronEnergyFraction(0),
		_motherParticleHistoryIndex(0),
		_hadronFormationLength(0),
		_history(0),
		_hardInteraction(false),
		_softInteraction(false),
		_outOfNucleus(false),
		_scatteringOnparticle(0),
		_numberOfCurrentGeneration(0),
		_indexInGeneration(0),
		_thetaHardping(0),
		_phiHardping(0),
		_thetaHardping2(0),
		_phiHardping2(0),
		_lastHard(false),
		_numberOfSoftCollisions(0),
		_numberOfHardCollisions(0),
		_preHadronFormationLength(0),
		_leftHadronFormationLength(0),
		_leftPreHadronFormationLength(0),
		_residualHadronFormationLength(0),
		_totalPathInNucleus(0),
		_energyLoss(0),
		_verbose(verbose),
		  _sinY1(0),
		  _cosY1(0),
		  _sinX2(0),
		  _cosX2(0)

	{
		_pythiaParticle = new Particle(0);
		_history = new vector <unsigned int>;
		_hadronNucleonCrossSection = 0;
		double preHadronCrossSectionNucleon = 10;
		double preHadronCrossSectionPiMeson = 7;
		double preHadronCrossSectionKMeson = 5;

		switch (this->id()) {
				case 2212:

					_hadronNucleonCrossSection = 25;
					_preHadronNucleonCrossSection = preHadronCrossSectionNucleon;

				break;

				case 2112:
					_hadronNucleonCrossSection = 25;
					_preHadronNucleonCrossSection = preHadronCrossSectionNucleon;
					//this->e(neutronMass);
				break;

				case -2212:

					_hadronNucleonCrossSection = 25;
					_preHadronNucleonCrossSection = preHadronCrossSectionNucleon;

				break;

				case -2112:
					_hadronNucleonCrossSection = 25;
					_preHadronNucleonCrossSection = preHadronCrossSectionNucleon;
					//this->e(neutronMass);
				break;

				case 211:
					_hadronNucleonCrossSection = 15;
					_preHadronNucleonCrossSection = preHadronCrossSectionPiMeson;

				break;

				case -211:
					_hadronNucleonCrossSection = 15;
					_preHadronNucleonCrossSection = preHadronCrossSectionPiMeson;
				break;
/*
				case 111:
					_hadronNucleonCrossSection = 15;
					_preHadronNucleonCrossSection = 7;
				break;
*/
				case 321:
					_hadronNucleonCrossSection = 10;
					_preHadronNucleonCrossSection = preHadronCrossSectionKMeson;
				break;

				case -321:
					_hadronNucleonCrossSection = 10;
					_preHadronNucleonCrossSection = preHadronCrossSectionKMeson;

				break;
/*				case 130:

					_hadronNucleonCrossSection = 10;
					_preHadronNucleonCrossSection = 5;

				break;
				case 3110:

					_hadronNucleonCrossSection = 10;
					_preHadronNucleonCrossSection = 5;

				break;
*/
				default:
					_hadronNucleonCrossSection =  15; //_hadronNucleonCrossSection =  0;
					_preHadronNucleonCrossSection = 10;//;/_preHadronNucleonCrossSection = 0;
				break;
		}



	}

	~hardpingParticle(){
	//	delete _pythiaParticle;
	//	delete _history;
	}


	//todo написать деструктор
	vector <unsigned int> * getHistory(void){
		return _history;
	}
	Particle* getPythiaParticle(void){

		_pythiaParticle->id(this->id());
		_pythiaParticle->p(this->p());
		_pythiaParticle->vProd(this->vProd());
		return _pythiaParticle;
	}
	void setPythiaParticleStatus(int status){
		_pythiaParticle->status(status);
	}
	bool isHard(void){	return _hardInteraction;}
	bool isSoft(void){	return _softInteraction;}
	bool isOut(void){	return _outOfNucleus;}
	void setSoft(bool b = true){ _softInteraction = b; }
	void setHard(bool b = true){ _hardInteraction = b; }
	void setOut (bool b = true){ _outOfNucleus = b; }
	void scatteringOnParticle(int id){ _scatteringOnparticle = id ;}
	int  getIdscatteringParticle(void){ return _scatteringOnparticle;}
	double getMaxTransverseMomentum(int type);
	double getNewPtInitialState(int type);
	void setThetaHardping(double theta){_thetaHardping = theta;}
	void setPhiHardping(double phi){_phiHardping = phi;}
	double getThetaHardping(){return _thetaHardping;}
	double getPhiHardping(){return _phiHardping;}

	void setThetaHardping2(double theta){_thetaHardping2 = theta;}
	void setPhiHardping2(double phi){_phiHardping2 = phi;}
	double getThetaHardping2(){return _thetaHardping2;}
	double getPhiHardping2(){return _phiHardping2;}

	void setInitialProjectileLabMomentum(double initialProjectileLabMomentum){
		this->pz(initialProjectileLabMomentum);
	}
	double getInitialProjectileLabMomentum(void){
		return this->pz();
	}
	double getRestMass(int id = 0){

		if(id == 0)id =this->id();
		id = abs(id);
		switch (id) {
		case 2212:

		  return protonMass;
		  break;
		case 2112:
			return neutronMass;
			//this->e(neutronMass);
		  break;
		case 13:
			return muonMass;
		  break;
		case 11:
			return electronMass;

		  break;
		case 211:
			return chargedPionMass;

		  break;

		case 111:
			return unchargedPionMass;

		  break;

		case 321:
			return chargedKaonMass;

		  break;
		case 311:
			return unchargedKaonMass;

		  break;

		default:
			return 0;
		  break;
		}
	}

	double getRestMass2(void){

		double restMass = 0;
		restMass = getRestMass();
		return restMass*restMass;

	}

	std::string getParticleName(){
		std::string	ElementName;

			switch (this->id()) {
			case -13:
				ElementName = "mu+";
			  break;
			case 13:
				ElementName = "mu-";
			  break;
			case 11:
				ElementName = "e-";
			  break;
			case -11:
				ElementName = "e+";
			  break;

			default:
				if(_verbose)cout<<"Error no found such element. "<<endl;
			  ElementName = "";
			  break;
			}
			return ElementName;
	}

	void setAngles(void){


		_thetaHardping = this->theta();
		_phiHardping = this->phi();
		return;
		char ch;

/*		_thetaHardping = acos(this->pz()/this->pAbs());
		if(this->pT() == 0){
			_phiHardping = 0;
		}else{
			_phiHardping   = acos(this->py()/this->pT());
		//	cout<<"_phiHardping = "<<_phiHardping<<endl;
		//	cout<<"asin phi"<<asin(this->px()/this->pT())+M_PIl<<endl;
		}
*/
		/*	cout.precision(18);
		cout<<"this->theta() = "<<this->theta()<<" this->phi() = "<<this->phi()<<endl;
		cout<<"thetaHardping = "<<_thetaHardping<<" phi = "<<_phiHardping<<endl;
		cout<<"sin1 "<<sin(this->phi())<<" sin2 "<<sin(_phiHardping+M_PIl/2.)<<endl;
		cout<<"cos1 "<<cos(this->phi())<<" cos2 "<<cos(_phiHardping+M_PIl/2.)<<endl;

		cout<<"tsin1 "<<sin(this->theta())<<" tsin2 "<<sin(_thetaHardping)<<endl;
		cout<<"tcos1 "<<cos(this->theta())<<" tcos2 "<<cos(_thetaHardping)<<endl;
	*/	//cin>>ch;
		//double sinY1, cosY1, sinX2, cosX2;
		//double px,py,pz, px1,py1,pz1,px2,py2,pz2;
		double x,y,z, x1,y1,z1,x2,y2,z2;

		double px,py,pz, px1,py1,pz1,px2,py2,pz2;
		double sinY1, cosY1, sinX2, cosX2;
		px = this->px();
		py = this->py();
		pz = this->pz();
		x = this->vProd().px();
		y = this->vProd().py();
		z = this->vProd().pz();
		sinY1 = px/sqrt(px*px+pz*pz);
		cosY1 = pz/sqrt(px*px+pz*pz);

		_phiHardping = asin(sinY1);
		double P,R;
		P = sqrt(pz*pz+py*py+px*px);
		cout.precision(12);
		sinX2 = -py/P;
		cosX2 = sqrt(px*px+pz*pz)/P;
		_thetaHardping = asin(sinX2);
		setTrigonometricFunctions(sinY1,cosY1,sinX2,cosX2);
		if(_verbose)cout<<"_phiHardping "<<_phiHardping<<" _thetaHardping "<<_thetaHardping<<endl;
		if(_verbose)cout<<"_phi "<<this->phi()<<" _theta "<<this->theta()<<endl;


		if(_verbose) cout<<"setangles sinY1 "<<sinY1<<" cosY1 "<<cosY1<<" sinX2 "<<sinX2<<" cosX2 "<<cosX2<<endl;
		if(_verbose) 		cout<<"px = "<<px<<" py = "<<py<<" pz = "<<pz<<endl;
		if(_verbose) 		cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<endl;
		 		px1 = px*cosY1 - pz*sinY1;
		 		py1 = py;
		 		pz1 = px*sinY1 + pz*cosY1;

		 		x1  =  x*cosY1  - z*sinY1;
		 		y1  =  y;
		 		z1  =  x*sinY1  + z*cosY1;

		if(_verbose)		cout<<"px1 = "<<px1<<" py1 = "<<py1<<" pz1 = "<<pz1<<endl;
 		if(_verbose)		cout<<"x1 = "<<x1<<" y1 = "<<y1<<" z1 = "<<z1<<endl;

		 		px2 = px1;
		 		py2 = py1*cosX2  + pz1*sinX2;
		 		pz2 = -py1*sinX2 + pz1*cosX2;

		 		x2  = x1;
		 		y2  = y1*cosX2  + z1*sinX2;
		 		z2  = -y1*sinX2 + z1*cosX2;
		  //       DX03=DX02
		  //       DY03=DY02*DCOSX2+DZ02*DSINX2
		  //       DZ03=-DY02*DSINX2+DZ02*DCOSX2
		 if(_verbose)cout<<"px2 = "<<px2<<" py2 = "<<py2<<" pz2 = "<<pz2<<endl;
		 if(_verbose)cout<<"x2 = "<<x2<<" y2 = "<<y2<<" z2 = "<<z2<<endl;

	}
	double getEcm(){
		double mProjectile = 0, mTarget = 0;
		double s = 0;
		mProjectile = this->getRestMass();
		mTarget = this->getRestMass(this->getIdscatteringParticle());
		s = mProjectile*mProjectile + 2*mTarget*this->e() + mTarget*mTarget;

		return sqrt(s);
	}
//	double getRandom(void){ return pythia->rndm.flat();}

	double recalculateXBjorken(nucleus* targetNucleus,Pythia * py, double projectileLabEnergy){
		char ch;
		double kEnergyLossCM = 0;
		double newEnergy = 0;
		double oldEnergy = 0;
		double x1New = 0;
		double x1Old = 0;
		double x2New = 0;
		double x2Old = 0;
		double energyMin = 0;
		double energyCM  = 0;
		double minHatMass = 2;// todo дописанно руками. переписать по нормальному.

		//aqrtkEnergyLoss = _kEnergyLoss;
		if(_verbose)cout<<" kEnergyLoss1 = "<<kEnergyLoss<<endl;
	//	newEnergy = projectileLabEnergy - kEnergyLoss;
	//	oldEnergy = projectileLabEnergy;


		oldEnergy = this->p().e();;
		//oldEnergy = _initialParticle;//this->p().e();
		if(_verbose)cout<<" newEnergy = "<<newEnergy<<" oldEnergy = "<<oldEnergy<<endl;
		kEnergyLossCM = sqrt(2*nucleonMass*nucleonMass + 2*projectileLabEnergy*nucleonMass) - sqrt(2*nucleonMass*nucleonMass + 2*(projectileLabEnergy - kEnergyLoss)*nucleonMass);

		newEnergy = this->p().e() - kEnergyLossCM;

		if(_verbose)cout<<" kEnergyLoss = "<<kEnergyLossCM<<endl;
//cin>>ch;

		double tau = this->tau();
		//double y   = this->y();
		if(_verbose)cout<<" x1 "<<this->getXBjorkenProjectile()<<endl;
		if(_verbose)cout<<" x2 "<<this->getXBjorkenTarget()<<endl;

		double y   = 0.5*log(this->getXBjorkenProjectile()/this->getXBjorkenTarget());
		if(_verbose)cout<<" p rcx = "<<this->p();
		cout.precision(12);
		if(_verbose)cout<<"tau = "<<tau<<" y = "<<y<<endl;
	//	y = 0.5*log((this->p().e()+this->p().pz())/(this->p().e()-this->p().pz()));
	//	cout<<"tau = "<<tau<<" y = "<<y<<endl;

		x1New = sqrt(tau)*exp(y);
		x2New = sqrt(tau)*exp(-y);
		if(abs(x1New) > 1 || abs(x2New) > 1){
			cout<<"something, something dark side "<<endl;
			cin>>ch;
		}
		//cout<<"cme "<<py->info.eCM()<<endl;
	//	cin>>ch;
		energyCM = sNN;//this->getEcm();
		//cout<<"tau "<<tau<<" y "<<y<<endl;
		//cout<<" p "<<this->p();
		if(_verbose)cout<<"x1New "<<x1New<<" x2New "<<x2New<<" newEnergy "<<newEnergy<<" oldEnergy "<<oldEnergy<<endl;
		energyMin =0.5*minHatMass*minHatMass/(energyCM*x2New);
		if(_verbose)cout<<" energyMin "<<energyMin<<endl;
		//cin>>ch;
		double rLMax = 0;
		double rLMin = 0;
		double ZMW0  = 0;
		double tempRandom = 0;
		double randomLenght = 0;
		double norm = 0;
		if(targetNucleus->A() == 9){
			//cin>>ch;
			rLMax = 6.0;
			ZMW0  = 0.622;
		}else{
			rLMax = 11.0;
			ZMW0  = 0.284;
		}
		norm = targetNucleus->ZMWLGaussIntegration5(rLMin,rLMax);
		if(_verbose)cout<<"norm "<<norm<<endl;
		double lenghtValue[200], lenghtProbability[200];
		for(int i = 0; i < 201; i++){
			lenghtValue[i] = rLMin + (rLMax - rLMin)*i/200.;
			lenghtProbability[i] = targetNucleus->ZMWLGaussIntegration5(rLMin,lenghtValue[i])/norm;

			if(_verbose)cout<<"lenghtValue[i] "<<lenghtValue[i]<<" lenghtProbability[i] "<<lenghtProbability[i]<<endl;
		}
		//cin>>ch;
		int jLow = 0;
		int jUp = 201;
		int jMedium = 0;
		int jIndex = 0;
		int loopCount = 0;


		tempRandom = getRandomFromFile();//py->rndm.flat();
		if(tempRandom <= ZMW0){
			newEnergy = x1New*oldEnergy;
			if(_verbose)cout<<"first case "<<endl;
			if(_verbose)cout<<"newEnergy "<<newEnergy<<" x1New "<<x1New<<" oldEnergy "<<oldEnergy<<endl;
		//	if(_verbose)cout<<"tempRandom "<<tempRandom<<" ZMW0 "<<ZMW0<<" oldEnergy "<<oldEnergy<<endl;

		}else{

			do{
				jLow = 0;
				jUp = 201;
				tempRandom = getRandomFromFile();//py->rndm.flat();
				if(_verbose)cout<<"jIndex1 "<<jIndex<<endl;
				if(loopCount > 50){
					newEnergy = energyMin;
					break;
				}
				while(jUp - jLow > 1){
					jMedium = (jLow + jUp)/2.;
					if(_verbose)cout<<"jMedium begin = "<<jMedium<<endl;

					//tempRandom = 0.5993830;
				//	cout<<"tempRandomFF "<<tempRandom<<endl;
				//	cout<<"lenghtProbability[200] = "<<lenghtProbability[200]<<" lenghtProbability[0] "<<lenghtProbability[0]<<endl;
					if(_verbose)cout<<"tempRandom =  "<<tempRandom<<" lenghtProbability[jMedium] "<<lenghtProbability[jMedium]<<endl;
					if(/*(lenghtProbability[200] > lenghtProbability[0]) == */tempRandom > lenghtProbability[jMedium]){
						jLow = jMedium;
					}else{
						jUp  = jMedium;
					}
				}
				jIndex = jLow;
				if(_verbose)cout<<"jIndex2 "<<jIndex<<endl;

				if(jIndex < 0)jIndex = 0;
				if(jIndex >= 200) jIndex = 199;
				randomLenght = (lenghtValue[jIndex] + lenghtValue[jIndex+1])/2.;

				if(_verbose)cout<<"randomLenght "<<randomLenght<<endl;
				newEnergy = x1New*oldEnergy -kEnergyLossCM*randomLenght;
				if(_verbose)cout<<" x1 "<<x1New<<" oldEnergy "<<oldEnergy<<" randomLenght "<<randomLenght<<" kEnergyLossTemp = "<<kEnergyLossCM<<endl;
				loopCount++;
				if(_verbose)cout<<" newEnergy = "<<newEnergy<<" energyMin = "<<energyMin<<endl;
			}while(newEnergy < energyMin);


		}

		//cout<<" newEnergy "<<newEnergy<<" oldEnergy "<<oldEnergy<<endl;
		x1New = newEnergy/oldEnergy;
		//cout<<x1New<<" "<<randomLenght<<endl;
		x1File<<x1New<<" "<<randomLenght<<endl;
		if(_verbose)cout<<"newEnergy2 = "<<newEnergy<<" x1New "<<x1New<<" oldEnergy "<<oldEnergy<<endl;
	//	cin>>ch;
		if(x1New > 1){
			cout<<"x1New > 1"<<endl;
			cin>>ch;
		}
		this->setXBjorkenProjectileRecalculated(x1New);
		//todo не пересчитанны tau и быстрота
		return x1New;

	}
	void setAngles(double sinPhi,double cosPhi,double sinTheta,double cosTheta){
		_sinPhi   = sinPhi;
		_cosPhi   = cosPhi;
		_sinTheta = sinTheta;
		_cosTheta = cosTheta;
	}
	void getAngles(double &sinPhi,double &cosPhi,double &sinTheta,double &cosTheta){
		sinPhi   = _sinPhi;
		cosPhi   = _cosPhi;
		sinTheta = _sinTheta;
		cosTheta = _cosTheta;
	}
	void setTrigonometricFunctions(double sinY1,double cosY1,double sinX2,double cosX2){
		_sinY1 = sinY1;
		_cosY1 = cosY1;
		_sinX2 = sinX2;
		_cosX2 = cosX2;
	}

	void getTrigonometricFunctions(double &sinY1,double &cosY1,double &sinX2,double &cosX2){
		sinY1   = _sinY1;
		cosY1   = _cosY1;
		sinX2   = _sinX2;
		cosX2   = _cosX2;
	}
	void setAngles2(void){


		//	_thetaHardping = this->theta();
		//	_phiHardping = this->phi();
			char ch;

			double px,py,pz, px1,py1,pz1,px2,py2,pz2;

			double cosPhi = 0;
			double sinPhi = 0;
			double sinTheta = 0;
			double cosTheta = 0;
			px = this->px();
			py = this->py();
			pz = this->pz();

			double absoluteNucleonMomentum = 0;
			absoluteNucleonMomentum = this->pAbs();

			if(px == 0 && py == 0){
				cosPhi = 1;
				sinPhi = 0;
			}else{
				 cosPhi = py/sqrt(px*px+py*py);
				 sinPhi = px/sqrt(px*px+py*py);
			}

			 sinTheta = sqrt(px*px+py*py)/absoluteNucleonMomentum;
			 cosTheta = pz/absoluteNucleonMomentum;

			 setAngles(sinPhi,cosPhi,sinTheta,cosTheta);

			 _phiHardping2 = asin(sinPhi);
			 _thetaHardping2 = asin(sinTheta);
	 if(_verbose)cout<<"_phiHardping2 "<<_phiHardping2<<" _thetaHardping2 "<<_thetaHardping2<<endl;
	 if(_verbose)cout<<"_phi2 "<<this->phi()<<" _theta2 "<<this->theta()<<endl;

		}

	void rotateHardping(void){
		//if(_verbose)cout<<"thetaHardping = "<<_thetaHardping<<" phi = "<<_phiHardping<<endl;
		//dispose momentum of particle along initial beam direction
		char ch;

		this->rot(_thetaHardping,0);
		this->rot(0,_phiHardping);


		return ;

		double sinY1, cosY1, sinX2, cosX2;
		double px,py,pz, px1,py1,pz1,px2,py2,pz2;
		double x,y,z, x1,y1,z1,x2,y2,z2;
 		px = this->px();
		py = this->py();
		pz = this->pz();

 		x = this->vProd().px();
 		y = this->vProd().py();
 		z = this->vProd().pz();

		sinY1 = px/sqrt(px*px+pz*pz);
		cosY1 = pz/sqrt(px*px+pz*pz);
		double P,R;
		P = sqrt(pz*pz+py*py+px*px);
		cout.precision(12);
		sinX2 = -py/P;
		cosX2 = sqrt(px*px+pz*pz)/P;


/*
		sinY1 = sin(_phiHardping);
		cosY1 = cos(_phiHardping);
		sinX2 = sin(_thetaHardping);
		cosX2 = cos(_thetaHardping);*/
		getTrigonometricFunctions(sinY1,cosY1, sinX2,cosX2);
		if(_verbose)cout<<"111sinY1 "<<sinY1<<" cosY1 "<<cosY1<<" sinX2 "<<sinX2<<" cosX2 "<<cosX2<<endl;
		if(_verbose)cout<<"px = "<<px<<" py = "<<py<<" pz = "<<pz<<endl;
		if(_verbose)cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<endl;
		//cin>>ch;
		px1 = px;
		py1 = py*cosX2 - pz*sinX2;
		pz1 = py*sinX2 + pz*cosX2;

		x1  = x;
		y1  = y*cosX2  - z*sinX2;
		z1  = y*sinX2  + z*cosX2;

		if(_verbose)cout<<"px1 = "<<px1<<" py1 = "<<py1<<" pz1 = "<<pz1<<endl;
		if(_verbose)cout<<"x1 = "<<x<<" y1 = "<<y1<<" z1 = "<<z1<<endl;

		px2 = px1*cosY1 + pz1*sinY1;
		py2 = py1;
		pz2 = -px1*sinY1 + pz1*cosY1;

		x2  = x1*cosY1 + z1*sinY1;
		y2  = y1;
		z2  = -x1*sinY1 + z1*cosY1;
		if(_verbose)cout<<"px2 = "<<px2<<" py2 = "<<py2<<" pz2 = "<<pz2<<endl;
		if(_verbose)cout<<"x2 = "<<x2<<" y2 = "<<y2<<" z2 = "<<z2<<endl;
		this->px(px2);
		this->py(py2);
		this->pz(pz2);
 //		this->vProd().px(x2);
// 		this->vProd().py(y2);
// 		this->vProd().pz(z2);

		Vec4 vector4(x2,y2,z2,0);
		this->vProd(vector4);

		if(_verbose)cout<<"x3 = "<<this->vProd().px()<<" y3 = "<<this->vProd().py()<<" z3 = "<<this->vProd().pz()<<endl;
//		cin>>ch;
	//	cout<<"_thetaHardping = "<<_thetaHardping<<endl;
	//	cout<<"_phiHardping = "<<_phiHardping<<endl;
		//this->p().rotHardpingTest(_thetaHardping,0);
		//this->p().rotHardpingTest(0,_phiHardping);
	//	cout<<"0 rotation "<<this->p();


	//	rotateAroundX(_thetaHardping);
	//	cout<<"1 rotation "<<this->p();
	//	rotateAroundZ(_phiHardping);


	//	cout<<"2 rotation "<<this->p();
	//	this->rot(0,_thetaHardping);
	//	this->rot(_phiHardping,0);
	}
	void rotateBackHardping(void){
		//if(_verbose)cout<<"theta = "<<_thetaHardping<<" phi = "<<_phiHardping<<endl;
		//set particle momentum along z' direction (moving along z' direction, px = py = 0)
		char ch;
		this->rot(0,-_phiHardping);

		this->rot(-_thetaHardping,0);
		return ;

		//cout<<"0 rotation "<<this->p();
	//	rotateAroundX(_thetaHardping);

		//cout<<"1 rotation "<<this->p();
	//	rotateAroundZ(-_phiHardping);
		//this->setAngles();
		//cout<<"2 rotation "<<this->p();
		double sinY1, cosY1, sinX2, cosX2;
		double px,py,pz, px1,py1,pz1,px2,py2,pz2;
		double x,y,z, x1,y1,z1,x2,y2,z2;

 		px = this->px();
		py = this->py();
		pz = this->pz();

 		x = this->vProd().px();
 		y = this->vProd().py();
 		z = this->vProd().pz();
		sinY1 = px/sqrt(px*px+pz*pz);
		cosY1 = pz/sqrt(px*px+pz*pz);
		double P,R;
		P = sqrt(pz*pz+py*py+px*px);
		cout.precision(12);
		sinX2 = -py/P;
		cosX2 = sqrt(px*px+pz*pz)/P;

/*
		sinY1 = sin(_phiHardping);
		cosY1 = cos(_phiHardping);
		sinX2 = sin(_thetaHardping);
		cosX2 = cos(_thetaHardping);
*/
		getTrigonometricFunctions(sinY1,cosY1, sinX2,cosX2);

		if(_verbose)cout<<"sfinY1 "<<sinY1<<" cosY1 "<<cosY1<<" sinX2 "<<sinX2<<" cosX2 "<<cosX2<<endl;
	//	cin>>ch;
		if(_verbose)cout<<"px = "<<px<<" py = "<<py<<" pz = "<<pz<<endl;
		if(_verbose)cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<endl;
	//	cin>>ch;
		px1 = px*cosY1 - pz*sinY1;
		py1 = py;
		pz1 = px*sinY1 + pz*cosY1;

		x1  =  x*cosY1  - z*sinY1;
		y1  =  y;
		z1  =  x*sinY1  + z*cosY1;

		if(_verbose)cout<<"px1 = "<<px1<<" py1 = "<<py1<<" pz1 = "<<pz1<<endl;
		if(_verbose)cout<<"x1 = "<<x1<<" y1 = "<<y1<<" z1 = "<<z1<<endl;

		px2 = px1;
		py2 = py1*cosX2  + pz1*sinX2;
		pz2 = -py1*sinX2 + pz1*cosX2;

		x2  = x1;
		y2  = y1*cosX2  + z1*sinX2;
		z2  = -y1*sinX2 + z1*cosX2;
 //       DX03=DX02
 //       DY03=DY02*DCOSX2+DZ02*DSINX2
 //       DZ03=-DY02*DSINX2+DZ02*DCOSX2
		if(_verbose)cout<<"px2 = "<<px2<<" py2 = "<<py2<<" pz2 = "<<pz2<<endl;
		if(_verbose)cout<<"x2 = "<<x2<<" y2 = "<<y2<<" z2 = "<<z2<<endl;
		this->px(px2);
		this->py(py2);
		this->pz(pz2);
		Vec4 vector4(x2,y2,z2,0);
		this->vProd(vector4);
 		//this->vProd().px(x2);
 		//this->vProd().py(y2);
 		//this->vProd().pz(z2);
		if(_verbose)cout<<"x3 = "<<this->vProd().px()<<" y3 = "<<this->vProd().py()<<" z3 = "<<this->vProd().pz()<<endl;
	//	cin>>ch;



	}



	void rotateAroundX(double thetaIn, double phiIn = 0){
		// phi assumed equal zero//
	  double xx = this->px();
	  double yy = this->py();
	  double zz = this->pz();
	  double tt = this->e();

//	  cout<<"xx = "<<xx<<" yy = "<<yy<<" zz = "<<zz<<endl;
	  double cthe = cos(thetaIn);
	  double sthe = sin(thetaIn);
	  double cphi = cos(phiIn);
	  double sphi = sin(phiIn);
	  double tmpx =   cphi* xx  +  cthe* sphi * yy + sthe * sphi * zz;
	  double tmpy =  -sphi* xx  +  cphi* cthe * yy + sthe * cphi * zz;
	  double tmpz =  -sthe* yy  +  cthe* zz;
	  xx          = tmpx;
	  yy          = tmpy;
	  zz          = tmpz;
	  Vec4 tempVector(xx, yy, zz, tt);
	  this->p(tempVector);
	//  cout<<"xx2 = "<<xx<<" yy2 = "<<yy<<" zz2 = "<<zz<<endl;
	}
	void rotateAroundZ(double phiIn, double thetaIn = 0){
		// theta assumed equal zero//
	  double xx = this->px();
	  double yy = this->py();
	  double zz = this->pz();
	  double tt = this->e();

	  //cout<<"xx = "<<xx<<" yy = "<<yy<<" zz = "<<zz<<endl;
	  double cthe = cos(thetaIn);
	  double sthe = sin(thetaIn);
	  double cphi = cos(phiIn);
	  double sphi = sin(phiIn);
	  double tmpx =   cphi* xx  +  cthe* sphi * yy + sthe * sphi * zz;
	  double tmpy =  -sphi* xx  +  cphi* cthe * yy + sthe * cphi * zz;
	  double tmpz =  -sthe* yy  +  cthe* zz;
	  xx          = tmpx;
	  yy          = tmpy;
	  zz          = tmpz;
	  Vec4 tempVector(xx, yy, zz, tt);
	  this->p(tempVector);
	  //cout<<"xx2 = "<<xx<<" yy2 = "<<yy<<" zz2 = "<<zz<<endl;
	}
	void setInexNumber(unsigned int i){
		_indexInGeneration = i;
	}
	unsigned int getIndexNumber(void){
		return _indexInGeneration;
	}
	void setNumberOfCurrentGeneration(unsigned int i){
		_numberOfCurrentGeneration = i;
	}
	unsigned int getNumberOfCurrentGeneration(void){
		return _numberOfCurrentGeneration;
	}
	double getXBjorkenProjectile(void){
		return _x1;
	}
	void setXBjorkenProjectile(double x1){
		 _x1 = x1;
	}
	double getXBjorkenProjectileRecalculated(void){
		return _x1Recalculated;
	}
	void setXBjorkenProjectileRecalculated(double x1){
		 _x1Recalculated = x1;
	}
	double getXBjorkenTarget(void){
		return _x2;
	}
	void setXBjorkenTarget(double x2){
		 _x2 = x2;
	}
	Vec4 getTransferred4Momentum(void){
		return _transferred4Momentum;
	}
	void setTransferred4Momentum(Vec4 transferred4Momentum){
		_transferred4Momentum = transferred4Momentum;
	}
	Vec4 getTransferredCM4Momentum(void){
		return _transferredCM4Momentum;
	}
	void setTransferredCM4Momentum(Vec4 transferred4Momentum){
		_transferredCM4Momentum = transferred4Momentum;
	}
	double getVirtualPhotonEnergy(void){
		return _virtualPhotonEnergy;
		//return _transferred4Momentum.e();
	}
	void setVirtualPhotonEnergy(double virtualPhotonEnergy){
		_virtualPhotonEnergy = virtualPhotonEnergy;
	}
	double getAbsQ2(void){
		//return abs(_transferred4Momentum*_transferred4Momentum);
		return abs(_transferred4Momentum.m2Calc());
	}
	double getQ2(void){
		return _transferred4Momentum.m2Calc();
	}
	unsigned int getMotherParticleHistoryIndex(void){
		return _motherParticleHistoryIndex;
	}
	void setMotherParticleHistoryIndex(unsigned int index){
		_motherParticleHistoryIndex = index;
	}

	double getHadronNucleonCrossSection(void){
		return _hadronNucleonCrossSection;
	}
	void setHadronParticleNucleonCrossSection(double crossSection){
		_hadronNucleonCrossSection = crossSection;
	}

	double getPreHadronNucleonCrossSection(void){
		return _preHadronNucleonCrossSection;
	}
	void setPreHadronNucleonCrossSection(double crossSection){
		_preHadronNucleonCrossSection = crossSection;
	}
	double getQuarkNucleonCrossSection(void){
		return _quarkNucleonCrossSection;
	}
	void setQuarkNucleonCrossSection(double crossSection){
		_quarkNucleonCrossSection = crossSection;
	}
	double getHadronEnergyFraction(void){
		return _hadronEnergyFraction;
	}
	void setHadronEnergyFraction(double fraction){
		_hadronEnergyFraction = fraction;
	}
	double getHadronFormationLength(void){
		return _hadronFormationLength;
	}
	void setHadronFormationLength(double formationLength){
		_hadronFormationLength = formationLength;
	}
	double getPreHadronFormationLength(void){
		return _preHadronFormationLength;
	}
	void setPreHadronFormationLength(double preHadronFormationLength){
		_preHadronFormationLength = preHadronFormationLength;
	}
	bool isPreviosCollisionHard(void){
		return _lastHard;
	}
	void setPreviosCollisionHard(bool lastHard = true){
		_lastHard = lastHard;
	}
	void increaseSoftCollisionNumber(void){
		_numberOfSoftCollisions++;
	}
	void increaseHardCollisionNumber(void){
		_numberOfHardCollisions++;
	}
	void decreaseSoftCollisionNumber(void){
		_numberOfSoftCollisions--;
	}
	void decreaseHardCollisionNumber(void){
		_numberOfHardCollisions--;
	}
	unsigned int getSoftCollisionNumber(void){
		return _numberOfSoftCollisions;
	}
	unsigned int getHardCollisionNumber(void){
		return _numberOfHardCollisions;
	}
	void setSoftCollisionNumber(unsigned int number){
		_numberOfSoftCollisions = number;
	}
	void setHardCollisionNumber(unsigned int number){
		_numberOfHardCollisions = number;
	}
	void setTotalPathInNucleus(double path){
		_totalPathInNucleus = path;
	}
	double getTotalPathInNucleus(void){
		return _totalPathInNucleus;
	}
	double getLeftHadronFormationLength(void){
		return _leftHadronFormationLength;
	}
	void setLeftHadronFormationLength(double formationLength){
		_leftHadronFormationLength = formationLength;
	}
	double getLeftPreHadronFormationLength(void){
		return _leftPreHadronFormationLength;
	}
	void setLeftPreHadronFormationLength(double formationLength){
		_leftPreHadronFormationLength = formationLength;
	}
	double getResidualHadronFormationLength(void){
		return _residualHadronFormationLength;
	}
	void setResidualHadronFormationLength(double formationLength){
		_residualHadronFormationLength = formationLength;
	}
	double getEnergyLoss(void){
		return _energyLoss;
	}
	void setEnergyLoss(double formationLength){
		_energyLoss = formationLength;
	}

	double getKEnergyLoss(void){
		return _kEnergyLoss;
	}
	void setKEnergyLoss(double energyLoss){
		_kEnergyLoss = energyLoss;
	}
	//todo написать функцию, которая по индексу частицы в истории возвращала бы указатель на эту частицу.
private:
	vector <unsigned int> * _history;
	bool _hardInteraction;
	bool _softInteraction;
	bool _outOfNucleus;
	int _scatteringOnparticle;
	double _thetaHardping;
	double _phiHardping;
	double _thetaHardping2;
	double _phiHardping2;
	unsigned int _indexInGeneration;// i_init - index number of particle in current generation
	unsigned int _numberOfCurrentGeneration;
	Vec4 _transferred4Momentum;
	Vec4 _transferredCM4Momentum;
	double _virtualPhotonEnergy;
	double _x1;
	double _x1Recalculated;
	double _x2;
	double _hadronNucleonCrossSection;
	double _preHadronNucleonCrossSection;
	double _hadronEnergyFraction;
	unsigned int _motherParticleHistoryIndex;
	double _hadronFormationLength;
	double _leftHadronFormationLength;
	double _preHadronFormationLength;
	double _leftPreHadronFormationLength;
	double _residualHadronFormationLength;
	bool   _lastHard;
	double _quarkNucleonCrossSection;
	unsigned int _numberOfSoftCollisions;
	unsigned int _numberOfHardCollisions;
	Particle* _pythiaParticle;
	double _totalPathInNucleus;
	double _energyLoss;
	double _cosTheta;
	double _sinTheta;
	double _cosPhi;
	double _sinPhi;
	int _verbose;
	double _sinY1, _cosY1, _sinX2, _cosX2;
	double _kEnergyLoss;
	//Rndm * _random;

};
class producedParticlesInformation{
public:
	producedParticlesInformation():nProducedParticles(NULL),nParticleInFinalState(NULL),nParticleForNewGeneration(NULL){
		nProducedParticles = new vector <int>;
		nParticleInFinalState = new vector <int>;
		nParticleForNewGeneration  = new vector <int>;

	}
	~producedParticlesInformation(){}
	int getNProducedParticles(int nGeneration){
		return nProducedParticles->at(nGeneration);
	}
	void setNProducedParticles(int nGeneration, int nParticles){
		nProducedParticles->at(nGeneration) = nParticles;
	}
	int getNParticleInFinalState(int nGeneration){
		return nParticleInFinalState->at(nGeneration);
	}
	void setNParticleInFinalState(int nGeneration, int nParticles){
		nParticleInFinalState->at(nGeneration) = nParticles;
	}
	int getNParticleForNewGeneration(int nGeneration){
		return nParticleForNewGeneration->at(nGeneration);
	}
	void setNParticleForNewGeneration(int nGeneration, int nParticles){
		nParticleForNewGeneration->at(nGeneration) = nParticles;
	}
	void increaseNParticleInFinalState(int nGeneration){nParticleInFinalState->at(nGeneration)++;}
	std::vector <int> * getVectorNProducedParticlesPointer(){return nProducedParticles;}
	std::vector <int> * getVectorNParticleInFinalStatePointer(){return nParticleInFinalState;}
	std::vector <int> * getVectorNParticleForNewGenerationPointer(){return nParticleForNewGeneration;}
protected:

	std::vector <int> *nProducedParticles;
	std::vector <int> *nParticleInFinalState;
	std::vector <int> *nParticleForNewGeneration;
};
class vectorHardping : protected Vec4{
public:
	/*vectorHardping(){
		cout<<"in constructor vectorHardping 0 "<<endl;
		Vec4 * pVec = new Vec4(0);
		setVector(*pVec);
		//setVector(pVec);
		_idParticle = 0 ;
		delete pVec;
		cout<<"here!!!"<<endl;
	}*/
	vectorHardping(int id = 0):_idParticle(id),_vector4Hardping(0),_theta(0),_phi(0),_impactParameterX(0),_impactParameterY(0){}
	vectorHardping(Vec4 & vector, int id){
		//cout<<"in constructor vectorHardping 1"<<endl;
		setVector (vector);
		_idParticle = id;
	}
	vectorHardping( Particle * particle):
		_idParticle(particle->id()),
		_phi(particle->phi()),
		_theta(particle->theta()),
		_vector4HardpingProduction(NULL)
	{
		//particle->
	//	cout<<"in constructor vectorHardping 2"<<endl;
		Vec4 * pVec = new Vec4(particle->p());
		setVector(*pVec);
		Vec4 * pVecProd = new Vec4(particle->vProd());
		setVectorProduction(*pVecProd);
		delete pVec;

	}
	~vectorHardping(){
		delete _vector4Hardping;
	}
	void setVector(Vec4 &v){
		//cout<<"in set vector"<<endl;
	//	delete _vectorHardping;
		_vector4Hardping = new Vec4(v);
	}
	void setVectorProduction(Vec4 &v){
		_vector4HardpingProduction = new Vec4(v);
	}
	void setHardpingVector(vectorHardping &v){
	//	cout<<"in setHardpingVector"<<endl;
	//	delete _vectorHardping;
		//v.getVector()
		_vector4Hardping = new Vec4(*(v.getVector()));
		_idParticle = v.getId();
		_impactParameterX = v.getImpactParameterX();
		_impactParameterY = v.getImpactParameterY();
		_theta = v.getTheta();
		_phi = v.getPhi();


	}
	Vec4 * getVector(void){
		return _vector4Hardping;
	}
	void setId(int id){_idParticle = 1;}
	int getId(void){return _idParticle;}
	void setTheta(double theta){_theta = theta;}
	double getTheta(void){return _theta;}
	void setPhi(double phi){_phi = phi;}
	double getPhi(void){return _phi;}
	void setImpactParameter(double bx, double by){
		_impactParameterX = bx;
		_impactParameterY = by;
	}
	double getAbsImpactParameter(void){
		return sqrt(_impactParameterX*_impactParameterX+_impactParameterY*_impactParameterY);
	}
	double getImpactParameterX(void){
		return _impactParameterX;
	}
	double getImpactParameterY(void){
		return _impactParameterY;
	}

//	vector->

private:
	Vec4 * _vector4Hardping;
	Vec4 * _vector4HardpingProduction;
	int _idParticle;
	double _theta;
	double _phi;
	double _impactParameterX;
	double _impactParameterY;

};

class generation{
	typedef std::vector< std::vector<hardpingParticle> > vectorConstruction;
public:
	generation(const int numberOfGeneration, int iSize, int jSize):
	_iSize(iSize),
	_jSize(jSize),
	//_matrix(0),
	_numberOfGeneration(numberOfGeneration)
	{
		/*
		_matrix = vectorConstruction(0);
		if(_matrix.size() == 0){
			cout<< "_matrix is not initialized. Initialization matrix ..."<<endl;
			_matrix = vectorConstruction(_iSize,std::vector<Vec4>(_jSize));
		}else{
			cout<< "matrix already initialized"<<endl;
		}
		cout<<"matrix initialized "<<endl;
		*/
		initializeNewMatrix();
		//cout<<

	};
	generation():
		_iSize(0),
		_jSize(0),
		//_matrix(0),
		_numberOfGeneration(0)
		{
			/*
			_matrix = vectorConstruction(0);
			if(_matrix.size() == 0){
				cout<< "_matrix is not initialized. Initialization matrix ..."<<endl;
				_matrix = vectorConstruction(_iSize,std::vector<Vec4>(_jSize));
			}else{
				cout<< "matrix already initialized"<<endl;
			}
			cout<<"matrix initialized "<<endl;
			*/
			initializeNewMatrix();
			//cout<<

		};
	~generation(){};

	vectorConstruction* initializeNewMatrix(){ //return pointer to the matrix
		//cout<<"matrix size "<<_matrix.size()<<endl;
		if(_matrix.size() == 0){
		//	cout<< "_matrix is not initialized. Initialization matrix ..."<<endl;
			_matrix = vectorConstruction(_iSize,std::vector<hardpingParticle>(_jSize));

		}else{
			cout<< "matrix already initialized"<<endl;
		}

		/*for (unsigned int i =0; i< _matrix.size(); i++){
			for (unsigned int j =0; j< _matrix[i].size(); j++){
				_matrix[i][j] = 0;
			}

		}*/
		//cout<<"matrix initialized "<<endl;

		return &_matrix;
	}
	vectorConstruction* getMatrix(){ //return pointer to the matrix
				return &_matrix;
	}
	int* initialize_int(){
		cout<<"int initialized "<<endl;

		return &_iSize;
	}
	/*
	void printMatrix(void){
		for(int i = 0; i < _matrix.size(); i++){
			for(int j=0;j< _matrix[i].size();j++){
				cout<<(_matrix[i][j])<<" ";
			}
			cout<<endl;
		}
	}
*/
	//vectorConstruction matrix(10,std::vector <int>(10));
	//std::vector< std::vector< int > >  m(1, std::vector<int>(1));
private:
	vectorConstruction _matrix;
	int _iSize;
	int _jSize;
	int _numberOfGeneration;
	Vec4 tempVec4;
};// end class generation

class Pythia;
class Hardping
{
	//Pythia8::Pythia * pythia;

public:

	Hardping(nucleus projectileNucleus, nucleus targetNucleus);
	Hardping(hardpingParticle incideniParticle, nucleus targetNucleus);

	~Hardping();
 
//	int    Z              () const { return _Z;                     }  ///< returns atomic number of nucleus
//	int    A              () const { return _A;                     }  ///< returns nucleon number of nucleus
//	double woodSaxonRadius() const { return 1.2 * pow(_A, 1. / 3.); }  ///< returns Wood-Saxon nuclear radius [fm] (Fermi model)
/*	double nuclearThicknessGauss5 (const double impactParameter, double zMin, double zMax) const;  ///< calculates nuclear thickness function for given distance b in impact parameter space (Eq. 4 in KN, PRC 60)
	double nuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;
	double renormalizedNuclearThicknessGauss12(double impactParameter, double zMin, double zMax)const;
	double renormalizedNuclearThicknessGauss5(double impactParameter, double zMin, double zMax)const;*/
	double probabilityOfNSoftScattering(double impactParameter, double zMin, double zMax, int numberOfCollisions)const;
	double integralProbabilityOfNSoftScattering(int numberOfCollisions) const;
	double rho0() const { return _rho0; }
	double getr0() const { return _r0;}
	int factorial(int i)const{return _factorials[i];}
	//void setPrecisionOfNuclearThicknessCalculation(double precision){_precisionOfNuclearThicknessCalculation = precision;}
	void randomGeneratorInitialize(int seed){
		_randomSeed = seed;
		//_random = new Rndm(_randomSeed);
		_random->init(seed);
	}
	double getRandom(void){ return pythia->rndm.flat();}
/*	double nuclearWoodSaxonDensity(const double r) const
	{ return 3./4.*_targetNucleus.A()/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/(1+M_PIl*M_PIl*woodSaxonSkinDepth()*woodSaxonSkinDepth()/getNuclearRadius()/getNuclearRadius())/(1. + exp((r - getNuclearRadius()) / woodSaxonSkinDepth())); } ///< Wood-Saxon nuclear density
	double nuclearGaussDensity(const double r) const
	{ 	//std::cout<< "i ma in nuclearGaussDensity"<<std::endl;
		return exp(-r*r/getNuclearRadius()/getNuclearRadius())*_targetNucleus.A()/M_PIl/getNuclearRadius()/getNuclearRadius()/getNuclearRadius()/2*M_2_SQRTPIl;}
	double getNuclearDensity(const double r) const
	{ return (_targetNucleus.A() < 10) ? nuclearGaussDensity(r) : nuclearWoodSaxonDensity(r);}*/
	double pathInNucleus(hardpingParticle * particleA);
	bool pathInNucleus2(hardpingParticle * particleA, double &zCoordinateOfCollision);
	Pythia8::Pythia * pythia;
	double getAbsPartonMomentum(void){
		_absPartonMomentum = sqrt(_xPartonMomentum*_xPartonMomentum+_yPartonMomentum*_yPartonMomentum+_zPartonMomentum*_zPartonMomentum);
		return _absPartonMomentum;
	}
	double getNewPtInitialState(hardpingParticle* particleA, int type);
//	Pythia* pythia;
//	Pythia* pythiaHard;
	Pythia8::generation * _generation;
	void initGeteration(int nOfGeneration,int nInitialParcicle){
		_generation = new generation(nOfGeneration,nInitialParcicle,0);
	}
	int pythiaInitialization(hardpingParticle * particleA ,hardpingParticle * particleB);
	void hardping(void);
	double getMaxTransverseMomentum(hardpingParticle* particleA,int type);
	void setEnergy(double Elab){
		_energyLab = Elab;
	}
	unsigned int calculateNumberOfParticlesAtPreviousGeneration(unsigned int numberOfGeneration){
		unsigned int firstWave = 0;
		unsigned int isumMax = 0;
		isumMax = (numberOfGeneration) ? _generations->at(numberOfGeneration-1).getMatrix()->size() : 0;
		if(_verbose)if(!_firstCall)cout<<" number of incident particle in previous generation "<<isumMax<<endl;
		////////////// calculating number of particles at previous generation/////////////////////////
		for(int isum = 0; isum < isumMax ;isum++){
			//cout<<" in summ !!"<<endl;
			if(_verbose)cout<<" in Generation "<<numberOfGeneration-1<<" at i = "<<isum<< " size equal = "<<_generations->at(numberOfGeneration-1).getMatrix()->at(isum).size()<<endl;
			firstWave +=  _generations->at(numberOfGeneration-1).getMatrix()->at(isum).size();

		}
		return firstWave;
	}
	void particleIndexation(unsigned int numberOfGeneration,unsigned int numberOfParticlesAtPreviousGeneration){

		unsigned int isumMax = (numberOfGeneration) ? _generations->at(numberOfGeneration-1).getMatrix()->size() : 0;
		_index->resize(numberOfParticlesAtPreviousGeneration);
		int	indexCount = 0;


		// indexation //////////////////////////////////////
		for(int i_index = 0; i_index < isumMax ;i_index++){
			for(unsigned int j_index = 0; j_index < _generations->at(numberOfGeneration-1).getMatrix()->at(i_index).size() ;j_index++){
				//cout<<" in summ !!"<<endl;

				_index->at(indexCount).i = i_index;
				_index->at(indexCount).j = j_index;
				indexCount++;
			}

		}
	}
	double getEnergy(void){return _energyLab;}
	void deleteParticleFromSoftCollision(void){
		if(_hardInteractionSummaryCount && _softToHard){


			//change 22.07.14
			//_softInteractionCount--;
			//change 22.07.14

			_softToHard = false;
			_vertexOfBadHardInteraction->resize(0);
			//todo delete particle from corresponding soft collisions
			if(_indexSoftToHard->size())do{
				if(_verbose)cout<<" _indexSoftToHard back=  "<<_indexSoftToHard->back()<<" _indexSoftToHard->size() "<<_indexSoftToHard->size()<<endl;
				_outOfNucleus->erase(_outOfNucleus->begin()+_indexSoftToHard->back());
				_indexSoftToHard->pop_back();
			}while(_indexSoftToHard->size()!= 0);



			if(_indexSoftToHardBadInit->size())do{
				if(_verbose)cout<<" _indexSoftToHard back=  "<<_indexSoftToHardBadInit->back()<<" _indexSoftToHardBadInit->size() "<<_indexSoftToHardBadInit->size()<<endl;
				_notInit->erase(_notInit->begin()+_indexSoftToHardBadInit->back());
				_indexSoftToHardBadInit->pop_back();
			}while(_indexSoftToHardBadInit->size()!= 0);

		}
	}
	std::vector <generation> *_generations;
	std::vector <index> * _index;
	std::vector <hardpingParticle> * _finalState;
	std::vector <hardpingParticle> * _notInit;
	std::vector <hardpingParticle> * _cutMassive;
	std::vector <hardpingParticle> * _outOfNucleus;
	std::vector <unsigned int>* _indexBadInitializations;
	std::vector <unsigned int>* _indexSoftToHard;
	std::vector <unsigned int>* _indexSoftToHardBadInit;
	std::vector <unsigned int>* _vertexOfInteraction;
	std::vector <unsigned int>* _vertexOfBadHardInteraction;
	double _impactParameterMax;
	double _impactParameterMin;
	hardpingParticle _initialParticle;
	void addGeneration(/*	Pyt hia8::generation * gen*/){
		//generation gener(1,1,1);
		if(_verbose)cout<<"i am in addGeneration"<<endl;
		if(_generation==NULL){
			cerr<<"pointer '_generation' is not initialized. nothing happened "<<endl;
			return;
		}
		_generations->push_back(*_generation);
		return;
	}
	int energyLoss(hardpingParticle* particleA, double zCoordinateOfCollisions){
		char ch;
		// procedure return 0 if particle is adsorbed and 1 if not
//		zCoordinateOfCollisions - this is a coordinate of particle up to which particle loses energy
		bool isNucleon = false;
		if(particleA->idAbs() == 2212 || particleA->idAbs() == 2112)isNucleon = true;
		unsigned int i_init = particleA->getIndexNumber();
		if(_verbose)cout<<" particleA p "<<particleA->p()<<endl;
		if(_verbose)cout<<" zCoordinateOfCollisions "<<zCoordinateOfCollisions<<endl;
		if(_verbose)cout<<" current zcoord "<<particleA->vProd().pz()<<endl;

		double passedHadronFormathionLenght = 0;

		if(_verbose)cout<<particleA->getLeftHadronFormationLength()<<endl;
		if(particleA->getLeftHadronFormationLength() > 0 && particleA->getSoftCollisionNumber() == 0){
			//cout<<"resid 1 "<<particleA->getResidualHadronFormationLength()<<endl;
			// если адрон прошел путь меньший чем его длина формирования, потери энергии не происходят, так как они уже учтены в pythia
			return 1;

		}else{

			if(_verbose)cout<<"resid 2 "<<particleA->getResidualHadronFormationLength()<<endl;
			if(_verbose)cout<<"left 2 "<<particleA->getLeftHadronFormationLength()<<endl;
			passedHadronFormathionLenght = particleA->getResidualHadronFormationLength();
			passedHadronFormathionLenght =0 ;
	    	//cin>>ch;
		}

		double	deltaPath = zCoordinateOfCollisions - particleA->vProd().pz() -passedHadronFormathionLenght;
		Vec4 momentum4;
		momentum4.p(particleA->p());
		if(particleA->vProd().pz() == -_maxZCoordinate) return 1;  //todo возможно второе условие дублирует это.
		if(particleA->getSoftCollisionNumber() == 0) return 1;
		if(_targetNucleus.getNuclearRadius()*_targetNucleus.getNuclearRadius() < particleA->vProd().pT2())return 1;

//if interaction occurred beyond boundary of nucleus hardron does not lose energy
		if(_verbose)cout<<deltaPath<<endl;

		if(deltaPath > 0 && (particleA->isHadron() || isNucleon)){

	if(zCoordinateOfCollisions == sqrt(_targetNucleus.getNuclearRadius()*_targetNucleus.getNuclearRadius() - particleA->vProd().pT2())){
		if(_verbose)cout<<"before "<<particleA->getTotalPathInNucleus()<<endl;
	//	particleA->setTotalPathInNucleus(particleA->getTotalPathInNucleus()+deltaPath);// учитывается путь от последней точки мягкого соударения до вылета из ядра
		//todo осмыслить: если раскоментировать строчку результат в 2а раза больше, если потери энергии происходят
		if(_verbose)cout<<"after "<<particleA->getTotalPathInNucleus()<<endl;
	//	cin>>ch;
	}
	if(_verbose)cout<<"totalPath "<<particleA->getTotalPathInNucleus()<<endl;
	//    cin>>ch;
			//pathInNucleiOutput<<deltaPath<<endl;
		double	deltaE = deltaPath*_kEnergyLoss;
		//suetin debug

	//	cin>>ch;
		//
			if(_verbose)cout<<"zCoordinateOfCollisions = "<<zCoordinateOfCollisions<<" zCoordinateOfCollisionsTemp = "<<zCoordinateOfCollisions<<" vecCoordinate.pz() = "<<particleA->vProd().pz()<<" deltaPath  = "<<deltaPath<<endl;
		//suetin debug
			if(_verbose)cout<<"particleA->e() - deltaE < particleA->m()"<<particleA->e() - deltaE <<" m "<<particleA->mCalc()<< endl;

			if(/*particleA->pz()*particleA->pz() + deltaE*deltaE - 2*deltaE*particleA->e() < 0 ||*/ particleA->e() - deltaE < particleA->mCalc() ){ // second condition for correct calculation particle momentum
				_indexBadInitializations->push_back(i_init);
				if(_verbose)cout << "particle "<<particleA->getHistory()->back()<<" is adsorbed "<<endl;
				//todo history is not initializing
				//if(_verbose)cout << "particle "<<particleA->getHistory()->back()<<" is adsorbed "<<endl;
				//continue;
				return 0;
				//todo think may be do not escape from cycle
			}
			particleA->setEnergyLoss(deltaE + particleA->getEnergyLoss());
			if(_verbose)cout<<" eL "<<particleA->getEnergyLoss()<<endl;
	//
			//   double	deltaP = sqrt(particleA->e()*particleA->e() - particleA->m2()) - sqrt((particleA->e() - deltaE)*(particleA->e() - deltaE)- particleA->m2());
			double	deltaP = 0, newAbsoluteMomentum = 0, newEnergy = 0, newPx = 0, newPy = 0, newPz = 0;
         //   abspi=sqrt(PATT(i,1)**2+PATT(i,2)**2+PATT(i,3)**2)
         //   PabsNew=sqrt(P4newmy**2-dmassi**2)
         //   PATT(I,1)=PATT(I,1)*PabsNew/abspi
         //   PATT(I,2)=PATT(I,2)*PabsNew/abspi
         //   PATT(I,3)=PATT(I,3)*PabsNew/abspi
         //   PATT(I,4)=P4newmy
			newEnergy =  particleA->e() - deltaE;
			newAbsoluteMomentum = sqrt(newEnergy*newEnergy - particleA->m2Calc());
			newPx = particleA->px()*newAbsoluteMomentum/particleA->pAbs();
			newPy = particleA->py()*newAbsoluteMomentum/particleA->pAbs();
			newPz = particleA->pz()*newAbsoluteMomentum/particleA->pAbs();

			if(_verbose)cout<<"m1 "<<particleA->mCalc()<<" p1 = "<<particleA->p()<<endl;

			particleA->e(newEnergy);
			particleA->px(newPx);
			particleA->py(newPy);
			particleA->pz(newPz);

			if(_verbose)cout<<"m "<<particleA->m2()<<" p = "<<particleA->p()<<endl;
		//	cin>>ch;
			//suetin debug
			/*
			deltaP = particleA->pz() - sqrt(particleA->pz()*particleA->pz() + deltaE*deltaE - 2*deltaE*particleA->e() );
	    	if(_verbose)cout<<particleA->p();
	    	if(_verbose)cout<<"deltaE = "<<deltaE<<" deltaP = "<<deltaP<<endl;
	    	//	cin>>ch;
			momentum4.e(particleA->e() - deltaE);
			momentum4.pz(particleA->pz() - deltaP);
			//momentum4.
			//particleA->p(momentum4);
			particleA->e(particleA->e() - deltaE);
			particleA->pz(particleA->pz() - deltaP);
			double restMass = 0;

			particleA->m(restMass);
			restMass = particleA->getRestMass();
			particleA->m(restMass);
			if(_verbose)cout<<"particleA->getRestMass() "<<particleA->getRestMass()<<endl;
			if(_verbose)cout<<"particleA p = "<<particleA->p();
			*/

		}else{

			if(_verbose)	cout<<"energy lose do not occurred, deltapath = "<<deltaPath<<" id particle = "<<particleA->id()<<endl;
			deltaPath = 0;
		}// end of deltaPath > 0 && particleA->isHadron()
		return 1;
	}// end energyLoss
	int getIdTargetNucleon(void){
		double probability = double(_targetNucleus.Z() - _scatteringOnProton)/(_targetNucleus.A() - _scatteringOnProton - _scatteringOnNeutron);
		if(_verbose)cout<< "probability = "<<probability<<endl;
		if(getRandom() < double(_targetNucleus.Z() - _scatteringOnProton)/(_targetNucleus.A() - _scatteringOnProton - _scatteringOnNeutron)){
			if(_verbose)cout<< " in proton collision "<<endl;
			//particleB->id(2212);
			_scatteringOnProton++;
			return 2212;

		}else{
			if(_verbose)cout<< " in neutron collision "<<endl;
			//particleB->id(2112);
			_scatteringOnNeutron++;
			return 2112;
		}
	}
	double getEnergyCut(void){return _energyCut;}
	void setVaribles(void);
	bool isBad(unsigned int index){
		for(unsigned int ibad = 0; ibad < _vertexOfBadHardInteraction->size(); ibad++){
			if(index == _vertexOfBadHardInteraction->at(ibad))return true;
		}
		return false;
	}

	void saveParticle4VectorsAfterEnergyLoss(hardpingParticle* particleA){
		char ch;
		//cout<<" in saveParticle4VectorsAfterEnergyLoss beagin"<<endl;
		//cin>>ch;
		if(_verbose)cout<<"saveParticle4VectorsAfterEnergyLoss1 = "<<particleA->p();
//todo suetin debug
	 	particleA->rotateHardping();

		if(_verbose)cout<<"saveParticle4VectorsAfterEnergyLoss2 = "<<particleA->p();
	//	cout.precision(12);
	//	cout<<"1particleA->theta() "<<particleA->theta()<<" 1particleA->phi() "<<particleA->phi()<<endl;
	//	particleA->setAngles();

		int numberOfGeneration = 0;
		numberOfGeneration = particleA->getNumberOfCurrentGeneration();
		if(numberOfGeneration)	{
		//	cout<<" particle momentum before energy loss "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).p();
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).p(particleA->p());
//		/	cout<<"saveParticle4VectorsAfterEnergyLoss3 "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).p()<<endl;
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).vProd(particleA->vProd());
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).id(particleA->id());
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).setPhiHardping(particleA->getPhiHardping());
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).setThetaHardping(particleA->getThetaHardping());
			if(_verbose)cout<<" if sca part "<<particleA->getIdscatteringParticle()<<endl;
			_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).scatteringOnParticle(particleA->getIdscatteringParticle());
		//	cout<<" particle momentum after  energy loss "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).p();
		}else{
			//todo extrapolate to the nucleus
		//	_generations->at(0).getMatrix()->at(0).at(0).scatteringOnParticle(particleA->getIdscatteringParticle());
		}
		if(numberOfGeneration&&0){
			cout.precision(12);
			//particleA->rotateBackHardping();
			cout<<"px "<<particleA->px()<<endl;
			cout<<"py "<<particleA->py()<<endl;
			cout<<"pz "<<particleA->pz()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss px "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).px()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss py "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).py()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss pz "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).pz()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss theta "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).theta()<<endl;
			cout<<"saveParticle4VectorsAfterEnergyLoss phi "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).phi()<<endl;
		}






	  	particleA->rotateBackHardping();// необходимо развернуть пучок вдоль z, чтобы корректно работала prepareNewGeneration
//		cout<<"saveParticle4VectorsAfterEnergyLoss3 = "<<particleA->p();
		//cout<<" in saveParticle4VectorsAfterEnergyLoss end"<<endl;
	//	cin>>ch;
	}
	bool findDrellYanPairs(int i_pyEv, hardpingParticle* particleA);


	bool addVertexOfInteraction(hardpingParticle* particleA){
		char ch;
		if(_generations->at(particleA->getNumberOfCurrentGeneration()).getMatrix()->at(particleA->getIndexNumber()).size() == 0){
			//this for situation when energy of all produced particles below energyCut
			if(_verbose)cout<<"getMatrix()->at(i_init).size() == 0"<<endl;
			_indexBadInitializations->push_back(particleA->getIndexNumber());
			return false;

		}else{
			// add vertex of interaction
			/*
			cout<<" history back particleA = "<<particleA->getHistory()->back()<<endl;


			cout<<"history particleA : ";
			for(unsigned int i = 0 ; i< particleA->getHistory()->size(); i++){
				cout<<particleA->getHistory()->at(i)<<" ";
			//	cin>>ch;
			}
			cout<<endl;
//			cin>>ch;

			 */
			_vertexOfInteraction->push_back(particleA->getHistory()->back());
			return true;
		}
	}


bool prepareNewGeneration(hardpingParticle* particleA,int i_pyEv);

void deleteBadInitializationParticlesFromNewGeneration(int numberOfGeneration){
	char ch;
	//cout<<" in deleteBadInitializationParticlesFromNewGeneration beagin"<<endl;
	//cin>>ch;
	if(_indexBadInitializations->size()){
	//	for(int i=0;i<_indexBadInitializations->size();i++){
	//
	//			cout<<"i = "<<i<<"  _indexBadInitializations = "<<_indexBadInitializations->at(i)<<endl;
	//	}

		do{
				//cout<<"  _indexBadInitializations = "<<_indexBadInitializations->back()<<endl;
	//		cout<<"bad size "<<_generations->at(numberOfGeneration).getMatrix()->at(_indexBadInitializations->back()).size()<<endl;//
	//		cout<<" val = "<<_generations->at(numberOfGeneration).getMatrix()->at(_indexBadInitializations->back()).at(0).p()<<endl;
			if(_verbose)cout<<" _indexBadInitializations =  "<<_indexBadInitializations->back()<<endl;
			_generations->at(numberOfGeneration).getMatrix()->erase(_generations->at(numberOfGeneration).getMatrix()->begin()+_indexBadInitializations->back());
			_indexBadInitializations->pop_back();
		}while(_indexBadInitializations->size()!= 0);
	}
//todo chash here
	for(unsigned int i_out = 0; i_out <_generations->at(numberOfGeneration).getMatrix()->size(); i_out++){

		for(unsigned int j_out = 0; j_out < _generations->at(numberOfGeneration).getMatrix()->at(i_out).size();j_out++){

			if(_verbose)cout<<" i = "<<i_out<<" j = "<<j_out<<" val = "<<_generations->at(numberOfGeneration).getMatrix()->at(i_out).at(j_out).p()<<" id = "<<_generations->at(numberOfGeneration).getMatrix()->at(i_out).at(j_out).id()<<" index = "<< _generations->at(numberOfGeneration).getMatrix()->at(i_out).at(j_out).getHistory()->back()<<endl;
		}
	}
	//cout<<" in deleteBadInitializationParticlesFromNewGeneration end"<<endl;
	//cin>>ch;
}
unsigned int getIndexHardOfHardCollision(void){

	unsigned int randomPossitionOfHardScattering = 0;
	char ch;
	//cout<<"_vertexOfInteraction size = "<<_vertexOfInteraction->size()<<endl;
	//cout<<"back = "<<_vertexOfInteraction->back()<<endl;

	//_softCollisionIterator = _vertexOfInteraction->begin();

	_indexOfHardCollosion = _vertexOfInteraction->back();

//	cout<<"_indexOfHardCollosion = "<<_indexOfHardCollosion<<endl;
//	cin>>ch;
	// we choose last of soft collisions. if hardping cant't initialize we get random position from _vertexOfInteraction

	if(isBad(_indexOfHardCollosion)){
		if(_vertexOfInteraction->size() == 1 || _vertexOfInteraction->size() == _vertexOfBadHardInteraction->size() +1){
				if(_verbose)cout<< " _vertexOfInteraction->size() == 1 "<<endl;
				_indexOfHardCollosion = 0;
				/// todo add condition for initial particle
			//	particleA->p(_initialParticle.p());
			//	particleA->vProd(_initialParticle.vProd());
				//todo may be add some new propertes

			}else{
				do{

					_softCollisionIterator = _vertexOfInteraction->begin();
					randomPossitionOfHardScattering = (unsigned int)(getRandom()*(_vertexOfInteraction->size() - 1)) +1;
					_softCollisionIterator = _softCollisionIterator + randomPossitionOfHardScattering;
					if(_verbose)cout<< "randomPossitionOfHardScattering = "<< randomPossitionOfHardScattering<<" index of hard scattering = "<< *_softCollisionIterator<<endl;
					_indexOfHardCollosion = *_softCollisionIterator;
			//		cout<<"_indexOfHardCollosion = "<<_indexOfHardCollosion<<endl;
			//		cin>>ch;
					if(_verbose)cout<< " is bad = "<< isBad(_indexOfHardCollosion)<<endl;
				}while(isBad(_indexOfHardCollosion));
			}//end of else _vertexOfInteraction->size() == 1
	}

	return _indexOfHardCollosion;
}
class pythia6Event{
public:

	pythia6Event(){
		pythia6Px = new vector<double>;
		pythia6Py = new vector<double>;
		pythia6Pz = new vector<double>;
		pythia6E = new vector<double>;
		pythia6Id = new vector<int>;
		leptonPosition = new vector<int>;
		pythia6Particle = new vector<hardpingParticle>;
	}
	~pythia6Event(){
		delete pythia6Px;
		delete pythia6Py;
		delete pythia6Pz;
		delete pythia6E;
		delete pythia6Id;
		delete pythia6Particle;
	}
	std::vector<double>* pythia6Px;
	std::vector<double>* pythia6Py;
	std::vector<double>* pythia6Pz;
	std::vector<double>* pythia6E;
	std::vector<  int >* pythia6Id;
	std::vector<  int >* leptonPosition;
	std::vector<hardpingParticle>* pythia6Particle;

};
double getQ2LeptonHadron(pythia6Event* py6Ev){

		double Q2 = 0;


		// in position 0 and 1 locate initial particles

 		for(int i = 3; i <= py6Ev->pythia6Particle->size(); i++ ){
			if(py6Ev->pythia6Particle->at(0).id() == py6Ev->pythia6Particle->at(i).id())py6Ev->leptonPosition->push_back(i);
		}
		Vec4 deltaLeptonP(0);

		//py6Ev->leptonPosition->at(py6Ev->leptonPosition->)
		deltaLeptonP = py6Ev->pythia6Particle->at(0).p() - py6Ev->pythia6Particle->at(py6Ev->leptonPosition->at(0)).p();

		double nu = 0; //energy of virtual photon
		nu = (deltaLeptonP * py6Ev->pythia6Particle->at(1).p())/py6Ev->pythia6Particle->at(1).mCalc();

		Q2 = deltaLeptonP*deltaLeptonP;

		double bl = 0;

		bl = nu/_kEnergyLoss;

		return Q2;
/*        IF(((ABS(IHNT2(5)).LT.11).OR.(ABS(IHNT2(5)).GT.18))) THEN
            z=2.0D0*sqrt(PATT(i,1)**2.0D0+PATT(i,2)**2.0D0)
&              /sqrt(PARI(33)*PARI(34)*PARI(12))
            DZA(I)=2.0D0*DSQRT(PATT(i,1)**2.0D0+PATT(i,2)**2.0D0)
&                   /DSQRT(DXFREL1*DXFR2*HINT1(1)*HINT1(1))
            z=DZA(I)
            IF(z.GT.1.0D0) THEN
                z=0.99D0
            ENDIF
            DZA(I)=z
            snu=dphadron/z
        ENDIF
        */
}
double getQ2HadronHadron(pythia6Event* py6Ev){

}
hardpingParticle* findHardParticle(/*hardpingParticle* particleA*/void){
	bool itIsOk = false;
	bool exit = false;
	double sinTheta = 0, cosTheta =0, sinPhi = 0, cosPhi = 0, anglePhi = 0, angleTheta = 0;
	char ch;
	if(_verbose)cout<<"_generations->size() = "<<_generations->size()<<endl;
	if(_verbose)cout<<"iam in = "<<endl;

	hardpingParticle* particleA;// = new hardpingParticle();
	//cout<<" p "<<particleA->p();
	//cin>>ch;
	for(unsigned int indexOfGeneration = 0; indexOfGeneration < _generations->size(); indexOfGeneration++ ){
		if (indexOfGeneration == 3){
			//cout<<"_generations->at(indexOfGeneration).getMatrix()->size() = "<<_generations->at(indexOfGeneration).getMatrix()->size()<<endl;
	//		cin>>ch;
		}
		//cout<<"indexOfGeneration = "<<indexOfGeneration<<" _generations->size() = "<<_generations->size()<<endl;
		if(exit)break;
		for(unsigned int initialParticleIndex = 0; initialParticleIndex < _generations->at(indexOfGeneration).getMatrix()->size(); initialParticleIndex++){
			//cout<<"at indexOfGeneration = "<<indexOfGeneration<< "_generations->at(indexOfGeneration).getMatrix()->size() = "<<_generations->at(indexOfGeneration).getMatrix()->size()<<endl;
			if(exit)break;
			for(unsigned int producedParticleIndex = 0; producedParticleIndex < _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).size(); producedParticleIndex++){
				if(exit)break;
				//cout<<"indexOfGeneration = "<<(int)indexOfGeneration<<" initialParticleIndex = "<<initialParticleIndex<<" producedParticleIndex = "<<producedParticleIndex<<endl;
				//cout<<"_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->size() = "<<_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->size()<<endl;
				//todo some particle havent history
				for(int i =0; i< _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->size(); i++){
				//	cout<<_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->at(i)<<endl;
				}
				if(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->back() ==
						_indexOfHardCollosion){
				//	cin>>ch;
			//todo строчка все косячит		//particleA = new hardpingParticle(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex));
					if(_verbose)	cout<<" diff = "<<_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).p();
					if(_verbose)pythia->event.list();
					pythia->event.append(*_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getPythiaParticle());
					if(_verbose)pythia->event.list();
					particleA = new hardpingParticle(pythia->event.back());
					particleA->e(sqrt(particleA->getRestMass()*particleA->getRestMass() + particleA->pAbs2()));
					particleA->m(particleA->getRestMass());
											//cout<<"incidentParticle "<<particleA->id()<<" p = "<<particleA->p();
					pythia->event.remove(pythia->event.size()-1,pythia->event.size()-1);
					if(_verbose)pythia->event.list();
					for(int i = 0; i < _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->size();i++){
							cout<<_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->at(i)<<" ";
					}
					cout<<endl;
					//cin>>ch;
					//particleA->
					//particleA->p(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).p());
					//particleA->vProd(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).vProd());
					particleA->scatteringOnParticle(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getIdscatteringParticle());
					particleA->setPhiHardping(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).phi());
					particleA->setThetaHardping(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).theta());
					particleA->getHistory()->clear();
					particleA->getHistory()->assign(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->begin(),_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getHistory()->end());

					_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getAngles(sinPhi, cosPhi, sinTheta, cosTheta);
					angleTheta = _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getThetaHardping();
					anglePhi = _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getPhiHardping();
					particleA->setPhiHardping(anglePhi);
					particleA->setThetaHardping(angleTheta);
					angleTheta = _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getThetaHardping2();
					anglePhi = _generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getPhiHardping2();
					particleA->setPhiHardping2(anglePhi);
					particleA->setThetaHardping2(angleTheta);
					particleA->setAngles(sinPhi, cosPhi, sinTheta, cosTheta);
					particleA->setSoftCollisionNumber(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getSoftCollisionNumber()-1);//минус 1, потому что одно мягкое стало жестким
					particleA->setEnergyLoss(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getEnergyLoss());

					//double sinY1 = 0, cosY1 = 0, sinX2 = 0, cosX2 = 0;
					//_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getTrigonometricFunctions(sinY1,cosY1,sinX2,cosX2);

					_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).getTrigonometricFunctions(sinPhi, cosPhi, sinTheta, cosTheta);
					particleA->setTrigonometricFunctions(sinPhi, cosPhi, sinTheta, cosTheta);
					cout<<"trig "<<sinPhi<<" "<<cosPhi<<" "<<sinTheta<<" "<<cosTheta<<endl;
				//	cin>>ch;
					cout<<"size "<<particleA->getHistory()->size()<<endl;
					for(int i = 0; i < particleA->getHistory()->size(); i++ ){
						cout<<particleA->getHistory()->at(i)<<" ";
					}
					cout<<endl;
				//	cin>>ch;
					if(_verbose)cout<<"particleA->isHadron() "<<particleA->isHadron()<<" id "<<particleA->id()<<endl;
					if(_verbose)cout<<"particleA->vProd()    "<<particleA->vProd();
					if(_verbose)cout<<"getIdscatteringParticle    "<<particleA->getIdscatteringParticle();
			//		cin>>ch;
					//particleA->id(_generations->at(indexOfGeneration).getMatrix()->at(initialParticleIndex).at(producedParticleIndex).id());
/*
					for(unsigned int i = 0 ; i < particleA->getHistory()->size(); i++){
						cout<<particleA->getHistory()->at(i)<<" ";
					}
					cout<<endl;
*/
					//todo may be add some new parameter
					if(_verbose)cout<<" particleA = "<<particleA->p()<<" idb "<<particleA->getIdscatteringParticle()<<endl;
					return particleA;
					//exit = true;
				//	cin>>ch;
					break;
				}// end of if getHistory()->back() ==_indexOfHardCollosion
			//	cout<<"1"<<endl;
			}// end producedParticleIndex cycle
		//	cout<<"2"<<endl;
		}// end initialParticleIndex cycle
	//	cout<<"3"<<endl;
	}// end indexOfGeneration cycle
//	cout<<"4"<<endl;
	return particleA;
}

bool ifNoHardCollisionHappened(unsigned int numberOfGeneration/*, hardpingParticle * particleA*/){
	char ch;
	//cout<<" in ifNoHardCollisionHappened beagin"<<endl;
	//cin>>ch;
	//Particle particleB
	bool noHardCollisionHappened = false;
	if(_generations->at(numberOfGeneration-1).getMatrix()->size() == 0 && _hardInteractionSummaryCount == 0)noHardCollisionHappened = true;
	if(noHardCollisionHappened == false)return noHardCollisionHappened;
	double sinTheta = 0, cosTheta =0, sinPhi = 0, cosPhi = 0, anglePhi = 0, angleTheta = 0;
	hardpingParticle* particleA ;//= new hardpingParticle();
	//particleA = new hardpingParticle();
	/*Particle newPythiaParticle;
	newPythiaParticle.pz(_energyLab);
	newPythiaParticle.id(2212);
	newPythiaParticle.status(0);
	//pythia->event.list();
	pythia->event.append(newPythiaParticle);
	//pythia->event.list();
	particleA = new hardpingParticle(pythia->event.back());
	particleA->e(sqrt(particleA->getRestMass()*particleA->getRestMass() + particleA->pAbs2()));
	particleA->m(particleA->getRestMass());

	//cout<<"incidentParticle "<<particleA->id()<<" p = "<<particleA->p();
	pythia->event.remove(pythia->event.size()-1,pythia->event.size()-1);
	for(int i = 0; i < particleA->getHistory()->size(); i++ ){
		cout<<particleA->getHistory()->at(i)<<" ";
	}

	cout<<"size beg "<<particleA->getHistory()->size()<<endl;
	cout<<"p "<<particleA->p()<<endl;
	cout<<"x "<<particleA->vProd().px()<<endl;
	cout<<"y "<<particleA->vProd().py()<<endl;
	cout<<"z "<<particleA->vProd().pz()<<endl;
	cout<<"size "<<particleA->vProd().e()<<endl;

	cin>>ch;
*/
	//pythia->event.clear();
	//pythia->event.list();
	//hardpingParticle* particleA;
	//particleA = new hardpingParticle(newPythiaParticle);
	//cout<<"particleA->isHadron() "<<particleA->isHadron()<<" id "<<particleA->id()<<endl;
	//cin>>ch;
	bool softToHard = false;
	//unsigned int randomPossitionOfHardScattering = 0 ;

//	cout<<"mat size _generations->at(numberOfGeneration-1).getMatrix()->size() = "<<_generations->at(numberOfGeneration-1).getMatrix()->size()<<" _hardInteractionCount = "<<_hardInteractionCount<<endl;
//	cin>>ch;
	if(_generations->at(numberOfGeneration-1).getMatrix()->size() == 0 && _hardInteractionSummaryCount == 0){


		if(_verbose)cout<<"_vertexOfInteraction->size() = "<<_vertexOfInteraction->size()<<" _vertexOfBadHardInteraction->size() = "<<_vertexOfBadHardInteraction->size()<<endl;
			if( _vertexOfInteraction->size() == _vertexOfBadHardInteraction->size())return softToHard;//condition for exit

				_softToHard = true;

				if(_verbose)cout<< "!!!!!!!!!!!!!!!!!!hard to soft part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
				if(_verbose)for(unsigned int iv = 0; iv < _vertexOfInteraction->size(); iv ++){
					cout<<_vertexOfInteraction->at(iv)<<" ";
				}
				if(_verbose)cout<<endl;
				if(_verbose)cout<<" vertexOfInteraction->size() "<<_vertexOfInteraction->size()<<" vertexOfBadHardInteraction->size() "<<_vertexOfBadHardInteraction->size()<<endl;


				if(getIndexHardOfHardCollision()){
					//when we get particle she is oriented along initial z direction
					particleA = findHardParticle();
					/*
					cout<<"size "<<particleA->getHistory()->size()<<endl;
					for(int i = 0; i < particleA->getHistory()->size(); i++ ){
						cout<<particleA->getHistory()->at(i)<<" ";
					}
					cout<<"size "<<particleA->getHistory()->size()<<endl;
					cout<<"p "<<particleA->p()<<endl;
					cout<<"x "<<particleA->vProd().px()<<endl;
					cout<<"y "<<particleA->vProd().py()<<endl;
					cout<<"z "<<particleA->vProd().pz()<<endl;
					cout<<"size "<<particleA->vProd().e()<<endl;
					*/
					//double t1,t2,t3,t4;
					//particleA->getAngles(t1,t2,t3,t4);

					//cout<<"t1234 "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;
					//particleA->getTrigonometricFunctions(t1,t2,t3,t4);
					//cout<<"t1234 2 "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;
				//	cin>>ch;

				}else{


						//Particle newPythiaParticle(*_initialParticle.getPythiaParticle());

					    if(_verbose)pythia->event.list();
						pythia->event.append(*_initialParticle.getPythiaParticle());
						if(_verbose)pythia->event.list();
						particleA = new hardpingParticle(pythia->event.back());
						particleA->e(sqrt(particleA->getRestMass()*particleA->getRestMass() + particleA->pAbs2()));
						particleA->m(particleA->getRestMass());
						//cout<<"incidentParticle "<<particleA->id()<<" p = "<<particleA->p();
						pythia->event.remove(pythia->event.size()-1,pythia->event.size()-1);
						if(_verbose)pythia->event.list();
						//particleA = new hardpingParticle(_initialParticle);
						if(_verbose)cout<<"particleA->isHadron() "<<particleA->isHadron()<<" id "<<particleA->id()<<endl;
						if(_verbose)cout<<"particleA->vProd()    "<<particleA->vProd();
						if(_verbose)cout<<"getIdscatteringParticle    "<<particleA->getIdscatteringParticle();

						particleA->scatteringOnParticle(_initialParticle.getIdscatteringParticle());
						particleA->getHistory()->clear();
						for(unsigned int i_history = 0; i_history < _initialParticle.getHistory()->size(); i_history++){
							particleA->getHistory()->push_back(_initialParticle.getHistory()->at(i_history));
						}
						_initialParticle.getAngles(sinPhi, cosPhi, sinTheta, cosTheta);
						angleTheta = _initialParticle.getThetaHardping();
						anglePhi = _initialParticle.getPhiHardping();
						particleA->setPhiHardping(anglePhi);
						particleA->setThetaHardping(angleTheta);
						angleTheta = _initialParticle.getThetaHardping2();
						anglePhi = _initialParticle.getPhiHardping2();
						particleA->setPhiHardping2(anglePhi);
						particleA->setThetaHardping2(angleTheta);
						particleA->setAngles(sinPhi, cosPhi, sinTheta, cosTheta);
						_initialParticle.getTrigonometricFunctions(sinPhi, cosPhi, sinTheta, cosTheta);
						particleA->setTrigonometricFunctions(sinPhi, cosPhi, sinTheta, cosTheta);
						particleA->setSoftCollisionNumber(_initialParticle.getSoftCollisionNumber());
						particleA->setEnergyLoss(_initialParticle.getEnergyLoss());
				}
				//cout<<"is h soft to hard "<<particleA->isHadron()<<endl;
			//	cin>>ch;
				findParticlesWithIndexOfHardCollisions();

				_generations->at(numberOfGeneration-1).getMatrix()->resize(1);

				_generations->at(numberOfGeneration-1).getMatrix()->at(0).push_back(*particleA);//at(0).push_back(*particleA);
				if(_verbose)cout<<" new generation = "<< _generations->at(numberOfGeneration-1).getMatrix()->at(0).at(0).p();
			}// end of if(_generations->at(numberOfGeneration-1).getMatrix()->size() == 0 && _hardInteractionCount == 0)
	delete particleA;
	return _softToHard;
	//cout<<" in ifNoHardCollisionHappened end"<<endl;
	//cin>>ch;
}
void findParticlesWithIndexOfHardCollisions(void){
	// find all particle which produced in soft collision with index _indexOfHardCollosion and put them into massive _indexSoftToHard
				for(unsigned int isoft = 0; isoft < _outOfNucleus->size(); isoft++){
					for(unsigned int jsoft = 0; jsoft <_outOfNucleus->at(isoft).getHistory()->size() ; jsoft++){
						if(_outOfNucleus->at(isoft).getHistory()->at(jsoft) == _indexOfHardCollosion) _indexSoftToHard->push_back(isoft);
					}
				}

				for(unsigned int isoftBadInit = 0; isoftBadInit < _notInit->size(); isoftBadInit++){
					for(unsigned int jsoftBadInit = 0; jsoftBadInit <_notInit->at(isoftBadInit).getHistory()->size() ; jsoftBadInit++){
						if(_notInit->at(isoftBadInit).getHistory()->at(jsoftBadInit) == _indexOfHardCollosion) _indexSoftToHardBadInit->push_back(isoftBadInit);
					}
				}
}
void getParticleFromPreviousGeneration(hardpingParticle* particleA){
	// get particle from previous generation
	char ch;
	unsigned int numberOfGeneration = 0;
	unsigned int i_init = 0;
	numberOfGeneration = particleA->getNumberOfCurrentGeneration();
	i_init = particleA->getIndexNumber();


	//particleA->getHistory()->clear();
	cout.precision(12);
	/*cout<<"hui "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).px()<<endl;
	cout<<"hui "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).py()<<endl;
	cout<<"hui "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(particleA->getIndexNumber()).i).at(_index->at(particleA->getIndexNumber()).j).pz()<<endl;





cout<<"px "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).px()<<endl;
cout<<"py "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).py()<<endl;
cout<<"pz "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).pz()<<endl;
cout<<"theta "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).theta()<<endl;
cout<<"phi "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).phi()<<endl;
*/
	*particleA = _generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j);
	//cout<<"is h getParticleFromPreviousGeneration "<<particleA->isHadron()<<endl;
	//cin>>ch;


	// set index number of particle in current generation
	particleA->setInexNumber(i_init);
	//and number of current generation
	particleA->setNumberOfCurrentGeneration(numberOfGeneration);


//	cout<<"theta "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).theta()<<" phi "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).phi()<<endl;
//	cout << "particleA->theta() 0 "<<particleA->theta()<<" particleA->phi() "<<particleA->phi()<<endl;
//	cout<<"p  getParticleFromPreviousGeneration= "<<particleA->p();
//	cout.precision(12);
//	cout << "particleA->theta() 0 "<<particleA->theta()<<" particleA->phi() "<<particleA->phi()<<endl;
	//particleA->rotateBackHardping();
	//cout<<"p  getParticleFromPreviousGeneration 2= "<<particleA->p();
	////////////////////////////////////////
//	cout<<"particleA id = "<<_generations->at(numberOfGeneration-1).getMatrix()->at(_index->at(i_init).i).at(_index->at(i_init).j).id()<<endl;
//	cout<<"particleA id = "<<particleA->id()<<endl;
	if(	!_softToHard){
		// todo may be we need remember angles after particle production immediately


		/// save angles relatively initial beam direction///
//		particleA->setAngles();
		//////////////////////////////////////////////////
	//	cout.precision(10);

		if(_verbose)cout<<" momentum  before  rotate   = "<<particleA->p()<<endl;
		if(_verbose)cout<<" coordinate before rotation = "<<particleA->vProd()<<endl;

		//dispose momentum of particle along z direction
	//	particleA->rotateBackHardping();
		//////////////////////////////////////////////////

		if(_verbose)cout<<" momentum  after rotate_    = "<<particleA->p()<<endl;
		if(_verbose)cout<<" coordinate after rotation_ = "<<particleA->vProd()<<endl;



		//////////////////////////////////////////////////////////
	}

}
void setInitinalImpactAndIndex(hardpingParticle* particleA){

	char ch;

	//particleA->vProd().px(gsl_ran_gaussian (gslRandomGenerator, 1));
	//particleA->vProd().py(gsl_ran_gaussian (gslRandomGenerator, 1));

	//cout<<"gsl_ran_gaussian "<<particleA->vProd();

//	cin>>ch;
		  // Function: double gsl_ran_gaussian_pdf (double x, double sigma).
		  //This function computes the probability density p(x) at x for a Gaussian distribution with standard deviation sigma, using the formula given above.


	////////// calculating impact parameter of incident particles//////////////////////////////////////////



	double xImpact = 0, yImpact = 0, zCoordinate = 0, impact = 0, phi = 0;
	double maxHalfPathInNucleus = 0;
	double R = 0;
	double temp1 = 0;

	if(particleA->isLepton()){


		R = _targetNucleus.getNuclearRadius();
//		xImpact = gsl_ran_gaussian (gslRandomGenerator, 1);
//		yImpact = gsl_ran_gaussian (gslRandomGenerator, 1);
//		zCoordinate = gsl_ran_gaussian (gslRandomGenerator, 1);

		//impact 	= getNewImpact();
		impact = getImpactParameter();
		phi = 2* M_PIl * getRandom();
//		temp1 = getRandomFromFile();
//		cout<<"temp 1 = "<<temp1<<endl;
//		phi = 2* M_PIl * temp1;
		 xImpact  =  impact*sin(phi);//getPointOfInteraction();//impact*sin(phi);
		yImpact	 =  impact*cos(phi);//getPointOfInteraction();//impact*cos(phi);
	//	DHALFMAXPATH=SQRT((R**(2.0D0))-((DIPA(ili))**(2.0D0)))
		maxHalfPathInNucleus = sqrt(R*R - impact*impact);

		zCoordinate = -maxHalfPathInNucleus + 2*maxHalfPathInNucleus*getRandom();//getPointOfInteraction();
		//temp1 = getRandomFromFile();
		//cout<<"temp 2 = "<<temp1<<endl;
		//zCoordinate = -maxHalfPathInNucleus + 2*maxHalfPathInNucleus*temp1;//getPointOfInteraction();
		cout.precision(12);
		if(_verbose)cout<<"impact "<<impact<<" phi "<<phi<<" z "<<zCoordinate<<endl;
	//	cin>>ch;
	//	xImpact = getPointOfInteraction();
	//	yImpact	= getPointOfInteraction();
	//	zCoordinate	= getPointOfInteraction();
	}else{

		getNucleusImpactParameter(xImpact,yImpact);
		zCoordinate = -_maxZCoordinate;
	}

	Vec4 vecCoordinate(0);
	//////////end of calculating impact parameter of incident particles//////////////////////////////////////////

	/////// set coordinate of incident particle////////////////////////////////////
//	coordinateFile>>xImpact>>yImpact>>zCoordinate;
 	//cout<<"coord23 "<<xImpact<<" "<<yImpact<<" "<<zCoordinate<<endl;
//  	cin>>ch;
	vecCoordinate.px(xImpact);
	vecCoordinate.py(yImpact);
	vecCoordinate.pz(zCoordinate); //todo suetin debug
/*
	vecCoordinate.px(0.5);
	vecCoordinate.py(0);
	vecCoordinate.pz(-2);
*/
	particleA->vProd(vecCoordinate);

	//////// end of set coordinate of incident particle ////////////////////

	// set sequentially number for each incident particle. It done for usable trace of particle in program

	particleA->getHistory()->clear();
	particleA->getHistory()->push_back(_indexParticle);
//	cout<<" initial index particle "<<particleA->getHistory()->back()<<endl;
	//cin>>ch;
	_indexParticle++;
	if(_verbose)cout<<"coordinate initial = "<<particleA->vProd();

	//cin>>ch;
	//////////////////////////////////////////////////////////////////////////////////////////////////////
}

double getNuclearDensity(double coordinate){


	if(_verbose)cout.precision(12);
	double nuclearDensity = 0;
	double r2 = _targetNucleus.getNuclearRadius()*_targetNucleus.getNuclearRadius();
	if(_verbose)cout<<"r2 = "<<r2<<endl;
	double r3 = _targetNucleus.getNuclearRadius()*r2;
	if(_verbose)cout<<"r3 = "<<r3<<endl;
	nuclearDensity =  exp(-coordinate*coordinate/r2)*_targetNucleus.A()/r3/M_PIl/sqrt(M_PIl);

	if(_verbose)cout<<"nuclearDensity = "<<nuclearDensity<<endl;
	return nuclearDensity;

}
double myGausLeft(double x){
	double sigma =  1.49316;
	double mean = -1.89812;
	double constant = 0.0862204;

	return 1.57*constant/sigma*exp(-(x - mean)*(x - mean)/2./sigma/sigma);


}
double myGausRight(double x){
	double sigma =  1.95555;
	double mean = 0.948823;
	double constant = 0.0434087;
	return 1.95*constant/sigma*exp(-(x - mean)*(x - mean)/2./sigma/sigma);


}
double getPointOfInteraction(void){

	double coordinate = 0;
	double myGaus = 0;
		//			suetin debug
	do{
		coordinate = -_targetNucleus.getNuclearRadius() + 2*getRandom()*_targetNucleus.getNuclearRadius();
	}while(getRandom() > getNuclearDensity(coordinate));// todo отнормировать все по y


/*
	do{
		coordinate = -_targetNucleus.getNuclearRadius() + 2*getRandom()*_targetNucleus.getNuclearRadius();
		(coordinate < 0)? myGaus = myGausLeft(coordinate) : myGaus = myGausRight(coordinate);  //fortran harping impactParameter for Kr
	}while(0.09*getRandom() > myGaus);
*/

	return coordinate;
}

double getNewImpact(void){

	double coordinate = 0;
	double probability = 0;

	double p0 = 0;
	double p1 = 0;
	double p2 = 0;
	double p3 = 0;
	double p4 = 0;
	double p5 = 0;
	double p6 = 0;
	double p7 = 0;
	double p8 = 0;
	double p9 = 0;
//	double random = 0;
	p0 = 0.00130852;
	p1 = 0.0114604;
	p2 = 0.0277417;
	p3 = -0.0140817;
	p4 = 0.0019539;
	p5 = 0.000102176;
	p6 = -5.81069*pow(10,-5);
	p7 = 6.73422*pow(10,-6);
	p8 = -3.47855*pow(10,-7);
	p9 = 6.98531*pow(10,-9);
	do{
		coordinate = 11.*getRandom();
		probability = p9*pow(coordinate,9) + p8*pow(coordinate,8) + p7*pow(coordinate,7) + p6*pow(coordinate,6) + p5*pow(coordinate,5) + p4*pow(coordinate,4) + p3*pow(coordinate,3) +p2*pow(coordinate,2) + p1*coordinate + p0;
		//coordinate = -_targetNucleus.getNuclearRadius() + 2*getRandom()*_targetNucleus.getNuclearRadius();
		//(coordinate < 0)? myGaus = myGausLeft(coordinate) : myGaus = myGausRight(coordinate);  //fortran harping impactParameter for Kr
	}while(0.06*getRandom() > probability);

/*	do{
		coordinate = -_targetNucleus.getNuclearRadius() + 2*getRandom()*_targetNucleus.getNuclearRadius();
		(coordinate < 0)? myGaus = myGausLeft(coordinate) : myGaus = myGausRight(coordinate);  //fortran harping impactParameter for Kr
	}while(0.09*getRandom() > myGaus);
*/

	return coordinate;
}
double getImpactParameter(void){
	double R = _targetNucleus.getNuclearRadius();
	double impact = 0;
	impact = R*getRandom();
	//impact = -R + 2*R*getRandomFromFile();
	//double temp = 0;
	//temp = getRandomFromFile();
	//cout<<"temp "<<temp<<endl;
	//impact = R*temp;
	return impact;
}
void softToHard(hardpingParticle* particleA,hardpingParticle* particleB){
	char ch;
	//cout<<" in softToHard beagin"<<endl;
	//cin>>ch;
	particleA->setHard();
	particleA->setSoft(false);
	//coordinateHardOutput<<particleA->vProd();
	if(particleA->getIdscatteringParticle() == 0){
		_scatteringOnNeutron = 0;
		_scatteringOnProton = 0;
		particleB->id(getIdTargetNucleon());
	}//end if(particleA->getIdscatteringParticle() == 0)
	else{
		if(_verbose)cout<< " else of (particleA->getIdscatteringParticle() == 0) "<<endl;
		particleB->id(particleA->getIdscatteringParticle());
		if(_verbose)cout<<"id particle B = "<<particleB->id()<<endl;
	} // else if(particleA->getIdscatteringParticle() == 0)
	//particleB->id(particleA->getIdscatteringParticle());
	//cout<<"here i am "<<particleB->id()<<endl;
//	cout<<" in softToHard end"<<endl;
	//cin>>ch;
}
bool checkEnergyCut(hardpingParticle* particleA){
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if after energy loss particle not have enough energy for initializing new collision.
	// if particle set out of nucleus put her into energy cut massive,
	// if not decrease her energy correspond to residual path in nucleus.
	// if energy became bellow zero or equal assume that particle is adsorbed if not add her into energy cut massive.
	double tempLenght = 0;//todo change a name
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//part for particles which can't initialize new event. for this particles - if z coordinate position more that nuclear radius assumed that particle set out of nucleus
	//	else we propose particle propagates in nuclear matter along the z direction and keeps losing energy until it reaches the boundary of nucleus (i.e. nuclear radius).
	bool isNotAdsorbed = true;
	if(particleA->e() < _energyCut){
		tempLenght = sqrt(_targetNucleus.getNuclearRadius()*_targetNucleus.getNuclearRadius() - particleA->vProd().pT2());
		//suetin debug
		isNotAdsorbed = energyLoss(particleA,tempLenght);
		//isNotAdsorbed = 1;
		// isNotAdsorbed = 0 - particle is adsorbed
		// isNotAdsorbed = 1 - all right
		if(isNotAdsorbed){
			particleA->setOut();
			//dispose momentum of particle along initial beam direction
			//suetin debug
			 particleA->rotateHardping();


			if(particleA->getVirtualPhotonEnergy() != 0 && particleA->isHadron()){

				particleA->setHadronEnergyFraction(particleA->e()/particleA->getVirtualPhotonEnergy());

				if(_verbose)cout<<" HadronEnergyFraction = "<<particleA->getHadronEnergyFraction()<<endl;
			}



			// put particle in massive
			_outOfNucleus->push_back(*particleA);
			// massive with index of particles which not initialize next wave
			_indexBadInitializations->push_back(particleA->getIndexNumber());
		}

	}//end of if (particleA->e() < _energyCut)
	return isNotAdsorbed;
}
void finalOutput(void){
	if(_verbose)for(unsigned int i=0;i < _finalState->size() ; i++){
		//_finalState->at(i).setHardpingVector(*particleA);
		cout<<"i = "<<i<<"  id = "<<_finalState->at(i).id()<<" energy = "<<_finalState->at(i).p();
		cout<< "index =  ";
		//for(unsigned int i = 0; i < _finalState->at(i).getHistory()->size(); i++){
		//	cout<<_finalState->at(i).getHistory()->at(i)<<" ";
		//}
	}


	if(_verbose)cout<<"final state size = "<<_finalState->size()<<endl;

	for(unsigned int i = 0 ; i < _cutMassive->size() ; i++ ){
			//_finalState->at(i).setHardpingVector(*particleA);
		if(_verbose)cout<<"i = "<<i<<"  id = "<<_cutMassive->at(i).id()<<" energy = "<<_cutMassive->at(i).e()<<endl;
		if(_verbose)cout<< "index =  ";
			for(unsigned int j = 0; j < _cutMassive->at(i).getHistory()->size(); j++){
					cout<<_cutMassive->at(i).getHistory()->at(j)<<" ";
			}
			if(_verbose)cout<<endl;
		}
	if(_verbose)cout<<"_cutMassive state size = "<<_cutMassive->size()<<endl;
	for(unsigned int ii = 0; ii<	_indexSoftToHard->size(); ii++){
		if(_verbose)cout<< " index of soft to hard  = "<<_indexSoftToHard->at(ii)<<" _indexSoftToHard->size() = "<<_indexSoftToHard->size()<<endl;
	}

	if(_verbose)for(unsigned int i = 0 ; i < _notInit->size() ; i++ ){
		//_finalState->at(i).setHardpingVector(*particleA);
		cout<<"i = "<<i<<"  id = "<<_notInit->at(i).id()<<" energy = "<<_notInit->at(i).p();
		cout<< "index =  ";
		for(unsigned int j = 0; j < _notInit->at(i).getHistory()->size(); j++){
				cout<<_notInit->at(i).getHistory()->at(j)<<" ";
		}
		cout<<endl;

	}
	if(_verbose)cout<<"_notInit state size = "<<_notInit->size()<<endl;

	if(_verbose)for(unsigned int i = 0 ; i < _outOfNucleus->size() ; i++ ){
		//_finalState->at(i).setHardpingVector(*particleA);
		cout<<"i = "<<i<<"  id = "<<_outOfNucleus->at(i).id()<<" energy = "<<_outOfNucleus->at(i).p();
		cout<< "index =  ";
		if(_verbose)for(unsigned int j = 0; j < _outOfNucleus->at(i).getHistory()->size(); j++){
				cout<<_outOfNucleus->at(i).getHistory()->at(j)<<" ";
		}
		if(_verbose)cout<<endl;
	}
	if(_verbose)cout<<"_outOfNucleus state size = "<<_outOfNucleus->size()<<endl;
	if(_verbose)cout<<"_vertexOfInteraction = "<<endl;
	if(_verbose)for(unsigned int iv = 0; iv < _vertexOfInteraction->size(); iv ++){
			cout<<_vertexOfInteraction->at(iv)<<" ";
	}
	if(_verbose)cout<< "number of hard scattering = "<<_hardInteractionSummaryCount<<endl;
	if(_verbose)cout<< "number of soft scattering = "<<_softInteractionSummaryCount<<endl;
}
	void notPythiaNext(hardpingParticle* particleA){
		if(_verbose)cout<<"in !pythiaNextFlag = "<<particleA->id()<<"  "<<particleA->p();

		int isNotAdsorbed = 0;
		if(_softToHard){
			_indexBadInitializations->push_back(particleA->getIndexNumber());
			_vertexOfBadHardInteraction->push_back(_indexOfHardCollosion);
			_indexSoftToHard->resize(0);
			return;
		}
		//todo may be fucked up here
		if(particleA->getNumberOfCurrentGeneration() && particleA->isHadron()){
			//////////////////////////////////////////////////////////////////////////////
			//if pythia can not initialize new event we assumed that projectile hardron //
			//lose energy during propagation in nucleus along the z direction			//
			//until it reaches the boundary of nucleus (i.e. nuclear radius).			//
			//////////////////////////////////////////////////////////////////////////////

			if(_verbose)cout<<"1before particleA->isHadron() = "<<particleA->p();
			//set particle momentum along z direction
			//particleA->setAngles();
			//particleA->rotateBackHardping();
			if(_verbose)cout<<"1after particleA->isHadron() = "<<particleA->p();

		//	cin>>lll;
			//suetin debug
			double tempLenght = sqrt(_targetNucleus.getNuclearRadius()*_targetNucleus.getNuclearRadius() - particleA->vProd().pT2());
			 isNotAdsorbed = energyLoss(particleA, tempLenght);
			//isNotAdsorbed = 1;
			// isNotAdsorbed = 0 - particle is adsorbed
			// isNotAdsorbed = 1 - all right
			if(isNotAdsorbed){
				particleA->setOut();
				//dispose momentum of particle along initial beam direction
				//suetin debug
				 particleA->rotateHardping();
			}else{return;}
		}// end of if(numberOfGeneration && particleA->isHadron())
		// put particle in massive
		_notInit->push_back(*particleA);
		// massive with index of particles which not initialize next wave
		_indexBadInitializations->push_back(particleA->getIndexNumber());
		return;

	}
	void useFortranMethod(bool value){
		_fortranHardping = value;
	}


	//double ::getRandomFromFile();
	void getNucleusImpactParameter(double &impactX,double &impactY){
	//	double impactParameterMax = min(3*_targetNucleus.getNuclearRadius(),_maxZCoordinate + 6.06); // from Fortran version. 6.06 fm - is proton characteristic

		//BB=DSQRT(BMIN**2+RAN(NSEED)*(BMAX**2-BMIN**2))



		double impactPar = sqrt(_impactParameterMin*_impactParameterMin + getRandom()*(_impactParameterMax*_impactParameterMax - _impactParameterMin*_impactParameterMin));
		//impactPar = getRandom()*impactParameterMax;


		double phiImpact = 2*M_PIl*getRandom();
		//cout<<" phi1 = "<<phiImpact<<endl;
		impactX = sin(phiImpact)*impactPar;
		impactY = cos(phiImpact)*impactPar;

		//suetindebug
		//impactX = 0;
		//impactY = 0;
		//suetindebug end
		return;

	}


	/*
	double getRandomFromFile(){
		double a;
		_randomFile>>a;
		return a;
	}

	*/

	int getNumberOfSoftInteraction(void){
		return _softInteractionSummaryCount;
	}


	bool   _firstCall;
	unsigned int _softInteractionSummaryCount;
	double getMaxZCoordinate(){
		return _maxZCoordinate;
	}
	void setVerbose(int verbose){
		_verbose = verbose;
	}
private:

	double woodSaxonSkinDepth() const { return 0.53;  }  ///< returns surface (0.53 fm for Au)

//	double getNuclearRadius  () const { return _r0 ; }

//	Pythia8::Pythia* _pythia;


	double _BQ;
	double _BP;
	double _BN;
	int    _Z;                      ///< atomic number of nucleus
	int    _A;                      ///< nucleon number of nucleus
	double _energyLab;
	int _randomSeed;
	Rndm * _random;
	double _xPartonCoord;
	double _yPartonCoord;
	double _zPartonCoord;
	double _xPartonMomentum;
	double _yPartonMomentum;
	double _zPartonMomentum;
	double _absPartonMomentum;
	Vec4 vek;
//	std::vector *a;

//	std::vector< std::vector< Vec4 > > * m(1, std::vector<Vec4>(1));
	double _energyCut;
	double _kEnergyLoss;
	double _r0;				// nuclear radius
	double _rho0;
	double _Q0;
	double _lastResultOfThicknessCalculation;
	double _lastErrorOfThicknessCalculation;
	unsigned int _scatteringOnProton;
	unsigned int _scatteringOnNeutron;
	//double _precisionOfNuclearThicknessCalculation;
	//bool   _firstCall;
	double _maxZCoordinate; //Maximum radial coordinator for target nucleons
	unsigned int _factorials[13];//={1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
	nucleus _targetNucleus;
	nucleus _projectileNucleus;
	bool 	_isScattering;
	hardpingParticle* _tempParticle;
	//hardpingParticle _initialParticle;
	hardpingParticle _incidentParticle;
	unsigned int _indexParticle;
	bool _softInteractionFlag;
	bool _hardInteractionFlag;
	bool _hardInteraction; // correspond to situation when hard collision do not occurred in regular regime of running
	unsigned int _hardInteractionSummaryCount;
	//unsigned int _softInteractionCount;
	unsigned int _indexOfHardCollosion;
	int _verbose;
	bool _softToHard;
	bool _fortranHardping;
	Timer _time;
	bool _cutMass;
	ifstream _randomFile;

	std::vector<unsigned int>::iterator _softCollisionIterator;
	std::vector<unsigned int> _indexOfDrellYanChain;
	TMacro _pythia6Event;
	ifstream* _pythia6File;

}; // end nucleus class

}// end of namespace pythia8
#endif  // NUCLEUS_H
