#include "../include/Pythia8/Pythia.h"
#include "../include/Pythia8/nucleus.h"

#include <fstream>
//#include "../Pythia8/timer.h"
#define HardQCD 1
#define SoftQCD 1
#define CMS 0
#define NumberOfEvents 10000000
#define verbose 2
#define nbin 200
typedef std::vector< std::vector<Pythia8::hardpingParticle> > vectorConstruction;
using namespace Pythia8;
std::string convertNumberOfEventsToString(unsigned int numEvent);
//std::string getElementName(int Z, int id = 0);
//std::string getElementName(nucleus A);
//std::string getElementName(hardpingParticle incidentParticle);
double getRandomFromFile();
double getFromWFile();
double getFromBeFile();
 Pythia* pythia8;
 ifstream _randomFile;
 ifstream W_File;
 ifstream Be_File;
 ofstream coordinateSoftOutput;
 ofstream coordinateHardOutput;
 ofstream pathInNucleiOutput;
 ofstream softCollisionsNumberOutput;
 ofstream probabilityOutput;
 ofstream impactParameterFileAfter;
 ofstream impactParameterFileBefore;
 ofstream deltaPtOutput;
 ofstream incidentParticlePt;
 const gsl_rng_type *gslRandomGeneratorType;
 gsl_rng *gslRandomGenerator;

int main() {
char ch;
unsigned int numEvent = 100000000;
Timer time2;
time2.start();
ofstream fileDrellYan;
std::string outputFilename,softCollisionsNumberFilename, randomGeneratorStateFilename,impactParameterAfterFilename,impactParameterBeforeFilename, pathInNucleusFilename, coordinateSoftFilename, coordinateHardFilename, targetElementName,projectileElementName, probabilityOutputFilename, deltaPtOutputFilename;
std::string incidentParticlePtFilename;
int Aproj = 0, Zproj = 0, Atarg = 0, Ztarg = 0, incidentParticleId = 0;
double initialProjectileLabMomentum = 27.6;

incidentParticleId = -11;//если ноль вызает инициализацию ядро-ядро, если нет -частица-ядро.

numEvent = 100000000;

Aproj = 1;
Zproj = 1;


Atarg = 184;
Ztarg = 74;
/*
TMacro m("/home/guest/programs/build/macros_/getEvent.C");
m.Exec();
cin>>ch;
*/
//const gsl_rng_type *gslRandomGeneratorType;
//gsl_rng *gslRandomGenerator;
gslRandomGeneratorType = gsl_rng_default;
gsl_rng_env_setup ();
gslRandomGenerator= gsl_rng_alloc (gslRandomGeneratorType);
_randomFile.open("/home/dsuetin/workspace/Pythia8180/Debug/randomNumbersFile.txt");
//W_File.open("/home/dsuetin/workspace/Hardping_201401/Debug/02.06.14/WnCol_fort.txt");
W_File.open("/home/dsuetin/workspace/Hardping_newold/Debug/08.10.14/W1nCol_fort.txt");
Be_File.open("/home/dsuetin/workspace/Hardping_201401/Debug/02.06.14/BenCol_fort.txt");


nucleus target    (Ztarg, Atarg);
nucleus projectile(Zproj, Aproj);

Hardping* hardping;
pythia8 = new Pythia("/home/guest/programs/pythia8180/xmldoc");


pythia8->rndm.init(19780503);
//pythia->rndm.init();

pythia8->readString("SoftQCD:all = off");
pythia8->readString("HardQCD:all = off");
//pythia8->rndm.readState(constCharStringRandomGeneratorFilename);
pythia8->init();
hardpingParticle* incidentParticle;
if(incidentParticleId){
	Particle newPythiaParticle;
	newPythiaParticle.pz(initialProjectileLabMomentum);
	newPythiaParticle.id(incidentParticleId);
	newPythiaParticle.status(0);
	pythia8->event.append(newPythiaParticle);
	incidentParticle = new hardpingParticle(pythia8->event.back());
	incidentParticle->e(sqrt(incidentParticle->getRestMass()*incidentParticle->getRestMass() + incidentParticle->pAbs2()));
	pythia8->event.clear();
}else{
	incidentParticle = new hardpingParticle();
}

projectile.setInitialProjectileLabMomentum(initialProjectileLabMomentum);

projectileElementName = (incidentParticle->id()==0 )? projectile.getNucleusName():incidentParticle->getParticleName();;
//cout<<" projectileElementName "<<projectileElementName<<endl;
//cin>>ch;
targetElementName = target.getNucleusName();


outputFilename = projectileElementName + "_" + targetElementName + "_"+ "800GeV_DY_" + convertNumberOfEventsToString(numEvent) + ".txt";
cout<< "outputFilename is "<<outputFilename<<endl;
//cin>>ch;

randomGeneratorStateFilename = outputFilename.substr(0,outputFilename.find_last_of("."));
pathInNucleusFilename	   	 = outputFilename.substr(0,outputFilename.find_last_of("."));
coordinateSoftFilename		 = outputFilename.substr(0,outputFilename.find_last_of("."));
coordinateHardFilename 		 = outputFilename.substr(0,outputFilename.find_last_of("."));
probabilityOutputFilename    = outputFilename.substr(0,outputFilename.find_last_of("."));
impactParameterAfterFilename    = outputFilename.substr(0,outputFilename.find_last_of("."));
impactParameterBeforeFilename   = outputFilename.substr(0,outputFilename.find_last_of("."));
softCollisionsNumberFilename	= outputFilename.substr(0,outputFilename.find_last_of("."));
deltaPtOutputFilename			= outputFilename.substr(0,outputFilename.find_last_of("."));
incidentParticlePtFilename		= outputFilename.substr(0,outputFilename.find_last_of("."));

//cout<<randomGeneratorStateFileName<<endl;

randomGeneratorStateFilename += "_randomGeneratorState";
pathInNucleusFilename        += "_pathInNucleus";
pathInNucleusFilename        += ".txt";
coordinateSoftFilename		 += "_coordSoftCollisions";
coordinateSoftFilename		 += ".txt";
coordinateHardFilename		 += "_coordHardCollisions";
coordinateHardFilename		 += ".txt";
probabilityOutputFilename    += "_probability";
probabilityOutputFilename    += ".txt";
impactParameterAfterFilename      += "_impactParameterAfter";
impactParameterAfterFilename      += ".txt";
impactParameterBeforeFilename     += "_impactParameterBefore";
impactParameterBeforeFilename     += ".txt";
softCollisionsNumberFilename	  += "_softCollisionsNumber";
softCollisionsNumberFilename	  += ".txt";
deltaPtOutputFilename			+= "_deltaPt";
deltaPtOutputFilename			+= ".txt";
incidentParticlePtFilename		+= "_incidentParticlePt";
incidentParticlePtFilename		+= ".txt";
const char * constCharStringFilename     		    =  outputFilename.c_str();
const char * constCharStringPathInNucleusFilename   = pathInNucleusFilename.c_str();
const char * constCharStringCoordinateSoftFilename  = coordinateSoftFilename.c_str();
const char * constCharStringCoordinateHardFilename  = coordinateHardFilename.c_str();
const char * constCharStringRandomGeneratorFilename = randomGeneratorStateFilename.c_str();
const char * constCharStringProbabilityOutputFilename = probabilityOutputFilename.c_str();
const char * constCharimpactParameterAfterFilename 		  = impactParameterAfterFilename.c_str();
const char * constCharimpactParameterBeforeFilename 	  = impactParameterBeforeFilename.c_str();
const char * constSoftCollisionsNumberFilename 	  = softCollisionsNumberFilename.c_str();
const char * constDeltaPtOutputFilename 		  = deltaPtOutputFilename.c_str();
const char * constIncidentParticlePtFilename 	  = incidentParticlePtFilename.c_str();
deltaPtOutput.open(constDeltaPtOutputFilename);
coordinateSoftOutput.open(constCharStringCoordinateSoftFilename);
coordinateHardOutput.open(constCharStringCoordinateHardFilename);
probabilityOutput.open(constCharStringProbabilityOutputFilename);
fileDrellYan.open(constCharStringFilename);
pathInNucleiOutput.open(constCharStringPathInNucleusFilename);
impactParameterFileAfter.open(constCharimpactParameterAfterFilename);
impactParameterFileBefore.open(constCharimpactParameterBeforeFilename);
softCollisionsNumberOutput.open(constSoftCollisionsNumberFilename);
incidentParticlePt.open(constIncidentParticlePtFilename);
if(fileDrellYan.is_open()){
	cout<<"all right file "<<constCharStringFilename<<" is opened "<<endl;
//	cin>>ch;
}else{
	cout<<"all not right file "<<constCharStringFilename<<" is not opened "<<endl;
//	cin>>ch;
}
if(pathInNucleiOutput.is_open()){
	cout<<"all right file "<<constCharStringPathInNucleusFilename<<" is opened "<<endl;
//	cin>>ch;
}else{
	cout<<"all not right file "<<constCharStringPathInNucleusFilename<<" is not opened "<<endl;
//	cin>>ch;
}
if(coordinateSoftOutput.is_open()){
	cout<<"all right file "<<constCharStringCoordinateSoftFilename<<" is opened "<<endl;
//	cin>>ch;
}else{
	cout<<"all not right file "<<constCharStringCoordinateSoftFilename<<" is not opened "<<endl;
//	cin>>ch;
}
if(coordinateHardOutput.is_open()){
	cout<<"all right file "<<constCharStringCoordinateHardFilename<<" is opened "<<endl;
//	cin>>ch;
}else{
	cout<<"all not right file "<<constCharStringCoordinateHardFilename<<" is not opened "<<endl;
//	cin>>ch;
}

	//cout<<convertNumberOfEventsToString(numEvent);
	//return 0;

	//pythia8->readString("SoftQCD:all = off");
	//pythia8->readString("HardQCD:all = off");
	pythia8->rndm.readState(constCharStringRandomGeneratorFilename);
//	pythia8->init();


	//настройки пифии в маине передаются в конструкторе хардпинг другой пифие
	//cin>>ch;
	time2.start();
	for(unsigned int iop = 0; iop < numEvent ; iop++){

		//cout<<" numEvent "<<iop<<endl;

		//hardping = new Hardping(projectile,target);
		hardping = (incidentParticle->id()/*incidentParticleId*/ == 0 )? new Hardping(projectile,target) : new Hardping(*incidentParticle,target);
		if(iop==218){
			cout<<"hui";
		}

		hardping->hardping();

		if(iop%100 == 0||1)cout<<"event number "<<iop<<endl;


		for(unsigned int ih = 0; ih < hardping->_finalState->size(); ih++){


			fileDrellYan<<hardping->_finalState->at(ih).id()<<" "<<hardping->_finalState->at(ih).p();
			cout<<"in main n = "<<iop<<" "<<"id = "<<hardping->_finalState->at(ih).id()<<" "<<" isHardron "<<hardping->_finalState->at(ih).isHadron()<<hardping->_finalState->at(ih).p();

		}

		delete hardping;
		//todo delete all dynamic memory

	}


}
/*std::string getElementName(nucleus A){
	std::string	ElementName;

	switch (A.Z()) {
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
				ElementName = "p";
			  break;
			case 2:
				ElementName = "D";
			  break;
			default:
			  cout<<"Error no found such element. "<<endl;
			  ElementName = "";
			  break;
			}
			return ElementName;

}
*/
std::string convertNumberOfEventsToString(unsigned int numEvent){
	std::string str = "10^";
	int countOfDigit = 0;
	char numOfDdigits;
	unsigned int temp = numEvent;
	do{
		temp  = temp/10.;

		countOfDigit++;

	}while(temp > 1);
//	numOfDdigits = (char)countOfDigit;
	numOfDdigits = (char)countOfDigit +48;
//	cout<<"countOfDigit = "<<countOfDigit<<endl;
//	cout<<"numOfDdigits = "<<numOfDdigits<<endl;
	str += numOfDdigits;

	return str;
}

double getRandomFromFile(/*std::ifstream file*/){
		double a;

		_randomFile>>a;
		//cout<<"a!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<a<<endl;
		return a;

}
double getFromWFile(/*std::ifstream file*/){
		double a;

		W_File>>a;
		//cout<<"a!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<a<<endl;
		return a;

}
double getFromBeFile(/*std::ifstream file*/){
		double a;

		Be_File>>a;
		//cout<<"a!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<a<<endl;
		return a;

}
