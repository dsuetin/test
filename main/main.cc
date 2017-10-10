#include "../include/Pythia8/Pythia.h"
#include "../include/Pythia8/nucleus.h"
#include <sstream>
#include <fstream>
//#include "../Pythia8/timer.h"
#define HardQCD 1
#define SoftQCD 1
#define CMS 0
#define NumberOfEvents 10000000
//#define verbose 0
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
double sNN;
 Pythia* pythia8;

 //huuuuuuuuuuuuuuui
 ifstream _randomFile;
 ifstream W_File;
 ifstream Be_File;
 ofstream coordinateSoftOutput;
 ofstream coordinateHardOutput;
 ofstream pathInNucleiOutput;
 ofstream softCollisionsNumberOutput;
 ofstream probabilityOutput;
 ofstream formationLenghtOutput;
 ofstream impactParameterFileBefore;
 ofstream deltaPtOutput;
 ofstream incidentParticlePt;
 ofstream numberOfSoftCollisions;
 ofstream totalPathOutput;
 ofstream pythiaEventFile;
 ofstream x1File;
 ofstream energyLossFile;
 ofstream energyLossFile2;
 //suetin debug
 ifstream pythia6File;
 ifstream pythia6Z0pFile;
 ifstream pythia6Z0nFile;
 ifstream pythia8ParticleFile;
 ifstream coordinateFile;
 ifstream softCollisionsNumberInput;
 int verbose = 0;
 double kEnergyLoss;
 double quarkNucleonCrossSection;
 bool quickVersion = false;
 int nSoft;
 //const gsl_rng_type *gslRandomGeneratorType;
// gsl_rng *gslRandomGenerator;

int main() {
char ch;
unsigned int numEvent = 10000000;
Timer time2;
time2.start();
std::ostringstream numberToStringConverter;
ofstream fileDrellYan;
ofstream transverseMomentumFile;
std::string transverseMomentumFilename;
std::string outputFilename,softCollisionsNumberFilename, randomGeneratorStateFilename,formationLenghtFilename,impactParameterBeforeFilename, pathInNucleusFilename, coordinateSoftFilename, coordinateHardFilename, targetElementName,projectileElementName, probabilityOutputFilename, deltaPtOutputFilename , x1OutputFilename,x1OutputFilename2;
std::string incidentParticlePtFilename;
std::string initialProjectileLabMomentumString;
std::string numberOfSoftCollisionsFilename;
std::string totalPathFilename;

int Aproj = 0, Zproj = 0, Atarg = 0, Ztarg = 0, incidentParticleId = 0;

double initialProjectileLabMomentum = 120;//27.6;

//incidentParticleId = -11;//если ноль вызает инициализацию ядро-ядро, если нет -частица-ядро.
quickVersion = true;
numEvent =  100000000;
   //numEvent =  11;//0;
   verbose = 0;
   kEnergyLoss = 3.5;//2.5;//2.5;//2.5;//0.1;//0.1;//2.5;
   quarkNucleonCrossSection = 10.;
Aproj = 1;
Zproj = 1;
//Kr(36,84)
//N(7,14)


Atarg = 2;//56;//184;//2;//40;//56;//184;//56;//9;//9;//184;//84;
Ztarg = 1;//74;//20;//26;//74;//26;//36;
//suetin debug

 TString filename = "/home/guest/workspace4/Hardping_newold/Debug/01.06.2015/pythia6event6.txt";
 TString filenameX1 = "/home/guest/workspace4/h/Debug/06.10.2015/x1.txt";
 TString filenameZ0p = "/home/guest/workspace/pythia6/Debug/06.10.2015/pythia6DYevent.txt";



         filenameZ0p = "/home/guest/workspace4/Pythia8/Debug/12.11.2015/p_p_800GeV_DY_ZO.txt";
         filenameZ0p = "/run/media/guest/DATA/E866 PythiaData/p_p_120GeV_DY_ZO.txt";
  TString filenameZ0n = "/home/guest/workspace4/Pythia8/Debug/12.11.2015/p_n_800GeV_DY_ZO.txt";
  	  	  filenameZ0n = "/run/media/guest/DATA/E866 PythiaData/p_n_120GeV_DY_ZO.txt";
 TString filenameQuick = "/home/guest/workspace4/Hardping_newold/Debug/01.06.2015/pythia6event6.txt";

 //filename = "/home/guest/workspace4/Hardping_newold/Debug/01.06.2015/pythia6event_c++.txt";
 TString filenamePythiaEvent = "/home/guest/workspace4/h/Debug/29.07.2015/pythia6event.txt";
pythia6File.open(filename,std::ifstream::binary);
pythia6Z0pFile.open(filenameZ0p,std::ifstream::binary);
pythia6Z0nFile.open(filenameZ0n,std::ifstream::binary);
pythiaEventFile.open(filenamePythiaEvent,std::ofstream::binary);
x1File.open(filenameX1,std::ofstream::binary);
pythia8ParticleFile.open(filenameQuick,std::ofstream::binary);
//pythiaEventFile.open(filename_initialParticlePythiaEvent);
if(pythia6File.is_open()){
	cout<<"ok"<<endl;
}else{
	cout<<"hui"<<endl;
}
if(pythia6Z0pFile.is_open()){
	cout<<"ok"<<endl;
}else{
	cout<<"hui"<<endl;
}
if(pythia6Z0nFile.is_open()){
	cout<<"ok"<<endl;
}else{
	cout<<"hui"<<endl;
}
if(pythia8ParticleFile.is_open()){
	cout<<"ok"<<endl;
}else{
	cout<<"hui, file pythia8ParticleFile do not open "<<endl;
}
int dummy = 0,idPythia6;
double pxPythia6,pyPythia6,pzPythia6,EPythia6;
double hadronFormLenght = 0, preHadronFormLenght = 0, virtualPhotonEnergy = 0, virtualPhotonEnergyOld = 0;
//pythia6File>>idPythia6>>pxPythia6>>pyPythia6>>pzPythia6>>EPythia6>>hadronFormLenght>>preHadronFormLenght>>virtualPhotonEnergy;
//cout<<"id = "<<idPythia6<<" p= "<<pxPythia6<<" "<<pyPythia6<<" "<<pzPythia6<<" "<<EPythia6<<endl;
//cin>>ch;
TString filenameCoordinate = "/home/guest/workspace4/Hardping_newold/Debug/20.04.2015/Kr_coordinateOfHardCollision.txt";
coordinateFile.open(filenameCoordinate,std::ifstream::binary);
if(coordinateFile.is_open()){
	cout<<"ok"<<endl;
}else{
	cout<<"hui"<<endl;
}
/*TString softCollisionsNumberI_initialParticlenputFilename = "/home/guest/workspace4/Hardping_newold/Debug/01.06.2015/numberOfSoftCollisions.txt";
softCollisionsNumberInput.open(softCollisionsNumberInputFilename,std::ifstream::binary);
if(softCollisionsNumberInput.is_open()){
	cout<<"ok"<<endl;
}else{
	cout<<"hui"<<endl;
}*/
//cin>>ch;
numberToStringConverter << initialProjectileLabMomentum;
initialProjectileLabMomentumString = numberToStringConverter.str();



/*
TMacro m("/home/guest/root-build/_macros/getPythia6Event.C");
for(int i =0; i< 10000; i++){
	cout<<"i = "<<i<<endl;
	m.Exec();
}

cin>>ch;
*/
//const gsl_rng_type *gslRandomGeneratorType;
//gsl_rng *gslRandomGenerator;
//gslRandomGeneratorType = gsl_rng_default;
//gsl_rng_env_setup ();
//gslRandomGenerator= gsl_rng_alloc (gslRandomGeneratorType);

//suetin debug

_randomFile.open("/home/guest/workspace4/h/randomNumbersFile.txt"/*"/home/guest/workspace4/Hardping_newold/randomSeq.txt"*/);

//

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
	incidentParticle->m(incidentParticle->getRestMass());
	cout<<"incidentParticle "<<incidentParticle->id()<<" p = "<<incidentParticle->p();
	pythia8->event.clear();
}else{
	incidentParticle = new hardpingParticle();
}

projectile.setInitialProjectileLabMomentum(initialProjectileLabMomentum);

projectileElementName = (incidentParticle->id()==0 )? projectile.getNucleusName():incidentParticle->getParticleName();;
//cout<<" projectileElementName "<<projectileElementName<<endl;
//cin>>ch;
targetElementName = target.getNucleusName();

//todo написать конвертацию энергии в строку + процесс
outputFilename = projectileElementName + "_" + targetElementName + "_"+ initialProjectileLabMomentumString+"GeV_" + convertNumberOfEventsToString(numEvent) + ".txt";
cout<< "outputFilename is "<<outputFilename<<endl;
//cin>>ch;

randomGeneratorStateFilename = outputFilename.substr(0,outputFilename.find_last_of("."));
pathInNucleusFilename	   	 = outputFilename.substr(0,outputFilename.find_last_of("."));
coordinateSoftFilename		 = outputFilename.substr(0,outputFilename.find_last_of("."));
coordinateHardFilename 		 = outputFilename.substr(0,outputFilename.find_last_of("."));
probabilityOutputFilename    = outputFilename.substr(0,outputFilename.find_last_of("."));
formationLenghtFilename      = outputFilename.substr(0,outputFilename.find_last_of("."));
impactParameterBeforeFilename   = outputFilename.substr(0,outputFilename.find_last_of("."));
softCollisionsNumberFilename	= outputFilename.substr(0,outputFilename.find_last_of("."));
deltaPtOutputFilename			= outputFilename.substr(0,outputFilename.find_last_of("."));
incidentParticlePtFilename		= outputFilename.substr(0,outputFilename.find_last_of("."));
numberOfSoftCollisionsFilename  = outputFilename.substr(0,outputFilename.find_last_of("."));
totalPathFilename				= outputFilename.substr(0,outputFilename.find_last_of("."));
transverseMomentumFilename      =  outputFilename.substr(0,outputFilename.find_last_of("."));
x1OutputFilename                =  outputFilename.substr(0,outputFilename.find_last_of("."));
x1OutputFilename2                =  outputFilename.substr(0,outputFilename.find_last_of("."));
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
formationLenghtFilename      += "_formationLenght";
formationLenghtFilename      += ".txt";
impactParameterBeforeFilename     += "_impactParameterBefore";
impactParameterBeforeFilename     += ".txt";
softCollisionsNumberFilename	  += "_summarySoftCollisionsNumber";
softCollisionsNumberFilename	  += ".txt";
deltaPtOutputFilename			+= "_deltaPt";
deltaPtOutputFilename			+= ".txt";
incidentParticlePtFilename		+= "_incidentParticlePt";
incidentParticlePtFilename		+= ".txt";
numberOfSoftCollisionsFilename  += "_numberOfSoftCollisions";
numberOfSoftCollisionsFilename  += ".txt";
totalPathFilename				+= "_totalPath";
totalPathFilename				+= ".txt";
transverseMomentumFilename      += "_transverseMomentum";
transverseMomentumFilename      += ".txt";
x1OutputFilename				+= "_energyLoss";
x1OutputFilename				+= ".txt";
x1OutputFilename2				+= "_energyLoss2";
x1OutputFilename2				+= ".txt";
const char * constCharStringTransverseMomentumFilename  =  transverseMomentumFilename.c_str();
const char * constCharStringFilename     		    =  outputFilename.c_str();
const char * constCharStringPathInNucleusFilename   = pathInNucleusFilename.c_str();
const char * constCharStringCoordinateSoftFilename  = coordinateSoftFilename.c_str();
const char * constCharStringCoordinateHardFilename  = coordinateHardFilename.c_str();
const char * constCharStringRandomGeneratorFilename = randomGeneratorStateFilename.c_str();
const char * constCharStringProbabilityOutputFilename = probabilityOutputFilename.c_str();
const char * constCharFormationLenghtFilename 		  = formationLenghtFilename.c_str();
const char * constCharimpactParameterBeforeFilename 	  = impactParameterBeforeFilename.c_str();
const char * constSoftCollisionsNumberFilename 	  = softCollisionsNumberFilename.c_str();
const char * constDeltaPtOutputFilename 		  = deltaPtOutputFilename.c_str();
const char * constIncidentParticlePtFilename 	  = incidentParticlePtFilename.c_str();
const char * constNumberOfSoftCollisionsFilename	  = numberOfSoftCollisionsFilename.c_str();
const char * constTotalPathFilename	  				  = totalPathFilename.c_str();
const char * constEnergyLossFilename	  			  = x1OutputFilename.c_str();
const char * constEnergyLossFilename2	  			  = x1OutputFilename2.c_str();
deltaPtOutput.open(constDeltaPtOutputFilename);
coordinateSoftOutput.open(constCharStringCoordinateSoftFilename);
coordinateHardOutput.open(constCharStringCoordinateHardFilename);
probabilityOutput.open(constCharStringProbabilityOutputFilename);
fileDrellYan.open(constCharStringFilename);
pathInNucleiOutput.open(constCharStringPathInNucleusFilename);
formationLenghtOutput.open(constCharFormationLenghtFilename);
impactParameterFileBefore.open(constCharimpactParameterBeforeFilename);
softCollisionsNumberOutput.open(constSoftCollisionsNumberFilename);
incidentParticlePt.open(constIncidentParticlePtFilename);
numberOfSoftCollisions.open(constNumberOfSoftCollisionsFilename);
totalPathOutput.open(constTotalPathFilename);
transverseMomentumFile.open(constCharStringTransverseMomentumFilename);
energyLossFile.open(constEnergyLossFilename);
energyLossFile2.open(constEnergyLossFilename2);
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
if(energyLossFile.is_open()){
	cout<<"all right file "<<constEnergyLossFilename<<" is opened "<<endl;
//	cin>>ch;
}else{
	cout<<"all not right file "<<constEnergyLossFilename<<" is not opened "<<endl;
//	cin>>ch;
}
if(energyLossFile2.is_open()){
	cout<<"all right file "<<constEnergyLossFilename2<<" is opened "<<endl;
//	cin>>ch;
}else{
	cout<<"all not right file "<<constEnergyLossFilename2<<" is not opened "<<endl;
//	cin>>ch;
}
double maxPathInNucleus = 0;
	//cout<<convertNumberOfEventsToString(numEvent);
	//return 0;

	//pythia8->readString("SoftQCD:all = off");
	//pythia8->readString("HardQCD:all = off");

	pythia8->rndm.readState(constCharStringRandomGeneratorFilename);
//	pythia8->init();


	//настройки пифии в маине передаются в конструкторе хардпинг другой пифие
	//cin>>ch;
	time2.start();
	hardpingParticle * particleA = new hardpingParticle();
	double zCoordinateOfCollision = 0;
	bool isScattering = false;
	int countCollision = 0;
	int intDummy = 0;
	Vec4 vec4;
	cout<<"vprod "<<particleA->vProd();
	particleA->id(2212);
	//particleA->setHadronFormationLength(1.0);
	//cin>>ch;

	for(unsigned int iop = 1; iop <= numEvent ; iop++){


//	  	_randomFile.open("/home/guest/workspace4/h/randomNumbersFile.txt"/*"/home/guest/workspace4/Hardping_newold/randomSeq.txt"*/);

		if(iop == 81){
			//cout<<"iop = "<<iop<<endl;
			cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<endl;
			//verbose = 1;
	//		cin>>ch;
		}
/*
  	 	cout<<"iop "<<iop<<endl;
		for(int iran = 0; iran < iop ; iran ++){
			if(iran == 0)continue;
			 getRandomFromFile();
			// cout<<"temp1 "<<getRandomFromFile()<<endl;
		}
*/
//		cin>>ch;

		if(verbose)cout<<"iop = "<<iop<<endl;
		if(verbose)cout<<"22222222222222222222222222222222222222222 "<<endl;
		//hardping = new Hardping(projectile,target);
		hardping = (incidentParticle->id()/*incidentParticleId*/ == 0 )? new Hardping(projectile,target) : new Hardping(*incidentParticle,target);


	/*	if(iop == 15){
			cout<<"stop";

		}
		countCollision = 0;
		vec4.pz(-hardping->getMaxZCoordinate());
		vec4.px(1.0);
		vec4.py(0.0);
		particleA->setHadronFormationLength(1.0);
		particleA->setLeftHadronFormationLength(1.0);
		particleA->vProd(vec4);
		particleA->setTotalPathInNucleus(0);
	//
		zCoordinateOfCollision = 0;
		cout<<" particleA "<<particleA->vProd()<<endl;
		cout<<" TotalPathInNucleus "<<particleA->getTotalPathInNucleus()<<endl;
		//cout << incidentParticle->getPythiaParticle()->p();
		//cin>>ch;*/
		//hardping->setVerbose(0);
		//if(iop == 4338)hardping->setVerbose(1);
		softCollisionsNumberInput>>intDummy>>nSoft;
		//cout<<"i = "<<intDummy<<" nSoft "<<nSoft<<endl;
	//	cin>>ch;
		hardping->hardping();
	/*	do{
			cout<<"coord b"<<particleA->vProd();
			isScattering = hardping->pathInNucleus2(particleA,zCoordinateOfCollision);
			if(isScattering){
				countCollision++;
				vec4.pz(zCoordinateOfCollision);
				particleA->vProd(vec4);
				cout<<"coord a"<<particleA->vProd();
			}
		}while(isScattering&&countCollision<(target.A()-1));
		numberOfSoftCollisions<<countCollision<<endl;
		cout<<"countCollision "<<countCollision<<endl;
		totalPathOutput<<particleA->getTotalPathInNucleus()<<endl;
	//	cin>>ch;

		if(iop%100 == 0||1)cout<<"event number "<<iop<<endl;
*/
		Vec4 summMomentum = 0;
		double totalEnergy =0;
		double absSystemMomentum =0;
		for(unsigned int ih = 0; ih < hardping->_finalState->size(); ih++){
			summMomentum += hardping->_finalState->at(ih).getSoftCollisionNumber();
		totalEnergy = summMomentum.e();
		absSystemMomentum = summMomentum.pAbs();
		}
		for(unsigned int ih = 0; ih < hardping->_finalState->size(); ih++){

			cout.precision(3);
			//if(1)cout<<"begin "<<endl;
			if(hardping->_finalState->at(ih).getPreHadronFormationLength() == 0 && hardping->_finalState->at(ih).isHadron()){
				cout<<hardping->_finalState->at(ih).getPreHadronFormationLength()<<endl;
			}
			//if(hardping->_finalState->at(ih).isHadron())fileDrellYan<<iop<<" "<<hardping->_finalState->at(ih).id()<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<" "<<hardping->_finalState->at(ih).getEnergyLoss()<<" "<<hardping->_finalState->at(ih).getPreHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getVirtualPhotonEnergy()<<" "<<hardping->_finalState->at(ih).getHadronEnergyFraction()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().px()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().py()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().pz()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().e()<<" "<<hardping->_finalState->at(ih).p();
			//cout.precision(14);
			fileDrellYan<<iop<<" "<<hardping->_finalState->at(ih).id()<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<" "<<hardping->_finalState->at(ih).getEnergyLoss()<<" "<<hardping->_finalState->at(ih).p();
			if(hardping->_finalState->at(ih).isHadron()||1)cout<<iop<<"huiii"<<hardping->_finalState->at(ih).id()<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<" "<<hardping->_finalState->at(ih).getEnergyLoss()<<" "<<hardping->_finalState->at(ih).getPreHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getVirtualPhotonEnergy()<<" "<<hardping->_finalState->at(ih).getHadronEnergyFraction()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().px()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().py()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().pz()<<" "<< hardping->_finalState->at(ih).getTransferred4Momentum().e()<<" "<<hardping->_finalState->at(ih).p();
			if(verbose)cout<<"tau "<<hardping->_finalState->at(ih).tau()<<endl;
			if(verbose)cout<<"y "<<hardping->_finalState->at(ih).y()<<endl;


		//	cin>>ch;
			//cout<< hardping->_finalState->at(ih).getTransferredCM4Momentum().px()<<endl;
			//cout<< hardping->_finalState->at(ih).getTransferredCM4Momentum().py()<<endl;
			//cout<< hardping->_finalState->at(ih).getTransferredCM4Momentum().pz()<<endl;
			//cout<< hardping->_finalState->at(ih).getTransferredCM4Momentum().e()<<endl;
			if(hardping->_finalState->at(ih).isHadron())cout<<iop<<" "<<hardping->_finalState->at(ih).id()<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<" "<<hardping->_finalState->at(ih).getPreHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getVirtualPhotonEnergy()<<" "<<hardping->_finalState->at(ih).getHadronEnergyFraction()<<" "<<hardping->_finalState->at(ih).p();
			//cout<<hardping->_finalState->at(ih).getTransferred4Momentum();
		//	if(hardping->_finalState->at(ih).isHadron())fileDrellYan<<hardping->_finalState->at(ih).vProd().px()<<" "<<hardping->_finalState->at(ih).vProd().py()<<" "<<hardping->_finalState->at(ih).vProd().pz()<<" "<< dummy<<" "<<'('<<dummy<<')'<<endl;
			if(hardping->_finalState->at(ih).isHadron())coordinateSoftOutput<<hardping->_finalState->at(ih).vProd().px()<<" "<<hardping->_finalState->at(ih).vProd().py()<<" "<<hardping->_finalState->at(ih).vProd().pz()<<" "<< dummy<<" "<<'('<<dummy<<')'<<endl;
			//cout.precision(17);
			//cout<<hardping->_finalState->at(ih).pT()<<endl;
		//	cin>>ch;
			if(hardping->_finalState->at(ih).isHadron())transverseMomentumFile<<iop<<" "<<hardping->_finalState->at(ih).pT()<<endl;
			if(hardping->_finalState->at(ih).isHadron())numberOfSoftCollisions<<iop<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<endl;

			//if(hardping->_finalState->at(ih).isHadron()||1)fileDrellYan<<iop<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<endl;
			if(hardping->_finalState->at(ih).isHadron())cout<<iop<<" "<<hardping->_finalState->at(ih).id()<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<" "<<hardping->_finalState->at(ih).getPreHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getHadronFormationLength()<<" "<<hardping->_finalState->at(ih).getVirtualPhotonEnergy()<<" "<<hardping->_finalState->at(ih).getHadronEnergyFraction()<<" "<<hardping->_finalState->at(ih).p();
			if(iop == 20){
				//cin>>ch;
				//cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<endl;
			}
		//	cout<<"iop "<<iop<<" nsc "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<endl;
	//		cout<<" hadr "<<hardping->_finalState->at(ih).isHadron()<<endl;
 	//	cin>>ch;
			if(hardping->_finalState->at(ih).getSoftCollisionNumber()){

			}
			if(hardping->_finalState->at(ih).isHadron()){

				if(verbose)cout<<"kopjj "<<sqrt(target.getNuclearRadius()*target.getNuclearRadius() - hardping->_finalState->at(ih).vProd().pT2()) - hardping->_finalState->at(ih).vProd().pz() - hardping->_finalState->at(ih).getTotalPathInNucleus()<<endl;
				maxPathInNucleus = (target.getNuclearRadius()*target.getNuclearRadius() > hardping->_finalState->at(ih).vProd().pT2() )? sqrt(target.getNuclearRadius()*target.getNuclearRadius() - hardping->_finalState->at(ih).vProd().pT2()) : 0;

				/*if(hardping->_finalState->at(ih).getHadronEnergyFraction() < 0.1){
					cout<<"i = "<<iop<<" z "<<hardping->_finalState->at(ih).getHadronEnergyFraction()<<endl;
					cin>>ch;
				}*/
				if(verbose)cout<<" imp "<<hardping->_finalState->at(ih).vProd().pT()<<"  z coorf "<<hardping->_finalState->at(ih).vProd().pz()<<" maxPathInNucleus "<<maxPathInNucleus<<endl;
				if(maxPathInNucleus)cout<<maxPathInNucleus - hardping->_finalState->at(ih).vProd().pz() - hardping->_finalState->at(ih).getTotalPathInNucleus()<<endl;

				//if(maxPathInNucleus)totalPathOutput<<maxPathInNucleus - hardping->_initialParticle.vProd().pz() /*- hardping->_finalState->at(ih).vProd().pz()/*- hardping->_finalState->at(ih).getTotalPathInNucleus()*/- hardping->_finalState->at(ih).getHadronFormationLength()- hardping->_finalState->at(ih).getPreHadronFormationLength()<<endl;
				//if(hardping->_finalState->at(ih).getEnergyLoss())totalPathOutput<<iop<<" "<<hardping->_finalState->at(ih).getEnergyLoss()<<endl;
				if(hardping->_finalState->at(ih).getEnergyLoss()||1)totalPathOutput<<hardping->_finalState->at(ih).getEnergyLoss()<<endl;
				if(maxPathInNucleus)cout<<"el "<<hardping->_finalState->at(ih).getEnergyLoss()<<endl;

				//cout.precision(12);
				if(maxPathInNucleus)cout<<"path1 "<<maxPathInNucleus - hardping->_initialParticle.vProd().pz() /*- hardping->_finalState->at(ih).vProd().pz()/*- hardping->_finalState->at(ih).getTotalPathInNucleus()*/- hardping->_finalState->at(ih).getHadronFormationLength()- hardping->_finalState->at(ih).getPreHadronFormationLength()<<endl;

				cout<<"_initialParticle.vProd().pz() "<<maxPathInNucleus - hardping->_initialParticle.vProd().pz()<<" _finalState->at(ih).vProd().pz() "<<maxPathInNucleus - hardping->_finalState->at(ih).vProd().pz()<<" getTotalPathInNucleus() "<<hardping->_finalState->at(ih).getTotalPathInNucleus()<<endl;

				if(maxPathInNucleus)cout<<"path2 " <<maxPathInNucleus - hardping->_finalState->at(ih).vProd().pz() +  hardping->_finalState->at(ih).getTotalPathInNucleus() - hardping->_finalState->at(ih).getHadronFormationLength()<<endl;
				if(maxPathInNucleus)cout<<"path3 " <<hardping->_finalState->at(ih).getEnergyLoss()<<endl;
			//	cin>>ch;
				if(verbose)cout<<hardping->_finalState->at(ih).vProd().pz()<<" "<<hardping->_initialParticle.vProd().pz()<<" phfl "<<hardping->_finalState->at(ih).getPreHadronFormationLength()<<endl;
				cout.precision(12);
				if(verbose)cout<<"hfl "<<hardping->_finalState->at(ih).getHadronFormationLength()<<" vfe "<<hardping->_finalState->at(ih).getVirtualPhotonEnergy()<<endl;
				if(hardping->_finalState->at(ih).getTotalPathInNucleus()){
					cout<<"energy loss "<<hardping->_finalState->at(ih).getTotalPathInNucleus() - hardping->_finalState->at(ih).getHadronFormationLength()<<" totpath "<<hardping->_finalState->at(ih).getTotalPathInNucleus()<< endl;
					cout<<" total path 2 "<<hardping->_finalState->at(ih).getTotalPathInNucleus()<<endl;
			//		cin>>ch;
				}
			//	cout<<"i = "<<iop<<" en lose "<<hardping->_finalState->at(ih).getEnergyLoss()<<endl;
			//	cin>>ch;
		if(hardping->_finalState->at(ih).getEnergyLoss()&&verbose){
			cout<<"i "<<iop<<endl;
			cout<<"getEnergyLoss "<<hardping->_finalState->at(ih).getEnergyLoss()<<endl;
			cout<<"p "<<hardping->_finalState->at(ih).p()<<" maxPathInNucleus "<<maxPathInNucleus<<" coord "<< hardping->_finalState->at(ih).vProd()<<endl;
			cout<<"getTotalPathInNucleus "<<hardping->_finalState->at(ih).getTotalPathInNucleus()<<endl;
		//	cin>>ch;
		}
			}
		//	cin>>ch;
			//if(hardping->_finalState->at(ih).isHadron())        cout<<iop<<" "<<hardping->_finalState->at(ih).id()<<" "<<hardping->_finalState->at(ih).getSoftCollisionNumber()<<" "<<hardping->_finalState->at(ih).getHardCollisionNumber()<<" "<<hardping->_finalState->at(ih).getVirtualPhotonEnergy()<<" "<<hardping->_finalState->at(ih).getHadronEnergyFraction()<<" "<<hardping->_finalState->at(ih).p();
			//cout<<"in main n = "<<iop<<" "<<"id = "<<hardping->_finalState->at(ih).id()<<" "<<" isHardron "<<hardping->_finalState->at(ih).isHadron()<<hardping->_finalState->at(ih).p();
			//if(1)cout<<"end "<<endl;

		}
		//suetin debug
//	  	_randomFile.close();
		delete hardping;
		//todo delete all dynamic memory

	}

	delete particleA;
	delete pythia8;
	delete incidentParticle;
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
