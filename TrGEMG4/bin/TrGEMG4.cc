
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4VModularPhysicsList.hh"
#include "G4RadioactiveDecayPhysics.hh"



#include "GEMDPG/TrGEMG4/interface/TrGEMDetectorConstruction.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMActionInitialization.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMAnalysis.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMPhysicsList.hh"

#include "Randomize.hh"

int main(int argc, char** argv) 
	{
  	// Set the Random engine
  	/*CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    	G4long seed = time(NULL);
    	CLHEP::HepRandom::setTheSeed(seed);
    	long seeds[2];
    	time_t systime = time(NULL);
    	seeds[0] = (long) systime;
    	seeds[1] = (long) (systime * G4UniformRand());
    	G4Random::setTheSeeds(seeds);*/

  	G4Random::setTheEngine(new CLHEP::RanecuEngine());
  	G4Random::setTheSeed(time(NULL)+38999008.);
  
  	G4RunManager* runManager = new G4RunManager;

  	runManager->SetUserInitialization(new TrGEMDetectorConstruction ) ;
  	G4VUserPhysicsList* physics = new FTFP_BERT_HP();
  	//G4VUserPhysicsList* physics = new QGSP_BERT();

  	runManager->SetUserInitialization(physics);
  	runManager->SetUserInitialization(new TrGEMActionInitialization(argv[1]));

  	// initialize G4 kernel
  	runManager->Initialize();

  	// get the pointer to the UI manager
  	G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
 		G4String temp = argv[2];
  	TrGEMAnalysis::GetInstance()->SetFileName(temp);
  	G4String command = "/run/beamOn ";
  	temp = argv[3];
  	UImanager->ApplyCommand("/run/initialize");
 		UImanager->ApplyCommand("/cuts/setLowEdge 10 eV");
  	UImanager->ApplyCommand(command+temp);

  	// job termination
  	delete runManager;
  	return 0;
	}
