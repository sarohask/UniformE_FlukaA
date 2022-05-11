
#include "TrGEMSteppingAction.hh"
#include "TrGEMAnalysis.hh"

#include "G4SteppingManager.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4HCtable.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include <string>
#include <sstream>
#include <iostream>

TrGEMSteppingAction::TrGEMSteppingAction() 
	: G4UserSteppingAction()
	{}

TrGEMSteppingAction::~TrGEMSteppingAction() 
	{}
	
void TrGEMSteppingAction::SetSteppingManagerPointer(G4SteppingManager* stManager)
	{
		stepManager = stManager;
		
	}

void TrGEMSteppingAction::UserSteppingAction(const G4Step*)
	{
		//G4cout<<"stacking ............. = "<<stepManager->GetfAlongStepDoItProcTriggered()<<"..process.................... = "<<stepManager->GetfCurrentProcess()->GetProcessName()<<G4endl;
	}


