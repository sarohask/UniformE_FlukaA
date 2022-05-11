
#include "TrGEMPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include <TMath.h>
#include <fstream>
#include "TRandom3.h"
#include "TH2.h"

TrGEMPrimaryGeneratorAction::TrGEMPrimaryGeneratorAction(char* partName_) 
	: partName(partName_)
	{
  	G4ParticleDefinition* particlen = G4ParticleTable::GetParticleTable()-> FindParticle(partName);
  	
  	std::string temp = partName;
  	if (temp == "neutron") eneRange = 15;
  	else if (temp == "proton") eneRange = 7;
  	else if (temp == "kaon+") eneRange = 7;
  	else if (temp == "kaon-") eneRange = 7;
  	else if (temp == "pi+") eneRange = 7;
  	else if (temp == "pi-") eneRange = 7;
  	else if (temp == "gamma") eneRange = 7;
  	else if (temp == "e-") eneRange = 6;
  	else if (temp == "e+") eneRange = 6;
  	else if (temp == "alpha") eneRange = 12;
  	else if (temp == "mu-") eneRange = 12;
  	else throw;
		
  	fParticleGun = new G4GeneralParticleSource();
  	fParticleGun->GetCurrentSource()->SetParticleDefinition(particlen);
  	
  	angularFlux = new TFile("/home/sukumar/MyDisk/Simulation_SuperChamber/Fluka_Inputs/New_AngularFlux.root") ;
  	inputFlux = new TFile("/home/sukumar/MyDisk/Simulation_SuperChamber/Fluka_Inputs/UniformYFlux.root") ;

   	TString particle = partName ;
   	h_Angle = (TH1D*)angularFlux->Get("h_Cosine_" + particle) ;
   	//h_Angle = (TH1D*)angularFlux->Get("h_RandomNumber") ;
   	h_InputFlux = (TH1D*)inputFlux->Get("h_RandomNumber") ;
   	
   	h_InputFlux->Fit("pol1","W") ;
   	fitFlux = h_InputFlux->GetFunction("pol1") ;
   	
   	//h_Angle->Fit("pol1", "W") ;
   	//fitAngle = h_Angle->GetFunction("pol1") ;
   	
   	h_Angle->Fit("pol5") ;
    fitAngle = h_Angle->GetFunction("pol5") ;
   	
   	/*if (particle == "neutron") 
   		{
   			h_Angle->Fit("pol5") ;
   			fitAngle = h_Angle->GetFunction("pol5") ;
   		}
   	else if (particle == "proton" || particle == "kaon+" || particle == "kaon+" || particle == "pi+" || particle == "pi-") 
   		{
   			h_Angle = (TH1D*)angularFlux->Get("h_Angle_neutron") ;
   			h_Angle->Fit("pol5") ;
   			fitAngle = h_Angle->GetFunction("pol5") ;
   		}
   	else 
   		{
   			h_Angle->Fit("pol4") ;
   			fitAngle = h_Angle->GetFunction("pol4") ;
   		}*/
	}

TrGEMPrimaryGeneratorAction::~TrGEMPrimaryGeneratorAction()
	{
  	delete fParticleGun;
	}

void TrGEMPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
	{
  	G4double x, y, xxx, yyy, lll, zzz, gX(0.), gY(0.), gZ(0.);
  	G4double cosx(0.), cosy(0.), cosz(0.);
  	G4double tt, kk;
 
 		gRandom->SetSeed(0);
 		
 		//gY = (2*G4UniformRand()-1)*800.0;
 		//gY = (2*G4UniformRand()-1)*641.5;
  	//gY = (2*G4UniformRand()-1)*608.580335699; 							//90% area
  	//gY = (2*G4UniformRand()-1)*672.810876101;							//110% area
 		gY = fitFlux->GetRandom(-641.5, 641.5);
 		xxx = (gY + 2230.67186252)/11.2637724419 ;
  	//xxx = (gY + 2116.20113939)/11.2637724419) ;				//90% area
  	//xxx = (gY + 2339.54838676)/11.2637724419) ;				//110% area
  	//cout<<fitFlux->GetRandom(-1., 1.)<<endl;
  	gX = (2*G4UniformRand()-1)*xxx;
 		
 		//y = G4UniformRand();
 		//if (y < 0.5) 
			{
				gZ = -39.55 ;													//cos (beta)
			}
		//if (y >= 0.5)
			{
				//gZ = -39.55 ;												//cos (beta)
			}
  																					// - for even    + for odd
  		
  	if (gZ < 0.)
  		{
  			gRandom->Uniform();
  	 		//cosz = cos(fitAngle->GetRandom(0., 2.*atan(1))) ;
  	 		//cosz = cos((fitAngle->GetRandom(2.*atan(1), 4.*atan(1)))-(2.*atan(1))) ;
  	 		cosz = fitAngle->GetRandom(0.,1.) ;

  	 		tt = G4UniformRand();
				if (tt < 0.5) 
					{
						cosy = abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;													//cos (beta)
					}
				if (tt >= 0.5)
					{
						cosy = -abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;												//cos (beta)
					}
				kk = G4UniformRand();
				if (kk < 0.5) 
					{
						cosx = abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}
				if (kk >= 0.5)
					{
						cosx = -abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}
			}
			
  	if (gZ > 0.)
  		{
  	 		//cosz = cos(fitAngle->GetRandom(2.*atan(1), 4.*atan(1))) ;
  	 		//cosz = cos((fitAngle->GetRandom(0., 2.*atan(1)))+(2.*atan(1))) ;
  	 		cosz = fitAngle->GetRandom(0.,-1.) ;
  	 		
  	 		tt = G4UniformRand();
				if (tt < 0.5) 
					{
						cosy = abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;													//cos (beta)
					}
				if (tt >= 0.5)
					{
						cosy = -abs(sqrt(1.0 - (cosz*cosz)) * G4UniformRand()) ;												//cos (beta)
					}
				kk = G4UniformRand();
				if (kk < 0.5) 
					{
						cosx = abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}
				if (kk >= 0.5)
					{
						cosx = -abs(sqrt(1.0 - (cosz*cosz) - (cosy*cosy))) ;																	//cos (alpha)
					}	
			}
		fParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(gX*mm, gY*mm, gZ*mm));
		//fParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0.*mm, 0.*mm, 0.*mm));
		fParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(cosx, cosy, cosz));
		//fParticleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0., 0.27973, 0.92175));

  	primaryEne = TMath::Power(10, (G4UniformRand()-1)*eneRange);
  	//primaryEne = 0.1;
  	//fParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(10. * primaryEne*GeV);
  	fParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(1*MeV);
 
  	//create vertex   
  	fParticleGun->GeneratePrimaryVertex(anEvent);
  	
  	primPDG = fParticleGun->GetParticleDefinition()->GetPDGEncoding() ;
  	primEnergy = fParticleGun->GetParticleEnergy() ;
  	//if (gZ < 0.) { primAngle = acos(fParticleGun->GetParticleMomentumDirection().z())*180./M_PI;	}
  	if (gZ < 0.) { primAngle = fParticleGun->GetParticleMomentumDirection().y();	}
  	if (gZ > 0.) { primAngle = acos((-1.)*fParticleGun->GetParticleMomentumDirection().z())*180./M_PI; }								//odd	
  	primMomX = fParticleGun->GetParticleMomentumDirection().x();
  	primMomY = fParticleGun->GetParticleMomentumDirection().y();
  	primMomZ = fParticleGun->GetParticleMomentumDirection().z();
  	primPosX = fParticleGun->GetParticlePosition().x() ;
  	primPosY = fParticleGun->GetParticlePosition().y() ;
  	primPosZ = fParticleGun->GetParticlePosition().z() ;
  	TrGEMAnalysis::GetInstance()->SavePrimary(primEnergy, primAngle, primPDG, primMomX, primMomY, primMomZ, primPosX, primPosY, primPosZ);
 }

