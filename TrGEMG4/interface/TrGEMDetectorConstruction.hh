
#ifndef TrGEMDetectorConstruction_hh
#define TrGEMDetectorConstruction_hh

#include <list>
#include <string>
#include <cmath>

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4Trd.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

class G4LogicalVolume ;
class G4PhysicalVolume ;

class TrGEMDetectorConstruction : public G4VUserDetectorConstruction
	{
  	public :
    	TrGEMDetectorConstruction() ;
    	virtual ~TrGEMDetectorConstruction() ;

    	void DefineMaterials() ;
    	G4VPhysicalVolume* Construct() ;
    	G4Trd* Trapezoid(G4String name, G4double width, G4double largeBase, G4double thinBase, G4double height) ;
    	void PlaceGeometry(G4int chamberNum, G4RotationMatrix *pRot, G4ThreeVector tlate, G4LogicalVolume *pMotherLogical) ; 
      void ConstructSDandField();

  	private :
  		static const G4int 	numChamber = 2 ;
  		static const G4int 	numSpacer = 4 ;
  		static const G4int 	numShortPad = 6 ;
  		static const G4int  numUnderPipePad = 6 ;
  		static const G4int  numFake = 4 ;
  		static const G4int  buildChamber = 2 ;
  		G4double theta = atan(1283./113.905);
  		
    	G4Material*        	fFR4Mat;
    	G4Material*        	fGasMat;
    	G4Material*        	fEmptyMat;
    	G4Material*        	fAirMat;
    	G4Material*        	fCuMat;
    	G4Material*        	fKaptonMat;
    	G4Material*        	fCoverMat;
    	G4Material*				 	fCoolMat;
    	G4Material*        	fChipMat;
    	G4Material*        	fPullMat;
    	G4double           	driftThinBase = 0.0; 
   	 	G4double           	driftLargeBase = 0.0;
    	G4double           	driftHeight = 0.0;
    	G4double           	foilThinBase = 0.0; 
   	 	G4double           	foilLargeBase = 0.0;
    	G4double           	foilHeight = 0.0;
    	G4double           	gapThinBase = 0.0; 
   	 	G4double           	gapLargeBase = 0.0;
    	G4double           	gapHeight = 0.0;
    	G4double           	spacerThinBase = 0.0; 
   	 	G4double           	spacerLargeBase = 0.0;
    	G4double           	spacerHeight = 0.0;
    	G4double						readoutThinBase = 0.0; 
  		G4double						readoutLargeBase = 0.0 ;		
  		G4double						readoutHeight = 0.0 ;	
  		G4double						zigThinBase = 0.0; 
  		G4double						zigLargeBase = 0.0 ;		
  		G4double						zigHeight = 0.0 ;	
  		G4double						frameInnerThinBase = 0.0; 
  		G4double						frameInnerLargeBase = 0.0 ;		
  		G4double						frameInnerHeight = 0.0 ;	
  		G4double						singleChamberThickness = 0.0 ;
  		G4double						gapInChambers	= 0.0 ;	
  		G4double						detectorThickness	= 0.0 ;		
    	
    	G4Box* 						 	worldBox ;
    	G4LogicalVolume* 	 	worldLog;
	G4Box* 						 	newBox ;
    	G4LogicalVolume* 	 	newLog;
    	
    	G4Trd* 							fake[numChamber][2];
    	G4LogicalVolume* 		fakeLog[numChamber][2] ;
    	
    	G4SubtractionSolid* lateralFake[numChamber];
    	G4LogicalVolume* 		lateralFakeLog[numChamber] ;

    	std::vector<G4Trd*>           trdCollection[numChamber] ;
    	std::vector<G4LogicalVolume*> trdLogCollection[numChamber] ;
	};

#endif

