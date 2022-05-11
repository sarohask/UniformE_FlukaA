
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

class G4LogicalVolume ;
class G4PhysicalVolume ;

class TrGEMDetectorConstruction : public G4VUserDetectorConstruction
	{
  	public :
    	TrGEMDetectorConstruction() ;
    	virtual ~TrGEMDetectorConstruction() ;

    	void DefineMaterials() ;
    	G4VPhysicalVolume* Construct() ;
    	G4Trd* Trapezoid(G4String name, G4double width) ;
    	void PlaceGeometry(G4RotationMatrix *pRot, G4ThreeVector tlate, G4LogicalVolume *pMotherLogical) ; 

  	private :
    	G4Material*        fFR4Mat;
    	G4Material*        fGasMat;
    	G4Material*        fEmptyMat;
    	G4Material*        fAirMat;
    	G4Material*        fCuMat;
    	G4Material*        fKaptonMat;
    	G4double           tripleGemThinBase; 
   	 	G4double           tripleGemLargeBase;
    	G4double           tripleGemHeight;

    	std::vector<G4Trd*>           trdCollection ;
    	std::vector<G4LogicalVolume*> trdLogCollection ;
	};

#endif

