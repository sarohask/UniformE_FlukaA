
#include "GEMDPG/TrGEMG4/interface/TrGEMDetectorConstruction.hh"
#include "GEMDPG/TrGEMG4/interface/TrGEMSensitiveDetector.hh"

#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4Para.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include "G4EqMagElectricField.hh"
#include "G4SystemOfUnits.hh"

#include "G4UniformElectricField.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"

TrGEMDetectorConstruction::TrGEMDetectorConstruction() 
	: fFR4Mat(0), fGasMat(0), fEmptyMat(0), fAirMat(0), fCuMat(0),fKaptonMat(0), fCoverMat(0), fCoolMat(0), fChipMat(0), fPullMat(0)    
	{
  	driftThinBase  = 282.174*mm ;					 						
  	driftLargeBase = 509.984*mm ;					
  	driftHeight    = 1283.0*mm ;																				
  	
  	foilThinBase  = 231.186*mm ;
  	foilLargeBase = 445.858*mm ;
  	foilHeight    = 1205.5*mm ;
  	
  	gapThinBase = foilThinBase + 2.*1.83030604881*mm ;
  	gapLargeBase = foilLargeBase + 2.*2.18542685941*mm ;
  	gapHeight = foilHeight + 4.*mm ;																									//2 mm each for large and thin
  	
  	spacerThinBase = gapThinBase + 2*4.57576512202*mm ;
  	spacerLargeBase = gapLargeBase + 2.*5.46356714852*mm ;
  	spacerHeight = gapHeight + 10.*mm ;																								//5 mm each for large and thin
  	
  	
  	readoutThinBase = spacerThinBase + (2.*9.151530244)*mm ;
  	readoutLargeBase = spacerLargeBase + (2.*10.927134297)*mm ;
  	readoutHeight = spacerHeight + 20.*mm ;																						//10 mm each for large and thin
  	
  	frameInnerThinBase = readoutThinBase - (2.*7.32122419523)*mm ;
  	frameInnerLargeBase = readoutLargeBase - (2.*8.74170743763)*mm ;
  	frameInnerHeight = readoutHeight - 16.*mm ;																				//8 mm each for large and thin
  	
  	zigThinBase = driftThinBase - (2.*9.151530244)*mm ;
  	zigLargeBase = driftLargeBase - (2.*10.927134297)*mm ;
  	zigHeight = driftHeight - 20.*mm ;																								//10 mm each for large and thin
	}

TrGEMDetectorConstruction::~TrGEMDetectorConstruction() 
	{}

void TrGEMDetectorConstruction::DefineMaterials() 
	{
  	G4NistManager* manager = G4NistManager::Instance() ;
  	// define Elements
  	G4Element* elC  = manager->FindOrBuildElement(6);
  	G4Element* elF  = manager->FindOrBuildElement(9);
  	G4Element* elSi = manager->FindOrBuildElement(14);
  	G4Element* elO  = manager->FindOrBuildElement(8);
  	G4Element* elH  = manager->FindOrBuildElement(1);

  	// define Materials
 	 	G4Material *Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si") ;
  	G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu") ;
  	G4Material *Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al") ;
 	 	G4Material *Kapton = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
 	 	G4Material *Water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
 	 	G4Material* Air  = manager->FindOrBuildMaterial("G4_AIR");
 	 	G4Material* Steel  = manager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  	
  	G4int nElement(0), nAtoms(0) ;
  	G4double density(0.), temperature(0.), pressure(0.), fractionMass(0.), mixtureDensity(0.0)  ;
  	G4String name, symbol ;
  	
  	//Quartz
  	G4Material* SiO2 =  new G4Material("quartz",density= 2.200*g/cm3, nElement=2);
  	SiO2->AddElement(elSi, nAtoms=1);
  	SiO2->AddElement(elO , nAtoms=2);

  	//from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
  	//Epoxy (for FR4 )
  	G4Material* Epoxy = new G4Material("Epoxy" , density = 1.2*g/cm3, nElement=2);
  	Epoxy->AddElement(elH, nAtoms=2);
  	Epoxy->AddElement(elC, nAtoms=2);
      
  	//FR4 (Glass + Epoxy)
  	G4Material* FR4 = new G4Material("FR4"  , density = 1.86*g/cm3, nElement=2);
  	FR4->AddMaterial(Epoxy, fractionMass=0.472);
  	FR4->AddMaterial(SiO2, fractionMass=0.528);
  	
  	// gases at STP conditions 
  	G4Material* Argon = manager->FindOrBuildMaterial("G4_Ar");
  	G4Material* CarbonDioxide = manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  	G4Material* Empty = manager->FindOrBuildMaterial("G4_Galactic");

  	// CF4 must be defined by hand
  	G4Material* CF4 = new G4Material(name="CF4", density=0.003884*g/cm3, nElement=2, kStateGas, temperature = 273.15*kelvin, pressure=1.0*atmosphere);
 	 	CF4->AddElement(elC, 1) ;
  	CF4->AddElement(elF, 4) ; 

  	// Ar:CO2 (70:30) @ STP conditions
  	mixtureDensity = (Argon->GetDensity() * 70/100.0 + CarbonDioxide->GetDensity() * 30/100.0) ;
  	G4Material *ArCO2 = new G4Material("Ar/CO2",mixtureDensity,2) ;
  	ArCO2->AddMaterial(Argon, 0.7) ;
  	ArCO2->AddMaterial(CarbonDioxide, 0.3) ;

  	// Ar:CO2:CF4 (45:15:40) @ STP conditions
  	mixtureDensity = (Argon->GetDensity() * 45/100.0 + CarbonDioxide->GetDensity() * 15/100.0 + CF4->GetDensity() * 40/100.0) ;
  	G4Material *ArCO2CF4 = new G4Material("Ar/CO2/CF4",mixtureDensity,3) ;
  	ArCO2CF4->AddMaterial(Argon, 0.45) ;
  	ArCO2CF4->AddMaterial(CarbonDioxide,0.15) ;
  	ArCO2CF4->AddMaterial(CF4,0.40) ;

  	// Choice of the gas
  	fEmptyMat = Empty ;
  	fCuMat = Cu;
  	fKaptonMat = Kapton;
  	fGasMat = ArCO2;
  	fAirMat = Air ;
  	fFR4Mat = FR4;
  	fCoverMat = Al ;
  	fChipMat = Si ;
  	fCoolMat = Water ;
  	fPullMat = Steel ;
	}

G4VPhysicalVolume* TrGEMDetectorConstruction::Construct() 
	{
  	G4GeometryManager::GetInstance()->OpenGeometry();
  	
  	// Define all materials and set global variables
  	DefineMaterials() ;
  	
  	// Visual attributes
  	G4VisAttributes *driftAttributes = new G4VisAttributes(G4Colour(0., 0.6, 0.1)) ; 										//Green
  	driftAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *gasAttributes = new G4VisAttributes(G4Colour(1., 0.1, 0.9)) ;										//Magenta
  	gasAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *spacerAttributes = new G4VisAttributes(G4Colour(1., 1., 0.)) ;										//Magenta
  	spacerAttributes->SetForceSolid(true) ;
  	G4VisAttributes *gemAttributes = new G4VisAttributes(G4Colour(0.8, 0.8, 0.1)) ;										//Light Brown
  	gemAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *kaptonAttributes = new G4VisAttributes(G4Colour(0.7, 0.5, 0.)) ;										//Dark Brown
  	kaptonAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *readoutAttributes = new G4VisAttributes(G4Colour(0., 0.6, 0.1)) ; 										//Green
  	readoutAttributes->SetForceWireframe(true) ;
  	G4VisAttributes *vFatAttributes = new G4VisAttributes(G4Colour(0., 0.6, 0.1)) ;
   	vFatAttributes->SetForceSolid(true) ;
   	G4VisAttributes *zigAttributes = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75)) ;
   	zigAttributes->SetForceSolid(true) ;
   	G4VisAttributes *pullOutAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;
   	pullOutAttributes->SetForceSolid(true) ;
   	G4VisAttributes *frameAttributes = new G4VisAttributes(G4Colour(0.3, 0.8, 0.8)) ;
   	frameAttributes->SetForceSolid(true) ;
   	G4VisAttributes *copperAttributes = new G4VisAttributes(G4Colour(0.8, 0.8, 0.1));
  	copperAttributes->SetForceWireframe(true);
  	G4VisAttributes *emptyAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0)) ;
   	emptyAttributes->SetForceWireframe(true) ;
   	G4VisAttributes *fakeAttributes = new G4VisAttributes(G4Colour(0., 0., 1.)) ;
  	fakeAttributes->SetForceSolid(true) ;
  	G4VisAttributes *coolingPadAttributes = new G4VisAttributes(G4Colour(0.8, 0.8, 0.1)) ;
   	coolingPadAttributes->SetForceSolid(true) ;
  	G4VisAttributes *coolingAttributes = new G4VisAttributes(G4Colour(0.9, 0.6, 0.)) ;
   	coolingAttributes->SetForceWireframe(true) ;
   	G4VisAttributes *coverAttributes = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75)) ;
  	coverAttributes->SetForceSolid(true) ;
  	G4VisAttributes *worldAttributes = new G4VisAttributes ;
  	worldAttributes->SetVisibility(false) ;
	
  	// rotation Matrix for layers
  	G4RotationMatrix* rotationPlacement = new G4RotationMatrix() ;
  	rotationPlacement->rotateY(M_PI / 2.0) ;
  	rotationPlacement->rotateX(M_PI / 2.0) ;
		
		//roatation matrix for fake lateral / external frame
  	G4RotationMatrix* rotationPlacement1 = new G4RotationMatrix() ;
   	rotationPlacement1->rotateY(M_PI / 2.0) ;
   	rotationPlacement1->rotateX(M_PI / 2.0) ;
   	
   	// Rotation Matrix for bend pipes
   	G4RotationMatrix* rotationPlacement2 = new G4RotationMatrix() ;
  	rotationPlacement2->rotateX(M_PI ) ;
   
   //rotaion matrix for staright pipes
   	G4RotationMatrix* rotationPlacement3 = new G4RotationMatrix() ;
   	rotationPlacement3->rotateY(0.) ;
   	rotationPlacement3->rotateX(M_PI/2.0 ) ;
   	rotationPlacement3->rotateZ(0.) ;
   	
   	//rotation matrix for vfat
   	G4RotationMatrix* rotationPlacement4 = new G4RotationMatrix() ;  	
   	rotationPlacement4->rotateZ(M_PI/2.0 ) ;
   	rotationPlacement4->rotateY(M_PI/2.0) ;
   	
   	//for frame lateral side
   	G4RotationMatrix* rotationPlacement5 = new G4RotationMatrix() ;  	
   	rotationPlacement5->rotateX(M_PI - atan(1283.0/113.905) ) ;
   	
   	//for frame lateral side
   	G4RotationMatrix* rotationPlacement6 = new G4RotationMatrix() ;  	
   	rotationPlacement6->rotateX(-M_PI + atan(1283.0/113.905) ) ;
   	
   	//for pull-outs lateral side
   	G4RotationMatrix* rotationPlacement7 = new G4RotationMatrix() ;  	
   	rotationPlacement7->rotateZ(M_PI - atan(1283.0/113.905) ) ;
   	
   	//for pull-outs lateral side
   	G4RotationMatrix* rotationPlacement8 = new G4RotationMatrix() ;  	
   	rotationPlacement8->rotateZ(-M_PI + atan(1283.0/113.905) ) ;
		
		G4RotationMatrix* rotationPlacement9 = new G4RotationMatrix() ;  	
   	rotationPlacement9->rotateZ(M_PI )  ;
   	
   	G4RotationMatrix* rotationPlacement10 = new G4RotationMatrix() ;  	
   	rotationPlacement10->rotateZ(M_PI+0.959931 )  ;
   	
   	G4RotationMatrix* rotationPlacement11 = new G4RotationMatrix() ;  	
   	rotationPlacement11->rotateZ(0.959931 )  ;
   	
   	G4RotationMatrix* rotationPlacement12 = new G4RotationMatrix() ;  
   	rotationPlacement12->rotateX(M_PI/2. )  ;	
   	rotationPlacement12->rotateY(M_PI/2. -0.60434096913)  ;
   	
   	G4RotationMatrix* rotationPlacement13 = new G4RotationMatrix() ;  
   	rotationPlacement13->rotateX(M_PI/2. )  ;	
   	rotationPlacement13->rotateY(M_PI/2. + 0.60434096913)  ;
   	
		//world
  	G4double worldSizeX = 1.0*m;
  	G4double worldSizeY = 1.0*m;
  	G4double worldSizeZ = 1.0*m;
  	
  	//world definition and placement
  	worldBox = new G4Box("WorldBox", worldSizeX, worldSizeY, worldSizeZ) ;
  	worldLog = new G4LogicalVolume( worldBox, fAirMat, "WorldLog") ;
  	worldLog->SetVisAttributes(worldAttributes) ;
		G4VPhysicalVolume* worldPhys = new G4PVPlacement(0, G4ThreeVector(), worldLog, "WorldSpace", 0, false, 0) ;
   
  	std::string layerNames[30]= { "DriftCopper1","DriftBoard","DriftCopper2",     			 	//Drift Board
    															"GasGap1",                                      				//Drift Gap
    															"Gem1Copper1","Gem1","Gem1Copper2",             			 	//GEM1
    															"GasGap2",                                      				//Transfer I Gap
    															"Gem2Copper1","Gem2","Gem2Copper2",             			 	//GEM2
    															"GasGap3",                                      				//Transfer II Gap
    															"Gem3Copper1","Gem3","Gem3Copper2",             			 	//GEM3
    															"GasGap4",                                      				//Induction Gap
    															"ReadCopper1","ReadoutBoard","ReadCopper2",     			 	//Readout Board
    															"GebCopper1",																						//GEB1
    															"GebInsulator1", 																				//GEB2
    															"GebRouting1", 																					//GEB3
    															"GebInsulator2", 																				//GEB4
    															"GebCopper2",  																					//GEB5
    															"GebInsulator3", 																				//GEB6
    															"GebCopper3",  																					//GEB7
    															"GebInsulator4", 																				//GEB8
    															"GebRouting2", 																					//GEB9
    															"GebInsulator5", 																				//GEB10
    															"GebCopper4"		 																				//GEB11																			 	
  															};
             
  	std::string layerNamesLog[30];
  	for(size_t i=1; i<30; i++) 
  		{ 
    		layerNamesLog[i]=layerNames[i]+"Log";
			}
		
  	G4Material* layerMat[30]= { fCuMat,fFR4Mat,fCuMat,      																//Drift Board
    														fGasMat,                    																//Drift Gap
    														fCuMat,fKaptonMat,fCuMat,   																//GEM1
    														fGasMat,                    																//Transfer I Gap
    														fCuMat,fKaptonMat,fCuMat,   																//GEM2
    														fGasMat,                    																//Transfer II Gap
    														fCuMat,fKaptonMat,fCuMat,   																//GEM3
    														fGasMat,                    																//Induction Gap
    													  fCuMat,fFR4Mat,fCuMat,      																//Readout Board
    													  fCuMat,																											//GEB1
    													  fFR4Mat, 																										//GEB2
    													  fFR4Mat,																										//GEB3
    													  fFR4Mat, 																										//GEB4
    													  fCuMat, 																										//GEB5
    													  fFR4Mat, 																										//GEB6
    													  fCuMat, 																										//GEB7
    													  fFR4Mat,																										//GEB8
    													  fFR4Mat, 																										//GEB9
    													  fFR4Mat,																										//GEB10
    													  fCuMat																											//GEB11
  														};

  	G4double thickness[30] = { 35.*um,3.2*mm,35.*um,      																	//Drift Board
    													 4.*mm,                     																	//Drift Gap
    													 5.*um,50*um,5.*um,         																	//gem1
    													 1.*mm,                     																	//Transfer I Gap
    													 5.*um,50*um,5.*um,         																	//gem2
    													 2.*mm,                     																	//Transfer II Gap
    													 5.*um,50.*um,5.*um,        																	//gem3
    													 1.*mm,                     																	//Induction Gap
    													 35.*um,3.2*mm,35.*um,      																	//Readout Board
    													 35.*um,																											//GEB1
    													 218*um, 																											//GEB2
    													 17.5*um, 																										//GEB3
    													 150*um, 																											//GEB4
    													 17.5*um, 																										//GEB5
    													 120*um,																											//GEB6
    													 17.5*um, 																										//GEB7
    													 150*um,																											//GEB8
    													 17.5*um, 																										//GEB9
    													 218*um,																											//GEB10
    													 35.*um 									 																		//GEB11			
  													 };
  																
  	G4double layerLargeBase[30] = { driftLargeBase,driftLargeBase,foilLargeBase,      			//Drift Board
    																gapLargeBase,                     											//Drift Gap
    																foilLargeBase,spacerLargeBase,foilLargeBase,         		//gem1
    																gapLargeBase,                     											//Transfer I Gap
    																foilLargeBase,spacerLargeBase,foilLargeBase,         		//gem2
    																gapLargeBase,                    												//Transfer II Gap
    																foilLargeBase,spacerLargeBase,foilLargeBase,        		//gem3
    																gapLargeBase,                     											//Induction Gap
    																foilLargeBase,readoutLargeBase,readoutLargeBase,     		//Readout Board
    																readoutLargeBase, 									 										//GEB
    																readoutLargeBase, readoutLargeBase,											//GEB
    																readoutLargeBase, 																			//GEB
    																readoutLargeBase, readoutLargeBase,											//GEB
    																readoutLargeBase,	readoutLargeBase,										  //GEB
    																readoutLargeBase,readoutLargeBase,										  //GEB
    																readoutLargeBase																				//GEB 
  																};
  	
  	G4double layerThinBase[30] = { driftThinBase,driftThinBase,foilThinBase,      					//Drift Board
    															 gapThinBase,                     												//Drift Gap
    															 foilThinBase,spacerThinBase,foilThinBase,         				//gem1
    															 gapThinBase,                     												//Transfer I Gap
    															 foilThinBase,spacerThinBase,foilThinBase,         				//gem2
    															 gapThinBase,                     												//Transfer II Gap
    															 foilThinBase,spacerThinBase,foilThinBase,        				//gem3
    															 gapThinBase,                     												//Induction Gap
    															 foilThinBase,readoutThinBase,readoutThinBase,     				//Readout Board
    															 readoutThinBase, 									 											//GEB
    															 readoutThinBase, readoutThinBase,												//GEB
    															 readoutThinBase, 																				//GEB
    															 readoutThinBase, readoutThinBase,												//GEB
    															 readoutThinBase, readoutThinBase,												//GEB
    															 readoutThinBase,readoutThinBase,													//GEB
    															 readoutThinBase																					//GEB 
  															 };
  	
  	G4double layerHeight[30] = { driftHeight,driftHeight,foilHeight,      									//Drift Board
    														 gapHeight,                     														//Drift Gap
    														 foilHeight,spacerHeight,foilHeight,         								//gem1
    														 gapHeight,                     														//Transfer I Gap
    														 foilHeight,spacerHeight,foilHeight,         								//gem2
    														 gapHeight,                     														//Transfer II Gap
    														 foilHeight,spacerHeight,foilHeight,        								//gem3
    														 gapHeight,                   															//Induction Gap
    														 foilHeight,readoutHeight,readoutHeight,     								//Readout Board
    														 readoutHeight, 									 													//GEB
    														 readoutHeight, readoutHeight,															//GEB
    														 readoutHeight, 																						//GEB
    														 readoutHeight, readoutHeight,															//GEB
    														 readoutHeight, readoutHeight,															//GEB
    														 readoutHeight,readoutHeight,																//GEB
    														 readoutHeight																							//GEB 
  														};
  	
  	G4VisAttributes *visAttributes[30] = { copperAttributes, driftAttributes, copperAttributes,					//Drift Board
  																		 		 gasAttributes,																								//Drift Gap
  																		 		 copperAttributes, gemAttributes, copperAttributes,						//gem1
  																		 		 gasAttributes,																								//Transfer I Gap
  																		 		 copperAttributes, gemAttributes, copperAttributes,						//gem2
  																		 		 gasAttributes,																								//Transfer II Gap
  																		 		 copperAttributes, gemAttributes, copperAttributes,						//gem3
  																		 		 gasAttributes,																								//Induction Gap
  																		 		 copperAttributes, readoutAttributes, copperAttributes,				//Readout Board
  																		 		 copperAttributes, 																						//GEB1
  																		 		 readoutAttributes, 																					//GEB2
  																		 		 readoutAttributes,																						//GEB3
  																		   	 readoutAttributes, 																					//GEB4
  																		 		 copperAttributes, 																						//GEB5
  																		 		 readoutAttributes,																						//GEB6
  																		 		 copperAttributes, 																						//GEB7
  																		 		 readoutAttributes,																						//GEB8
  																		 	 	 readoutAttributes, 																					//GEB9
  																		 	 	 readoutAttributes,																						//GEB10 
  																		 		 copperAttributes																							//GEB11
  																			};
  																			
  	// Chamber values
  	double chamberZPos[numChamber] = { -36.55, 1.85 };
  	
  	// Fake values
  	double fakeThickness = 0.1*mm ;
  	double fakeZPos[numChamber][numFake] = {{-36.6*mm, -1.8*mm}, {1.8*mm, 36.6*mm}} ;
  	
  	// Lateral fake values
  	double lateralFakeThick = 34.9*mm ;
		double lateralFakeZPos[numChamber] = {-19.2*mm, 19.2*mm} ;
  	
  	// Spacers values			
  	double spacerThickness[numSpacer] = {3.*mm, 1.*mm, 2.*mm, 1.*mm} ;
  	double spacerPosition[numChamber][numSpacer] = {{-31.815*mm, -29.755*mm, -28.195*mm, -26.635*mm}, {6.585*mm, 8.645*mm, 10.205*mm, 11.765*mm}} ;          	//symmetric
  	
  	// Pull-outs values
  	G4double pullOutX = 19.*mm;
		G4double pullOutY = 4.*mm ;
		G4double pullOutZ = 7.*mm ;
																																								
  	double pulloutXPosition[numChamber][64] = {	{-196.*mm, -147.*mm, -98.*mm, -49.*mm, 0., 49.*mm, 98*mm, 147.*mm, 196.*mm, -98.*mm, -49.*mm, 0., 49.*mm, 98.*mm, 
  																							(180.296366454-588.*cos(theta))*mm, (180.296366454-539.*cos(theta))*mm, (180.296366454-490.*cos(theta))*mm,
  																							(180.296366454-441.*cos(theta))*mm, (180.296366454-392.*cos(theta))*mm, (180.296366454-343.*cos(theta))*mm, 
  																							(180.296366454-294.*cos(theta))*mm, (180.296366454-245.*cos(theta))*mm, (180.296366454-196.*cos(theta))*mm, 
  																							(180.296366454-147.*cos(theta))*mm, (180.296366454-98.*cos(theta))*mm, (180.296366454-49.*cos(theta))*mm, 
  																							(180.296366454)*mm, 
  																							(180.296366454+49.*cos(theta))*mm, (180.296366454+98.*cos(theta))*mm, (180.296366454+147.*cos(theta))*mm, 
  																							(180.296366454+196.*cos(theta))*mm, (180.296366454+245.*cos(theta))*mm, (180.296366454+294.*cos(theta))*mm, 
  																							(180.296366454+343.*cos(theta))*mm, (180.296366454+392.*cos(theta))*mm, (180.296366454+441.*cos(theta))*mm, 	
  																							(180.296366454+490.*cos(theta))*mm, (180.296366454+539.*cos(theta))*mm, (180.296366454+588.*cos(theta))*mm, 
  																							-(180.296366454-588.*cos(theta))*mm, -(180.296366454-539.*cos(theta))*mm, -(180.296366454-490.*cos(theta))*mm,
  																							-(180.296366454-441.*cos(theta))*mm, -(180.296366454-392.*cos(theta))*mm, -(180.296366454-343.*cos(theta))*mm,
  																							-(180.296366454-294.*cos(theta))*mm, -(180.296366454-245.*cos(theta))*mm, -(180.296366454-196.*cos(theta))*mm,
  																							-(180.296366454-147.*cos(theta))*mm, -(180.296366454-98.*cos(theta))*mm, -(180.296366454-49.*cos(theta))*mm, 
  																							-180.296366454*mm, 
  																							-(180.296366454+49.*cos(theta))*mm, -(180.296366454+98.*cos(theta))*mm, -(180.296366454+147.*cos(theta))*mm,
  																							-(180.296366454+196.*cos(theta))*mm, -(180.296366454+245.*cos(theta))*mm, -(180.296366454+294.*cos(theta))*mm, 
  																							-(180.296366454+343.*cos(theta))*mm, -(180.296366454+392.*cos(theta))*mm, -(180.296366454+441.*cos(theta))*mm, 
  																							-(180.296366454+490.*cos(theta))*mm, -(180.296366454+539.*cos(theta))*mm, -(180.296366454+588.*cos(theta))*mm},
  																							{-196.*mm, -147.*mm, -98.*mm, -49.*mm, 0., 49.*mm, 98*mm, 147.*mm, 196.*mm, -98.*mm, -49.*mm, 0., 49.*mm, 98.*mm, 
  																							(180.296366454-588.*cos(theta))*mm, (180.296366454-539.*cos(theta))*mm, (180.296366454-490.*cos(theta))*mm,
  																							(180.296366454-441.*cos(theta))*mm, (180.296366454-392.*cos(theta))*mm, (180.296366454-343.*cos(theta))*mm, 
  																							(180.296366454-294.*cos(theta))*mm, (180.296366454-245.*cos(theta))*mm, (180.296366454-196.*cos(theta))*mm, 
  																							(180.296366454-147.*cos(theta))*mm, (180.296366454-98.*cos(theta))*mm, (180.296366454-49.*cos(theta))*mm,
  																							180.296366454*mm, 
  																							(180.296366454+49.*cos(theta))*mm, (180.296366454+98.*cos(theta))*mm, (180.296366454+147.*cos(theta))*mm, 
  																							(180.296366454+196.*cos(theta))*mm, (180.296366454+245.*cos(theta))*mm, (180.296366454+294.*cos(theta))*mm, 
  																							(180.296366454+343.*cos(theta))*mm, (180.296366454+392.*cos(theta))*mm, (180.296366454+441.*cos(theta))*mm, 
  																							(180.296366454+490.*cos(theta))*mm, (180.296366454+539.*cos(theta))*mm, (180.296366454+588.*cos(theta))*mm, 
  																							-(180.296366454-588.*cos(theta))*mm, -(180.296366454-539.*cos(theta))*mm, -(180.296366454-490.*cos(theta))*mm,
  																							-(180.296366454-441.*cos(theta))*mm, -(180.296366454-392.*cos(theta))*mm, -(180.296366454-343.*cos(theta))*mm,
  																							-(180.296366454-294.*cos(theta))*mm, -(180.296366454-245.*cos(theta))*mm, -(180.296366454-196.*cos(theta))*mm, 
  																							-(180.296366454-147.*cos(theta))*mm, -(180.296366454-98.*cos(theta))*mm, -(180.296366454-49.*cos(theta))*mm, 
  																							-180.296366454*mm, 
  																							-(180.296366454+49.*cos(theta))*mm, -(180.296366454+98.*cos(theta))*mm, -(180.296366454+147.*cos(theta))*mm,
  																							-(180.296366454+196.*cos(theta))*mm, -(180.296366454+245.*cos(theta))*mm, -(180.296366454+294.*cos(theta))*mm, 
  																							-(180.296366454+343.*cos(theta))*mm, -(180.296366454+392.*cos(theta))*mm, -(180.296366454+441.*cos(theta))*mm,
  																							-(180.296366454+490.*cos(theta))*mm, -(180.296366454+539.*cos(theta))*mm, -(180.296366454+588.*cos(theta))*mm}
  																						};
  	double pulloutYPosition[numChamber][64] = {	{613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 
  																							-613.75*mm, -613.75*mm, -613.75*mm, -613.75*mm, -613.75*mm, 
  																							-(588.*sin(theta))*mm, -(539.*sin(theta))*mm, -(490.*sin(theta))*mm, -(441.*sin(theta))*mm, -(392.*sin(theta))*mm,
  																							-(343.*sin(theta))*mm, -(294.*sin(theta))*mm, -(245.*sin(theta))*mm, -(196.*sin(theta))*mm, -(147.*sin(theta))*mm, 
  																							-(98.*sin(theta))*mm, -(49.*sin(theta))*mm, 
  																							0.,
  																							(49.*sin(theta))*mm, (98.*sin(theta))*mm, (147.*sin(theta))*mm, (196.*sin(theta))*mm, (245.*sin(theta))*mm,
  																							(294.*sin(theta))*mm, (343.*sin(theta))*mm, (392.*sin(theta))*mm, (441.*sin(theta))*mm, (490.*sin(theta))*mm, 
  																							(539.*sin(theta))*mm, (588.*sin(theta))*mm, 
  																							-(588.*sin(theta))*mm, -(539.*sin(theta))*mm, -(490.*sin(theta))*mm, -(441.*sin(theta))*mm, -(392.*sin(theta))*mm,
  																							-(343.*sin(theta))*mm, -(294.*sin(theta))*mm, -(245.*sin(theta))*mm, -(196.*sin(theta))*mm, -(147.*sin(theta))*mm,
  																							-(98.*sin(theta))*mm, -(49.*sin(theta))*mm, 
  																							0., 
 																								(49.*sin(theta))*mm, (98.*sin(theta))*mm, (147.*sin(theta))*mm, (196.*sin(theta))*mm, (245.*sin(theta))*mm, 
 																								(294.*sin(theta))*mm, (343.*sin(theta))*mm, (392.*sin(theta))*mm, (441.*sin(theta))*mm, (490.*sin(theta))*mm, 
 																								(539.*sin(theta))*mm, (588.*sin(theta))*mm}, 
  																							{613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 613.75*mm, 
  																							-613.75*mm, -613.75*mm, -613.75*mm, -613.75*mm, -613.75*mm, 
  																							-(588.*sin(theta))*mm, -(539.*sin(theta))*mm, -(490.*sin(theta))*mm, -(441.*sin(theta))*mm, -(392.*sin(theta))*mm,
  																							-(343.*sin(theta))*mm, -(294.*sin(theta))*mm, -(245.*sin(theta))*mm, -(196.*sin(theta))*mm, -(147.*sin(theta))*mm, 
  																							-(98.*sin(theta))*mm, -(49.*sin(theta))*mm, 
  																							0.,
  																							(49.*sin(theta))*mm, (98.*sin(theta))*mm, (147.*sin(theta))*mm, (196.*sin(theta))*mm, (245.*sin(theta))*mm, 
  																							(294.*sin(theta))*mm, (343.*sin(theta))*mm, (392.*sin(theta))*mm, (441.*sin(theta))*mm, (490.*sin(theta))*mm, 
  																							(539.*sin(theta))*mm, (588.*sin(theta))*mm, 
  																							-(588.*sin(theta))*mm, -(539.*sin(theta))*mm, -(490.*sin(theta))*mm, -(441.*sin(theta))*mm, -(392.*sin(theta))*mm,
  																							-(343.*sin(theta))*mm, -(294.*sin(theta))*mm, -(245.*sin(theta))*mm, -(196.*sin(theta))*mm, -(147.*sin(theta))*mm,
  																							-(98.*sin(theta))*mm, -(49.*sin(theta))*mm, 
  																							0., 
  																							(49.*sin(theta))*mm, (98.*sin(theta))*mm, (147.*sin(theta))*mm, (196.*sin(theta))*mm, (245.*sin(theta))*mm,
  																							(294.*sin(theta))*mm, (343.*sin(theta))*mm, (392.*sin(theta))*mm, (441.*sin(theta))*mm, (490.*sin(theta))*mm,  										
  																							(539.*sin(theta))*mm, (588.*sin(theta))*mm}
  																						};
  	double pulloutZPosition[numChamber] = {-29.815*mm, 8.585*mm} ;									//symmetric
  																	
  	G4RotationMatrix* rotationPullOut[numChamber][64] = { {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8},
  																												{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, rotationPlacement7, 
  																												rotationPlacement7, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, rotationPlacement8, 
  																												rotationPlacement8}
  																											};
  																											
  	// VFats values
  	double vFatX = 46*mm ;
	 	double vFatY = 43*mm ;
	 	double vFatThick = 1.*mm ;
	 	double copperThick = 17.5*um ;
	 	double insulatorThick = 310*um ;
	 	G4double chipX = 7.*mm ;
	 	G4double chipY = 9.*mm ;
	 	G4double chipThick = 350*um ;
	 			
  	double vFatZPos[numChamber] = {-15.834*mm, 22.566*mm} ;							//symmetric			// Adding 38.4 mm in first position
  	
  	// OptoHybrid values
  	G4double optoHybridThickness = 1.6*mm ;
	 	G4double optoHybridWidth = 220.13*mm ;
	 	G4double optoHybridHeight = 140.60*mm ;
  	double optoHybridZPos[numChamber] = {-10.034*mm, 28.366*mm} ;				//symmetric			// Adding 38.4 mm in first position
  	
  	// Cooling copper pads values
  	double coolPadThickness = 1.*mm ;
  	G4double coolBigHeight = 281.493*mm ;
	 	G4double coolBigThinBase = 144.*mm ;
	 	G4double coolBigLargeBase = coolBigThinBase + (2.*coolBigHeight*113.905)/1283.;
	 	G4double coolShortHeight = 51.*mm ;		
  	G4double optoCoolThick = 1.*mm ;
	 	G4double optoCoolWidth = optoHybridWidth + 24.*mm ;
	 	G4double optoCoolHeight = optoHybridHeight + 8.*mm ;
  	G4double underPipePadWidth = 28*mm;
	 	double underPipePadYPos[numUnderPipePad +1] = {-360.007*mm, -255.49*mm, -120.425*mm, 41.1356*mm, 203.112*mm, 354.237*mm, 462.909*mm} ;
  	double coolPadsZPos[numChamber] = {-14.834*mm, 23.566*mm} ;                // Adding 38.4 mm in first position
  	double coolShortPadYPos[numShortPad] = {-255.49*mm, -120.425*mm, 41.1356*mm, 203.112*mm, 354.237*mm, 462.909*mm} ;					
  	double coolShortPadWidth[numShortPad] = {208.012*mm, 231.994*mm, 260.681*mm, 289.442*mm, 316.276*mm, 335.576*mm} ;
  	double optoCoolPadsZPos[numChamber] = {-5.834*mm, 32.566*mm} ;
  	
  	// Cooling pipe valeus
  	G4double inRadius = 6./2.*mm ; // 3. mm
	 	G4double outRadius = 8./2.*mm ; // 4. mm
	 	G4double shift = 90.*mm ;
	 	G4double longPipeLength = 490.*mm;
	 	G4double shortPipeLength = optoCoolHeight/2.;
	 	G4double slantPipeLength = 42.2520922264*mm;
	 	G4double aperture = 40.*mm ;
	 	double largeSideBendXPos[4] = {-108.065*mm, -50.*mm, 108.065*mm, 50.*mm} ;		
	 	double largeSideBendYPos[4] = {492.90*mm, 428.5*mm, 492.90*mm, 428.5*mm} ;		
	 	G4RotationMatrix* rotationLargeSideBend[4] = { rotationPlacement9, 0, rotationPlacement11, rotationPlacement10 };
	 			
	 	double largeSideSlantXPos[2] = {-79.0325000002*mm, 79.0325000002*mm} ;		
	 	double largeSideSlantYPos[2] = {460.7*mm, 460.7*mm} ;		
	 	G4RotationMatrix* rotationLargeSideSlant[2] = { rotationPlacement12, rotationPlacement13 };	
	 	
	 	double straightLongXPos[2] = {aperture, -aperture} ;		
	 	double straightLongYPos[2] = {(-driftHeight/2.+shift+longPipeLength), (-driftHeight/2.+shift+longPipeLength)} ;		
	 	G4RotationMatrix* rotationStraightLong[2] = { rotationPlacement3, rotationPlacement };
	 	
	 	double straightShortXPos[2] = {(optoCoolWidth/2.-outRadius), (-optoCoolWidth/2.+outRadius)} ;		
	 	double straightShortYPos[2] = {(driftHeight/2.-shortPipeLength), (driftHeight/2.-shortPipeLength)} ;		
	 	G4RotationMatrix* rotationStraightShort[2] = { rotationPlacement3, rotationPlacement };				
  	
  	double coolPipeZPos[numChamber] = {-10.334*mm, 28.066*mm} ;
 		
 		// External frame values
  	double frameThickness = 7.25*mm ;
  	double frameZPos[numChamber] = {-29.69*mm, 8.71*mm} ;                //symmetric
  	
  	// External zig values
  	G4double zigThickness = 11.481*mm ;
  	double zigZPos[numChamber] = {-27.5745*mm, 10.8255*mm} ;
  	
  	// Kapton cover values
  	G4double kaptonCoverThickness = 0.050*mm ;
  	double attachedKaptonCoverZPos[numChamber] = {-21.809*mm, 16.591*mm} ;
  	double bulgeKaptonCoverZPos[numChamber] = {-2.875*mm, 35.525*mm} ;
  	double sideKaptonCoverXPos[2] = {137.778211389*mm, -137.778211389*mm} ;
  	
  	//Aluminium cover values
  	G4double coverThickness = 1.*mm ;
  	double attachedCoverZPos[numChamber] = {-21.284*mm, 17.116*mm} ;
  	double bulgeCoverZPos[numChamber] = {-2.35*mm, 36.05*mm} ;
  	double sideCoverXPos[2] = {138.305472996*mm, -138.305472996*mm} ;
  	//double sideCoverXPos[2] = {139.04343323*mm, -139.04343323*mm} ;
  	G4RotationMatrix* rotationSideCover[numChamber][2] = {{rotationPlacement7, rotationPlacement8}, {rotationPlacement7, rotationPlacement8}};
  																						
  	// For loop --- 2 Chambers														 	
  	for(G4int chamber = 0; chamber < buildChamber ; chamber++)
  		{				
				G4Trd* strato[numChamber];
  			G4LogicalVolume* logicStrato[numChamber];
  
  			for(G4int lyr=0; lyr<30; lyr++)
  				{
    				strato[chamber] = Trapezoid(layerNames[lyr], thickness[lyr],layerLargeBase[lyr],layerThinBase[lyr], layerHeight[lyr]) ;
    				logicStrato[chamber] = new G4LogicalVolume (strato[chamber], layerMat[lyr],layerNamesLog[lyr]);   
    				logicStrato[chamber]->SetVisAttributes(visAttributes[lyr]) ;
    				trdCollection[chamber].push_back(strato[chamber]) ;
    				trdLogCollection[chamber].push_back(logicStrato[chamber]) ;
  				} 
  		
				PlaceGeometry(chamber, rotationPlacement,G4ThreeVector(0.,0.,chamberZPos[chamber]),worldLog) ;
				
				for(G4int oFake = 0; oFake < 2; oFake++)
  				{				
						fake[chamber][oFake] = new G4Trd("Fake"+to_string(oFake)+to_string(chamber),  fakeThickness/2., fakeThickness/2., driftThinBase/2., driftLargeBase/2., driftHeight/2.);
  					fakeLog[chamber][oFake] = new G4LogicalVolume(fake[chamber][oFake], fAirMat, "FakeLog");
						fakeLog[chamber][oFake] ->SetVisAttributes(new G4VisAttributes(*fakeAttributes)) ;
			
						new G4PVPlacement ( rotationPlacement,
									  						G4ThreeVector(0., 0., fakeZPos[chamber][oFake]),
									  						fakeLog[chamber][oFake] ,
									  						"Fake"+to_string(oFake)+"Phy"+to_string(chamber),
									  						worldLog,
									 		 					false,
									  						0);	
					}	
				
				G4Trd* outerFake = new G4Trd("OuterFake",  lateralFakeThick/2., lateralFakeThick/2., (driftThinBase + 2.*0.09151530244*mm)/2, (driftLargeBase + 2.*0.10927134296*mm)/2,
		 															(driftHeight + 0.2*mm)/2);
    		G4Trd* innerFake = new G4Trd("InnerFake",  (lateralFakeThick+0.1*mm)/2., (lateralFakeThick+0.1*mm)/2., (driftThinBase)/2, (driftLargeBase)/2, (driftHeight)/2);		
		
				lateralFake[chamber] = new G4SubtractionSolid("LateralFake"+to_string(chamber), outerFake, innerFake);
				lateralFakeLog[chamber] = new G4LogicalVolume(lateralFake[chamber], fAirMat, "LateralFakeLog");
				lateralFakeLog[chamber]->SetVisAttributes(new G4VisAttributes(*fakeAttributes)) ;
						
				new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0., 0., lateralFakeZPos[chamber]),
									  				lateralFakeLog[chamber],
									  				"LateralFakePhy"+to_string(chamber),
									  				worldLog, 
									  				false,
									  				0);																		
  			
  			// other stuff behind the readout starts
  			
  			// Placement of Spacers
				for(G4int spacer = 0; spacer < numSpacer; spacer++)					  				
					{
    				G4Trd* outSpacerTrd = new G4Trd("OutSpacerTrd",  spacerThickness[spacer]/2.0, spacerThickness[spacer]/2.0, (spacerThinBase)/2, (spacerLargeBase)/2,
    																				(spacerHeight)/2);													
						G4Trd* inSpacerTrd = new G4Trd("InSpacerTrd", (spacerThickness[spacer]+0.1)/2.0, (spacerThickness[spacer]+0.1)/2.0, (gapThinBase)/2.0, (gapLargeBase)/2, 
																						(gapHeight)/2);
						
						G4SubtractionSolid* spacerTrd = new G4SubtractionSolid("SpacerTrd", outSpacerTrd, inSpacerTrd);
						G4LogicalVolume* spacerTrdLog = new G4LogicalVolume(spacerTrd, fFR4Mat, "SpacerTrdLog");
						spacerTrdLog->SetVisAttributes(new G4VisAttributes(*spacerAttributes)) ;
						new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0, 0.,spacerPosition[chamber][spacer]),
									  				spacerTrdLog,
									  				"SpacerTrd"+to_string(spacer)+"Phy"+to_string(chamber),
									  				worldLog,
									  				false,
									  				0);
					}
  			
  			// Placement of Pull-Outs Starts
				for(G4int pull = 0; pull < 64; pull++)		
					{		  				
						G4Box* pullOut = new G4Box("PullOut",pullOutX /2., pullOutY/2., pullOutZ/2.) ;
						G4LogicalVolume* pullOutLog = new G4LogicalVolume(pullOut, fPullMat, "PullOutLog");
						pullOutLog->SetVisAttributes(new G4VisAttributes(*pullOutAttributes)) ;
						new G4PVPlacement ( rotationPullOut[chamber][pull],
									  						G4ThreeVector(pulloutXPosition[chamber][pull], pulloutYPosition[chamber][pull],pulloutZPosition[chamber]),
									  						pullOutLog,
									  						"PullOut"+to_string(pull)+"Phy"+to_string(chamber),
									  						worldLog,
									  						false,
									  						0);
					}	
  			
  			// Placement of External Frame
  			G4Box* pOut = new G4Box("POut", 7.255*mm /2., 19.*mm/2., 4.10*mm/2.) ;	
    		G4Trd* outFrame = new G4Trd("OutFrame", (frameThickness)/2.0, (frameThickness)/2.0, (readoutThinBase)/2, (readoutLargeBase)/2, (readoutHeight)/2.);
				G4Trd* inFrame = new G4Trd("InFrame", (frameThickness+0.2)/2.0, (frameThickness+0.2)/2.0, (frameInnerThinBase)/2.0, (frameInnerLargeBase)/2., (frameInnerHeight)/2.);
				
				// Created a Frame		
				G4SubtractionSolid* frame = new G4SubtractionSolid("Frame", outFrame, inFrame);
				
				// Subtracting Pull-Outs sizes from Thin & Large sides  				
				G4SubtractionSolid* f1 = new G4SubtractionSolid("Frame", frame, pOut, 0, G4ThreeVector(0., -196.*mm, 613.75*mm));
				G4SubtractionSolid* f2 = new G4SubtractionSolid("Frame", f1, pOut, 0, G4ThreeVector(0., -147.*mm, 613.75*mm));			  				
				G4SubtractionSolid* f3 = new G4SubtractionSolid("Frame", f2, pOut, 0, G4ThreeVector(0., -98.*mm, 613.75*mm));
				G4SubtractionSolid* f4 = new G4SubtractionSolid("Frame", f3, pOut, 0, G4ThreeVector(0., -49.*mm, 613.75*mm));					  				
				G4SubtractionSolid* f5 = new G4SubtractionSolid("Frame", f4, pOut, 0, G4ThreeVector(0., 0., 613.75*mm));  	
				G4SubtractionSolid* f6 = new G4SubtractionSolid("Frame", f5, pOut, 0, G4ThreeVector(0., 49.*mm, 613.75*mm));
				G4SubtractionSolid* f7 = new G4SubtractionSolid("Frame", f6, pOut, 0, G4ThreeVector(0., 98.*mm, 613.75*mm));
				G4SubtractionSolid* f8 = new G4SubtractionSolid("Frame", f7, pOut, 0, G4ThreeVector(0., 147.*mm, 613.75*mm));
				G4SubtractionSolid* f9 = new G4SubtractionSolid("Frame", f8, pOut, 0, G4ThreeVector(0., 196.*mm, 613.75*mm));					  				
				G4SubtractionSolid* f10 = new G4SubtractionSolid("Frame", f9, pOut, 0, G4ThreeVector(0., -98.*mm, -613.75*mm));
				G4SubtractionSolid* f11 = new G4SubtractionSolid("Frame", f10, pOut, 0, G4ThreeVector(0., -49.*mm, -613.75*mm));					  				
				G4SubtractionSolid* f12 = new G4SubtractionSolid("Frame", f11, pOut, 0, G4ThreeVector(0., 0., -613.75*mm));  	
				G4SubtractionSolid* f13 = new G4SubtractionSolid("Frame", f12, pOut, 0, G4ThreeVector(0., 49.*mm, -613.75*mm));
				G4SubtractionSolid* f14 = new G4SubtractionSolid("Frame", f13, pOut, 0, G4ThreeVector(0., 98.*mm, -613.75*mm));
				
				// Subtracting Pull-Outs sizes from Non parallel side -- Right one from bottom view  
				G4SubtractionSolid* f15 = new G4SubtractionSolid("Frame", f14, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-588.*cos(theta))*mm, -(588.*sin(theta))*mm));
				G4SubtractionSolid* f16 = new G4SubtractionSolid("Frame", f15, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-539.*cos(theta))*mm, -(539.*sin(theta))*mm));
				G4SubtractionSolid* f17 = new G4SubtractionSolid("Frame", f16, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-490.*cos(theta))*mm, -(490.*sin(theta))*mm));
				G4SubtractionSolid* f18 = new G4SubtractionSolid("Frame", f17, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-441.*cos(theta))*mm, -(441.*sin(theta))*mm));
				G4SubtractionSolid* f19 = new G4SubtractionSolid("Frame", f18, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-392.*cos(theta))*mm, -(392.*sin(theta))*mm));
				G4SubtractionSolid* f20 = new G4SubtractionSolid("Frame", f19, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-343.*cos(theta))*mm, -(343.*sin(theta))*mm));
				G4SubtractionSolid* f21 = new G4SubtractionSolid("Frame", f20, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-294.*cos(theta))*mm, -(294.*sin(theta))*mm));
				G4SubtractionSolid* f22 = new G4SubtractionSolid("Frame", f21, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-245.*cos(theta))*mm, -(245.*sin(theta))*mm));
				G4SubtractionSolid* f23 = new G4SubtractionSolid("Frame", f22, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-196.*cos(theta))*mm, -(196.*sin(theta))*mm));
				G4SubtractionSolid* f24 = new G4SubtractionSolid("Frame", f23, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-147.*cos(theta))*mm, -(147.*sin(theta))*mm));
				G4SubtractionSolid* f25 = new G4SubtractionSolid("Frame", f24, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-98.*cos(theta))*mm, -(98.*sin(theta))*mm));
				G4SubtractionSolid* f26 = new G4SubtractionSolid("Frame", f25, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454-49.*cos(theta))*mm, -(49.*sin(theta))*mm));
				G4SubtractionSolid* f27 = new G4SubtractionSolid("Frame", f26, pOut, rotationPlacement5, G4ThreeVector(0., 180.296366454*mm, 0.));			
				G4SubtractionSolid* f28 = new G4SubtractionSolid("Frame", f27, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+49.*cos(theta))*mm, (49.*sin(theta))*mm));
				G4SubtractionSolid* f29 = new G4SubtractionSolid("Frame", f28, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+98.*cos(theta))*mm, (98.*sin(theta))*mm));
				G4SubtractionSolid* f30 = new G4SubtractionSolid("Frame", f29, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+147.*cos(theta))*mm, (147.*sin(theta))*mm));
				G4SubtractionSolid* f31 = new G4SubtractionSolid("Frame", f30, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+196.*cos(theta))*mm, (196.*sin(theta))*mm));
				G4SubtractionSolid* f32 = new G4SubtractionSolid("Frame", f31, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+245.*cos(theta))*mm, (245.*sin(theta))*mm));
				G4SubtractionSolid* f33 = new G4SubtractionSolid("Frame", f32, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+294.*cos(theta))*mm, (294.*sin(theta))*mm));
				G4SubtractionSolid* f34 = new G4SubtractionSolid("Frame", f33, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+343.*cos(theta))*mm, (343.*sin(theta))*mm));
				G4SubtractionSolid* f35 = new G4SubtractionSolid("Frame", f34, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+392.*cos(theta))*mm, (392.*sin(theta))*mm));
				G4SubtractionSolid* f36 = new G4SubtractionSolid("Frame", f35, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+441.*cos(theta))*mm, (441.*sin(theta))*mm));
				G4SubtractionSolid* f37 = new G4SubtractionSolid("Frame", f36, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+490.*cos(theta))*mm, (490.*sin(theta))*mm));
				G4SubtractionSolid* f38 = new G4SubtractionSolid("Frame", f37, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+539.*cos(theta))*mm, (539.*sin(theta))*mm));
				G4SubtractionSolid* f39 = new G4SubtractionSolid("Frame", f38, pOut, rotationPlacement5, G4ThreeVector(0., (180.296366454+588.*cos(theta))*mm, (588.*sin(theta))*mm));
				
				// Subtracting Pull-Outs sizes from Non parallel side -- Left one from bottom view  
				G4SubtractionSolid* f40 = new G4SubtractionSolid("Frame", f39, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-588.*cos(theta))*mm, -(588.*sin(theta))*mm));
				G4SubtractionSolid* f41 = new G4SubtractionSolid("Frame", f40, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-539.*cos(theta))*mm, -(539.*sin(theta))*mm));
				G4SubtractionSolid* f42 = new G4SubtractionSolid("Frame", f41, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-490.*cos(theta))*mm, -(490.*sin(theta))*mm));
				G4SubtractionSolid* f43 = new G4SubtractionSolid("Frame", f42, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-441.*cos(theta))*mm, -(441.*sin(theta))*mm));
				G4SubtractionSolid* f44 = new G4SubtractionSolid("Frame", f43, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-392.*cos(theta))*mm, -(392.*sin(theta))*mm));
				G4SubtractionSolid* f45 = new G4SubtractionSolid("Frame", f44, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-343.*cos(theta))*mm, -(343.*sin(theta))*mm));
				G4SubtractionSolid* f46 = new G4SubtractionSolid("Frame", f45, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-294.*cos(theta))*mm, -(294.*sin(theta))*mm));
				G4SubtractionSolid* f47 = new G4SubtractionSolid("Frame", f46, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-245.*cos(theta))*mm, -(245.*sin(theta))*mm));
				G4SubtractionSolid* f48 = new G4SubtractionSolid("Frame", f47, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-196.*cos(theta))*mm, -(196.*sin(theta))*mm));
				G4SubtractionSolid* f49 = new G4SubtractionSolid("Frame", f48, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-147.*cos(theta))*mm, -(147.*sin(theta))*mm));
				G4SubtractionSolid* f50 = new G4SubtractionSolid("Frame", f49, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-98.*cos(theta))*mm, -(98.*sin(theta))*mm));
				G4SubtractionSolid* f51 = new G4SubtractionSolid("Frame", f50, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454-49.*cos(theta))*mm, -(49.*sin(theta))*mm));
				G4SubtractionSolid* f52 = new G4SubtractionSolid("Frame", f51, pOut, rotationPlacement6, G4ThreeVector(0., -180.296366454*mm, 0.));
				G4SubtractionSolid* f53 = new G4SubtractionSolid("Frame", f52, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+49.*cos(theta))*mm, (49.*sin(theta))*mm));
				G4SubtractionSolid* f54 = new G4SubtractionSolid("Frame", f53, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+98.*cos(theta))*mm, (98.*sin(theta))*mm));
				G4SubtractionSolid* f55 = new G4SubtractionSolid("Frame", f54, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+147.*cos(theta))*mm, (147.*sin(theta))*mm));
				G4SubtractionSolid* f56 = new G4SubtractionSolid("Frame", f55, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+196.*cos(theta))*mm, (196.*sin(theta))*mm));
				G4SubtractionSolid* f57 = new G4SubtractionSolid("Frame", f56, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+245.*cos(theta))*mm, (245.*sin(theta))*mm));
				G4SubtractionSolid* f58 = new G4SubtractionSolid("Frame", f57, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+294.*cos(theta))*mm, (294.*sin(theta))*mm));
				G4SubtractionSolid* f59 = new G4SubtractionSolid("Frame", f58, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+343.*cos(theta))*mm, (343.*sin(theta))*mm));
				G4SubtractionSolid* f60 = new G4SubtractionSolid("Frame", f59, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+392.*cos(theta))*mm, (392.*sin(theta))*mm));
				G4SubtractionSolid* f61 = new G4SubtractionSolid("Frame", f60, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+441.*cos(theta))*mm, (441.*sin(theta))*mm));
				G4SubtractionSolid* f62 = new G4SubtractionSolid("Frame", f61, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+490.*cos(theta))*mm, (490.*sin(theta))*mm));
				G4SubtractionSolid* f63 = new G4SubtractionSolid("Frame", f62, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+539.*cos(theta))*mm, (539.*sin(theta))*mm));
				G4SubtractionSolid* externalFrame = new G4SubtractionSolid("ExternalFrame", f63, pOut, rotationPlacement6, G4ThreeVector(0., -(180.296366454+588.*cos(theta))*mm,
																																		(588.*sin(theta))*mm));
				
				G4LogicalVolume* externalFrameLog = new G4LogicalVolume(externalFrame, fFR4Mat, "ExternalFrameLogical");	
				externalFrameLog->SetVisAttributes(new G4VisAttributes(*frameAttributes)) ;
				
				new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0., 0., frameZPos[chamber]),
									  				externalFrameLog,
									  				"ExternalFramePhy"+to_string(chamber),
									  				worldLog,
									  				false,
									  				0);
  			
  			//vfats starts
   			G4String coordX ;
   			G4String coordY ;
   			std::vector<G4double> x ;
   			std::vector<G4double> y ;
   			std::ifstream coordFile ("GEMDPG/TrGEMG4/data/VFatCoord.txt", std::ios::in) ;
   			if(coordFile) 
   				{
      			while(!coordFile.eof()) 
      				{
	 							coordFile >> coordX ; 
	 							coordFile >> coordY ;
	 							if(coordX == "" || coordY == "") continue ;
	 							x.push_back(std::strtod(coordX, NULL)) ; 
	 							y.push_back(std::strtod(coordY, NULL)) ; 
      				}
      			coordFile.close() ;
   				}
   	
   			// Placement of VFats -- Here lies the positionning of every single VFat2 module
	 			for(G4int vfat = 0; vfat < 24; vfat++) 
	 				{
	 					G4Box* vFatModule = new G4Box("VFatModule", vFatThick/2., vFatX/2., vFatY/2.) ; // G4Box wants half sizes
	 					G4LogicalVolume* vFatModuleLog = new G4LogicalVolume(vFatModule, fAirMat, "VFatModuleLog") ;
	 					vFatModuleLog->SetVisAttributes(new G4VisAttributes(*vFatAttributes)) ;
	 			
	 					G4Box* copperVFat = new G4Box("CopperVFat", copperThick/2., vFatX/2., vFatY/2.) ;
	    			G4LogicalVolume* copperVFatLog = new G4LogicalVolume(copperVFat, fCuMat, "CopperVFatLog") ;
	    			
	    			G4Box* insulatorVFat = new G4Box("InsulatorVFat", insulatorThick/2., vFatX/2., vFatY/2.) ;
	    			G4LogicalVolume* insulatorVFatLog = new G4LogicalVolume(insulatorVFat, fFR4Mat, "InsulatorVFatLog") ;
	    			
	    			G4Box* chipVFat = new G4Box("ChipVFat", chipThick/2., chipX/2., chipY/2.) ;
	    			G4LogicalVolume* chipVFatLog = new G4LogicalVolume(chipVFat, fChipMat, "ChipVFatLog") ;
	 					
	    			// VFat vonfiguration
	    			G4int copyNo = 0. ;
	    			
	    			new G4PVPlacement(0, G4ThreeVector(vFatThick/2.+chipThick/2., 0., 0.), chipVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), vFatModuleLog, false, 
	    												copyNo++);
	    			new G4PVPlacement(0, G4ThreeVector(vFatThick/2.-copperThick/2., 0., 0.), copperVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), vFatModuleLog, false, 
	    												copyNo++) ; 
	   	 			new G4PVPlacement(0, G4ThreeVector(vFatThick/2.-copperThick-insulatorThick/2., 0., 0.), insulatorVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), 
	   	 												vFatModuleLog, false, copyNo++) ; 
	    			new G4PVPlacement(0, G4ThreeVector(insulatorThick/2.+copperThick/2., 0., 0.), copperVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), vFatModuleLog, 
	    												false, copyNo++) ; 
	    			new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), insulatorVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), vFatModuleLog, false, copyNo++) ; 
	    			new G4PVPlacement(0, G4ThreeVector(-insulatorThick/2.-copperThick/2., 0., 0.), copperVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), vFatModuleLog, 
	    												false, copyNo++) ; 
	    			new G4PVPlacement(0, G4ThreeVector(-vFatThick/2.+copperThick+insulatorThick/2., 0., 0.), insulatorVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), 
	    												vFatModuleLog, false, copyNo++) ; 
	    			new G4PVPlacement(0, G4ThreeVector(-vFatThick/2.+copperThick/2., 0., 0.), copperVFatLog, "VFatModule"+to_string(vfat)+"Phy"+to_string(chamber), vFatModuleLog, false, 
	    												copyNo++) ; 
	    												
	    			// Placement of a unit of VFat
	    			new G4PVPlacement	( rotationPlacement4, 
	    													G4ThreeVector( y[vfat], -x[vfat] + driftHeight/2.0, vFatZPos[chamber]), 
	    													vFatModuleLog, 
	    													"VFat"+to_string(vfat)+"Phy"+to_string(chamber), 
	    													worldLog, 
	    													false, 
	    													vfat) ;
	    			G4cout<<"VFat Positions.. = "<<G4ThreeVector( y[vfat], -x[vfat] + driftHeight/2.0, vFatZPos[chamber])<<G4endl;
	 				}
	 			
	 			// Placement of Optohybrid
	 			G4Box* optoHybrid = new G4Box("OptoHybrid", optoHybridThickness/2., optoHybridWidth/2.*mm, optoHybridHeight/2.*mm) ;
	 			G4LogicalVolume* optoHybridLog = new G4LogicalVolume(optoHybrid, fFR4Mat, "OptoHybridLog") ;
	 			optoHybridLog->SetVisAttributes(new G4VisAttributes(*vFatAttributes)) ;
	 			new G4PVPlacement	( rotationPlacement4, 
	 												G4ThreeVector(0., -4.*mm+(-optoHybridHeight+driftHeight)/2., optoHybridZPos[chamber]), 
	 												optoHybridLog, 
	 												"OptoHybridPhy"+to_string(chamber), 
	 												worldLog, 
	 												false, 
	 												0) ;
				
				// Placement of Optohybrid Cooling Pad
	 			G4Box* optoCoolPad = new G4Box("OptoCoolPad", optoCoolThick/2., optoCoolWidth/2., optoCoolHeight/2.) ; 
	 			G4LogicalVolume* optoCoolPadLog = new G4LogicalVolume(optoCoolPad, fCuMat, "OptoCoolPadLog") ;
    		optoCoolPadLog->SetVisAttributes(new G4VisAttributes(*coolingPadAttributes)) ;
	 					
	 			new G4PVPlacement ( rotationPlacement4, 
	 											 		G4ThreeVector(0., -4.*mm+(-optoHybridHeight+driftHeight)/2., optoCoolPadsZPos[chamber]), 
	 											 		optoCoolPadLog, 
	 											 		"OptoCoolPadPhys"+to_string(chamber), 
	 											 		worldLog, 
	 											 		false, 
	 											 		0) ;
	 				
	 			// Placement of Cooling Pads
	 			// Bigger Pad
	 			G4Trd* bigCoolPad = new G4Trd("BigCoolPad", coolPadThickness/2., coolPadThickness/2., coolBigThinBase/2., coolBigLargeBase/2., coolBigHeight/2.) ; 
	 			G4LogicalVolume* bigCoolPadLog = new G4LogicalVolume(bigCoolPad, fCuMat, "BigCoolPadLog") ;
    		bigCoolPadLog->SetVisAttributes(new G4VisAttributes(*coolingPadAttributes)) ;
	 					
	 			new G4PVPlacement ( rotationPlacement4, 
	 											 		G4ThreeVector(0., (-driftHeight + coolBigHeight)/2., coolPadsZPos[chamber]), 
	 											 		bigCoolPadLog, 
	 											 		"BigCoolPadPhys"+to_string(chamber), 
	 											 		worldLog, 
	 											 		false, 
	 											 		0) ;
	 			
	 			// Smaller Pads	
	 			G4Box* smallCoolPad[numShortPad] ;
	 			G4LogicalVolume* smallCoolPadLog[numShortPad] ;		
	 			
	 			for(G4int sPad = 0; sPad < 6; sPad++) 
	 				{	
	 					smallCoolPad[sPad] = new G4Box("SmallCoolPad", coolPadThickness/2., coolShortPadWidth[sPad]/2., coolShortHeight/2.) ; 
	 					smallCoolPadLog[sPad] = new G4LogicalVolume(smallCoolPad[sPad], fCuMat, "SmallCoolPadLog") ;
    				smallCoolPadLog[sPad]->SetVisAttributes(new G4VisAttributes(*coolingPadAttributes)) ;
    			
	 					new G4PVPlacement ( rotationPlacement4, 
	 															G4ThreeVector(0, coolShortPadYPos[sPad], coolPadsZPos[chamber]), 
	 															smallCoolPadLog[sPad], 
	 															"SmallCoolPad"+to_string(sPad)+"Phys"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 				}
	 			
	 			// Under Pipe Pads
	 			G4Box* underPipePad[numUnderPipePad] ;
	 			G4LogicalVolume* underPipePadLog[numUnderPipePad] ;	
	 			
	 			for(G4int uPad = 0; uPad < 6; uPad++) 
	 				{	
	 					underPipePad[uPad] = new G4Box("UnderPipeCoolPad", underPipePadWidth/2.,  (underPipePadYPos[uPad+1]-51.*mm/2. - underPipePadYPos[uPad])/2., coolPadThickness/2.) ; 
	 					underPipePadLog[uPad] = new G4LogicalVolume(underPipePad[uPad], fCuMat, "UnderPipeCoolPadLog") ;
    				underPipePadLog[uPad]->SetVisAttributes(new G4VisAttributes(*coolingPadAttributes)) ;
	 	
	 					new G4PVPlacement ( 0, 
	 															G4ThreeVector(aperture, (underPipePadYPos[uPad+1]-51.*mm/2. + underPipePadYPos[uPad])/2, coolPadsZPos[chamber]), 
	 															underPipePadLog[uPad], 
	 															"UnderPipeCoolPadR"+to_string(uPad)+"Phys"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 															
	 					new G4PVPlacement ( 0, 
	 															G4ThreeVector(-aperture, (underPipePadYPos[uPad+1]-51.*mm/2. + underPipePadYPos[uPad])/2, coolPadsZPos[chamber]), 
	 															underPipePadLog[uPad], 
	 															"UnderPipeCoolPadL"+to_string(uPad)+"Phys"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 				}
	 													
	 			// Placement of Cooling Pipe
	 			G4Torus* thinSideBendPipe = new G4Torus("ThinSideBendPipe", inRadius, outRadius, aperture, 0.*degree, 180.*degree) ; 
	 			G4LogicalVolume* thinSideBendPipeLog = new G4LogicalVolume(thinSideBendPipe, fCuMat, "ThinSideBendPipeLog") ; 
    		thinSideBendPipeLog->SetVisAttributes(new G4VisAttributes(*coolingAttributes)) ;
         		
	 			G4Torus* thinSideBendWater = new G4Torus("ThinSideBendWater", 0., inRadius, aperture, 0.*degree, 180.*degree) ; 
	 			G4LogicalVolume* thinSideBendWaterLog = new G4LogicalVolume(thinSideBendWater, fCoolMat, "ThinSideBendWaterLog") ; 
	 			
	 			G4Torus* largeSideBendPipe = new G4Torus("LargeSideBendPipe", inRadius, outRadius, 10.*mm, 0.*degree, 55.*degree) ; 
	 			G4LogicalVolume* largeSideBendPipeLog = new G4LogicalVolume(largeSideBendPipe, fCuMat, "LargeSideBendPipeLog") ; 
    		largeSideBendPipeLog->SetVisAttributes(new G4VisAttributes(*coolingAttributes)) ;
         		
	 			G4Torus* largeSideBendWater = new G4Torus("LargeSideBendWater", 0., inRadius, 10*mm, 0.*degree, 55.*degree) ; 
	 			G4LogicalVolume* largeSideBendWaterLog = new G4LogicalVolume(largeSideBendWater, fCoolMat, "LargeSideBendWaterLog") ; 
	 					
	 			G4Tubs* straightLongPipe = new G4Tubs("StraightLongPipe", inRadius, outRadius, longPipeLength, 0., 2*M_PI) ;	
	 			G4LogicalVolume* straightLongPipeLog = new G4LogicalVolume(straightLongPipe, fCuMat, "straightPipeLog") ;
    		straightLongPipeLog->SetVisAttributes(new G4VisAttributes(*coolingAttributes)) ;
         		
    		G4Tubs* straightLongWater = new G4Tubs("StraightLongWater", 0., inRadius, longPipeLength, 0., 2*M_PI) ;
	 			G4LogicalVolume* straightLongWaterLog = new G4LogicalVolume(straightLongWater, fCoolMat, "StraightWaterLog") ;
	 			
	 			G4Tubs* straightShortPipe = new G4Tubs("StraightShortPipe", inRadius, outRadius, shortPipeLength, 0., 2*M_PI) ;	
	 			G4LogicalVolume* straightShortPipeLog = new G4LogicalVolume(straightShortPipe, fCuMat, "StraightShortPipeLog") ;
    		straightShortPipeLog->SetVisAttributes(new G4VisAttributes(*coolingAttributes)) ;
         		
    		G4Tubs* straightShortWater = new G4Tubs("StraightShortWater", 0., inRadius, shortPipeLength, 0., 2*M_PI) ;
	 			G4LogicalVolume* straightShortWaterLog = new G4LogicalVolume(straightShortWater, fCoolMat, "StraightShortWaterLog") ;
	 			
	 			G4Tubs* slantShortPipe = new G4Tubs("SlantShortPipe", inRadius, outRadius, slantPipeLength, 0., 2*M_PI) ;	
	 			G4LogicalVolume* slantShortPipeLog = new G4LogicalVolume(slantShortPipe, fCuMat, "SlantShortPipeLog") ;
    		slantShortPipeLog->SetVisAttributes(new G4VisAttributes(*coolingAttributes)) ;
         		
    		G4Tubs* slantShortWater = new G4Tubs("SlantShortWater", 0., inRadius, slantPipeLength, 0., 2*M_PI) ;
	 			G4LogicalVolume* slantShortWaterLog = new G4LogicalVolume(slantShortWater, fCoolMat, "SlantShortWaterLog") ;
	 			
	 			new G4PVPlacement ( rotationPlacement2, 
	 													G4ThreeVector(0.,-driftHeight/2.+shift, coolPipeZPos[chamber]), 
	 													thinSideBendPipeLog, 
	 													"ThinSideBendPipePhy"+to_string(chamber), 
	 													worldLog, 
	 													false, 
	 													0) ;
	 			new G4PVPlacement ( rotationPlacement2, 
	 													G4ThreeVector(0.,-driftHeight/2.+shift, coolPipeZPos[chamber]), 
	 													thinSideBendWaterLog, 
	 													"ThinSideBendWaterPhy"+to_string(chamber),
	 													worldLog, 
	 													false, 
	 													0) ;
	 			
	 			// Bends on lrge base side
	 			for(G4int largeSideBend = 0; largeSideBend < 4; largeSideBend++) 
	 				{
	 					new G4PVPlacement ( rotationLargeSideBend[largeSideBend], 
	 															G4ThreeVector(largeSideBendXPos[largeSideBend], largeSideBendYPos[largeSideBend], coolPipeZPos[chamber]), 
	 															largeSideBendPipeLog, 
	 															"LargeSideBendPipe"+to_string(largeSideBend)+"Phy"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 					new G4PVPlacement ( rotationLargeSideBend[largeSideBend], 
	 															G4ThreeVector(largeSideBendXPos[largeSideBend], largeSideBendYPos[largeSideBend], coolPipeZPos[chamber]),
	 															largeSideBendWaterLog, 
	 															"LargeSideBendWater"+to_string(largeSideBend)+"Phy"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;	
	 				}						
	 				
	 			for(G4int largeSideSlant = 0; largeSideSlant < 2; largeSideSlant++) 
	 				{
	 					new G4PVPlacement ( rotationLargeSideSlant[largeSideSlant], 
	 															G4ThreeVector(largeSideSlantXPos[largeSideSlant], largeSideSlantYPos[largeSideSlant], coolPipeZPos[chamber]), 
	 															slantShortPipeLog, 
	 															"SlantShortPipe"+to_string(largeSideSlant)+"Phy"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 					new G4PVPlacement ( rotationLargeSideSlant[largeSideSlant], 
	 															G4ThreeVector(largeSideSlantXPos[largeSideSlant], largeSideSlantYPos[largeSideSlant], coolPipeZPos[chamber]), 
	 															slantShortWaterLog, 
	 															"SlantShortWater"+to_string(largeSideSlant)+"Phy"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 				}
	 				
	 			for(G4int straightLong = 0; straightLong < 2; straightLong++) 
	 				{
	 					new G4PVPlacement ( rotationStraightLong[straightLong], 
	 															G4ThreeVector(straightLongXPos[straightLong], straightLongYPos[straightLong], coolPipeZPos[chamber]),
	 															straightLongPipeLog, 
	 															"StraightLongPipe"+to_string(straightLong)+"Phy"+to_string(chamber),
	 															worldLog, 
	 															false, 
	 															0) ;
	 					new G4PVPlacement ( rotationStraightLong[straightLong], 
	 															G4ThreeVector(straightLongXPos[straightLong], straightLongYPos[straightLong], coolPipeZPos[chamber]),
	 															straightLongWaterLog, 
	 															"StraightLongWater"+to_string(straightLong)+"Phy"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 				}
	 							
	 			for(G4int straightShort = 0; straightShort < 2; straightShort++) 
	 				{
	 					new G4PVPlacement ( rotationStraightShort[straightShort], 
	 															G4ThreeVector(straightShortXPos[straightShort], straightShortYPos[straightShort], coolPipeZPos[chamber]),
	 															straightShortPipeLog, 
	 															"StraightShortPipe"+to_string(straightShort)+"Phy"+to_string(chamber),
	 															worldLog, 
	 															false, 
	 															0) ;
	 					new G4PVPlacement ( rotationStraightShort[straightShort], 
	 															G4ThreeVector(straightShortXPos[straightShort], straightShortYPos[straightShort], coolPipeZPos[chamber]),
	 															straightShortWaterLog, 
	 															"StraightShortWater"+to_string(straightShort)+"Phy"+to_string(chamber), 
	 															worldLog, 
	 															false, 
	 															0) ;
	 				}								
	 											
	 				
  			// Placement of Aluminium Zig Starts
  			G4Trd* outZig = new G4Trd("OutZig",  (zigThickness)/2.0, (zigThickness)/2.0, (driftThinBase)/2, (driftLargeBase)/2, (driftHeight)/2);													
				G4Trd* inZig = new G4Trd("InCover", (zigThickness+0.1*mm)/2.0, (zigThickness+0.1*mm)/2.0, (zigThinBase)/2.0, (zigLargeBase)/2, (zigHeight)/2);
  			
  			G4SubtractionSolid* zigFrame = new G4SubtractionSolid("ZigFrame", outZig, inZig);
				G4LogicalVolume* zigFrameLog = new G4LogicalVolume(zigFrame, fCoverMat, "ZigFrameLog");
				zigFrameLog->SetVisAttributes(new G4VisAttributes(*zigAttributes)) ;
				new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0, 0., zigZPos[chamber]),
									  				zigFrameLog ,
									  				"ZigFramePhy"+to_string(chamber),
									  				worldLog,
									  				false,
									  				0);
									  				
				// Kapton Cover Starts
  			G4Trd* outKaptonCover = new G4Trd("OutKaptonCover",  (kaptonCoverThickness)/2.0, (kaptonCoverThickness)/2.0, (driftThinBase)/2, (driftLargeBase)/2, (driftHeight)/2);			
				G4Trd* inKaptonCover = new G4Trd("InKaptonCover", (kaptonCoverThickness+0.1*mm)/2.0, (kaptonCoverThickness+0.1*mm)/2.0, (driftThinBase-120.572380569*mm)/2.0, 	
																					(driftLargeBase-120.572380569*mm)/2, (driftHeight+1.*mm)/2);
  			
  			G4SubtractionSolid* attachedKaptonCover = new G4SubtractionSolid("AttachedKaptonCover", outKaptonCover, inKaptonCover);
				G4LogicalVolume* attachedKaptonCoverLog = new G4LogicalVolume(attachedKaptonCover, fKaptonMat, "AttachedKaptonCoverLog");
				attachedKaptonCoverLog->SetVisAttributes(new G4VisAttributes(*gemAttributes)) ;
				new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0., 0., attachedKaptonCoverZPos[chamber]),
									  				attachedKaptonCoverLog,
									  				"AttachedKaptonCoverPhy"+to_string(chamber),
									  				worldLog,
									  				false,
									  				0);
				// Bulge Kapton Cover				  		
				G4Trd* bulgeKaptonCover = new G4Trd("BulgeKaptonCover", (kaptonCoverThickness)/2.0, (kaptonCoverThickness)/2.0, (driftThinBase-120.472380569*mm)/2.0, 	
																						(driftLargeBase-120.472380569*mm)/2, (driftHeight)/2);
				G4LogicalVolume* bulgeKaptonCoverLog = new G4LogicalVolume(bulgeKaptonCover, fKaptonMat, "BulgeKaptonCoverLog");
				bulgeKaptonCoverLog->SetVisAttributes(new G4VisAttributes(*gemAttributes)) ;
				new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0., 0., bulgeKaptonCoverZPos[chamber]),
									  				bulgeKaptonCoverLog,
									  				"BulgeKaptonCoverPhy"+to_string(chamber),
									  				worldLog,
									  				false,
									  				0);
				// Side Kapton Cover	
				for(G4int sideK = 0; sideK < 2; sideK++) 
	 				{				  				
						//G4Para* sideKaptonCover = new G4Para("SideKaptonCover", 1283.*mm/2., 1.00393322705*mm/2., (bulgeKaptonCoverZPos[chamber]-attachedKaptonCoverZPos[chamber])/2., theta, 0.,0.) ;
						G4Box* sideKaptonCover = new G4Box("SideKaptonCover", 1288.04633031*mm/2.,  0.050*mm /2., (bulgeKaptonCoverZPos[chamber]-attachedKaptonCoverZPos[chamber])/2.) ;
						G4LogicalVolume* sideKaptonCoverLog = new G4LogicalVolume(sideKaptonCover, fKaptonMat, "SideKaptonCoverLog");
						sideKaptonCoverLog->SetVisAttributes(new G4VisAttributes(*gemAttributes)) ;
						new G4PVPlacement ( rotationSideCover[chamber][sideK],
									  						G4ThreeVector(sideKaptonCoverXPos[sideK], 0., (attachedKaptonCoverZPos[chamber]+bulgeKaptonCoverZPos[chamber])/2.),
									  						sideKaptonCoverLog,
									  						"SideKaptonCover"+to_string(sideK)+"Phy"+to_string(chamber),
									  						worldLog,
									  						false,
									  						0);	
					}	
  			
  			// Aluminium Cover Starts
  			G4Trd* outCover = new G4Trd("OutCover",  (coverThickness)/2.0, (coverThickness)/2.0, (driftThinBase)/2, (driftLargeBase)/2, (driftHeight)/2);													
				G4Trd* inCover = new G4Trd("InCover", (coverThickness+0.1*mm)/2.0, (coverThickness+0.1*mm)/2.0, (driftThinBase-120.4719872476*mm)/2.0, 
																		(driftLargeBase-120.4719872476*mm)/2, (driftHeight+1.*mm)/2);
  			
  			G4SubtractionSolid* attachedCover = new G4SubtractionSolid("AttachedCover", outCover, inCover);
				G4LogicalVolume* attachedCoverLog = new G4LogicalVolume(attachedCover, fCoverMat, "AttachedCoverLog");
				attachedCoverLog->SetVisAttributes(new G4VisAttributes(*coverAttributes)) ;
				new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0., 0., attachedCoverZPos[chamber]),
									  				attachedCoverLog,
									  				"AttachedCoverPhy"+to_string(chamber),
									  				worldLog,
									  				false,
									  				0);
				// Bulge Cover				  		
				G4Trd* bulgeCover = new G4Trd("BulgeCover", (coverThickness)/2.0, (coverThickness)/2.0, (driftThinBase-118.464120792*mm)/2.0, (driftLargeBase-118.464120792*mm)/2, 	
																			(driftHeight)/2);											
  			G4LogicalVolume* bulgeCoverLog = new G4LogicalVolume(bulgeCover, fCoverMat, "BulgeCoverLog");
				//bulgeCoverLog->SetVisAttributes(new G4VisAttributes(*coverAttributes)) ;
				new G4PVPlacement ( rotationPlacement,
									  				G4ThreeVector(0., 0., bulgeCoverZPos[chamber]),
									  				bulgeCoverLog,
									  				"BulgeCoverPhy"+to_string(chamber),
									  				worldLog,
									  				false,
									  				0);
				// Side Cover	
				for(G4int side = 0; side < 2; side++) 
	 				{	
	 					//G4Para* sideCover = new G4Para("SideCover", 1283.*mm/2., 1.00393322705*mm/2., (bulgeCoverZPos[chamber]-attachedCoverZPos[chamber])/2., theta, 0.,0.) ;		  				
						G4Box* sideCover = new G4Box("SideCover", 1288.04633031*mm/2.,  1.*mm /2., (bulgeCoverZPos[chamber]-attachedCoverZPos[chamber])/2.) ;
						G4LogicalVolume* sideCoverLog = new G4LogicalVolume(sideCover, fCoverMat, "SideCoverLog");
						sideCoverLog->SetVisAttributes(new G4VisAttributes(*coverAttributes)) ;
						new G4PVPlacement ( rotationSideCover[chamber][side],
									  						G4ThreeVector(sideCoverXPos[side], 0., (attachedCoverZPos[chamber]+bulgeCoverZPos[chamber])/2.),
									  						sideCoverLog,
									  						"SideCover"+to_string(side)+"Phy"+to_string(chamber),
									  						worldLog,
									  						false,
									  						0);	
					}	
			}
  	return worldPhys ;
	}

G4Trd* TrGEMDetectorConstruction::Trapezoid(G4String name, G4double width, G4double largeBase, G4double thinBase, G4double height) 
	{
  	G4Trd* shape = new G4Trd( name,
                           	 	width/2, width/2,
                           		thinBase/2,
                           		largeBase/2,
                           		height/2) ;
  	return shape ;
	}

void TrGEMDetectorConstruction::PlaceGeometry(G4int chamberNum, G4RotationMatrix *pRot, G4ThreeVector tlate, G4LogicalVolume* pMotherLogical) 
	{
  	G4double XTranslation = 0. ;
   	
   	for(size_t i=0 ; i<trdCollection[chamberNum].size() ; i++) 
   		{
      	// i counts as the copyNo
      	G4String layerName = trdCollection[chamberNum].at(i)->GetName() ;
      	XTranslation += trdCollection[chamberNum].at(i)->GetXHalfLength1() ;
      	G4ThreeVector position = tlate + G4ThreeVector(XTranslation,0,0).transform(G4RotationMatrix(*pRot).inverse()) ;
      	G4cout <<"For Chamber "<<chamberNum<< ".......Volume (" << i << ") " << layerName << " the position is " << G4BestUnit(XTranslation,"Length") << G4endl ;
  
				new G4PVPlacement(pRot, position, trdLogCollection[chamberNum].at(i), layerName+"Phy"+to_string(chamberNum), pMotherLogical, false, i) ;

   			XTranslation += trdCollection[chamberNum].at(i)->GetXHalfLength1() ;
			}
	}

void TrGEMDetectorConstruction::ConstructSDandField()
	{
		G4SDManager* sdman = G4SDManager::GetSDMpointer() ;
		
		TrGEMSensitiveDetector* sensitive = new TrGEMSensitiveDetector("/GasGap") ;
  	sdman->AddNewDetector(sensitive) ;
  	
  	for(G4int chamber = 0; chamber < buildChamber ; chamber++)
  		{
  			fakeLog[chamber][0]->SetSensitiveDetector(sensitive) ;
  			fakeLog[chamber][1]->SetSensitiveDetector(sensitive) ;
  			lateralFakeLog[chamber]->SetSensitiveDetector(sensitive) ;
				trdLogCollection[chamber][3]->SetSensitiveDetector(sensitive) ;									//GasGap1Log	
  			trdLogCollection[chamber][7]->SetSensitiveDetector(sensitive) ;									//GasGap2Log	
  			trdLogCollection[chamber][11]->SetSensitiveDetector(sensitive) ;								//GasGap3Log
  			trdLogCollection[chamber][15]->SetSensitiveDetector(sensitive) ;								//GasGap4Log		
 		}
	}
