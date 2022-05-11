
#ifndef TrGEMActionInitialization_h
#define TrGEMActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class TrGEMActionInitialization : public G4VUserActionInitialization
	{
  	public:
    	TrGEMActionInitialization(char*);
    	virtual ~TrGEMActionInitialization();

    	virtual void Build() const;
    	virtual void BuildForMaster() const;
  	private:
    	char* partName;
    };
    
#endif 
