// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch
//

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>



class G4Event;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

};


#endif
