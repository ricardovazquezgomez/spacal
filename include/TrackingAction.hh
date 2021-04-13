// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch
#ifndef TrackingAction_H
#define TrackingAction_H 1

#include "G4UserTrackingAction.hh"
#include "TrackInformation.hh"



class TrackingAction : public G4UserTrackingAction
{
public:
  TrackingAction();
  ~TrackingAction();

public:
  void PreUserTrackingAction(const G4Track* aTrack);
  void PostUserTrackingAction(const G4Track* aTrack);
};

#endif
