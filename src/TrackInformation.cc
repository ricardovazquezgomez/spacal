#include "TrackInformation.hh"

using namespace CLHEP;



G4Allocator<TrackInformation> aTrackInformationAllocator;

TrackInformation::TrackInformation()
{
    particleDefinition = 0;
    particleName = "";
    particleTrackID = 0;
    particlePosition = G4ThreeVector(0.,0.,0.);
    particleMomentum = G4ThreeVector(0.,0.,0.);
    particleEnergy = 0.;
    particleProdTime = 0.;
    particleTime = 0.;

    parentDefinition = 0;
    parentName = "";
    parentTrackID = 0;
    parentPosition = G4ThreeVector(0.,0.,0.);
    parentMomentum = G4ThreeVector(0.,0.,0.);
    parentEnergy = 0.;
    parentTime = 0.;

    primaryID = 0;
    primaryPDGID = 0;
    primaryEnergy = 0;
}



TrackInformation::TrackInformation(const G4Track* aTrack)
{
    particleDefinition = aTrack->GetDefinition();
    particleName = aTrack->GetDefinition()->GetParticleName();
    particleTrackID = aTrack->GetTrackID();
    particlePosition = aTrack->GetPosition();
    particleMomentum = aTrack->GetMomentum();
    particleEnergy = aTrack->GetTotalEnergy();
    particleTime = aTrack->GetGlobalTime();

    parentDefinition = aTrack->GetDefinition();
    parentName = aTrack->GetDefinition()->GetParticleName();
    parentTrackID = aTrack->GetTrackID();
    parentPosition = aTrack->GetPosition();
    parentMomentum = aTrack->GetMomentum();
    parentEnergy = aTrack->GetTotalEnergy();
    parentTime = aTrack->GetGlobalTime();

    if(aTrack->GetParentID() ==0) { //is primary
        primaryID = aTrack->GetTrackID();
        primaryPDGID = (G4int) aTrack->GetParticleDefinition()->GetPDGEncoding();
        primaryEnergy = aTrack->GetTotalEnergy();
        //G4cout<<"Primary particle: "<<primaryID<<" with energy: "<<primaryEnergy<<G4endl;
    }
    else{
        primaryID = -1;
        primaryPDGID = -1;
        primaryEnergy = -1.;
    }
}



TrackInformation::TrackInformation(const TrackInformation* aTrackInfo)
{
    particleDefinition = aTrackInfo->particleDefinition;
    particleName = aTrackInfo->particleName;
    particleTrackID = aTrackInfo->particleTrackID;
    particlePosition = aTrackInfo->particlePosition;
    particleMomentum = aTrackInfo->particleMomentum;
    particleEnergy = aTrackInfo->particleEnergy;
    particleProdTime = aTrackInfo->particleProdTime;
    particleTime = aTrackInfo->particleTime;

    parentDefinition = aTrackInfo->parentDefinition;
    parentName = aTrackInfo->parentName;
    parentTrackID = aTrackInfo->parentTrackID;
    parentPosition = aTrackInfo->parentPosition;
    parentMomentum = aTrackInfo->parentMomentum;
    parentEnergy = aTrackInfo->parentEnergy;
    parentTime = aTrackInfo->parentTime;

    primaryID = aTrackInfo->primaryID;
    primaryPDGID = aTrackInfo->primaryPDGID;
    primaryEnergy = aTrackInfo->primaryEnergy;
}


//At each step, update the particle information
void TrackInformation::SetParticleInformation(const TrackInformation* aTrackInfo)
{
    particleDefinition = aTrackInfo->particleDefinition;
    particleName = aTrackInfo->particleName;
    particleTrackID = aTrackInfo->particleTrackID;
    particlePosition = aTrackInfo->particlePosition;
    particleMomentum = aTrackInfo->particleMomentum;
    particleEnergy = aTrackInfo->particleEnergy;
    particleTime = aTrackInfo->particleTime;


}



//When produced save the information of the parent
void TrackInformation::SetParentInformation(const TrackInformation* aTrackInfo)
{
    parentDefinition = aTrackInfo->particleDefinition;
    parentName = aTrackInfo->particleName;
    parentTrackID = aTrackInfo->particleTrackID;
    parentPosition = aTrackInfo->particlePosition;
    parentMomentum = aTrackInfo->particleMomentum;
    parentEnergy = aTrackInfo->particleEnergy;
    parentTime = aTrackInfo->particleTime;

    primaryID = aTrackInfo->primaryID;
    primaryPDGID = aTrackInfo->primaryPDGID;
    primaryEnergy = aTrackInfo->primaryEnergy;
}



TrackInformation::~TrackInformation()
{}



void TrackInformation::Print() const
{
    G4cout << ">>>>>> TrackInformation::Print()::track ID " << particleTrackID
        << " (" << particleName << ")"
        << " at " << particlePosition
        << " with " << particleEnergy/GeV
        << " GeV" << G4endl;

    G4cout << ">>>>>> TrackInformation::Print()::parent track ID " << parentTrackID
        << " (" << parentName << ")"
        << " at " << parentPosition
        << " with " << parentEnergy/GeV
        << " GeV" << G4endl;
}
