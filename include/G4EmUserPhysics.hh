// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch
//
// ********************************************************************
// * License and Disclaimer *
// * *
// * The Geant4 software is copyright of the Copyright Holders of *
// * the Geant4 Collaboration. It is provided under the terms and *
// * conditions of the Geant4 Software License, included in the file *
// * LICENSE and available at http://cern.ch/geant4/license . These *
// * include a list of copyright holders. *
// * *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work make any representation or warranty, express or implied, *
// * regarding this software system or assume any liability for its *
// * use. Please see the license in the file LICENSE and URL above *
// * for the full disclaimer and the limitation of liability. *
// * *
// * This code implementation is the result of the scientific and *
// * technical work of the GEANT4 collaboration. *
// * By using, copying, modifying or distributing the software (or *
// * any work based on the software) you agree to acknowledge its *
// * use in resulting scientific publications, and indicate your *
// * acceptance of all terms of the Geant4 Software license. *
// ********************************************************************
//
/// \file hadronic/Hadr01/include/G4EmUserPhysics.hh
/// \brief Definition of the G4EmUserPhysics class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName: G4EmUserPhysics
//
// Author: V.Ivanchenko 11.07.2012
//
// Modified:
//
//----------------------------------------------------------------------------
//
// This class shows how extra EM options can be defined on top of
// any reference Physics List
//

#ifndef G4EmUserPhysics_h
#define G4EmUserPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmUserPhysics : public G4VPhysicsConstructor
{
public:

  G4EmUserPhysics(const G4int& scint, const G4int& cher);

  virtual ~G4EmUserPhysics();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
  // virtual void SetCuts();

private:
  G4int switchOnScintillation;
  G4int switchOnCerenkov;

  G4Cerenkov * theCerenkovProcess;
  G4Scintillation * theScintillationProcess;
  G4OpAbsorption * theAbsorptionProcess;
  G4OpRayleigh * theRayleighScatteringProcess;
  G4OpMieHG * theMieHGScatteringProcess;
  G4OpBoundaryProcess * theBoundaryProcess;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
