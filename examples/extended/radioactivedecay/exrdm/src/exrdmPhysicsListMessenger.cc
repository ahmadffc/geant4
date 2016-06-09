//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "exrdmPhysicsListMessenger.hh"

#include "exrdmPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysicsListMessenger::exrdmPhysicsListMessenger(exrdmPhysicsList* pPhys)
:pPhysicsList(pPhys)
{   
  physDir = new G4UIdirectory("/exrdm/phys/");
  physDir->SetGuidance("physics control.");

  gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/phys/setGCut",this);  
  gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut",false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>0.0");
  gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/phys/setECut",this);  
  electCutCmd->SetGuidance("Set electron cut.");
  electCutCmd->SetParameterName("Ecut",false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>0.0");
  electCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  protoCutCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/phys/setPCut",this);  
  protoCutCmd->SetGuidance("Set positron cut.");
  protoCutCmd->SetParameterName("Pcut",false);
  protoCutCmd->SetUnitCategory("Length");
  protoCutCmd->SetRange("Pcut>0.0");
  protoCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  allCutCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/phys/setCuts",this);
  allCutCmd->SetGuidance("Set cut for all.");
  allCutCmd->SetParameterName("cut",false);
  allCutCmd->SetUnitCategory("Length");
  allCutCmd->SetRange("cut>0.0");
  allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  eCutCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/phys/TargetCuts",this);
  eCutCmd->SetGuidance("Set cuts for the target");
  eCutCmd->SetParameterName("Ecut",false);
  eCutCmd->SetUnitCategory("Length");
  eCutCmd->SetRange("Ecut>0.0");
  eCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mCutCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/phys/DetectorCuts",this);
  mCutCmd->SetGuidance("Set cuts for the Detector");
  mCutCmd->SetParameterName("Ecut",false);
  mCutCmd->SetUnitCategory("Length");
  mCutCmd->SetRange("Ecut>0.0");
  mCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  pListCmd = new G4UIcmdWithAString("/exrdm/phys/SelectPhysics",this);
  pListCmd->SetGuidance("Select modula physics list.");
  pListCmd->SetParameterName("PList",false);
  pListCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysicsListMessenger::~exrdmPhysicsListMessenger()
{
  delete physDir;
  delete gammaCutCmd;
  delete electCutCmd;
  delete protoCutCmd;
  delete allCutCmd;
  delete pListCmd;
  delete eCutCmd;
  delete mCutCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if( command == gammaCutCmd )
   { pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));}

  if( command == electCutCmd )
   { pPhysicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));}

  if( command == protoCutCmd )
   { pPhysicsList->SetCutForPositron(protoCutCmd->GetNewDoubleValue(newValue));}

  if( command == allCutCmd )
    {
      G4double cut = allCutCmd->GetNewDoubleValue(newValue);
      pPhysicsList->SetCutForGamma(cut);
      pPhysicsList->SetCutForElectron(cut);
      pPhysicsList->SetCutForPositron(cut);
    }

  if( command == eCutCmd )
   { pPhysicsList->SetTargetCut(eCutCmd->GetNewDoubleValue(newValue));}

  if( command == mCutCmd )
   { pPhysicsList->SetDetectorCut(mCutCmd->GetNewDoubleValue(newValue));}

  if( command == pListCmd )
   { pPhysicsList->SelectPhysicsList(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......