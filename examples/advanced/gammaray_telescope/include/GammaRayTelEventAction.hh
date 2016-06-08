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
//
// $Id: GammaRayTelEventAction.hh,v 1.9 2001/12/04 11:40:28 flongo Exp $
// GEANT4 tag $Name: geant4-04-00 $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelEventAction  ------
//           by R.Giannitrapani, F. Longo & G.Santin (13 nov 2000)
//
//- inclusion of Digits by F.Longo & R.Giannitrapani (24 oct 2001)
//  20.11.01 G.Santin: new analysis management, modified according to GammaRayTelAnalysis
// ************************************************************


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelEventAction_h
#define GammaRayTelEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelEventAction : public G4UserEventAction
{
public:

  GammaRayTelEventAction();
  virtual ~GammaRayTelEventAction();
  
public:
  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  
  void SetDrawFlag   (G4String val)  {drawFlag = val;};
  
private:

  G4int       trackerCollID;                
  G4int       calorimeterCollID;                
  G4int       anticoincidenceCollID;                
  G4String    drawFlag;
};

#endif

    




