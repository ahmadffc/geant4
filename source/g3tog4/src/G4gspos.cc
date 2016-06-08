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
// $Id: G4gspos.cc,v 1.12 2001/07/11 09:59:02 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// by I.Hrivnacova, 13.10.99

#include "G3G4Interface.hh"
#include "G3VolTable.hh"
#include "G3toG4.hh"
#include "G3Pos.hh"
#include "globals.hh"

void G4CreateCloneVTE(G3VolTableEntry* vte, G3VolTableEntry* mvte,
              G4double pars[], G4int npar, G4int num,
              G4double x, G4double y, G4double z, G4int irot, G4String vonly);

void PG4gspos(G4String tokens[])
{
        // fill the parameter containers
    G3fillParams(tokens,PTgspos);
  
        // interpret the parameters
    G4String name = Spar[0];
    G4String moth = Spar[1];
    G4String only = Spar[2];
    G4int num = Ipar[0];
    G4int irot = Ipar[1];
        // all parameters are passed to G4gsxxx methods
        // in G3 default units         
    //G4double x = Rpar[0]*cm;
    //G4double y = Rpar[1]*cm;
    //G4double z = Rpar[2]*cm;
    G4double x = Rpar[0];
    G4double y = Rpar[1];
    G4double z = Rpar[2];
  
    G4gspos(name, num, moth, x, y, z, irot, only);
}

void G4gspos(G4String vname, G4int num, G4String vmoth, G4double x,
             G4double y, G4double z, G4int irot, G4String vonly)
{
  // find VTEs
  G3VolTableEntry* vte = G3Vol.GetVTE(vname);
  G3VolTableEntry* mvte = G3Vol.GetVTE(vmoth);

  if (vte == 0) {
    G4Exception("G4gspos: '" + vname + "' has no VolTableEntry");
  } 
  else if (mvte == 0) {
    G4Exception("G4gspos: '" + vmoth + "' has no VolTableEntry");
  } 
  else { 
    if (!vte->HasNegPars()) {
      // position vector
      G4ThreeVector* offset = new G4ThreeVector(x*cm, y*cm, z*cm);

      // create a G3Pos object and add it to the vte
      G3Pos* aG3Pos = new G3Pos(vmoth, num, offset, irot, vonly);              
      vte->AddG3Pos(aG3Pos);

      // loop over all mothers
      for (G4int i=0; i<mvte->GetNoClones(); i++) {
                       // (mvte is retrieved from its "master" name
                       //  -> there is no need to call GetMasterClone()
        G3VolTableEntry* mvteClone = mvte->GetClone(i);
        vte->AddMother(mvteClone);
        mvteClone->AddDaughter(vte);
      }     
    } 
    else {
      // if vte has neg parameters 
      // a new vte clone copy is created for each mother (clone copy)  
      // and its parameters are derived from it if possible

      G4CreateCloneVTE(vte, mvte, vte->GetRpar(), vte->GetNpar(), num,
                     x, y, z, irot, vonly);
    }    
  }      
}