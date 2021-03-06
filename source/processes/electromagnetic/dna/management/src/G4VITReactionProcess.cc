//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4VITReactionProcess.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4VITReactionProcess.hh"

G4VITReactionProcess::G4VITReactionProcess() : fpReactionTable(0), fpChanges(0)
{
    //ctor
}

G4VITReactionProcess::~G4VITReactionProcess()
{
    //dtor
}

G4VITReactionProcess::G4VITReactionProcess(const G4VITReactionProcess& /*other*/)
{
    //copy ctor
    fpChanges = 0;
    fpReactionTable = 0;
}

G4VITReactionProcess& G4VITReactionProcess::operator=(const G4VITReactionProcess& rhs)
{
    //assignment operator
    if (this == &rhs) return *this; // handle self assignment
//    fApplicableType1 = rhs.fApplicableType1;
//    fApplicableType2 = rhs.fApplicableType2;

    fName       = rhs.fName;
    return *this;
}

//void G4VITReactionProcess::GetApplicableTypes(G4ITType& type1, G4ITType& type2) const
//{
//    type1 = fApplicableType1;
//    type2 = fApplicableType2;
//}
