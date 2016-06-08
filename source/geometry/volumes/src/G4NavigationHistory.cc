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
// $Id: G4NavigationHistory.cc,v 1.4 2001/07/11 10:00:32 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// 
// G4NavigationHistory Implementation P.Kent August 96

#include "G4NavigationHistory.hh"
#include "G4ios.hh"

G4std::ostream&
operator << (G4std::ostream& os, const G4NavigationHistory& nav)
{
  G4cout << "History depth="<<nav.GetDepth()<< G4endl;
  for (G4int i=0;i<=nav.GetDepth();i++)
  {
      os << "Level=["<<i<<"]: " ;
      if( nav.GetVolume(i) != 0 ) {
	 os   << "Phys Name=["<< nav.GetVolume(i)->GetName()
	      << "] Type=[";
	 switch(nav.GetVolumeType(i))
	   {
	   case kNormal:
	     os <<"N";
	     break;
	   case kReplica:
	     os <<"R" << nav.GetReplicaNo(i);
	     break;
	   case kParameterised:
	     os <<"P" << nav.GetReplicaNo(i);
	     break;
	   }
	 os << "]";
      }else{
	 os << "Phys = <Null>";
      }
      os << G4endl;
  }
  return os;
}
