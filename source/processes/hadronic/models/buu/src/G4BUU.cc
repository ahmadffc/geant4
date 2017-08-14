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
#include "globals.hh"

//#include <cmath>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4BUU.hh"
#include "G4GenericIon.hh"
#include "G4INCLCascade.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4String.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4CompetitiveFission.hh"

#include "FakeModel.hh"

G4BUU::G4BUU(G4VPreCompoundModel * const aPreCompound) :
  G4VIntraNuclearTransportModel("BUU"),
  thePreCompoundModel(aPreCompound),
  complainedAboutPreCompound(false),
  theIonTable(G4IonTable::GetIonTable())
{
  if(!thePreCompoundModel)
    {
      G4HadronicInteraction* p =
  	G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
      thePreCompoundModel = static_cast<G4VPreCompoundModel*>(p);
      if(!thePreCompoundModel) { thePreCompoundModel = new G4PreCompoundModel; }
    }
  
  // Use the environment variable G4INCLXX_NO_DE_EXCITATION to disable de-excitation
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  theDeExcitation = static_cast<G4VPreCompoundModel*>(p);
  if(!theDeExcitation) { theDeExcitation = new G4PreCompoundModel; }

  // set the fission parameters for G4ExcitationHandler
  G4VEvaporationChannel * const theFissionChannel =
    theDeExcitation->GetExcitationHandler()->GetEvaporation()->GetFissionChannel();
  G4CompetitiveFission * const theFissionChannelCast = dynamic_cast<G4CompetitiveFission *>(theFissionChannel);

  // theBackupModel = new G4BinaryLightIonReaction;
  // theBackupModelNucleon = new G4BinaryCascade;

  theModel = new FakeModel("blob-hot.root");
  currentFragment = 0;
}

G4BUU::~G4BUU()
{
  delete theModel;
}

G4HadFinalState* G4BUU::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  // G4cout<<"G4BUU::ApplyYourself"<<G4endl;
  
  G4HadProjectile const *aProjectileTrack = &aTrack;
  
  G4ParticleDefinition const * const trackDefinition = aTrack.GetDefinition();
  const G4bool isIonTrack = trackDefinition->GetParticleType()==G4GenericIon::GenericIon()->GetParticleType();
  const G4int trackA = trackDefinition->GetAtomicMass();
  const G4int trackZ = (G4int) trackDefinition->GetPDGCharge();
  const G4int nucleusA = theNucleus.GetA_asInt();
  const G4int nucleusZ = theNucleus.GetZ_asInt();

  // Calculate the total four-momentum in the entrance channel
  const G4double trackKinE = aTrack.GetKineticEnergy();
  const G4double theNucleusMass = theIonTable->GetIonMass(nucleusZ, nucleusA);
  const G4double theTrackMass = trackDefinition->GetPDGMass();
  const G4double theTrackEnergy = trackKinE + theTrackMass;
  const G4double theTrackMomentumAbs2 = theTrackEnergy*theTrackEnergy - theTrackMass*theTrackMass;
  const G4double theTrackMomentumAbs = ((theTrackMomentumAbs2>0.0) ? std::sqrt(theTrackMomentumAbs2) : 0.0);
  const G4ThreeVector theTrackMomentum = aTrack.Get4Momentum().getV().unit() * theTrackMomentumAbs;

  G4LorentzVector goodTrack4Momentum(theTrackMomentum, theTrackEnergy);
  G4LorentzVector fourMomentumIn;
  fourMomentumIn.setE(theTrackEnergy + theNucleusMass);
  fourMomentumIn.setVect(theTrackMomentum);

  // INCL assumes the projectile particle is going in the direction of
  // the Z-axis. Here we construct proper rotation to convert the
  // momentum vectors of the outcoming particles to the original
  // coordinate system.
  G4LorentzVector projectileMomentum = aProjectileTrack->Get4Momentum();

  // INCL++ assumes that the projectile is going in the direction of
  // the z-axis. In principle, if the coordinate system used by G4
  // hadronic framework is defined differently we need a rotation to
  // transform the INCL++ reaction products to the appropriate
  // frame. Please note that it isn't necessary to apply this
  // transform to the projectile because when creating the INCL++
  // projectile particle; \see{toINCLKineticEnergy} needs to use only the
  // projectile energy (direction is simply assumed to be along z-axis).
  G4RotationMatrix toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4RotationMatrix toLabFrame3 = toZ.inverse();
  G4LorentzRotation toLabFrame(toLabFrame3);
  // However, it turns out that the projectile given to us by G4
  // hadronic framework is already going in the direction of the
  // z-axis so this rotation is actually unnecessary. Both toZ and
  // toLabFrame turn out to be unit matrices as can be seen by
  // uncommenting the folowing two lines:
  // G4cout <<"toZ = " << toZ << G4endl;
  // G4cout <<"toLabFrame = " << toLabFrame << G4endl;

  G4bool theEventIsGood = false;
  const G4int initialA = trackA + nucleusA;
  const G4int initialZ = trackZ + nucleusZ;
  do
    {
      G4int totalEventA = 0, totalEventZ = 0;
      
      // Make sure the output data structure is clean.
      const G4int nSecondaries = theResult.GetNumberOfSecondaries();
      for(G4int j=0; j<nSecondaries; ++j)
	delete theResult.GetSecondary(j)->GetParticle();      
      theResult.Clear(); 
      theResult.SetStatusChange(stopAndKill);

  std::list<G4Fragment> remnants;
  
  // If the collision was run in inverse kinematics, prepare to boost back
  // all the reaction products
  // if(inverseKinematics) {
  //   toDirectKinematics = new G4LorentzRotation(toInverseKinematics->inverse());
  // }
  
  G4LorentzVector fourMomentumOut;
  G4int firstOfThisEvent = theModel->GetFirstFragmentIdInEvent();
  G4int lastOfThisEvent = theModel->GetLastFragmentIdInEvent(); 
  for (currentFragment = firstOfThisEvent;
       currentFragment<= lastOfThisEvent;
       currentFragment++)
    {
      // G4cout<<"G4BUU::ApplyYourself currentFragment: "<<currentFragment<<G4endl;
      
      theModel->GetEntry(currentFragment);

      G4int A_frag = theModel->A;
      G4int Z_frag = theModel->Z;
      G4double Ek_frag = theModel->EK;
      G4double px_frag = theModel->px;
      G4double py_frag = theModel->py;
      G4double pz_frag = theModel->pz;      
      G4ThreeVector spin(theModel->spinx, theModel->spiny, theModel->spinz);
      
      G4double excitationEnergy = theModel->Eecc;

      totalEventA += A_frag;
      totalEventZ += Z_frag;      
      
      if(excitationEnergy == 0)
	{      
	  // G4int PDGCode = 0;
	  // G4DynamicParticle *p = toG4Particle(theModel->A, theModel->Z, theModel->EK, theModel->px, theModel->py, theModel->pz);
	  G4DynamicParticle *p = toG4Particle(A_frag, Z_frag, Ek_frag, px_frag, py_frag, pz_frag);
	  if(p != 0)
	    {
	      G4LorentzVector momentum = p->Get4Momentum();
	      // Apply the toLabFrame rotation
	      momentum *= toLabFrame;
	      // // Apply the inverse-kinematics boost
	      // if(inverseKinematics)
	      // 	{
	      // 	  momentum *= *toDirectKinematics;
	      // 	  momentum.setVect(-momentum.vect());
	      // 	}
	      // Set the four-momentum of the reaction products
	      p->Set4Momentum(momentum);
	      fourMomentumOut += momentum;
	      theResult.AddSecondary(p);
	      
	    }
	  else
	    {
	      G4String message = "the model produced a particle that couldn't be converted to Geant4 particle.";
	      G4cout<<message<<G4endl;
	    }
	}
      else //it has to be de-excited
	{
	  // const G4double nuclearMass = G4NucleiProperties::GetNuclearMass(theModel->A, theModel->Z) + theModel->Eecc;
	  const G4double nuclearMass = G4NucleiProperties::GetNuclearMass(A_frag, Z_frag) + excitationEnergy;

	  // G4LorentzVector fourMomentum(theModel->px, theModel->py, theModel->pz,
	  // 			       nuclearMass + theModel->EK);
	  G4LorentzVector fourMomentum(px_frag, py_frag, pz_frag,
				       nuclearMass + Ek_frag);

	  // Apply the toLabFrame rotation
	  fourMomentum *= toLabFrame;
	  spin *= toLabFrame3;
	  
	  G4Fragment remnant(A_frag, Z_frag, fourMomentum );	  
	  remnant.SetAngularMomentum(spin);
	  G4ReactionProductVector *deExcitationResult = theDeExcitation->DeExcite(remnant);
	  for(G4ReactionProductVector::iterator fragment = deExcitationResult->begin();
	      fragment != deExcitationResult->end(); ++fragment)
	    {
	      const G4ParticleDefinition *def = (*fragment)->GetDefinition();
	      if(def != 0)
		{
		  G4DynamicParticle *theFragment = new G4DynamicParticle(def, (*fragment)->GetMomentum());
		  theResult.AddSecondary(theFragment);
		}
	    }
	  for(G4ReactionProductVector::iterator fragment = deExcitationResult->begin(); fragment != deExcitationResult->end(); ++fragment)
	    {
	      
	      delete (*fragment);
	    }
	  deExcitationResult->clear();
	  delete deExcitationResult;
	}
      
    }//loop on all the products of the model

  if( (totalEventA == initialA) && (totalEventZ == initialZ) )
    {
      theEventIsGood = true;
    }
  else
    {
      G4int epReportLevel = getenv("G4Hadronic_epReportLevel") ?
	strtol(getenv("G4Hadronic_epReportLevel"),0,10) : 0;
      G4ExceptionDescription desc;
      desc << "Warning: non-conservation detected of mass or charge in BUU, will "
	   << (epReportLevel<0 ? "abort the event" :  "re-sample the interaction") << G4endl
	   << " initial A = " << initialA << " total fragments A = " << totalEventA << G4endl
	   << " initial Z = " << initialZ << " total fragments Z = " << totalEventZ << G4endl;

      // const G4int nSecondaries = theResult.GetNumberOfSecondaries();
      // for(G4int j=0; j<nSecondaries; ++j)
      // 	{
      // 	  theResult.GetSecondary(j)->GetParticle();
      // 	  desc<<
      // 	}
	  
      G4Exception("G4BUU::ApplyYourself", "had???", 
		  epReportLevel<0 ? EventMustBeAborted : JustWarning, desc);
    }
  
    }while(!theEventIsGood);

  return &theResult;
}

G4ReactionProductVector* G4BUU::Propagate(G4KineticTrackVector* , G4V3DNucleus* ) {
  return 0;
}

G4ParticleDefinition *G4BUU::toG4ParticleDefinition(G4int A, G4int Z) const {
  if     (A == 1 && Z == 1)  return G4Proton::Proton();
  else if(A == 1 && Z == 0)  return G4Neutron::Neutron();
  else if(A == 0 && Z == 1)  return G4PionPlus::PionPlus();
  else if(A == 0 && Z == -1) return G4PionMinus::PionMinus();
  else if(A == 0 && Z == 0)  {
    // if (PDGCode == 221) return G4Eta::Eta();
    // else if (PDGCode == 22)
    return G4Gamma::Gamma();
    // else return G4PionZero::PionZero();
  }
  else if(A == 2 && Z == 1)  return G4Deuteron::Deuteron();
  else if(A == 3 && Z == 1)  return G4Triton::Triton();
  else if(A == 3 && Z == 2)  return G4He3::He3();
  else if(A == 4 && Z == 2)  return G4Alpha::Alpha();
  else if(A > 0 && Z > 0 && A > Z) { // Returns ground state ion definition
    return theIonTable->GetIon(Z, A, 0);
  } else { // Error, unrecognized particle
    return 0;
  }
}

G4DynamicParticle *G4BUU::toG4Particle(G4int A, G4int Z, //G4int PDGCode,
						 G4double kinE,
						 G4double px,
                                                 G4double py, G4double pz) const {
  const G4ParticleDefinition *def = toG4ParticleDefinition(A, Z);//, PDGCode);
  if(def == 0) { // Check if we have a valid particle definition
    return 0;
  }
  const G4double energy = kinE * MeV;
  const G4ThreeVector momentum(px, py, pz);
  const G4ThreeVector momentumDirection = momentum.unit();
  G4DynamicParticle *p = new G4DynamicParticle(def, momentumDirection, energy);
  return p;
}

void G4BUU::ModelDescription(std::ostream& outFile) const {
   outFile
     << "A fake model\n";
}

G4String const &G4BUU::GetDeExcitationModelName() const {
  return theDeExcitation->GetModelName();
}
