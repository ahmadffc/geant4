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

#ifndef G4BUU_hh
#define G4BUU_hh 1

#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4HadronicInteraction.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "globals.hh"

// Geant4 de-excitation
#include "G4ExcitationHandler.hh"

// Binary cascade
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"

// PreCompound
#include "G4VPreCompoundModel.hh"
#include "G4PreCompoundModel.hh"

// G4IonTable
#include "G4IonTable.hh"

// fission
#include "G4VLevelDensityParameter.hh"
#include "G4FissionProbability.hh"

// #include <fstream>
// #include <iostream>

/** \brief INCL++ intra-nuclear cascade
 *
 * Interface for INCL++. This interface handles basic hadron bullet particles
 * (protons, neutrons, pions), as well as light ions.
 *
 * Example usage in case of protons:
 * @code
 * G4INCLXXInterface* inclModel = new G4INCLXXInterface;
 * inclModel -> SetMinEnergy(0.0 * MeV); // Set the energy limits
 * inclModel -> SetMaxEnergy(3.0 * GeV);
 *
 * G4ProtonInelasticProcess* protonInelasticProcess = new G4ProtonInelasticProcess();
 * G4ProtonInelasticCrossSection* protonInelasticCrossSection =  new G4ProtonInelasticCrossSection();
 *
 * protonInelasticProcess -> RegisterMe(inclModel);
 * protonInelasticProcess -> AddDataSet(protonInelasticCrossSection);
 *
 * particle = G4Proton::Proton();
 * processManager = particle -> GetProcessManager();
 * processManager -> AddDiscreteProcess(protonInelasticProcess);
 * @endcode
 * The same setup procedure is needed for neutron, pion and generic-ion
 * inelastic processes as well.
 */
class FakeModel;

class G4BUU : public G4VIntraNuclearTransportModel {
public:
  G4BUU(G4VPreCompoundModel * const aPreCompound = 0);
  ~G4BUU(); // Destructor

  G4int operator==(G4BUU& right) {
    return (this == &right);
  }

  G4int operator!=(G4BUU& right) {
    return (this != &right);
  }

  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus); // Idle

  /**
   * Main method to apply the INCL physics model.
   * @param aTrack the projectile particle
   * @param theNucleus target nucleus
   * @return the output of the INCL physics model
   */
  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,  G4Nucleus& theNucleus);

  using G4VIntraNuclearTransportModel::SetDeExcitation;

  virtual void ModelDescription(std::ostream& outFile) const;

  G4String const &GetDeExcitationModelName() const;

private:
  /// \brief Dummy copy constructor to shut up Coverity warnings
  G4BUU(const G4BUU &rhs);

  /// \brief Dummy assignment operator to shut up Coverity warnings
  G4BUU &operator=(G4BUU const &rhs);

    /// \brief Convert an INCL particle to a G4DynamicParticle
  G4DynamicParticle *toG4Particle(G4int A, G4int Z, G4double kinE, G4double px, G4double py, G4double pz) const;
  
  /// \brief Convert A and Z to a G4ParticleDefinition
  G4ParticleDefinition *toG4ParticleDefinition (G4int A, G4int Z) const;

  G4VPreCompoundModel *thePreCompoundModel;

  G4HadFinalState theResult;

  // G4HadronicInteraction *theBackupModel;
  // G4HadronicInteraction *theBackupModelNucleon;


  G4bool complainedAboutPreCompound;

  G4IonTable * const theIonTable;

  G4bool dumpRemnantInfo;

  FakeModel* theModel;
  G4int currentFragment;
};

#endif
