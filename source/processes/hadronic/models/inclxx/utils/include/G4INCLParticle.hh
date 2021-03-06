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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * G4INCLParticle.hh
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef PARTICLE_HH_
#define PARTICLE_HH_

#include "G4INCLThreeVector.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLParticleType.hh"
#include "G4INCLParticleSpecies.hh"
#include "G4INCLLogger.hh"
#include "G4INCLUnorderedVector.hh"
#include "G4INCLAllocationPool.hh"
#include <sstream>
#include <string>

namespace G4INCL {

  class Particle;

  class ParticleList : public UnorderedVector<Particle*> {
    public:
      void rotatePositionAndMomentum(const G4double angle, const ThreeVector &axis) const;
      void rotatePosition(const G4double angle, const ThreeVector &axis) const;
      void rotateMomentum(const G4double angle, const ThreeVector &axis) const;
      void boost(const ThreeVector &b) const;
  };

  typedef ParticleList::const_iterator ParticleIter;
  typedef ParticleList::iterator       ParticleMutableIter;

  class Particle {
  public:
    Particle();
    Particle(ParticleType t, G4double energy, ThreeVector const &momentum, ThreeVector const &position);
    Particle(ParticleType t, ThreeVector const &momentum, ThreeVector const &position);
    virtual ~Particle() {}

    /** \brief Copy constructor
     *
     * Does not copy the particle ID.
     */
    Particle(const Particle &rhs) :
      theZ(rhs.theZ),
      theA(rhs.theA),
      theParticipantType(rhs.theParticipantType),
      theType(rhs.theType),
      theEnergy(rhs.theEnergy),
      theFrozenEnergy(rhs.theFrozenEnergy),
      theMomentum(rhs.theMomentum),
      theFrozenMomentum(rhs.theFrozenMomentum),
      thePosition(rhs.thePosition),
      nCollisions(rhs.nCollisions),
      nDecays(rhs.nDecays),
      thePotentialEnergy(rhs.thePotentialEnergy),
      rpCorrelated(rhs.rpCorrelated),
      uncorrelatedMomentum(rhs.uncorrelatedMomentum),
      theHelicity(rhs.theHelicity),
      emissionTime(rhs.emissionTime),
      outOfWell(rhs.outOfWell),
      theMass(rhs.theMass)
      {
        if(rhs.thePropagationEnergy == &(rhs.theFrozenEnergy))
          thePropagationEnergy = &theFrozenEnergy;
        else
          thePropagationEnergy = &theEnergy;
        if(rhs.thePropagationMomentum == &(rhs.theFrozenMomentum))
          thePropagationMomentum = &theFrozenMomentum;
        else
          thePropagationMomentum = &theMomentum;
        // ID intentionally not copied
        ID = nextID++;
      }

  protected:
    /// \brief Helper method for the assignment operator
    void swap(Particle &rhs) {
      std::swap(theZ, rhs.theZ);
      std::swap(theA, rhs.theA);
      std::swap(theParticipantType, rhs.theParticipantType);
      std::swap(theType, rhs.theType);
      if(rhs.thePropagationEnergy == &(rhs.theFrozenEnergy))
        thePropagationEnergy = &theFrozenEnergy;
      else
        thePropagationEnergy = &theEnergy;
      std::swap(theEnergy, rhs.theEnergy);
      std::swap(theFrozenEnergy, rhs.theFrozenEnergy);
      if(rhs.thePropagationMomentum == &(rhs.theFrozenMomentum))
        thePropagationMomentum = &theFrozenMomentum;
      else
        thePropagationMomentum = &theMomentum;
      std::swap(theMomentum, rhs.theMomentum);
      std::swap(theFrozenMomentum, rhs.theFrozenMomentum);
      std::swap(thePosition, rhs.thePosition);
      std::swap(nCollisions, rhs.nCollisions);
      std::swap(nDecays, rhs.nDecays);
      std::swap(thePotentialEnergy, rhs.thePotentialEnergy);
      // ID intentionally not swapped

      std::swap(theHelicity, rhs.theHelicity);
      std::swap(emissionTime, rhs.emissionTime);
      std::swap(outOfWell, rhs.outOfWell);

      std::swap(theMass, rhs.theMass);
      std::swap(rpCorrelated, rhs.rpCorrelated);
      std::swap(uncorrelatedMomentum, rhs.uncorrelatedMomentum);
    }

  public:

    /** \brief Assignment operator
     *
     * Does not copy the particle ID.
     */
    Particle &operator=(const Particle &rhs) {
      Particle temporaryParticle(rhs);
      swap(temporaryParticle);
      return *this;
    }

    /**
     * Get the particle type.
     * @see G4INCL::ParticleType
     */
    G4INCL::ParticleType getType() const {
      return theType;
    };

    /// \brief Get the particle species
    virtual G4INCL::ParticleSpecies getSpecies() const {
      return ParticleSpecies(theType);
    };

    void setType(ParticleType t) {
      theType = t;
      switch(theType)
      {
        case DeltaPlusPlus:
          theA = 1;
          theZ = 2;
          break;
        case Proton:
        case DeltaPlus:
          theA = 1;
          theZ = 1;
          break;
        case Neutron:
        case DeltaZero:
          theA = 1;
          theZ = 0;
          break;
        case DeltaMinus:
          theA = 1;
          theZ = -1;
          break;
        case PiPlus:
          theA = 0;
          theZ = 1;
          break;
        case PiZero:
		case Eta:
		case Omega:
		case EtaPrime:
		case Photon:
          theA = 0;
          theZ = 0;
          break;
        case PiMinus:
          theA = 0;
          theZ = -1;
          break;
        case Composite:
         // INCL_ERROR("Trying to set particle type to Composite! Construct a Cluster object instead" << '\n');
          theA = 0;
          theZ = 0;
          break;
        case UnknownParticle:
          theA = 0;
          theZ = 0;
          INCL_ERROR("Trying to set particle type to Unknown!" << '\n');
          break;
      }

      if( !isResonance() && t!=Composite )
        setINCLMass();
    }

    /**
     * Is this a nucleon?
     */
    G4bool isNucleon() const {
      if(theType == G4INCL::Proton || theType == G4INCL::Neutron)
	return true;
      else
	return false;
    };

    ParticipantType getParticipantType() const {
      return theParticipantType;
    }

    void setParticipantType(ParticipantType const p) {
      theParticipantType = p;
    }

    G4bool isParticipant() const {
      return (theParticipantType==Participant);
    }

    G4bool isTargetSpectator() const {
      return (theParticipantType==TargetSpectator);
    }

    G4bool isProjectileSpectator() const {
      return (theParticipantType==ProjectileSpectator);
    }

    virtual void makeParticipant() {
      theParticipantType = Participant;
    }

    virtual void makeTargetSpectator() {
      theParticipantType = TargetSpectator;
    }

    virtual void makeProjectileSpectator() {
      theParticipantType = ProjectileSpectator;
    }

    /** \brief Is this a pion? */
    G4bool isPion() const { return (theType == PiPlus || theType == PiZero || theType == PiMinus); }

    /** \brief Is this a eta? */
    G4bool isEta() const { return (theType == Eta); }

    /** \brief Is this a omega? */
    G4bool isOmega() const { return (theType == Omega); }

    /** \brief Is this a etaprime? */
    G4bool isEtaPrime() const { return (theType == EtaPrime); }

    /** \brief Is it a resonance? */
    inline G4bool isResonance() const { return isDelta(); }

    /** \brief Is it a Delta? */
    inline G4bool isDelta() const {
      return (theType==DeltaPlusPlus || theType==DeltaPlus ||
          theType==DeltaZero || theType==DeltaMinus);
    }

    /** \brief Returns the baryon number. */
    G4int getA() const { return theA; }

    /** \brief Returns the charge number. */
    G4int getZ() const { return theZ; }

    G4double getBeta() const {
      const G4double P = theMomentum.mag();
      return P/theEnergy;
    }

    /**
     * Returns a three vector we can give to the boost() -method.
     *
     * In order to go to the particle rest frame you need to multiply
     * the boost vector by -1.0.
     */
    ThreeVector boostVector() const {
      return theMomentum / theEnergy;
    }

    /**
     * Boost the particle using a boost vector.
     *
     * Example (go to the particle rest frame):
     * particle->boost(particle->boostVector());
     */
    void boost(const ThreeVector &aBoostVector) {
      const G4double beta2 = aBoostVector.mag2();
      const G4double gamma = 1.0 / std::sqrt(1.0 - beta2);
      const G4double bp = theMomentum.dot(aBoostVector);
      const G4double alpha = (gamma*gamma)/(1.0 + gamma);

      theMomentum = theMomentum + aBoostVector * (alpha * bp - gamma * theEnergy);
      theEnergy = gamma * (theEnergy - bp);
    }

    /** \brief Lorentz-contract the particle position around some center
     *
     * Apply Lorentz contraction to the position component along the
     * direction of the boost vector.
     *
     * \param aBoostVector the boost vector (velocity) [c]
     * \param refPos the reference position
     */
    void lorentzContract(const ThreeVector &aBoostVector, const ThreeVector &refPos) {
      const G4double beta2 = aBoostVector.mag2();
      const G4double gamma = 1.0 / std::sqrt(1.0 - beta2);
      const ThreeVector theRelativePosition = thePosition - refPos;
      const ThreeVector transversePosition = theRelativePosition - aBoostVector * (theRelativePosition.dot(aBoostVector) / aBoostVector.mag2());
      const ThreeVector longitudinalPosition = theRelativePosition - transversePosition;

      thePosition = refPos + transversePosition + longitudinalPosition / gamma;
    }

    /** \brief Get the cached particle mass. */
    inline G4double getMass() const { return theMass; }

    /** \brief Get the INCL particle mass. */
    inline G4double getINCLMass() const {
      switch(theType) {
        case Proton:
        case Neutron:
        case PiPlus:
        case PiMinus:
        case PiZero:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:
          return ParticleTable::getINCLMass(theType);
          break;

        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return theMass;
          break;

        case Composite:
          return ParticleTable::getINCLMass(theA,theZ);
          break;

        default:
          INCL_ERROR("Particle::getINCLMass: Unknown particle type." << '\n');
          return 0.0;
          break;
      }
    }

    /** \brief Get the tabulated particle mass. */
    inline virtual G4double getTableMass() const {
      switch(theType) {
        case Proton:
        case Neutron:
        case PiPlus:
        case PiMinus:
        case PiZero:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:
          return ParticleTable::getTableParticleMass(theType);
          break;

        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return theMass;
          break;

        case Composite:
          return ParticleTable::getTableMass(theA,theZ);
          break;

        default:
          INCL_ERROR("Particle::getTableMass: Unknown particle type." << '\n');
          return 0.0;
          break;
      }
    }

    /** \brief Get the real particle mass. */
    inline G4double getRealMass() const {
      switch(theType) {
        case Proton:
        case Neutron:
        case PiPlus:
        case PiMinus:
        case PiZero:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:
          return ParticleTable::getRealMass(theType);
          break;

        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return theMass;
          break;

        case Composite:
          return ParticleTable::getRealMass(theA,theZ);
          break;

        default:
          INCL_ERROR("Particle::getRealMass: Unknown particle type." << '\n');
          return 0.0;
          break;
      }
    }

    /// \brief Set the mass of the Particle to its real mass
    void setRealMass() { setMass(getRealMass()); }

    /// \brief Set the mass of the Particle to its table mass
    void setTableMass() { setMass(getTableMass()); }

    /// \brief Set the mass of the Particle to its table mass
    void setINCLMass() { setMass(getINCLMass()); }

    /**\brief Computes correction on the emission Q-value
     *
     * Computes the correction that must be applied to INCL particles in
     * order to obtain the correct Q-value for particle emission from a given
     * nucleus. For absorption, the correction is obviously equal to minus
     * the value returned by this function.
     *
     * \param AParent the mass number of the emitting nucleus
     * \param ZParent the charge number of the emitting nucleus
     * \return the correction
     */
    G4double getEmissionQValueCorrection(const G4int AParent, const G4int ZParent) const {
      const G4int ADaughter = AParent - theA;
      const G4int ZDaughter = ZParent - theZ;

      // Note the minus sign here
      G4double theQValue;
      if(isCluster())
        theQValue = -ParticleTable::getTableQValue(theA, theZ, ADaughter, ZDaughter);
      else {
        const G4double massTableParent = ParticleTable::getTableMass(AParent,ZParent);
        const G4double massTableDaughter = ParticleTable::getTableMass(ADaughter,ZDaughter);
        const G4double massTableParticle = getTableMass();
        theQValue = massTableParent - massTableDaughter - massTableParticle;
      }

      const G4double massINCLParent = ParticleTable::getINCLMass(AParent,ZParent);
      const G4double massINCLDaughter = ParticleTable::getINCLMass(ADaughter,ZDaughter);
      const G4double massINCLParticle = getINCLMass();

      // The rhs corresponds to the INCL Q-value
      return theQValue - (massINCLParent-massINCLDaughter-massINCLParticle);
    }

    /**\brief Computes correction on the transfer Q-value
     *
     * Computes the correction that must be applied to INCL particles in
     * order to obtain the correct Q-value for particle transfer from a given
     * nucleus to another.
     *
     * Assumes that the receving nucleus is INCL's target nucleus, with the
     * INCL separation energy.
     *
     * \param AFrom the mass number of the donating nucleus
     * \param ZFrom the charge number of the donating nucleus
     * \param ATo the mass number of the receiving nucleus
     * \param ZTo the charge number of the receiving nucleus
     * \return the correction
     */
    G4double getTransferQValueCorrection(const G4int AFrom, const G4int ZFrom, const G4int ATo, const G4int ZTo) const {
      const G4int AFromDaughter = AFrom - theA;
      const G4int ZFromDaughter = ZFrom - theZ;
      const G4int AToDaughter = ATo + theA;
      const G4int ZToDaughter = ZTo + theZ;
      const G4double theQValue = ParticleTable::getTableQValue(AToDaughter,ZToDaughter,AFromDaughter,ZFromDaughter,AFrom,ZFrom);

      const G4double massINCLTo = ParticleTable::getINCLMass(ATo,ZTo);
      const G4double massINCLToDaughter = ParticleTable::getINCLMass(AToDaughter,ZToDaughter);
      /* Note that here we have to use the table mass in the INCL Q-value. We
       * cannot use theMass, because at this stage the particle is probably
       * still off-shell; and we cannot use getINCLMass(), because it leads to
       * violations of global energy conservation.
       */
      const G4double massINCLParticle = getTableMass();

      // The rhs corresponds to the INCL Q-value for particle absorption
      return theQValue - (massINCLToDaughter-massINCLTo-massINCLParticle);
    }

    /** \brief Get the the particle invariant mass.
     *
     * Uses the relativistic invariant
     * \f[ m = \sqrt{E^2 - {\vec p}^2}\f]
     **/
    G4double getInvariantMass() const {
      const G4double mass = std::pow(theEnergy, 2) - theMomentum.dot(theMomentum);
      if(mass < 0.0) {
        INCL_ERROR("E*E - p*p is negative." << '\n');
        return 0.0;
      } else {
        return std::sqrt(mass);
      }
    };

    /// \brief Get the particle kinetic energy.
    inline G4double getKineticEnergy() const { return theEnergy - theMass; }

    /// \brief Get the particle potential energy.
    inline G4double getPotentialEnergy() const { return thePotentialEnergy; }

    /// \brief Set the particle potential energy.
    inline void setPotentialEnergy(G4double v) { thePotentialEnergy = v; }

    /**
     * Get the energy of the particle in MeV.
     */
    G4double getEnergy() const
    {
      return theEnergy;
    };

    /**
     * Set the mass of the particle in MeV/c^2.
     */
    void setMass(G4double mass)
    {
      this->theMass = mass;
    }

    /**
     * Set the energy of the particle in MeV.
     */
    void setEnergy(G4double energy)
    {
      this->theEnergy = energy;
    };

    /**
     * Get the momentum vector.
     */
    const G4INCL::ThreeVector &getMomentum() const
    {
      return theMomentum;
    };

    /** Get the angular momentum w.r.t. the origin */
    virtual G4INCL::ThreeVector getAngularMomentum() const
    {
      return thePosition.vector(theMomentum);
    };

    /**
     * Set the momentum vector.
     */
    virtual void setMomentum(const G4INCL::ThreeVector &momentum)
    {
      this->theMomentum = momentum;
    };

    /**
     * Set the position vector.
     */
    const G4INCL::ThreeVector &getPosition() const
    {
      return thePosition;
    };

    virtual void setPosition(const G4INCL::ThreeVector &position)
    {
      this->thePosition = position;
    };

    G4double getHelicity() { return theHelicity; };
    void setHelicity(G4double h) { theHelicity = h; };

    void propagate(G4double step) {
      thePosition += ((*thePropagationMomentum)*(step/(*thePropagationEnergy)));
    };

    /** \brief Return the number of collisions undergone by the particle. **/
    G4int getNumberOfCollisions() const { return nCollisions; }

    /** \brief Set the number of collisions undergone by the particle. **/
    void setNumberOfCollisions(G4int n) { nCollisions = n; }

    /** \brief Increment the number of collisions undergone by the particle. **/
    void incrementNumberOfCollisions() { nCollisions++; }

    /** \brief Return the number of decays undergone by the particle. **/
    G4int getNumberOfDecays() const { return nDecays; }

    /** \brief Set the number of decays undergone by the particle. **/
    void setNumberOfDecays(G4int n) { nDecays = n; }

    /** \brief Increment the number of decays undergone by the particle. **/
    void incrementNumberOfDecays() { nDecays++; }

    /** \brief Mark the particle as out of its potential well
     *
     * This flag is used to control pions created outside their potential well
     * in delta decay. The pion potential checks it and returns zero if it is
     * true (necessary in order to correctly enforce energy conservation). The
     * Nucleus::applyFinalState() method uses it to determine whether new
     * avatars should be generated for the particle.
     */
    void setOutOfWell() { outOfWell = true; }

    /// \brief Check if the particle is out of its potential well
    G4bool isOutOfWell() const { return outOfWell; }

    void setEmissionTime(G4double t) { emissionTime = t; }
    G4double getEmissionTime() { return emissionTime; };

    /** \brief Transverse component of the position w.r.t. the momentum. */
    ThreeVector getTransversePosition() const {
      return thePosition - getLongitudinalPosition();
    }

    /** \brief Longitudinal component of the position w.r.t. the momentum. */
    ThreeVector getLongitudinalPosition() const {
      return *thePropagationMomentum * (thePosition.dot(*thePropagationMomentum)/thePropagationMomentum->mag2());
    }

    /** \brief Rescale the momentum to match the total energy. */
    const ThreeVector &adjustMomentumFromEnergy();

    /** \brief Recompute the energy to match the momentum. */
    G4double adjustEnergyFromMomentum();

    G4bool isCluster() const {
      return (theType == Composite);
    }

    /// \brief Set the frozen particle momentum
    void setFrozenMomentum(const ThreeVector &momentum) { theFrozenMomentum = momentum; }

    /// \brief Set the frozen particle momentum
    void setFrozenEnergy(const G4double energy) { theFrozenEnergy = energy; }

    /// \brief Get the frozen particle momentum
    ThreeVector getFrozenMomentum() const { return theFrozenMomentum; }

    /// \brief Get the frozen particle momentum
    G4double getFrozenEnergy() const { return theFrozenEnergy; }

    /// \brief Get the propagation velocity of the particle
    ThreeVector getPropagationVelocity() const { return (*thePropagationMomentum)/(*thePropagationEnergy); }

    /** \brief Freeze particle propagation
     *
     * Make the particle use theFrozenMomentum and theFrozenEnergy for
     * propagation. The normal state can be restored by calling the
     * thawPropagation() method.
     */
    void freezePropagation() {
      thePropagationMomentum = &theFrozenMomentum;
      thePropagationEnergy = &theFrozenEnergy;
    }

    /** \brief Unfreeze particle propagation
     *
     * Make the particle use theMomentum and theEnergy for propagation. Call
     * this method to restore the normal propagation if the
     * freezePropagation() method has been called.
     */
    void thawPropagation() {
      thePropagationMomentum = &theMomentum;
      thePropagationEnergy = &theEnergy;
    }

    /** \brief Rotate the particle position and momentum
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotatePositionAndMomentum(const G4double angle, const ThreeVector &axis) {
      rotatePosition(angle, axis);
      rotateMomentum(angle, axis);
    }

    /** \brief Rotate the particle position
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotatePosition(const G4double angle, const ThreeVector &axis) {
      thePosition.rotate(angle, axis);
    }

    /** \brief Rotate the particle momentum
     *
     * \param angle the rotation angle
     * \param axis a unit vector representing the rotation axis
     */
    virtual void rotateMomentum(const G4double angle, const ThreeVector &axis) {
      theMomentum.rotate(angle, axis);
      theFrozenMomentum.rotate(angle, axis);
    }

    std::string print() const {
      std::stringstream ss;
      ss << "Particle (ID = " << ID << ") type = ";
      ss << ParticleTable::getName(theType);
      ss << '\n'
        << "   energy = " << theEnergy << '\n'
        << "   momentum = "
        << theMomentum.print()
        << '\n'
        << "   position = "
        << thePosition.print()
        << '\n';
      return ss.str();
    };

    std::string dump() const {
      std::stringstream ss;
      ss << "(particle " << ID << " ";
      ss << ParticleTable::getName(theType);
      ss << '\n'
        << thePosition.dump()
        << '\n'
        << theMomentum.dump()
        << '\n'
        << theEnergy << ")" << '\n';
      return ss.str();
    };

    long getID() const { return ID; };

    /**
     * Return a NULL pointer
     */
    ParticleList const *getParticles() const {
      INCL_WARN("Particle::getParticles() method was called on a Particle object" << '\n');
      return 0;
    }

    /** \brief Return the reflection momentum
     *
     * The reflection momentum is used by calls to getSurfaceRadius to compute
     * the radius of the sphere where the nucleon moves. It is necessary to
     * introduce fuzzy r-p correlations.
     */
    G4double getReflectionMomentum() const {
      if(rpCorrelated)
        return theMomentum.mag();
      else
        return uncorrelatedMomentum;
    }

    /// \brief Set the uncorrelated momentum
    void setUncorrelatedMomentum(const G4double p) { uncorrelatedMomentum = p; }

    /// \brief Make the particle follow a strict r-p correlation
    void rpCorrelate() { rpCorrelated = true; }

    /// \brief Make the particle not follow a strict r-p correlation
    void rpDecorrelate() { rpCorrelated = false; }

    /// \brief Get the cosine of the angle between position and momentum
    G4double getCosRPAngle() const {
      const G4double norm = thePosition.mag2()*thePropagationMomentum->mag2();
      if(norm>0.)
        return thePosition.dot(*thePropagationMomentum) / std::sqrt(norm);
      else
        return 1.;
    }

  protected:
    G4int theZ, theA;
    ParticipantType theParticipantType;
    G4INCL::ParticleType theType;
    G4double theEnergy;
    G4double *thePropagationEnergy;
    G4double theFrozenEnergy;
    G4INCL::ThreeVector theMomentum;
    G4INCL::ThreeVector *thePropagationMomentum;
    G4INCL::ThreeVector theFrozenMomentum;
    G4INCL::ThreeVector thePosition;
    G4int nCollisions;
    G4int nDecays;
    G4double thePotentialEnergy;
    long ID;

    G4bool rpCorrelated;
    G4double uncorrelatedMomentum;

  private:
    G4double theHelicity;
    G4double emissionTime;
    G4bool outOfWell;

    G4double theMass;
    static G4ThreadLocal long nextID;

    INCL_DECLARE_ALLOCATION_POOL(Particle)
  };
}

#endif /* PARTICLE_HH_ */
