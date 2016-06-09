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
// $Id: G4VLowEnergyDiscretePhotonProcess.hh,v 1.2 2005/05/20 15:19:18 pia Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// --------------------------------------------------------------
//
// File name:     G4VLowEnergyDiscretePhotonProcess.hh
//
// Author:        Capra Riccardo
//
// Creation date: May 2005
//
// History:
// -----------
// 02 May 2005  R. Capra         1st implementation
// 20 May 2005  MGP              Modified protected member functions to return a const pointer
//
//----------------------------------------------------------------



#ifndef   G4VLowEnergyDiscretePhotonProcess_hh
#define  G4VLowEnergyDiscretePhotonProcess_hh
 
// Base class
#include "G4VDiscreteProcess.hh"

// Forward declaration
class G4String;
class G4ParticleDefinition;
class G4DynamicParticle;
class G4Track;
class G4VCrossSectionHandler;
class G4VEMDataSet;
class G4VDataSetAlgorithm;
 
// G4VLowEnergyDiscretePhotonProcess
// A common class for Rayleigh and Compton processes
class G4VLowEnergyDiscretePhotonProcess : public G4VDiscreteProcess
{
public:
  //   Class constructor
  //
  //   Creates crossSectionHandler and scatterFunctionData
  //   scatterFunctionData are loaded from the scatterFile file, and are interpolated with scatterInterpolation algorithm
  //   processName The name of the process
  //   aCrossSectionFileName The name of the cross-section data file
  //   aScatterFileName The name of the scatter function data
  //   aScatterInterpolation The interpolation algorithm
  //   aLowEnergyLimit The lower energy limit of validity against data
  //   aHighEnergyLimit The higher energy limit of validity against data
  G4VLowEnergyDiscretePhotonProcess(const G4String& processName, 
				    const G4String& aCrossSectionFileName, 
				    const G4String& aScatterFileName, 
				    G4VDataSetAlgorithm* aScatterInterpolation, 
				    G4double aLowEnergyLimit, 
				    G4double aHighEnergyLimit);

  //   Class destructor
  //
  //        Deletes crossSectionHandler, meanFreePathTable and scatterFunctionData
  virtual ~G4VLowEnergyDiscretePhotonProcess(void);



  //   Checks if the process is applicable to the particle
  //
  //   For processes inheriting from this class the only applicable particle is the photon
  //   particleDefinition Is the particle to be checked for the applicability of the process
  //   true only if the particle is a photon
  virtual G4bool IsApplicable(const G4ParticleDefinition& particleDefinition);

  //   Updates the crossSectionHandler and builds the meanFreePathTable from it
  //
  //        crossSectionHandler data is loaded from crossSectionFile file
  //   photon particle is always a photon
  virtual void BuildPhysicsTable(const G4ParticleDefinition& photon);
   
  //   Calls GetMeanFreePath (for testing purpose only)
  //   aTrack the particle momentum for which the mean free path must be evaluates (a photon)
  //   previousStepSize the size of the prevous step (not used by this implementation)
  //   condition the consition to be updated (not used by this implementation)
  //   The mean free path evaluated by GetMeanFreePath
  inline G4double DumpMeanFreePath(const G4Track& aTrack, 
				   G4double previousStepSize, 
				   G4ForceCondition* condition) 
  {return GetMeanFreePath(aTrack, previousStepSize, condition);}
   
   
   
  //   Returns the low energy limit of the process
  //   The low energy limit of the process
  inline G4double GetLowEnergyLimit(void) const {return lowEnergyLimit;}
   
  //   Returns the high energy limit of the process
  //   The high energy limit of the process
  inline G4double GetHighEnergyLimit(void) const {return highEnergyLimit;}



protected:
  //   Evaluates the process mean free path
  //
  //    Mean free path evaluation is based on meanFreePathTable generated by the crossSectionHandler with data taken from the crossSectionFile.
  //   The method uses lowEnergyLimit and highEnergyLimit
  //   aTrack the particle momentum for which the mean free path must be evaluates (a photon)
  //   previousStepSize the size of the prevous step (not used by this implementation)
  //   condition the consition to be updated (not used by this implementation)
  //   The mean free path for the process
  virtual G4double GetMeanFreePath(const G4Track &aTrack, 
				   G4double previousStepSize, 
				   G4ForceCondition* condition);
   
  //   Returns the cross-section handler
  //   The cross-section handler
  inline const G4VCrossSectionHandler* GetCrossSectionHandler(void) const {return crossSectionHandler;}
   
  //   Returns the mean free path table
  //   The mean free path table
  inline const G4VEMDataSet* GetMeanFreePathTable(void) const {return meanFreePathTable;}
   
  //   Returns the scatter function data
  //   The scatter function data
  inline const G4VEMDataSet* GetScatterFunctionData(void) const {return scatterFunctionData;}



  //   Verifies if the polarization vector is orthonormal to the photon direction
  //
  //        This method is used by polarized processes. It returns always a vector orthonormal to the photon direction.
  //        When G4DynamicParticle::GetPolarization() is well defined the returned vector matches this vector.
  //        When G4DynamicParticle::GetPolarization() is sensibly different from a null-vector, the returned vector is
  //        orthonormal to the photon direction and on the plane made with G4DynamicParticle::GetPolarization()
  //        When G4DynamicParticle::GetPolarization() is almost a null-vector, the returned vector is a random
  //        vector orthonormal to the photon direction
  //   photon The incoming photon
  //  s Always a vector ortogonal to the photon direction.
  G4ThreeVector GetPhotonPolarization(const G4DynamicParticle& photon);



private:

  //   Hides copy constructor
  G4VLowEnergyDiscretePhotonProcess(const G4VLowEnergyDiscretePhotonProcess& );

  //   Hides assignment operator
  G4VLowEnergyDiscretePhotonProcess& operator=(const G4VLowEnergyDiscretePhotonProcess& );

  //   Low energy limit
  G4double lowEnergyLimit;

  //   High energy limit
  G4double highEnergyLimit;

  //   Cross-section file name
  G4String crossSectionFileName;

  //   Cross-section handler
  G4VCrossSectionHandler* crossSectionHandler;

  //   The mean free path table based on cross-section data
  G4VEMDataSet* meanFreePathTable;

  //   Scatter function table
  G4VEMDataSet* scatterFunctionData;
};

#endif /* G4VLowEnergyDiscretePhotonProcess_hh */