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
// $Id: G4Hyperbola.cc,v 1.7 2001/07/11 09:59:45 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Hyperbola.cc
//
// ----------------------------------------------------------------------

#include "G4Hyperbola.hh"
#include "G4CurvePoint.hh"

G4Hyperbola::G4Hyperbola()
{
}

G4Hyperbola::~G4Hyperbola()
{
}

G4Hyperbola::G4Hyperbola(const G4Hyperbola& right)
  : Focus1(right.Focus1), Focus2(right.Focus2),
    ProjFocus1(right.ProjFocus1), ProjFocus2(right.ProjFocus2),
    semiAxis(right.semiAxis), semiImagAxis(right.semiImagAxis),
    ratioAxisImagAxis(right.ratioAxisImagAxis),
    toUnitHyperbola(right.toUnitHyperbola), forTangent(right.forTangent)
{
  pShift    = right.pShift;
  position  = right.position;
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;
}

G4Hyperbola& G4Hyperbola::operator=(const G4Hyperbola& right)
{
  if (&right == this) return *this;

  Focus1 = right.Focus1;
  Focus2 = right.Focus2;
  ProjFocus1   = right.ProjFocus1;
  ProjFocus2   = right.ProjFocus2;
  semiAxis     = right.semiAxis;
  semiImagAxis = right.semiImagAxis;
  ratioAxisImagAxis = right.ratioAxisImagAxis;
  toUnitHyperbola   = right.toUnitHyperbola;
  forTangent        = right.forTangent;
  pShift    = right.pShift;
  position  = right.position;
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;

  return *this;
}

G4Curve* G4Hyperbola::Project(const G4Transform3D& tr)
{
  G4Exception("G4Hyperbola::Project");

  G4Point3D newLocation= tr*position.GetLocation();
  newLocation.setZ(0);
  G4double axisZ= (tr*position.GetPZ()).unit().z();

  if (abs(axisZ)<kAngTolerance) 
  {
    return 0;
  }

  G4Vector3D newAxis(0, 0, axisZ>0? +1: -1);

  // get the parameter of an endpoint of an axis
  // (this is a point the distance of which from the center is extreme)
  G4Vector3D xPrime = tr*position.GetPX();
  xPrime.setZ(0);
  
  G4Vector3D yPrime = tr*position.GetPY();
  yPrime.setZ(0);
  
  G4Vector3D a = G4Vector3D( semiAxis*xPrime );
  G4Vector3D b = G4Vector3D( semiImagAxis*yPrime );

  G4double xval = -2*a*b/(a.mag2()+b.mag2());
  
  G4double u= (0.5*log((1+xval)/(1-xval)))/2;  // atanh(xval)/2

  // get the coordinate axis directions and the semiaxis lengths
  G4Vector3D sAxis= G4Vector3D( a*cosh(u)+b*sinh(u) );
  
  //!!!!!!!!!!!!
  G4Vector3D sImagAxis= G4Vector3D( a*cosh(u+pi/2)+b*sinh(u+pi/2) );
  
  //!!!!!!!!!!!!
  G4double   newSemiAxis     = sAxis.mag();
  G4double   newSemiImagAxis = sImagAxis.mag();
  G4Vector3D newRefDirection = sAxis;

  // create the new hyperbola
  G4Axis2Placement3D newPosition;
  newPosition.Init(newRefDirection, newAxis, newLocation);

  G4Hyperbola* r= new G4Hyperbola;
  r->Init(newPosition, newSemiAxis, newSemiImagAxis);

  // introduce the shift in the parametrization
  // maybe the Sign must be changed?
  r->SetPShift(u);
  
  // set the bounds when necessary
  if (IsBounded()) 
    r->SetBounds(GetPStart(), GetPEnd());
  
  return r;
}


void G4Hyperbola::InitBounded()
{
  // the bbox must include the start and endpoints as well as the
  // extreme points if they lie on the curve
  bBox.Init(GetStart(), GetEnd());

  // the parameter values
  // belonging to the points with an extreme x, y and z coordinate
  for (G4int i=0; i<3; i++) 
  {
    G4double x_i= position.GetPX()(i);
    
    if (abs(x_i) <= kAngTolerance) 
    {
      G4double tanhu= - (semiImagAxis*position.GetPY()(i)) / (semiAxis*x_i);
      
      if (abs(tanhu)<=1) 
      {
        G4double u= 0.5*log((1+tanhu)/(1-tanhu));  // atanh(tanhu)
	if (IsPOn(u))
	  bBox.Extend(GetPoint(u));
      }
    }
  }
}

G4bool G4Hyperbola::Tangent(G4CurvePoint& cp, G4Vector3D& v)
{
  // The tangent is computed from the 3D point representation
  // for all conics. An alternaive implementation (based on
  // the parametric point) might be worthwhile adding
  // for efficiency.
  
  const G4Axis2Placement3D& pos= *(GetPosition());
  G4Point3D p= pos.GetToPlacementCoordinates() * cp.GetPoint();
  
  v= forTangent*p.y()*pos.GetPX() + p.x()*pos.GetPY();
  
  return true;
}