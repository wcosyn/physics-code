#include "FourVector.h"
#include <TLorentzVector.h>


//_____________________________________________________________________
FourVector<double> ToFourVector(const TLorentzVector& v)
{
  return FourVector<double>(v.E(),v.X(),v.Y(),v.Z());
}

//_____________________________________________________________________
TLorentzVector ToLorentzVector(const FourVector<double>& v)
{
  return TLorentzVector(v[1],v[2],v[3],v[0]);
}
