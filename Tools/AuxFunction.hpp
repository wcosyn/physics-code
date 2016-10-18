#ifndef GUAuxFunctions_H
#define GUAuxFunctions_H
 
// add the sincos function on MAC because sincos is not part of math.h
#ifdef __APPLE__ // possibly other conditions
inline void sincos(double x, double *s, double *c){
  __sincos(x,s,c);
}
#endif
 
#endif