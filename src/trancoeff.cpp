#include <iostream>
#include <iomanip>
#include "eos.h"
#include "trancoeff.h"
#include "inc.h"

TransportCoeff::TransportCoeff(double _etaS, double _zetaS, int _zetaSparam, EoS *_eos, int _etaSparam, double _ah, double _al, double _aRho, double _T0, double _etaSMin, double _eEtaSMin) {
 etaS0 = _etaS;
 zetaS0 = _zetaS;
 zetaSparam = _zetaSparam; //0 - basic ,1 ,2 - arxiv:1910.12930, 3 arxiv:2103.09848
 eos = _eos;
 etaSparam = _etaSparam;
 ah = _ah;
 al = _al;
 aRho =_aRho;
 etaSMin = _etaSMin;
 eEtaSMin = _eEtaSMin;
 T0=_T0;
}

void TransportCoeff::printZetaT()
{
 std::cout << "------zeta(T):\n";
 for(double e=0.1; e<3.0; e+=0.1){
  double T, mub, muq, mus, p;
  eos->eos(e, 0., 0., 0., T, mub, muq, mus, p);
  const double s = eos->s(e, 0, 0, 0);
  std::cout << std::setw(14) << T << std::setw(14) << zeta(e, s, T)/s << std::endl;
 }
 std::cout << "---------------:\n";
}

double TransportCoeff::zeta(double e, double s, double T)
{
 double T_p=0.180;

 double T_peak=0.165;
 double T_width=0.010;
 double B_norm=0.24;
 double B_width=1.5;
 double T_peak2=0.160;
 double B_norm2=0.13;
 double B1=0.01;
 double B2=0.12;

 if(zetaSparam==0)
    return (zetaS0 * (1. / 3. - eos->cs2(e)) / (exp((0.16 - T) / 0.001) + 1.)) * s;
 else if(zetaSparam==1)
 {
    if(T<0.180)
       return (0.03+(0.08*exp(((T/T_p)-1.)/(0.0025)))+(0.22*exp(((T/T_p)-1)/(0.0022)))) * s;
    else if(T>=0.180 && T<0.200)
       return (27.55*(T/T_p)-13.45-(13.77*(T/T_p)*(T/T_p))) * s;
    else if(T>=0.200)
       return (0.001+(0.9*exp(-((T/T_p)-1.)/(0.0025)))+(0.25*exp(-((T/T_p)-1.)/(0.13)))) * s;
 }
 else if(zetaSparam==2)
 {
    if(T>T_peak)
       return (B_norm*((B_width*B_width)/((((T/T_peak)-1.)*((T/T_peak)-1.))+(B_width*B_width)))) * s;
    else if(T<=T_peak)
       return (B_norm*(exp(-((T-T_peak)/T_width)*((T-T_peak)/T_width)))) * s;
 }
 else if(zetaSparam==3)
 {
    if(T<T_peak2)
       return (B_norm2*exp(-((T-T_peak2)*(T-T_peak2)/(B1*B1)))) * s;
    else if(T>=T_peak2)
       return (B_norm2*exp(-((T-T_peak2)*(T-T_peak2)/(B2*B2)))) * s;
 }
}

void TransportCoeff::getEta(double e, double p, double nb, double s, double T, double &_eta, double &_zeta) {
 _eta = eta(e, p, nb, s, T);
 _zeta = zeta(e,s,T);
}

double TransportCoeff::eta(double e, double p, double nb, double s, double T)
{
  if (etaSparam == 0){
      return etaS0*s;
  }
  else if (etaSparam == 1){
      return (etaSMin +  ((T>T0) ? ah*(T-T0) :  al*(T-T0))) * s;
  }
  else if (etaSparam == 2){
      return (std::max(0.0, etaSMin + ((e>eEtaSMin) ? ( (ah*(e-eEtaSMin)+aRho*nb) ): al*(e-eEtaSMin)+aRho*nb))) * s;
  }
}

void TransportCoeff::getTau(double e, double p, double nb, double s, double T, double &_taupi, double &_tauPi) {
 if (T > 0.) {
  _taupi = 5. / 5.068 * eta(e, p, nb, s, T) / (s * T);
  _tauPi = 6.0 / 5.068 * zeta(e, s, T) / (s * T);
 } else {
  _taupi = _tauPi = 0.;
 }
}
