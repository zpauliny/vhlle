#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <math.h>

#include "fld.h"
#include "particle.h"
#include "icPartSMASH.h"

using namespace std;

Particle::Particle(){
   B = 0;
   Q = 0;
   S = 0;

   tau = 0;
   x = 0;
   y = 0;
   eta = 0;
   mt = 0.0;
   px = 0.0;
   py = 0.0;
   rap = 0.0; 
   R = 0.0;

   x_traveled = 0.0;
   e_init = 0.0;
   isDissolved = false;

   deT = 0.0;
   deX = 0.0;
   deY = 0.0;
   deZ = 0.0;
}

Particle::Particle(Fluid *f, double _R, int _B, int _Q, int _S, double _tau, double _x, 
      double _y, double _eta, double _mt, double _px, double _py, double _rap){
    B = _B;
    Q = _Q;
    S = _S;

    tau = _tau;
    x = _x;
    y = _y;
    eta = _eta;
    mt = _mt;
    px = _px;
    py = _py;
    rap = _rap; 

    R = _R;
    x_traveled = 0.0;
    e_init = mt * cosh(rap);
    isDissolved = false;

    deT = 0.0;
    deX = 0.0;
    deY = 0.0;
    deZ = 0.0;
    
    const double xmin = f->getX(0);
    const double xmax = f->getX(f->getNX() - 1);
    const double ymin = f->getY(0);
    const double ymax = f->getY(f->getNY() - 1);
    const double zmin = f->getZ(0);
    const double zmax = f->getZ(f->getNZ() - 1);
    const double dx = f->getDx();
    const double dy = f->getDy();
    const double dz = f->getDz();
        
    ixc = round((x - xmin) / dx);
    iyc = round((y - ymin) / dy);
    izc = round((eta - zmin) / dz);
    
      
    const double range = 2.0;
    
    nsmoothx = static_cast<int>(range * R / dx);
    nsmoothy = static_cast<int>(range * R / dy);
    nsmoothz = static_cast<int>(range * R / dz);
    
    gauss_norm = calculateNorm(f,R);
    //gauss_norm = 1.0;
}

void Particle::setR(Fluid *f, double _R) {
   R = _R;

   const double xmin = f->getX(0);
   const double xmax = f->getX(f->getNX() - 1);
   const double ymin = f->getY(0);
   const double ymax = f->getY(f->getNY() - 1);
   const double zmin = f->getZ(0);
   const double zmax = f->getZ(f->getNZ() - 1);
   const double dx = f->getDx();
   const double dy = f->getDy();
   const double dz = f->getDz();
       
   ixc = round((x - xmin) / dx);
   iyc = round((y - ymin) / dy);
   izc = round((eta - zmin) / dz);
   
     
   const double range = 2.0;
   
   nsmoothx = static_cast<int>(range * R / dx);
   nsmoothy = static_cast<int>(range * R / dy);
   nsmoothz = static_cast<int>(range * R / dz);
   gauss_norm = calculateNorm(f,R);
}

double Particle::calculateNorm(Fluid *f, double R) {
// ATTENTION
// now in Milne coordinates
   int ix, iy, iz;
   const int nx = f->getNX();
   const int ny = f->getNY();
   const int nz = f->getNZ();
   const double xmin = f->getX(0);
   const double xmax = f->getX(nx - 1);
   const double ymin = f->getY(0);
   const double ymax = f->getY(ny - 1);
   const double zmin = f->getZ(0);
   const double zmax = f->getZ(nz - 1);
   const double dx = f->getDx();
   const double dy = f->getDy();
   const double dz = f->getDz();

   double norm;
   const double gammaz = cosh(rap - eta);
   // hard-coded here for 200 GeV, TODO
   const double Rgx = 1.0;
   const double Rgy = 1.0;
   const double Rgz = 1.1;
   
   for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
    for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
     for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
      if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
       
       const double xdiff = x - (xmin + ix * dx);
       const double ydiff = y - (ymin + iy * dy);
       const double zdiff = eta - (zmin + iz * dz);
       spatialVector rdiff {xdiff, ydiff, zdiff};
       velocityVector velocity = velocityHyperbolic(mt, px, py, rap, eta, zdiff, tau);
       //norm += smoothingKernelInvariant(rdiff, velocity, R, tau);
       norm += smoothingKernel(rdiff, gammaz, tau,  Rgx, Rgy, Rgz);
       
     }
    return norm;
}

double Particle::getWeight(double _xdiff, double _ydiff, double _zdiff) {

   const double gammaz = cosh(rap - eta);
   // hard-coded here for 200 GeV, TODO
   const double Rgx = 1.0;
   const double Rgy = 1.0;
   const double Rgz = 1.1;

   spatialVector rdiff {_xdiff, _ydiff, _zdiff};
   velocityVector velocity {velocityHyperbolic(mt, px, py, rap, eta, _zdiff, tau)};
   //double weight = 1. / gauss_norm * smoothingKernelInvariant(rdiff, velocity, R, tau);
   double weight = 1. / gauss_norm * smoothingKernel(rdiff, gammaz, tau, Rgx, Rgy, Rgz);

   return weight;
}

double Particle::getWeight(int ix, int iy, int iz, Fluid *f, 
        double R) {

   const int nx = f->getNX();
   const int ny = f->getNY();
   const int nz = f->getNZ();
   const double xmin = f->getX(0);
   const double xmax = f->getX(nx - 1);
   const double ymin = f->getY(0);
   const double ymax = f->getY(ny - 1);
   const double zmin = f->getZ(0);
   const double zmax = f->getZ(nz - 1);
   const double dx = f->getDx();
   const double dy = f->getDy();
   const double dz = f->getDz();

   const double xdiff = x - (xmin + ix * dx);
   const double ydiff = y - (ymin + iy * dy);
   const double zdiff = eta - (zmin + iz * dz);

   double weight = getWeight(xdiff, ydiff, zdiff);
   return weight;
}

double Particle::getE(int iz, Fluid *f) {

   const int nz = f->getNZ();
   const double zmin = f->getZ(0);
   const double dz = f->getDz();

   const double zdiff {eta - (zmin + iz * dz)};
   double E = mt * cosh(rap - eta + zdiff);
   return E;
}

double Particle::getE() {
   
   return (mt * cosh(rap - eta));
}

double Particle::getPEta(int iz, Fluid *f) {
   
   const int nz {f->getNZ()};
   const double zmin {f->getZ(0)};
   const double dz {f->getDz()};

   const double zdiff {eta - (zmin + iz * dz)};
   double pEta = mt * sinh(rap - eta + zdiff);
   return pEta;
}

double Particle::getPEta(){
   
   return (mt * sinh(rap - eta));
}

// update particle at the end of the timestep
void Particle::updateParticle(double _dE, double _dtau) {

   // store old values
   double px_old = px;
   double py_old = py; 

   // change in energy + dE distribution
   double e = mt * cosh(rap);
   double e_new = e +_dE;
   
   double pz = mt * sinh(rap);
   double p = sqrt(px * px + py * py + pz * pz);
   double vx = px / p;
   double vy = py / p;
   double vz = pz / p;
   
   deT = -1.0 * _dE;
   deX = -1.0 * vx * _dE;
   deY = -1.0 * vy * _dE;
   deZ = -1.0 * vy * _dE;

   // keep in mind, massless: pt = mt
   double alpha = abs(mt / pz);
   pz = sqrt(e_new * e_new / (alpha * alpha + 1));
   mt = alpha * pz;
   double beta = abs(px / py);
   py = sqrt(mt * mt / (beta * beta + 1));
   px = beta * py;
   rap = 0.5 * log((e_new + pz) / (e_new - pz));

   // change in position
   // use Cartesian coordinates then transform
   p = sqrt(px * px + py * py + pz * pz);
   vx = px / p;
   vy = py / p;
   vz = pz / p;
   
   double dt = _dtau * cosh(eta);
   double t = tau * cosh(eta) + dt;
   double z = tau * sinh(eta);
   x += vx * dt;
   y += vy * dt;
   z += vz * dt;
   tau += _dtau;
   eta = 0.5 * log((t + z) / (t - z));
}

double Particle::getDE(double _dtau, double _deta) {
   // formulae acc. to Pablos 2202.03414
   
   // choosing constant x_stop
   double x_stop = 3.0;

   // logic here: velocity is c, dx = dt = 
   //double dx = _dtau * cosh(eta) + tau * sinh(eta) * _deta;
   double dx = _dtau * cosh(eta);
   x_traveled += dx;
   double xratio = x_traveled * x_traveled / x_stop / x_stop;
   double de_dx = -4.0 / M_PI * e_init * xratio / sqrt(1.0 - xratio);
   double dE = de_dx * dx;
   
   // correct for the last deposition od energy
   double e = mt * cosh(rap);
   double e_new = e + dE;
   if (e_new <= 0.) {
      dE = -1.0 * e;
      isDissolved = true;
      cout << "Jet dissolved in hydro \n";
   }
   return dE;
}

void Particle::setMassless() {
   // particle's mass is set to 0
   // momentum is conserved
   double m = sqrt(mt * mt - px * px - py * py);
   double pz = mt * sinh(rap);
   double e_new = sqrt(px * px + py * py + pz * pz);
   mt = e_new / cosh(rap);
   double alpha = abs(px / py);
   py = sqrt(mt * mt / (alpha * alpha + 1));
   px = alpha * py;
}

double Particle::getDEtau(double _eta) {
   
   double dEtau = cosh(_eta) * deT - sinh(_eta) * deZ;
   return dEtau;
}

double Particle::getDEeta(double _eta) {

   double dEeta = - 1.0 / tau * sinh(_eta) * deT + 1.0 / tau * cosh(_eta) * deZ;
   return dEeta;
}