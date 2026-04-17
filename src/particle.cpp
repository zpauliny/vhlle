#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "fld.h"
#include "particle.h"
#include "cll.h"
#include "eos.h"
#include "rmn.h"

using namespace std;

Particle::Particle(){
   B = 0;
   Q = 0;
   S = 0;

   t = 0;
   x = 0;
   y = 0;
   z = 0;
   e = 0.0;
   px = 0.0;
   py = 0.0;
   pz = 0.0;

   eventNo = 0;
   R = 0.0;
}

Particle::Particle(Fluid *f, double _R, int _B, int _Q, int _S, double _t, double _x,
      double _y, double _z, double _e, double _px, double _py, double _pz, int _pdg, int _eventNo){
    B = _B;
    Q = _Q;
    S = _S;
    pdg = _pdg;

    t = _t;
    x = _x;
    y = _y;
    z = _z;
    e = _e;
    px = _px;
    py = _py;
    pz = _pz;

    eventNo = _eventNo;
    R = _R;
    scale = 1.0;

    const double xmin = f->getX(0);
    const double ymin = f->getY(0);
    const double zmin = f->getZ(0);
    const double dx = f->getDx();
    const double dy = f->getDy();
    const double dz = f->getDz();

    ixc = round((x - xmin) / dx);
    iyc = round((y - ymin) / dy);
    izc = round((z - zmin) / dz);

    const double range = 2.0;

    nsmoothx = static_cast<int>(range * R / dx);
    nsmoothy = static_cast<int>(range * R / dy);
    nsmoothz = static_cast<int>(range * R / dz);

    gauss_norm = calculateNorm(f,R);
}



double Particle::calculateNorm(Fluid *f, double R) {
   const int nx = f->getNX();
   const int ny = f->getNY();
   const int nz = f->getNZ();
   const double xmin = f->getX(0);
   const double ymin = f->getY(0);
   const double zmin = f->getZ(0);
   const double dx = f->getDx();
   const double dy = f->getDy();
   const double dz = f->getDz();

   const double m = sqrt(e * e - px * px - py * py - pz * pz);

   const double ux = px / m;
   const double uy = py / m;
   const double uz = pz / m;

   double norm = 0;

   for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
    for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
     for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
      if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {

       const double xdiff = x - (xmin + ix * dx);
       const double ydiff = y - (ymin + iy * dy);
       const double zdiff = z - (zmin + iz * dz);

       const double rr = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
       const double ru = xdiff * ux + ydiff * uy + zdiff * uz;
       norm +=
           exp((-rr - ru * ru) / R / R);  // this term may become big, >800 !
     }
   return norm;
}

double Particle::getWeight(double _xdiff, double _ydiff, double _zdiff) {
   const double m = sqrt(e * e - px * px - py * py - pz * pz);
   const double ux = px / m;
   const double uy = py / m;
   const double uz = pz / m;

   const double rr = _xdiff * _xdiff + _ydiff * _ydiff + _zdiff * _zdiff;
   const double ru = _xdiff * ux + _ydiff * uy + _zdiff * uz;

   return 1. / gauss_norm *
                    exp((-rr - ru * ru) / R / R);
}

double Particle::getWeight(int ix, int iy, int iz, Fluid *f,
        double R) {
   const double xmin = f->getX(0);
   const double ymin = f->getY(0);
   const double zmin = f->getZ(0);
   const double dx = f->getDx();
   const double dy = f->getDy();
   const double dz = f->getDz();

   const double xdiff = x - (xmin + ix * dx);
   const double ydiff = y - (ymin + iy * dy);
   const double zdiff = z - (zmin + iz * dz);

   return getWeight(xdiff, ydiff, zdiff);
}

void Particle::energyLoss(double energyLoss0, double dt, Fluid* f, EoS* eos, double* dp){
   double p = sqrt(px*px + py*py + pz*pz);
   double tau = t;
   double m2 = max(0.0, e*e - px*px - py*py - tau*tau * pz*pz);

   Cell* myCell = f->getCell(ixc, iyc, izc);
   double Q[7];
   myCell->getQ(Q);
   // cell velocity: vx, vy, vz
   double e_cell, p_cell, nb, nq, ns, vx, vy, vz;
   transformPV(eos, Q, e_cell, p_cell, nb, nq, ns, vx, vy, vz);

   double p_eta = pz;
   double p_tau = e; 
   double vtau = 1.0 / sqrt(1.0 - vx*vx - vy*vy - vz*vz);
   double E_lrf = max( p_tau*vtau - px*vx - py*vy - tau*tau * p_eta*vz, 0.0);
   double p_lrf = sqrt(max((E_lrf*E_lrf - m2 ), 1.0e-12));
   
   dp[0] = -energyLoss0 * dt * (p_tau - E_lrf*vtau) / p_lrf;
   dp[1] = -energyLoss0 * dt * (px - E_lrf*vx) / p_lrf;
   dp[2] = -energyLoss0 * dt * (py - E_lrf*vy) / p_lrf;
   dp[3] = -energyLoss0 * dt * (p_eta - E_lrf*vz) / p_lrf;

   double mass = sqrt(m2);
   updateMomentum(dp, tau, mass);
   updatePosition(dt);
}

void Particle::updateMomentum(double* dp, double tau, double mass){
   e += dp[0];
   px += dp[1];
   py += dp[2];
   pz += dp[3];

   double spatial = px*px + py*py + tau*tau * pz*pz;
   e = sqrt(mass * mass + spatial);
}

void Particle::updatePosition(double dt){
   t += dt;
   x += dt * px / e;
   y += dt * py / e;
   z += dt * pz / e;
}