#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "fld.h"
#include "particle.h"

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
void Particle::energyLoss(double energyLoss0, double dt){
   double p = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
   const double m = pow(e,2) - pow(m,2);

   double phi = std::asin(py/px);
   double theta = std::acos(pz/p);

   double v  = sqrt(p/e);
   double dx =  v * cos(phi) * sin(theta) * dt;
   double dy =  v * sin(phi) * sin(theta) * dt;
   double dz =  v * cos(theta) * dt;
   x += dx;
   y += dy;
   z += dz;
   double dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
   double de = energyLoss0 * dr;
   e -= de;
   p = sqrt(pow(e,2) - pow(m,2));
}

