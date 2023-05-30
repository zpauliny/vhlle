class Fluid;

class Particle {
    private:
     int B, Q, S;               // baryon number, charge, strangeness
     double t, x, y, z;      // particle position
     double e, px, py, pz;   // particle momentum
     int ixc, iyc, izc;      // cell coordinates of the particle
     double R;               // sigma of the gaussian
     int nsmoothx, nsmoothy,
      nsmoothz;              // smoothly distribute to +- this many cells
     double gauss_norm;      // normalization of the smearing kernel
     double calculateNorm(Fluid *f, double R);
                             // function to calculate the norm of the kernel
         
    public:
     Particle();
     Particle(Fluid *f, double _R, int _B, int _Q, int _S, double _t, double _x, 
      double _y, double _z, double _e, double _px, double _py, double _pz);
     ~Particle(){};
     double getWeight(int ix, int iy, int iz, Fluid *f, double R);
     double getWeight(double _xdiff, double _ydiff, double _zdiff);


     inline int getIxc(void) { return ixc; }
     inline int getIyc(void) { return iyc; }
     inline int getIzc(void) { return izc; }

     inline int getNsmoothX(void) { return nsmoothx; }
     inline int getNsmoothY(void) { return nsmoothy; }
     inline int getNsmoothZ(void) { return nsmoothz; }

     inline int getB(void) { return B; }
     inline int getQ(void) { return Q; }
     inline int getS(void) { return S; }

     inline double getT(void) { return t; }
     inline double getX(void) { return x; }
     inline double getY(void) { return y; }
     inline double getZ(void) { return z; }
     
     inline double getE(void) { return e; }
     inline double getPx(void) { return px; }
     inline double getPy(void) { return py; }
     inline double getPz(void) { return pz; }
};