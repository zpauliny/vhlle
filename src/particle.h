class Fluid;

class Particle {
    private:
     int B, Q, S;               // baryon number, charge, strangeness
     double tau, x, y, eta;      // particle position
     double mt, px, py, rap;   // particle momentum
     int ixc, iyc, izc;      // cell coordinates of the particle
     double R;               // sigma of the gaussian
     int nsmoothx, nsmoothy,
      nsmoothz;              // smoothly distribute to +- this many cells
     double gauss_norm;      // normalization of the smearing kernel
     double calculateNorm(Fluid *f, double R);
                             // function to calculate the norm of the kernel
     double x_traveled;      // to keep track how far has the "jet" gotten
     double e_init;          // initial energy of the particle
     bool isDissolved;       // flag when all energy has been pumped into hydro
     double deT, deX, deY, deZ;
         
    public:
     Particle();
     Particle(Fluid *f, double _R, int _B, int _Q, int _S, double _tau, double _x, 
      double _y, double _eta, double _mt, double _px, double _py, double _rap);
     ~Particle(){};
     double getWeight(int ix, int iy, int iz, Fluid *f, double R);
     double getWeight(double _xdiff, double _ydiff, double _zdiff);
     double getE(int iz, Fluid *f);
     double getE();
     double getPEta(int iz, Fluid *f);
     double getPEta();
     double getDE(double _dtau);
     void updateParticle(double _dE, double _dtau);
     void setMassless();
     void setR(Fluid *f, double _R); // to set R to a value appropriate for jets
     double getDEtau(double _eta);
     double getDEeta(double _eta);

     inline int getIxc(void) { return ixc; }
     inline int getIyc(void) { return iyc; }
     inline int getIzc(void) { return izc; }

     inline int getNsmoothX(void) { return nsmoothx; }
     inline int getNsmoothY(void) { return nsmoothy; }
     inline int getNsmoothZ(void) { return nsmoothz; }

     inline int getB(void) { return B; }
     inline int getQ(void) { return Q; }
     inline int getS(void) { return S; }

     inline double getT(void) { return tau; }
     inline double getX(void) { return x; }
     inline double getY(void) { return y; }
     inline double getZ(void) { return eta; }
     
     inline double getMt(void) { return mt; }
     inline double getPx(void) { return px; }
     inline double getPy(void) { return py; }
     inline double getRap(void) { return rap; }
     
     inline bool isFluidized(void) { return isDissolved; }

     inline double getDEt(void) { return deT; }
     inline double getDEx(void) { return deX; }
     inline double getDEy(void) { return deY; }
     inline double getDEz(void) { return deZ; }
};