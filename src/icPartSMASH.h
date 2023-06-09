#include <vector>
#include <queue>

class Fluid;
class EoS;
class Particle;

class IcPartSMASH {
private:
 int nx, ny, nz, nevents;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***T00, ***T0x, ***T0y, ***T0z, ***QB, ***QE, ***QS;
 // auxiliary particle values for reading from file
 double Tau_val, X_val, Y_val, Eta_val, Mt_val, Px_val, Py_val, Rap_val;
 int Id_val, Charge_val, Baryon_val, Strangeness_val;
  // auxiliary particle arrays
 std::vector<double> Tau, X, Y, Eta, Mt, Px, Py, Rap;
 std::vector<int> Id, Charge;

 double tau0;
 double Rgx, Rgy, Rgz;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;
 void makeSmoothTable(int npart);
 bool isKernelInvariant;

 // To properly add baryon contributions to net-baryon number
 // PDG IDs corresponding to baryons included in SMASH-1.8
 int PDG_Codes_Baryons[392] =
              {1112, 1114, 1116, 1118, 1212, 1214, 1216, 1218, 2112, 2114, 2116,
               2118, 2122, 2124, 2126, 2128, 2212, 2214, 2216, 2218, 2222, 2224,
               2226, 2228, 3112, 3114, 3116, 3118, 3122, 3124, 3126, 3128, 3212,
               3214, 3216, 3218, 3222, 3224, 3226, 3228, 3312, 3314, 3322, 3324,
               3334, 4112, 4114, 4122, 4132, 4212, 4214, 4222, 4224, 4232, 4312,
               4314, 4322, 4324, 4332, 4334, 4412, 4414, 4422, 4424, 4432, 4434,
               4444, 5112, 5114, 5122, 5132, 5142, 5212, 5214, 5222, 5224, 5232,
               5242, 5312, 5314, 5322, 5324, 5332, 5334, 5342, 5412, 5414, 5422,
               5424, 5432, 5434, 5442, 5444, 5512, 5514, 5522, 5524, 5532, 5534,
               5542, 5544, 5554, 11112, 11114, 11116, 11212, 11216, 12112, 12114,
               12116, 12122, 12126, 12212, 12214, 12216, 12222, 12224, 12226,
               13112, 13114, 13116, 13122, 13124, 13126, 13212, 13214, 13216,
               13222, 13224, 13226, 13314, 13324, 14122, 21112, 21114, 21212,
               21214, 22112, 22114, 22122, 22124, 22212, 22214, 22222, 22224,
               23112, 23114, 23122, 23124, 23126, 23212, 23214, 23222, 23224,
               31214, 32112, 32124, 32212, 33122, 42112, 42212, 43122, 53122,
               103316, 103326, 203312, 203316, 203322, 203326, 203338, 9902114,
               9902118, 9902214, 9902218, 9903118, 19903129, 9903218, 9903228,
               9912114, 9912214, 9922114, 9922116, 19922119, 9922214, 9922216,
               19922219, 9932114, 19932119, 9932214, 19932219, 9952112, 9952212,
               9962112, 9962212, 9972112, 9972212, -1112, -1114, -1116, -1118,
               -1212, -1214, -1216, -1218, -2112, -2114, -2116, -2118, -2122,
               -2124, -2126, -2128, -2212, -2214, -2216, -2218, -2222, -2224,
               -2226, -2228, -3112, -3114, -3116, -3118, -3122, -3124, -3126,
               -3128, -3212, -3214, -3216, -3218, -3222, -3224, -3226, -3228,
               -3312, -3314, -3322, -3324, -3334, -4112, -4114, -4122, -4132,
               -4212, -4214, -4222, -4224, -4232, -4312, -4314, -4322, -4324,
               -4332, -4334, -4412, -4414, -4422, -4424, -4432, -4434, -4444,
               -5112, -5114, -5122, -5132, -5142, -5212, -5214, -5222, -5224,
               -5232, -5242, -5312, -5314, -5322, -5324, -5332, -5334, -5342,
               -5412, -5414, -5422, -5424, -5432, -5434, -5442, -5444, -5512,
               -5514, -5522, -5524, -5532, -5534, -5542, -5544, -5554, -11112,
               -11114, -11116, -11212, -11216, -12112, -12114, -12116, -12122,
               -12126, -12212, -12214, -12216, -12222, -12224, -12226, -13112,
               -13114, -13116, -13122, -13124, -13126, -13212, -13214, -13216,
               -13222, -13224, -13226, -13314, -13324, -14122, -21112, -21114,
               -21212, -21214, -22112, -22114, -22122, -22124, -22212, -22214,
               -22222, -22224, -23112, -23114, -23122, -23124, -23126, -23212,
               -23214, -23222, -23224, -31214, -32112, -32124, -32212, -33122,
               -42112, -42212, -43122, -53122, -103316, -103326, -203312,
               -203316, -203322, -203326, -203338, -9902114, -9902118, -9902214,
               -9902218, -9903118, -19903129, -9903218, -9903228, -9912114,
               -9912214, -9922114, -9922116, -19922119, -9922214, -9922216,
               -19922219, -9932114, -19932119, -9932214, -19932219, -9952112,
               -9952212, -9962112, -9962212, -9972112, -9972212};

public:
 IcPartSMASH(Fluid *f, const char *filename, double _Rgt, double _Rgz, double tau0, int _smoothingType);
 IcPartSMASH(Fluid *f, const char *filename, double _Rgt, double _Rgz, std::queue<Particle>* particles);
 ~IcPartSMASH();
 void setIC(Fluid *f, EoS *eos);
 void setIC(Fluid *f, EoS *eos, std::queue<Particle>* particles, double* ctime);
};

struct spatialVector {
    double xdiff;
    double ydiff;
    double zdiff;
};

struct velocityVector {
    double vx;
    double vy;
    double vz;
};

velocityVector velocityHyperbolic(double _mt, double _px, double _py, double _y, 
  double _eta, double _etaDiff, double _tau);

double smoothingKernel(spatialVector _r, double _gammaz, double _tau, 
  double _Rgx, double _Rgy, double _Rgz);

double smoothingKernelInvariant(spatialVector _r, velocityVector _v, double _R, 
  double _tau);