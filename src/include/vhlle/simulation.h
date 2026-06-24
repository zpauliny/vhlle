#pragma once

#include <ctime>
#include <deque>
#include <string>

class EoS;
class Fluid;
class Hydro;
class TransportCoeff;
class VtkOutput;
struct Particle;

class Simulation {
public:
    explicit Simulation(int argc, char** argv);
    ~Simulation();

    void readParameters(const char* parFile);
    void printParameters();
    void setup();
    bool step();
    void finalize();

private:
    void readCommandLine(int argc, char** argv);
    void checkGridDimension(int n, char axis);
    void checkGridBorders(double min, double max, const std::string& axis);
    bool parse_bool(const std::string& value);
    void expandGrid2x();

    EoS*                  eos      = nullptr;
    EoS*                  eosH     = nullptr;
    TransportCoeff*       trcoeff  = nullptr;
    Fluid*                f        = nullptr;
    Hydro*                h        = nullptr;
    VtkOutput*            vtk_out  = nullptr;
    std::deque<Particle>* particles = nullptr;

    int    nx {100}, ny {100}, nz {100};
    int    eosType {1}, eosTypeHadron {0};
    int    etaSparam {0}, zetaSparam {0};
    int    icModel {1}, glauberVariable {1};
    int    smoothingType {0}, minParticlesFO {15};
    bool   freezeoutOnly {false}, freezeoutExtend {false};
    bool   vorticityOn {false}, cartesian {false};
    double xmin {-5.0}, xmax {5.0}, ymin {-5.0}, ymax {5.0};
    double etamin {-5.0}, etamax {5.0};
    double tau0 {1.0}, tauMax {20.0}, tauResize {4.0}, dtau {0.05};
    double etaS {0.08}, zetaS {0.0}, eCrit {0.5};
    double etaSEpsilonMin {5.}, etaSMin {0.08};
    double etaSShiftMuB {0.}, etaSScaleMuB {0.};
    double al {0.}, ah {0.}, aRho {0.}, T0 {0.15};
    double zetaSPeakEpsilon {5.}, zetaSScaleBeta {0.103};
    double zetaSSigmaMinus {0.1}, zetaSSigmaPlus {0.1};
    double epsilon0, Rgt {1.0}, Rgz {1.0};
    double impactPar, s0ScaleFactor, gaussian_sigma {0.0};
    std::string collSystem, outputDir {"data"}, isInputFile, vtk_values {""};

    double ctime    {0.0};
    bool   resized  {false};
    int    nelements {1};
    time_t tstart   {0};
    time_t tinit    {0};
    double timeInitFO {0.0};
};