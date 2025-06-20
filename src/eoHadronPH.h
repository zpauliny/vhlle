class EoS;

// auxiliary EoS table class
class EoSauxPH {
public:
    EoSauxPH(const std::string filename_tab, const std::string filename_base);
    double emin() const { return e_min; }
    double nbmin() const { return nb_min; }
    double e_therm_max() const { return e_min + d_e * (e_steps - 1); }
    double emin(double nb); // MIN energy density at given baryon density that the table covers
    double emax(double nb); // MAX energy density at given baryon density that the table covers
    double nbmax() const { return nb_min + d_nb * (nb_steps - 1); }
    double p(double e, double nb);
    void eos(double e, double nb, double nq, double ns, double& T, double& mub, double& muq, double& mus, double& p);
private:
    double e_min;  // minimum energy density
    double d_e;    // energy density step
    int e_steps;   // number of energy density steps
    double nb_min;  // minimum baryon density
    double d_nb;    // baryon density step
    int nb_steps;   // number of baryon density steps
    std::vector<double> etab;  // steps in energy density
    std::vector<double> nbtab;  // steps in baryon density
    // 2D arrays holding the thermodynamic quantities
    std::vector<std::vector<double>> ptab;  // pressure table
    std::vector<std::vector<double>> Ttab;   // temperature table
    std::vector<std::vector<double>> stab;   // entropy table
    std::vector<std::vector<double>> mubtab;  // baryon chemical potential table
    std::vector<std::vector<double>> muqtab;  // electric chemical potential table
    std::vector<std::vector<double>> mustab;  // strangeness chemical potential table
    // ---------- baseline table
    double nb_base_min;  // minimum baryon density for the baseline table
    double d_nb_base;    // baryon density step for the baseline table
    int nb_base_steps;   // number of baryon density steps for the baseline table
    std::vector<double> e_base;  // Fermi energy density
    std::vector<double> nb_base;  // steps in baryon density for the baseline table
};


class EoShadronPH : public EoS {
public:
    EoShadronPH();
    ~EoShadronPH();
    //virtual double emin(double nb) { return eosSmall->emin(nb); } // minimum energy density for given nb where the EoS is robust
    //virtual double emax(double nb) { return eosBig->emax(nb); } // maximum tabulated energy density for given nb
    virtual double p(double e, double nb, double nq, double ns);
    virtual void eos(double e, double nb, double nq, double ns, double& T, double& mub, double& muq, double& mus, double& p);
private:
    EoSauxPH *eosSmall, *eosBig;  // pointer to the auxiliary EoS objects holding the big and small EoS tables
};
