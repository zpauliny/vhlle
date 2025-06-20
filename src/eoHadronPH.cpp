#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "eos.h"
#include "eoHadronPH.h"

EoSauxPH::EoSauxPH(const std::string filename_tab, const std::string filename_base)
{
  // ##### reading the tables  ######
  // Open the file
  std::ifstream file_tab(filename_tab);
  if (!file_tab.is_open()) {
    std::cerr << "Error: Could not open file " << filename_tab << std::endl;
    return;
  }
  // read the dimensions of the data
  std::string line;
  std::getline(file_tab, line);
  std::istringstream iss(line);
  // Read the header line and extract the number of steps and the ranges
  // The first line contains the ranges and steps for e and nb
  iss >> nb_min >> d_nb >> nb_steps;
  // Read the second line for e_min, d_e, e_steps
  std::getline(file_tab, line);
  iss.clear();
  iss.str(line);
  iss >> e_min >> d_e >> e_steps;
  etab.resize(e_steps);
  nbtab.resize(nb_steps);
  // 2D arrays holding the thermodynamic quantities
  ptab.resize(nb_steps, std::vector<double>(e_steps, 0.0));  // is it correct?
  Ttab.resize(nb_steps, std::vector<double>(e_steps, 0.0));
  stab.resize(nb_steps, std::vector<double>(e_steps, 0.0));
  mubtab.resize(nb_steps, std::vector<double>(e_steps, 0.0));
  muqtab.resize(nb_steps, std::vector<double>(e_steps, 0.0));
  mustab.resize(nb_steps, std::vector<double>(e_steps, 0.0));
  // Read the data from the file
  int in = 0, ie = 0;
  while(std::getline(file_tab, line)) {
      std::replace(line.begin(), line.end(), 'D', 'E');
      std::istringstream iss(line);
      // there are some blank lines in the file, therefore we need to check if the reading was successful
      if(iss >> nbtab[in] >> etab[ie] >> ptab[in][ie] >> stab[in][ie] >> Ttab[in][ie]
           >> mubtab[in][ie] >> mustab[in][ie] >> muqtab[in][ie]) {
      ie++;
      if(ie == e_steps) {
        ie = 0;
        in++;
        if(in == nb_steps) break; // Stop reading if we have filled the arrays
      }
      } else {
        //std::cerr << "Error: Incomplete data at line " << in * e_steps + ie + 1 << std::endl;
      }
  }
  // Close the file
  file_tab.close();
  std::cout << "EoSauxPH: e range is " << etab[0] << " to " << etab[e_steps - 1] << std::endl;
  std::cout << "EoSauxPH: nb range is " << nbtab[0] << " to " << nbtab[nb_steps - 1] << std::endl;
  // ##### reading the baseline file  ######
  std::ifstream file_base(filename_base);
  if (!file_base.is_open()) {
    std::cerr << "Error: Could not open file " << filename_base << std::endl;
    return;
  }
  // Read the baseline data
  // The first line contains the ranges and steps for nb
  std::getline(file_base, line);
  std::istringstream iss_base(line);
  iss_base >> nb_base_min >> d_nb_base >> nb_base_steps;
  nb_base.resize(nb_base_steps);
  e_base.resize(nb_base_steps);
  // Read the data from the file
  int ibase = 0;
  while (std::getline(file_base, line)) {
    std::replace(line.begin(), line.end(), 'D', 'E');
    std::istringstream iss(line);
    if (iss >> nb_base[ibase] >> e_base[ibase]) {
      ibase++;
    }
    if(ibase == nb_base_steps+1) {
      std::cerr << "Error: Too many lines in the baseline file." << std::endl;
      break; // Stop reading if we have filled the arrays
    }
  }
  std::cout << "EoSauxPH: nb_base range is " << nb_base[0] << " to " << nb_base[nb_base_steps - 1] << std::endl;
  // Close the file
  file_base.close();
}


double EoSauxPH::emin(double nb)
{
  // Check if the input value is within the valid range
  // this part I may need to delete since the EoS is supposed to always return a value
 // if (nb < nbtab[0] || nb > nbtab[nb_steps - 1]) {
 //   std::cerr << "Error: Input values out of range." << std::endl;
 //   return -1.0; // Return an error value
 // }
  // Compute indices for interpolation (assuming equidistant grids), in baryon density first
  double fn = (nb - nb_min) / d_nb;
  int in = static_cast<int>(fn);
  // Clamp indices to valid range for interpolation
  if (in < 0) in = 0;
  if (in > nb_steps - 2) in = nb_steps - 2;
  double dn = fn - in;
  double ebas = e_base[in] + dn * (e_base[in + 1] - e_base[in]);
  return ebas + etab[1]; // Fermi energy density at given baryon density
}


double EoSauxPH::emax(double nb)
{
  // Check if the input value is within the valid range
  // this part I may need to delete since the EoS is supposed to always return a value
 // if (nb < nbtab[0] || nb > nbtab[nb_steps - 1]) {
 //   std::cerr << "Error: Input values out of range." << std::endl;
 //   return -1.0; // Return an error value
 // }
  // Compute indices for interpolation (assuming equidistant grids), in baryon density first
  double fn = (nb - nb_min) / d_nb;
  int in = static_cast<int>(fn);
  // Clamp indices to valid range for interpolation
  if (in < 0) in = 0;
  if (in > nb_steps - 2) in = nb_steps - 2;
  double dn = fn - in;
  double ebas = e_base[in] + dn * (e_base[in + 1] - e_base[in]);
  return ebas + etab[e_steps - 1]; // Fermi energy density at given baryon density
}


double EoSauxPH::p(double e, double nb)
{
  // Check if the input values are within the valid range
  // this part I may need to delete since the EoS is supposed to always return a value
 // if (nb < nbtab[0] || nb > nbtab[nb_steps - 1]) {
 //   std::cerr << "Error: Input values out of range." << std::endl;
 //   return -1.0; // Return an error value
 // }
  // Compute indices for interpolation (assuming equidistant grids), in baryon density first
  double fn = (nb - nb_min) / d_nb;
  int in = static_cast<int>(fn);
  // Clamp indices to valid range for interpolation
  if (in < 0) in = 0;
  if (in > nb_steps - 2) in = nb_steps - 2;
  double dn = fn - in;
  double ebas = e_base[in] + dn * (e_base[in + 1] - e_base[in]);
  e = e - ebas;  // subtract the Fermi energy density for given baryon density
  double fe = (e - e_min) / d_e;
  int ie = static_cast<int>(fe);
  // Clamp indices to valid range for interpolation
  if (ie < 0) ie = 0;
  if (ie > e_steps - 2) ie = e_steps - 2;
  double de = fe - ie;
  // Bilinear interpolation
  double wn [2] = {1.0 - dn, dn};
  double we [2] = {1.0 - de, de};
  double p = 0.0;
  for (int jn = 0; jn < 2; jn++) {
    for (int je = 0; je < 2; je++) {
      p += wn[jn] * we[je] * ptab[in + jn][ie + je];
    }
  }
  return std::max(0.0, p);  // Ensure pressure is non-negative
}


void EoSauxPH::eos(double e, double nb, double nq, double ns, double& T, double& mub, double& muq, double& mus, double& p)
{
  // Check if the input values are within the valid range
  //if (nb < nbtab[0] || nb > nbtab[nb_steps - 1]) {
  //  std::cerr << "Error: Input values out of range." << std::endl;
  //  return;
  //}
  // Compute indices for interpolation (assuming equidistant grids), in baryon density first
  double fn = (nb - nb_min) / d_nb;
  int in = static_cast<int>(fn);
  // Clamp indices to valid range for interpolation
  if (in < 0) {  // negative baryon density is unphysical
    T = mub = muq = mus = 0.0;
    p = 0.0;
    return;
  }
  if (in > nb_steps - 2) in = nb_steps - 2;
  double dn = fn - in;
  double ebas = e_base[in] + dn * (e_base[in + 1] - e_base[in]);
  e = e - ebas;  // subtract the Fermi energy density for given baryon density
  double fe = (e - e_min) / d_e;
  int ie = static_cast<int>(fe);
  // Clamp indices to valid range for interpolation
  if (ie < 0) {  // negative e' is clearly unphysical
    T = mub = muq = mus = 0.0;
    p = 0.0;
    return;
  }
  if (ie > e_steps - 2) ie = e_steps - 2;
  double de = fe - ie;
  // Bilinear interpolation
  double wn [2] = {1.0 - dn, dn};
  double we [2] = {1.0 - de, de};
  p = 0.0; T = 0.0; mub = 0.0; muq = 0.0; mus = 0.0;
  for (int jn = 0; jn < 2; jn++) {
    for (int je = 0; je < 2; je++) {
      p += wn[jn] * we[je] * ptab[in + jn][ie + je];
      T += wn[jn] * we[je] * Ttab[in + jn][ie + je];
      mub += wn[jn] * we[je] * mubtab[in + jn][ie + je];
      muq += wn[jn] * we[je] * muqtab[in + jn][ie + je];
      mus += wn[jn] * we[je] * mustab[in + jn][ie + je];
    }
  }
  //p = std::max(0.0, p);  // Ensure pressure is non-negative
  T = std::max(0.0, T);  // Ensure temperature is non-negative
}


EoShadronPH::EoShadronPH()
{
  // Initialize the EoSauxPH objects with the appropriate filenames
  eosSmall = new EoSauxPH("eos/EoSTable1.dat", "eos/Baseline1.dat");
  eosBig = new EoSauxPH("eos/EoSTable2.dat", "eos/Baseline2.dat");
}

EoShadronPH::~EoShadronPH()
{
  // Clean up the dynamically allocated EoSauxPH objects
  delete eosSmall;
  delete eosBig;
}

double EoShadronPH::p(double e, double nb, double nq, double ns)
{
  // Call the appropriate EoSauxPH object based on baryon and energy density
  // note: nq and ns are ignored since the original tables are only in (e,nb)
  if (nb < eosSmall->nbmax() && e < eosSmall->emax(nb)) {
    return eosSmall->p(e, nb);
  } else {
    return eosBig->p(e, nb);
  }
}

void EoShadronPH::eos(double e, double nb, double nq, double ns, double& T, double& mub, double& muq, double& mus, double& p)
{
  // Call the appropriate EoSauxPH object based on the baryon and energy density
  if (e < eosSmall->emin(nb)) {
    T = mub = muq = mus = 0.0;
    p = 0.0;
    return;
  } else if (nb < eosSmall->nbmax() && e < eosSmall->emax(nb)) {
    eosSmall->eos(e, nb, nq, ns, T, mub, muq, mus, p);
  } else {
    eosBig->eos(e, nb, nq, ns, T, mub, muq, mus, p);
  }
}
