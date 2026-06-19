/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016   arXiv:1312.4160                   *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it.                   *
*                                                                             *
*******************************************************************************/

#include <iostream>

#include <vhlle/simulation.h>

int main(int argc, char** argv) {
    try {
        Simulation simulation(argc, argv);
        simulation.printParameters();
        simulation.setup();
        while (simulation.step());
        simulation.finalize();
    } catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "FATAL: unknown exception\n";
        return 1;
    }
    return 0;
}