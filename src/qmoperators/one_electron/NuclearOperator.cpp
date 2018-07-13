#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NuclearOperator.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

NuclearOperator::~NuclearOperator() { }

double NuclearOperator::trace(const Nuclei &nucs_b) const {
    const Nuclei &nucs_a = this->nuclei;
    double E_nuc = 0.0;
    for (int a = 0; a < nucs_a.size(); a++) {
        const Nucleus &nuc_a = nucs_a[a];
        double Z_a = nuc_a.getCharge();
        const double *R_a = nuc_a.getCoord();
        for (int b = 0; b < nucs_b.size(); b++) {
            const Nucleus &nuc_b = nucs_b[b];
            double Z_b = nuc_b.getCharge();
            const double *R_b = nuc_b.getCoord();
            E_nuc += Z_a*Z_b/math_utils::calc_distance(R_a, R_b);
        }
    }
    return E_nuc;
}

} //namespace mrchem
