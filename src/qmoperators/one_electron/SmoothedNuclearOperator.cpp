#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "SmoothedNuclearOperator.h"
#include "Nucleus.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

void SmoothedNuclearOperator::setupU_2(double prec) {
    int oldprec = Printer::setPrecision(5);
    Printer::printHeader(0, "Setting up smoothed nuclear potential");
    println(0, " Nr  Element         Charge        Precision     Smoothing ");
    Printer::printSeparator(0, '-');

    const Nuclei &nucs = this->nuclei;
    std::vector<DoubleFunction> funcs;
    for (int k = 0; k < nucs.size(); k++) {
        double Z_k = nucs[k].getCharge();
        double Z_5 = std::pow(Z_k, 5.0);
        double c = std::pow(0.00435*prec/Z_5, 1.0/3.0);

        std::stringstream symbol;
        symbol << nucs[k].getElement().getSymbol();
        symbol << "  ";
        printout(0, std::setw(3) << k+1 << "     ");
        printout(0, symbol.str()[0] << symbol.str()[1]);
        printout(0, std::setw(22) << Z_k);
        printout(0, std::setw(14) << prec);
        printout(0, std::setw(14) << c << std::endl);

        // u(r) ~ 1/r
        auto u = [] (double r) -> double {
            return - std::erf(r)/r
                   - 1.0/(3.0*MATHCONST::sqrt_pi)*(std::exp(-r*r)
                   + 16.0*std::exp(-4.0*r*r));
        };

        // f_k(r) = Z_k*u(r/c)/c ~ Z_k/r
        const Nucleus &nuc_k = nucs[k];
        auto f_k = [nuc_k, u, c] (const double *r) -> double {
            const double Z_k = nuc_k.getCharge();
            const double *R_k = nuc_k.getCoord();
            const double R = math_utils::calc_distance(R_k, r);
            return Z_k*u(R/c)/c;
        };
        funcs.push_back(f_k);
    }

    // f(r) = \sum_k f_k(r)
    auto f = [funcs] (const double *r) -> double {
        double result = 0.0;
        for (int k = 0; k < funcs.size(); k++) {
            result += funcs[k](r);
        }
        return result;
    };
    this->U_2.setReal(f);

    Printer::printSeparator(0, '=', 2);
    Printer::setPrecision(oldprec);
}

} //namespace mrchem
