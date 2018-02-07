#pragma once

#include <complex>

#include "MRCPP/Printer"
#include "MRCPP/MWFunctions"

using namespace mrcpp;

class QMFunction {
public:
 QMFunction(mrcpp::FunctionTree<3> *r = 0, mrcpp::FunctionTree<3> *i = 0) : re(r), im(i) { }
    QMFunction(const QMFunction &func) : re(0), im(0) { }
    QMFunction &operator=(const QMFunction &func) { NOT_IMPLEMENTED_ABORT; }
    virtual ~QMFunction() { clearReal(true); clearImag(true); }

    void deepCopy(QMFunction &inp);
    void shallowCopy(const QMFunction &inp);

    bool hasReal() const { if (this->re == 0) return false; return true; }
    bool hasImag() const { if (this->im == 0) return false; return true; }

    void allocReal();
    void allocImag();

    void clearReal(bool free);
    void clearImag(bool free);

    void setReal(mrcpp::FunctionTree<3> *real) { this->re = real; }
    void setImag(mrcpp::FunctionTree<3> *imag) { this->im = imag; }

    mrcpp::FunctionTree<3> &real() { return *this->re; }
    mrcpp::FunctionTree<3> &imag() { return *this->im; }

    const mrcpp::FunctionTree<3> &real() const { return *this->re; }
    const mrcpp::FunctionTree<3> &imag() const { return *this->im; }

    int getNNodes(int type = Total) const;
    double getSquareNorm(int type = Total) const;
    std::complex<double> dot(QMFunction &ket);

    void normalize();
    void operator*=(double c);

    enum NumType { Total, Real, Imag };

protected:
    mrcpp::FunctionTree<3> *re;
    mrcpp::FunctionTree<3> *im;
};

