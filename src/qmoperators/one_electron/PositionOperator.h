#pragma once

#include "RankOneTensorOperator.h"
#include "AnalyticPotential.h"

namespace mrchem {

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const double *o = 0) {
        double o_x = (o != 0) ? o[0] : 0.0;
        double o_y = (o != 0) ? o[1] : 0.0;
        double o_z = (o != 0) ? o[2] : 0.0;

        auto f_x = [o_x] (const double *r) -> double { return r[0] - o_x; };
        auto f_y = [o_y] (const double *r) -> double { return r[1] - o_y; };
        auto f_z = [o_z] (const double *r) -> double { return r[2] - o_z; };

        this->r_x.setReal(f_x);
        this->r_y.setReal(f_y);
        this->r_z.setReal(f_z);

        RankOneTensorOperator &r = (*this);
        r[0] = r_x;
        r[1] = r_y;
        r[2] = r_z;
    }

protected:
    AnalyticPotential r_x;
    AnalyticPotential r_y;
    AnalyticPotential r_z;
};

} //namespace mrchem
