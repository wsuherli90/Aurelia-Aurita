#ifndef SYZYGY_BASIS_H
#define SYZYGY_BASIS_H

#include <vector>
#include <cmath>

namespace Aurelia {
namespace Crypto {
namespace Invariant {


class SyzygyBasis {
public:
    struct BasisProperties {
        double density;
        double algebraic_rank;
        bool has_linear_syzygies;
    };

    static BasisProperties extractFromDirectors(const std::vector<double>& eigenvalues) {
        BasisProperties props;
        

        double lambda1 = eigenvalues[0];
        

        
        double normalized_order = 1.5 * (lambda1 - 0.33333333);
        if (normalized_order < 0) normalized_order = 0;
        
        props.algebraic_rank = 1.0 - normalized_order; 
        
        props.has_linear_syzygies = (props.algebraic_rank < 0.2);
        
        return props;
    }
};

} // namespace Invariant
} // namespace Crypto
} // namespace Aurelia

#endif
