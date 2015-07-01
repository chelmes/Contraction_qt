/*! \file CrossOperator.h
    \brief Class to build the open ended scissor diagrams of 4 pt meson functions

    Each crossdiagram for mesons can be divided into several scissors. This class
    initializes 2 scissors per Instance and implements swapping and constructing them.
    Functions: -
*/

#ifndef CROSSOPERATOR_H_
#define CROSSOPERATOR_H_
//!
//! \class CrossOperator
//! \brief Class to build the open ended scissor diagrams of 4 pt meson functions
//!
#include <algorithm>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <typeinfo>
#include <vector>

#include "global_data.h"
#include "BasicOperator.h"
#include "typedefs.h"
#include "VdaggerV.h"

namespace LapH {

class CrossOperator{

public:
  //! \brief Empty Ctor

  //! Empty Constructor inheriting from X as private variable X is a 5
  //! vector of 5d arrays of Eigen-matrices. Each element has dimensions number
  //! of elements, index of quark line sink, index of quark line source and 3
  //! random vectors.
  CrossOperator() : X() {};
  //! \brief Filled Ctor, initialized with <number> of Xes

  //! Constructor initializing X with <number> empty scissors. Ensurance of 
  //! enough memory.
  CrossOperator(const size_t number);
  //! \brief Default Dtor
  ~CrossOperator() {};
  //! \brief Constructs one cross operator as an vector entry in X.

  //! construct builds a scissor object out of a basic operator comprising the
  //! according quark lines and randomvectors and a vdaggerv object taking care
  //! eigensystems. The argument nb describes the position in the vector and
  //! thus the scissor to build.
  //! t_source is the source time index of the diagram and t_sink is the sink
  //! time index of the diagram. type refers to if the second sink index should
  //! lie before (type = 1) or after (typ[e = 0) the argument sink time index
  void construct(const BasicOperator& basic, const VdaggerV& vdaggerv, 
                 const size_t nb, const int t_source, const int t_sink,
                 const size_t type);

  void compute_X(const BasicOperator& basic, const size_t id_si, 
                 const Eigen::MatrixXcd& Q2, const Eigen::MatrixXcd& VdaggerV, 
                 Eigen::MatrixXcd& X);

  void swap(const size_t nb1, const size_t nb2);

  inline const Eigen::MatrixXcd& operator()(const size_t nb, 
                                            const size_t id_so,
                                            const size_t id_si, 
                                            const size_t rnd1, 
                                            const size_t rnd2, 
                                            const size_t rnd3) const {
    return X[nb][id_si][id_so][rnd1][rnd2][rnd3];
  }

private:
  std::vector<array_Xcd_d5_eigen> X;

};

}// end of namespace

#endif // CROSSOPERATOR_H_
