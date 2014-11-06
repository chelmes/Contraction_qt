/*
 * BasicOperator.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <vector>

#include "GlobalData.h"
#include "Perambulator.h"
#include "propagator_io.h"
#include "quark.h"
#include "typedefs.h"
#include "VdaggerV.h"

// struct for Look-up table in create_gamma and get_operator. To read as
// "in column i the row[i]-element is non-zero and its value is value[i]"
// As Gamma matrices are 4x4 matrices, row and value are 4-vectors

struct lookup {
  int row[4];
  std::complex<double> value[4];
  };

class BasicOperator {

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  BasicOperator();
//  ~BasicOperator ();

  void init_operator_u (const int particle_no, const int t_source, 
                        const int t_sink, const char dilution, const int displ);
  void init_operator_d(const int particle_no, const int t_source, 
                       const int t_sink, const char dilution, const int displ);
  void get_operator_charged(array_Xcd_d2_eigen& op_1, const int particle_no, 
                            const int t_sink, 
                            const int dirac, const int p) const;
  void get_operator_g5(vec_Xcd_eigen& op_1, const int particle_no, 
                       const int dirac, const int p) const;
  void get_operator_uncharged(vec_Xcd_eigen& op_1, const int particle_no, 
                              const int dirac, const int p) const;

  void read_rnd_vectors_from_file (const int config_i);

  inline void set_basic(const size_t config){
    peram.read_perambulators_from_file(config);
    read_rnd_vectors_from_file(config);
    vdaggerv.build_vdaggerv(config);
    vdaggerv.build_rvdaggervr(config, rnd_vec);
  }

protected:
  LapH::Perambulator peram;
  std::vector<LapH::RandomVector> rnd_vec;
  LapH::VdaggerV vdaggerv;
  array_Xcd_d4_eigen contraction_dagger;
  array_Xcd_d5_eigen contraction;
  std::vector<struct lookup>  gamma;

};

#endif /* BASICOPERATOR_H_ */
