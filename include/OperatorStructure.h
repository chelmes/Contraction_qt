#ifndef OPERATORSTRUCTURE_H_
#define OPERATORSTRUCTURE_H_ 

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "CrossOperator.h" 
#include "GlobalData.h"
#include "BasicOperator.h"
#include "Perambulator.h"
#include "typedefs.h"

#include "omp.h"

namespace LapH {

  struct pdg{
    size_t id;
    std::array<int,3> p;
    std::array<int,3> dis;
    std::array<size_t,4> gamma;
  };

  struct pdg_C2 {
    size_t p_sq;
    size_t dg_so;
    size_t dg_si;
    std::list<std::pair<size_t, size_t> > index;
  };

//  struct pdg_C4 {
//    size_t p_sq_so;
//    size_t p_sq_si;
//    size_t dg_so;
//    size_t dg_si;
//    std::list<std::tuple<size_t, size_t, size_t, size_t> > index;
//  };


  void init_from_infile(std::vector<pdg>& op_Corr, 
                        std::vector<pdg_C2>& op_C2);

  void set_default(std::vector<pdg>& op);
  void set_default(std::vector<pdg>& op, std::vector<pdg_C2>& op_C2);
 
  void create_momenta();

}// end of namespace

#endif // OPERATORSTRUCTURE_H_