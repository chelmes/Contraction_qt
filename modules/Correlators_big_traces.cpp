#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/****************** Compute second cross diagram ******************************/
/******************************************************************************/
void LapH::Correlators::build_C4_4(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();

  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  // TODO: must be changed in GlobalData {
  // TODO: }

  std::cout << "\n\tcomputing the traces of 2 pi_+/-:\r";
  clock_t time = clock();
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of 2 pi_+/-: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t_source = 0; t_source < Lt; ++t_source){
      const int t_source_1 = (t_source + 1) % Lt;
      // As a Test construct X from scratch
      //if(t_source != 0){
      //  if(t_source%2 == 0){
      //    //X.swap(1, 0);
      //    X.construct(basic, vdaggerv, 0, t_source, t_sink, 1);
      //    X.construct(basic, vdaggerv, 1, t_source_1, t_sink_1, 0);
      //  }
      //  else{
      //    //X.swap(0, 1);
      //    X.construct(basic, vdaggerv, 0, t_source, t_sink_1, 0);
      //    X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
      //  }
      //}
      //else{
      //    X.construct(basic, vdaggerv, 0, t_source, t_sink, 1);
      //    X.construct(basic, vdaggerv, 1, t_source_1, t_sink_1, 0);
      //}

      if(t_source != 0){
        if(t_source%2 == 0){
          //X.swap(1, 0);
          X.construct(basic, vdaggerv, 0, t_source_1, t_sink, 0);
          X.construct(basic, vdaggerv, 1, t_source, t_sink, 1);
        }
        else{
          //X.swap(0, 1);
          X.construct(basic, vdaggerv, 0, t_source_1, t_sink, 0);
          X.construct(basic, vdaggerv, 1, t_source, t_sink, 1);
        }
      }
      else{
          X.construct(basic, vdaggerv, 0, t_source_1, t_sink, 0);
          X.construct(basic, vdaggerv, 1, t_source, t_sink, 1);
      }
      // Might be more efficient to start loop at t_sink+1 and completely skip
      // t_source == t_sink
      if(t_source == t_sink)
        continue;

      // The parallelisation is not done with #pragma omp for because it is 
      // incompatible with auto loops
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        // different quantum number combinations denoted by the index i may be 
        // handled by different tasks. Their results are summed up in the end.
        // TODO: Further speedup by different parallelisation might be possible
//        #pragma omp task shared(op.id, i, i_0, i_1, i_2, i_3)
        for(const auto& i : op.index_pt){

          size_t i_0 = op_C4[i].index_Q2[0];
          size_t i_1 = op_C4[i].index_Corr[0];
          size_t i_2 = op_C4[i].index_Q2[1];
          size_t i_3 = op_C4[i].index_Corr[1];
          // complete diagramm. combine X and Y to four-trace
          // C4_mes = tr(D_u^-1(t_source     | t_sink      ) Gamma 
          //             D_d^-1(t_sink       | t_source + 1) Gamma 
          //             D_u^-1(t_source + 1 | t_sink + 1  ) Gamma 
          //             D_d^-1(t_sink + 1   | t_source    ) Gamma)
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          //for(const auto& rnd_it : rnd_vec_index) {

          //  if(t_source%2 == 0)

          //    priv_C4 += (X(1, i_0, i_1, rnd_it[3], rnd_it[0], rnd_it[1]) *
          //                X(0, i_2, i_3, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          //  else
          //    priv_C4 += std::conj(
          //               (X(1, i_0, i_1, rnd_it[3], rnd_it[0], rnd_it[1]) *
          //                X(0, i_2, i_3, rnd_it[1], rnd_it[2], rnd_it[3])).trace());
          //}
          for(const auto& rnd_it : rnd_vec_index) {
            // rnd_it[0], rnd_it[2] correspond to 1st quark
            // rnd_it[1], rnd_it[3] correspond to 2nd quark 
            if(t_source%2 == 0)

              priv_C4 += (X(1, i_0, i_1, rnd_it[0], rnd_it[1], rnd_it[2]) *
                          X(0, i_2, i_3, rnd_it[2], rnd_it[3], rnd_it[0])).trace();
              
            else
              priv_C4 += std::conj(
                         (X(1, i_0, i_1, rnd_it[0], rnd_it[1], rnd_it[2]) *
                          X(0, i_2, i_3, rnd_it[2], rnd_it[3], rnd_it[0])).trace());
          }
          #pragma omp critical
          {
            C4_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region
    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the traces of 2 k: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//TODO: Call that build_C4_3() ?
void LapH::Correlators::compute_meson_4pt_cross_trace(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();

  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  // TODO: must be changed in GlobalData {
  // TODO: }

  std::cout << "\n\tcomputing the traces of 2 k:\r";
  clock_t time = clock();
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of 2 k: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t_source = 0; t_source < Lt; ++t_source){
      const int t_source_1 = (t_source + 1) % Lt;
      // As a Test construct X from scratch
      if(t_source != 0){
        if(t_source%2 == 0){
          //X.swap(1, 0);
          X.construct(basic, vdaggerv, 0, t_source, t_sink, 0);
          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
        }
        else{
          //X.swap(0, 1);
          X.construct(basic, vdaggerv, 0, t_source, t_sink, 0);
          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
        }
      }
      else{
          X.construct(basic, vdaggerv, 0, t_source,   t_sink, 0);
          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
      }

      // Might be more efficient to start loop at t_sink+1 and completely skip
      // t_source == t_sink
      if(t_source == t_sink)
        continue;

      // The parallelisation is not done with #pragma omp for because it is 
      // incompatible with auto loops
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        // different quantum number combinations denoted by the index i may be 
        // handled by different tasks. Their results are summed up in the end.
        // TODO: Further speedup by different parallelisation might be possible
//        #pragma omp task shared(op.id, i, i_0, i_1, i_2, i_3)
        for(const auto& i : op.index_pt){

          size_t i_0 = op_C4[i].index_Q2[0];
          size_t i_1 = op_C4[i].index_Corr[0];
          size_t i_2 = op_C4[i].index_Q2[1];
          size_t i_3 = op_C4[i].index_Corr[1];
          // complete diagramm. combine X and Y to four-trace
          // C4_mes = tr(D_u^-1(t_source     | t_sink      ) Gamma 
          //             D_d^-1(t_sink       | t_source + 1) Gamma 
          //             D_u^-1(t_source + 1 | t_sink + 1  ) Gamma 
          //             D_d^-1(t_sink + 1   | t_source    ) Gamma)
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          //for(const auto& rnd_it : rnd_vec_index) {

          //  if(t_source%2 == 0)

          //    priv_C4 += (X(1, i_0, i_1, rnd_it[3], rnd_it[0], rnd_it[1]) *
          //                X(0, i_2, i_3, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          //  else
          //    priv_C4 += std::conj(
          //               (X(1, i_0, i_1, rnd_it[3], rnd_it[0], rnd_it[1]) *
          //                X(0, i_2, i_3, rnd_it[1], rnd_it[2], rnd_it[3])).trace());
          //}
          for(const auto& rnd_it : rnd_vec_index) {

            if(t_source%2 == 0)

              //priv_C4 += (X(1, i_0, i_1, rnd_it[0], rnd_it[1], rnd_it[2]) *
              //            X(0, i_2, i_3, rnd_it[2], rnd_it[3], rnd_it[0])).trace();
              priv_C4 += (X(1, i_0, i_1, rnd_it[0], rnd_it[1], rnd_it[2]) *
                          X(0, i_2, i_3, rnd_it[2], rnd_it[3], rnd_it[0])).trace();
              
            else
              //priv_C4 += std::conj(
              //           (X(1, i_0, i_1, rnd_it[0], rnd_it[1], rnd_it[2]) *
              //            X(0, i_2, i_3, rnd_it[2], rnd_it[3], rnd_it[0])).trace());
              priv_C4 += std::conj(
                         (X(1, i_0, i_1, rnd_it[0], rnd_it[1], rnd_it[2]) *
                          X(0, i_2, i_3, rnd_it[2], rnd_it[3], rnd_it[0])).trace());
          }
          #pragma omp critical
          {
            C4_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region
    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the traces of 2 pi_+/-: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//TODO: is that still necessary?
void LapH::Correlators::write_C4_3(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();

  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  const size_t norm1 = Lt*rnd_vec_index.size();

  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();

  if(op_C4_IO.size() == 0)
    return;

  // normalisation
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile - filename
  // C4_mes  - boost structure containing all correlators

  sprintf(outfile, "%s/C4_3_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I2+_3", C4_mes);

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//TODO: is that still necessary?
void LapH::Correlators::write_C4_4(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();

  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  const size_t norm1 = Lt*rnd_vec_index.size();

  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();

  if(op_C4_IO.size() == 0)
    return;

  // normalisation
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile - filename
  // C4_mes  - boost structure containing all correlators

  sprintf(outfile, "%s/C4_4_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I2+_4", C4_mes);

}
