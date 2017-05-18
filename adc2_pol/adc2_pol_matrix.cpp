#include "adc2_pol_matrix.hpp"
#include "analysis/adc_analyzer.hpp"
#include "integral_table.hpp"
#include "./scf_data/scf_data_reader.hpp"

#include <iostream>
#include <sstream>


using namespace std;

#if defined INTEL
#define INIT_DATA                adc2_pol_adapter_mp_init_data_
#define SELECT_ATOM_IS           select_fano_mp_select_atom_is_
#define SELECT_ATOM_D            select_fano_mp_select_atom_d_
#define GET_DIAG_ADC2_DIRECT     get_matrix_mp_get_diag_adc2_direct_
#define GET_OFFDIAG_ADC2_DIRECT  get_matrix_mp_get_offdiag_adc2_direct_
#define GET_PHPH_ADC2            get_matrix_mp_get_phph_adc2_
#define GET_PH_2P2H              get_matrix_mp_get_ph_2p2h_
#define GET_2P2H2P2H_DG2S        get_matrix_mp_get_2p2h2p2h_dg2s_
#elif defined PGI
#define INIT_DATA                adc2_pol_adapter_init_data_
#define SELECT_ATOM_IS           select_fano_select_atom_is_
#define SELECT_ATOM_D            select_fano_select_atom_d_
#define GET_DIAG_ADC2_DIRECT     get_matrix_get_diag_adc2_direct_
#define GET_OFFDIAG_ADC2_DIRECT  get_matrix_get_offdiag_adc2_direct_
#define GET_PHPH_ADC2            get_matrix_get_phph_adc2_
#define GET_PH_2P2H              get_matrix_get_ph_2p2h_
#define GET_2P2H2P2H_DG2S        get_matrix_get_2p2h2p2h_dg2s_
#elif defined GCC
#define INIT_DATA                __adc2_pol_adapter_MOD_init_data
#define SELECT_ATOM_IS           __select_fano_MOD_select_atom_is
#define SELECT_ATOM_D            __select_fano_MOD_select_atom_d
#define GET_DIAG_ADC2_DIRECT     __get_matrix_MOD_get_diag_adc2_direct
#define GET_OFFDIAG_ADC2_DIRECT  __get_matrix_MOD_get_offdiag_adc2_direct
#define GET_PHPH_ADC2            __get_matrix_MOD_get_phph_adc2
#define GET_PH_2P2H              __get_matrix_MOD_get_ph_2p2h
#define GET_2P2H2P2H_DG2S        __get_matrix_MOD_get_2p2h2p2h_dg2s
#endif



extern "C" {
  void INIT_DATA(int*, int*, int*, int*, double*, int*, int*, int*, int*);
  void SELECT_ATOM_IS(int*);
  void SELECT_ATOM_D(int*, int*);
  void GET_DIAG_ADC2_DIRECT(int*, int*, int*, double*);
  void GET_OFFDIAG_ADC2_DIRECT(int*,int*,double*);
  void GET_PHPH_ADC2(int*,int*,double*); 
  void GET_PH_2P2H(int*,int*,int*,int*,double*);
  void GET_2P2H2P2H_DG2S(int*,int*,double*);
}


extern Integral_table* integral_table;

int ADC2_pol_matrix::accept_analyzer( ADC_analyzer& adc_an)
{
  return adc_an.analyze(*this);
}


ADC2_pol_matrix::ADC2_pol_matrix(SCF_data_reader& phis, unsigned int sym, ListOfOrbs holes) 
  : phis_(&phis), sym_(sym), count_mult(0)
{
  int nBas = phis_->number_orbitals();
  int* orbs = new int[nBas];
  int* syms = new int[nBas];
 
  int hc_size = holes.size();
  if (!hc_size)
     for(int i = 0; i< phis_->number_occupied(); i++)
         holes.insert(i+1);
  hc_size = phis_->number_occupied(); 
  int* hc = new int[hc_size];
  set<unsigned>::iterator it;
  int i = 0;
  for (it = holes.begin(); it != holes.end(); i++,it++)
    hc[i] = *it;

  
  double* en = new double[nBas];
  
  for(int i = 0; i < nBas; i++) {
    orbs[i] = i+1;
    syms[i] = phis_->irrep(orbs[i]-1)+1;
    en[i] = phis_->energy(orbs[i]-1);
  }

  int nirreps = phis_->number_irreps();
  int nOcc = phis_->number_occupied();

 

  int dummy_sym = sym+1;
  INIT_DATA(&nBas, &nOcc, orbs, 
				 syms, en, &nirreps, &dummy_sym, &hc_size, hc);


  kpq = new int [(nBas*nBas*4*nOcc*nOcc+1)*7];

  int flag =-1;
  SELECT_ATOM_IS(kpq);
  SELECT_ATOM_D(kpq, &flag);


  delete [] orbs;
  delete [] syms;
  delete [] en;

  main_size = kpq[0];
  sat_size = kpq[1]+kpq[2]+kpq[3]+2*kpq[4];
 
}





ADC2_pol_matrix::~ADC2_pol_matrix()
{
   delete [] kpq;
}


// Builds the whole matrix
void ADC2_pol_matrix::build_matrix(Blas_matrix& mat)
{
  int size = this->size();
  mat.allocate(size,size);
  Blas_matrix diag(size);
  mat = 0.;
  diag = 0.;

  GET_OFFDIAG_ADC2_DIRECT(&size,kpq,&mat(0,0));
  GET_DIAG_ADC2_DIRECT(&main_size, &sat_size, kpq, &diag(0));


  mat.add_diag(diag);

 
}

std::string ADC2_pol_matrix::get_conf(unsigned int i) const{
  ostringstream oss;
  i++;
  if (i <= main_size)
    oss << '<' << kpq[7*i+2] << ',' << kpq[7*i+4] << '|';
  else if (i <= size())
    oss << '<' << kpq[7*i+2] << ',' << kpq[7*i+4] << ',' << kpq[7*i+6] << ',' << kpq[7*i+1] << '|';
  
  return oss.str();
}



int ADC2_pol_matrix::operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count) 
{
  
  // Emulate Lanczos subspace iterations, see adc_analyzer.cpp
  if (count_mult < 2*main_size) {
    // main part
    multiply_main_block(block_out, block_in, count);
    //main/sat block
    multiply_main_sat_block(block_out, block_in, count);
  }
  // sat blocks
  multiply_sat_block(block_out, block_in, count);

  count_mult+=count;
  return count;
  
}



void ADC2_pol_matrix::multiply_main_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count) 
{

  Blas_matrix main_block(main_size, main_size);
  main_block = 0.;

  
  GET_PHPH_ADC2(&main_size,kpq,&main_block(0,0)); 

  for (unsigned i = 0; i < main_block.rows(); i++)
    for (unsigned j = 0; j < i; j++)
      main_block(i,j) = main_block(j,i);

//   for(int i = 0; i < main_block.rows(); i++)
//     for(int j = 0; j < main_block.cols(); j++)
//       if (main_block(i,j) < 10e-6) main_block(i,j) = 0.;
      

  Blas_matrix** vec_in = block_in; 
  Blas_matrix** vec_out = block_out;
  
  for (int i = 0; i < count; i++) {
    (**vec_out)(0, main_size, 0, 1)
      .dgemv('N', main_block, (**vec_in)(0, main_size, 0, 1));
    vec_in++; vec_out++;
  }

  
}

void ADC2_pol_matrix::multiply_main_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, 
					     int count) 
{
  Blas_matrix sat_block(sat_size, main_size);
  sat_block = 0.;
  int i1 = 1;
  GET_PH_2P2H(&sat_size,&i1,&main_size,kpq,&sat_block(0,0));

//   for(int i = 0; i < sat_block.rows(); i++)
//     for(int j = 0; j < sat_block.cols(); j++)
//       if (sat_block(i,j) < 10e-6) sat_block(i,j) = 0.;
      


  Blas_matrix** vec_in = block_in; 
  Blas_matrix** vec_out = block_out;
  for (int i = 0; i < count; i++) {
    (**vec_out)(main_size, sat_size, 0, 1)
      .dgemv('N', sat_block, (**vec_in)(0, main_size, 0, 1));
    (**vec_out)(0, main_size, 0, 1)
      .dgemv('T', sat_block, (**vec_in)(main_size, sat_size, 0, 1));
    vec_out++; vec_in++;
  }


}

void ADC2_pol_matrix::multiply_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, 
					int count) 
{

  Blas_matrix sat_diag(sat_size);
  sat_diag = 0.;

  GET_2P2H2P2H_DG2S(&sat_size,kpq,&sat_diag(0));

#pragma omp parallel for
  for(int col = main_size; col < size(); col++) {
    Blas_matrix** vec_in = block_in, **vec_out = block_out;
    //if (sat_diag(col-main_size) < 10e-6) continue;
    for(int vec = 0; vec < count; vec++) {
      (**vec_out)(col) += sat_diag(col-main_size) * (**vec_in)(col);
      vec_out++, vec_in++;
    }
  }

}


