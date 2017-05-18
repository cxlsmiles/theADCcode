#include "nd_adc3_cap_matrix.hpp"
#include "analysis/adc_analyzer.hpp"
#include "integral_table.hpp"
#include "integral_blocks.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "ndadc3_ip/self_energy.hpp"

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>

extern "C" {
#include "my_globals.h"
} 

using namespace std;

// The file contains the implementation of the nd_adc3_matrix class methods. 
// These are based on Joerg's nd_adc3 functions called here to build
// parts of ND ADC3 matrix. 

//Integral_table* ND_ADC3_CAP_matrix::int_tab = 0;
SCF_data_reader* ND_ADC3_CAP_matrix::phis_ = 0;
Triangular_matrix<double>* ND_ADC3_CAP_matrix::dens_mat = 0;
extern Integral_table* integral_table;

int ND_ADC3_CAP_matrix::accept_analyzer( ADC_analyzer& adc_an)
{
  return adc_an.analyze(*this);
}


ND_ADC3_CAP_matrix::ND_ADC3_CAP_matrix(SCF_data_reader& phis, Triangular_matrix<double>* mat, 
				       unsigned int sym, unsigned int mode) 
  : sym_(sym), count_mult(0) 
{

  cap_mat = new Triangular_matrix<double>(mat->rows());
  for(int i = 0; i < cap_mat->rows(); i++)
    for(int j = i; j < cap_mat->rows(); j++)
      (*cap_mat)(i,j) = (*mat)(i,j);
  
  phis_ = &phis;


  // Allocate and initialize the structures required by Joerg's code
  inp = new inp_t;
  scf = new scf_t;
  symtab = new symtab_t;

  //Initialize some of the inp fields
  inp->sym = sym_+1;
  inp->debug = 0;
  inp->affinity_mode = 0;
  
  scf_data(scf);
  make_symtab(inp->sym, 0, 0, 0, scf, symtab);
  
  dim_1h = symtab->nOcc[inp->sym];
  dim_2h1p = calc_dim_2h1p(inp,symtab);

  // Load the integrals
  //  if (!int_tab) int_tab = new Integral_table(*phis_, OOOO|OOVO|VOVO|OOVV|VOVV|VVVV);
  

  // Load the density
  if (!dens_mat) {
    Blas_matrix rho(phis_->number_orbitals(), phis_->number_orbitals());
    rho = 0.;
    Self_energy(*integral_table,phis).density_3plus(rho);
    //Self_energy(*(adapt.int_tab),phis).density_3(rho);
    dens_mat = new Triangular_matrix<double>(phis_->number_orbitals());
    for (int i = 0; i < phis_->number_orbitals(); i++)
      for (int j = 0; j <= i; j++)
	(*dens_mat)(i,j) = rho(i,j);
  }
//   cout << "densmat\n";
//   for (int i= 0; i < dens_mat->rows(); i++)
//     for(int j = 0; j <= i; j++)
//       cout << i << " " << j << " "<<(*dens_mat)(i,j) << endl;



  //Initialize the dimensions of the matrix
  add_configs();
  
  // Intialize the expectation value
  d_null = d_null_null = 0.;
  double* cap_elements = &cap_mat->operator()(0,0);
  double* rho_0        = &dens_mat->operator()(0,0);
  calc_d_null_null(symtab,cap_elements,&d_null_null,inp->debug);
  d_null_null *= 2.; // Get the correct spin free value
  cap_calc_d_null(scf, rho_0, cap_elements,&d_null, inp->debug);
  d_null *= 2.;// D(2) must not include D(0) if the density is full(i.e. +0th order)!!!
  //cout << "d_null: " << d_null << endl;

}


void ND_ADC3_CAP_matrix::add_configs()
{
  // Initialize the main configurations
  
  ostringstream oss;
  for(int i= 0; i < dim_1h; i++) {
    oss << '<' << symtab->occ[inp->sym][i] << '|';
    configs_.push_back(oss.str());
    oss.str("");
    oss.clear();
  }

  // Initialize the satellite configurations
  int k_sym,l_sym,a_sym,kI_sym; // variables needed 
  int *k,*l,*a;                 // by the FOR_ALL_2H1P_AKL macro, see adc_macros.h
  
  int i = dim_1h;
  FOR_ALL_2H1P_AKL({
      oss << '<' << *k << ',' << *l <<  ',' << *a << ",I|";
      configs_.push_back(oss.str()); oss.str(""); oss.clear();
      i++;
      oss << '<' << *k << ',' << *l <<  ',' << *a << ",II|";
      configs_.push_back(oss.str()); oss.str(""); oss.clear();
      i++;
    },
    {
      oss << '<' << *k << ',' << *k << ',' << *a << '|';
      configs_.push_back(oss.str()); oss.str(""); oss.clear();
      i++;
    });

  
//    for(i = 0; i < configs_.size(); i++)
//      cout << i << ':' << configs_[i] << endl;
  

}

// Initialize the scf structure as if Joerg has done it,
// see read_scf_data.c in the original code
void ND_ADC3_CAP_matrix::scf_data(scf_t* scf)
{
  scf->nSym = phis_->number_irreps();
  scf->nBas = phis_->number_orbitals();

  /* read list of active orbitals */
  int n=scf->nBas;
  scf->loa=(int*)malloc((n+1)*sizeof(int));

  for(int i = 1; i <= n; i++)
    scf->loa[i] = i;

  /* read energies */
  n=scf->nBas;
  scf->epsi=(double*)malloc((n+1)*sizeof(double));

  for(int i = 1; i <= n; i++)
    scf->epsi[i] = phis_->energy(i-1);

  scf->occ=(double*)malloc((n+1)*sizeof(double));
  for(int i = 1; i <= n; i++)
    if (i <= phis_->number_occupied())
      scf->occ[i] = 2.;
    else
      scf->occ[i] = 0.;

  /* read symmetry information */
  n=scf->nBas;
  scf->sym=(int*)malloc((n+1)*sizeof(int));
  
  for(int i = 1; i <= n; i++)
    scf->sym[i] = phis_->irrep(i-1)+1;
  
}


ND_ADC3_CAP_matrix::~ND_ADC3_CAP_matrix()
{
  delete inp;
  for(int i=1;i<=symtab->nSym;i++) {
    free(symtab->occ[i]);
    free(symtab->vir[i]);
  }
  free(scf->loa); free(scf->epsi); free(scf->occ); free(scf->sym);
  delete scf;
  delete symtab;
  delete cap_mat;
}


// Builds the CAP matrix
void ND_ADC3_CAP_matrix::build_matrix(Blas_matrix& mat)
{

  mat.allocate(dim_1h + dim_2h1p, dim_1h + dim_2h1p);
  mat = 0.;

  build_main_block(mat);
  build_main_sat_block(mat);
  build_sat_block(mat);
  
}



void ND_ADC3_CAP_matrix::build_main_block(Blas_matrix& mat)
{
  double *d11=0;
  double* rho_0 = &dens_mat->operator()(0,0);
  double* cap_elements = &cap_mat->operator()(0,0);
  d11 = new double [dim_1h*(dim_1h+1)/2];
  //make_matrix(dim_1h,dim_1h,'p','z',&d11);
  for (int i = 0; i < dim_1h*(dim_1h+1)/2; i++) d11[i] = 0.;
  cap_calc_d11(inp,scf,symtab,&d_null,d11,rho_0,cap_elements);
  
  for(unsigned col = 0; col < dim_1h; col++) {
    Rectangular_proxymatrix main_row(d11+col*(col+1)/2, 1, col+1);
    mat(col,1,0,col+1).daxpy(1., main_row);
  }
  
  delete []d11; 

}

void ND_ADC3_CAP_matrix::build_main_sat_block(Blas_matrix& mat)
{
  double *d12=0;
  double* cap_elements = &cap_mat->operator()(0,0);
  
  // Build the upper-triangular main/sat block,
  d12 = new double [dim_1h*dim_2h1p];
  //make_matrix(dim_1h,dim_2h1p,'g','z',&d12);
  for (int i = 0; i < dim_1h*dim_2h1p; i++) d12[i] = 0.;
  cap_calc_d12(inp,scf,symtab,d12, cap_elements);//my_calc_d12.;l
  
  // but we need the lower triangle, which is transposed
  Transposed_proxymatrix trans(d12,  dim_1h, dim_2h1p);
  mat(dim_1h,dim_2h1p,0,dim_1h).daxpy(1.,trans);
  delete [] d12;
}

void ND_ADC3_CAP_matrix::build_sat_block(Blas_matrix& mat)
{
  
  double* cap_elements = &cap_mat->operator()(0,0);

  // Build the off diagonal elements
  double* d_col0=new double[dim_2h1p+1];
  double* d_col1=new double[dim_2h1p+1];
  
  int k_sym,l_sym,a_sym,kI_sym; // variables needed by
  int *k,*l,*a;                 // the FOR_ALL_2H1P_AKL macro
  
  int row = 0;
  // Note here that the diagonal elements are not included
  FOR_ALL_2H1P_AKL({
      /* akl part of 2h1p loop (k!=l) */
      cap_calc_d22_off_cols(inp,scf,symtab,k,l,a,row,d_col0,d_col1,cap_elements);
      // Add type I elements
      Rectangular_proxymatrix sat_row1(d_col0,1,row);

      mat(dim_1h+row,1,dim_1h,row).daxpy(1., sat_row1);
      row++;
      // Add type II elements
      Rectangular_proxymatrix sat_row2(d_col1,1,row);
      mat(dim_1h+row,1,dim_1h,row).daxpy(1., sat_row2);
      row++;
    },
    {
      /* akk part of 2h1p loop */
      cap_calc_d22_off_cols(inp,scf,symtab,k,l,a,row,d_col0,d_col1,cap_elements);
      Rectangular_proxymatrix sat_row(d_col0,1,row);
      mat(dim_1h+row,1,dim_1h,row).daxpy(1., sat_row);
      row++;
    });
  
  delete [] d_col0;
  delete [] d_col1;

  // Build the diagonal elements
  Blas_matrix d22_diag(dim_2h1p);
  d22_diag = 0.;

  cap_calc_d22_diag(inp,scf,symtab,&d_null_null,&d22_diag(0,0), cap_elements);
  mat(dim_1h,dim_2h1p,dim_1h,dim_2h1p).add_diag(d22_diag);
}


