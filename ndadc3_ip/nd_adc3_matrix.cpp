#include "nd_adc3_matrix.hpp"
#include "analysis/adc_analyzer.hpp"
#include "integral_table.hpp"
#include "integral_blocks.hpp"
#include "scf_data/scf_data_reader.hpp"

#include "id_self_energy.hpp"

#include <cstdlib>
#include <sstream>

#include <iostream>
#include <ctime>

extern "C" {
#include "globals.h"
} 

using namespace std;

//******************************************
//#define THIRD
//#undef THIRD
/* Introduced value 2 for mode - not quite meaningful,
   but the easiest way to use the same class for an
   adc2ip extended calculation.
*/
//******************************************

// The file contains the implementation of the nd_adc3_matrix class methods. 
// These are based on Joerg's nd_adc3 functions called here to build
// parts of ND ADC3 matrix. 

// Integral_table* ND_ADC3_matrix::int_tab = 0;
//SCF_data_reader* ND_ADC3_matrix::phis_ = 0;
Blas_matrix* ND_ADC3_matrix::sigma = 0;
extern Integral_table* integral_table;
 
int ND_ADC3_matrix::accept_analyzer( ADC_analyzer& adc_an)
{
  return adc_an.analyze(*this);
} 


ND_ADC3_matrix::ND_ADC3_matrix(SCF_data_reader& phis, unsigned int sym, unsigned spin, unsigned int mode) 
  : phis_(&phis), sym_(sym), count_mult(0), dim_1h(0), dim_2h1p(0)
{
  if (mode == 2) third_order = 0;
  else third_order = 1;

  // Allocate and initialize the structures required by Joerg's code
  inp = new inp_t;
  scf = new scf_t;
  symtab = new symtab_t;

  //Initialize some of the inp fields
  inp->sym = sym_+1;
  inp->debug = 0;
  if (mode != 2) inp->affinity_mode = mode; // AFFINITY MODE HASN'T BEEN TESTED!
  else inp->affinity_mode = 0;

  scf_data(scf);
  make_symtab(inp->sym, inp->affinity_mode, 0, 0, scf, symtab);

  dim_1h = symtab->nOcc[inp->sym];
  dim_2h1p = calc_dim_2h1p(inp,symtab);

//   // Load the integrals
//   if (!int_tab) int_tab = new Integral_table(*phis_, OOOO|OOVO|VOVO|OOVV|VOVV|VVVV);

  //Initialize the dimensions of the matrix
  add_configs();

  if (third_order) {
//#ifdef THIRD
      if (!sigma) {
        sigma = new Blas_matrix(phis_->number_orbitals(), phis_->number_orbitals()); 

    
        *sigma = 0.;

        time_t t1, t2;
        time(&t1);
        cout << "Computing static self-energy, 4+: ";
        flush(cout);
        Self_energy(*integral_table, *phis_).static_selfenergy_4plus(*sigma);
        //idSelf_energy(*(adapt.int_tab), *phis_).static_selfenergy_4plus(*sigma);
        //Self_energy(*integral_table, *phis_).static_selfenergy_3(*sigma);
        time(&t2);
        cout << int (t2 - t1)<< "sec" << endl;
        //*sigma *=27.21139;sigma->print();
        //throw;
  
      }
// #endif
 }

  sigma_ = new Blas_matrix(dim_1h,dim_1h);
  //  cout << "dim1h:"<< dim_1h << endl;
  *sigma_ = 0.;
  if (third_order) {
// #ifdef THIRD 
      for(int i_= 0; i_ < dim_1h; i_++)
        for(int j_= 0; j_ < dim_1h; j_++) {
          unsigned i = symtab->occ[inp->sym][i_] - 1;
          unsigned j = symtab->occ[inp->sym][j_] - 1;
          (*sigma_)(i_,j_) = (*sigma)(i,j);
        }
//#endif
  }
  if (inp->affinity_mode)
    *sigma_ *= -1.;


}


void ND_ADC3_matrix::add_configs()
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
void ND_ADC3_matrix::scf_data(scf_t* scf)
{
  scf->nSym = phis_->number_irreps();
  scf->nBas = phis_->number_orbitals();

  /* read list of active orbitals */
  int n=scf->nBas;
  scf->loa= new int [n+1];

  for(int i = 1; i <= n; i++)
    scf->loa[i] = i;

  /* read energies */
  n=scf->nBas;
  scf->epsi= new double [n+1];

  for(int i = 1; i <= n; i++)
    scf->epsi[i] = phis_->energy(i-1);

  scf->occ= new double [n+1];
  for(int i = 1; i <= n; i++)
    if (i <= phis_->number_occupied())
      scf->occ[i] = 2.;
    else
      scf->occ[i] = 0.;

  /* read symmetry information */
  n=scf->nBas;
  scf->sym=new int [n+1];
  
  for(int i = 1; i <= n; i++)
    scf->sym[i] = phis_->irrep(i-1)+1;
  
}


ND_ADC3_matrix::~ND_ADC3_matrix()
{
  delete inp;
  for(int i=1;i<=symtab->nSym;i++) { 
    delete [] symtab->occ[i];
    delete [] symtab->vir[i];
  }
  delete [] scf->loa; delete [] scf->epsi; 
  delete [] scf->occ; delete [] scf->sym;
  delete scf;
  delete symtab;
  delete sigma_;
}


// Builds the whole matrix
void ND_ADC3_matrix::build_matrix(Blas_matrix& mat)
{
  mat.allocate(dim_1h + dim_2h1p, dim_1h + dim_2h1p);
  mat = 0.;

  build_main_block(mat);
  build_main_sat_block(mat);
  build_sat_block(mat);

}



void ND_ADC3_matrix::build_main_block(Blas_matrix& mat)
{
  
  double *c11=0,*f_matrix=0;
  c11 = new double[dim_1h*(dim_1h+1)/2];
  f_matrix = new double [dim_1h*dim_1h];
  for (int i = 0; i < dim_1h*(dim_1h+1)/2; i++) c11[i] = 0.;
  for (int i = 0; i < dim_1h*dim_1h; i++) f_matrix[i] = 0.;
  //make_matrix(dim_1h,dim_1h,'p','z',&c11);
  //make_matrix(dim_1h,dim_1h,'g','z',&f_matrix);  
  calc_k1(inp,scf,symtab,c11,f_matrix);
  calc_c11_2(inp,scf,symtab,c11,f_matrix);
  if (third_order) {
//#ifdef THIRD
     calc_c11_3(inp,scf,symtab,c11,f_matrix);
//#endif
  }
  for(unsigned col = 0; col < dim_1h; col++) {
    Rectangular_proxymatrix main_row(c11+col*(col+1)/2, 1, col+1);
    mat(col,1,0,col+1).daxpy(1., main_row);
  }

  mat(0, dim_1h, 0, dim_1h).daxpy(-1., *sigma_);

  delete [] c11; delete [] f_matrix;
}

void ND_ADC3_matrix::build_main_sat_block(Blas_matrix& mat)
{
  double *c12=0;
  // Build the upper-triangular main/sat block,
  c12 = new double [dim_1h*dim_2h1p];
  for (int i = 0; i < dim_1h*dim_2h1p; i++) c12[i] = 0.;
  //make_matrix(dim_1h,dim_2h1p,'g','z',&c12);
  calc_c12_1(inp,scf,symtab,c12);
  if (third_order) { 
//#ifdef THIRD
      calc_c12_2(inp,scf,symtab,c12);
//#endif
  }
  // but we need the lower triangle, which is transposed
  Transposed_proxymatrix trans(c12,  dim_1h, dim_2h1p);
  // a sign error in Joerg's code
  mat(dim_1h,dim_2h1p,0,dim_1h).daxpy(-1.,trans);

  //cout << "susu --> the matrix: only main and coupling" << endl;
  //mat.print();



  delete [] c12;
}

void ND_ADC3_matrix::build_sat_block(Blas_matrix& mat)
{
  // Build the off diagonal elements
  double* c_col0=new double[dim_2h1p+1];
  double* c_col1=new double[dim_2h1p+1];
  
  int k_sym,l_sym,a_sym,kI_sym; // variables needed by
  int *k,*l,*a;                 // the FOR_ALL_2H1P_AKL macro
  
  int row = 0;
  // Note here that the diagonal elements are not included
  FOR_ALL_2H1P_AKL({
      /* akl part of 2h1p loop (k!=l) */
      calc_c22_1_cols(inp,scf,symtab,k,l,a,row,c_col0,c_col1);
      // Add type I elements
      Rectangular_proxymatrix sat_row1(c_col0,1,row);
      mat(dim_1h+row,1,dim_1h,row).daxpy(1., sat_row1);
      row++;
      // Add type II elements
      Rectangular_proxymatrix sat_row2(c_col1,1,row);
      mat(dim_1h+row,1,dim_1h,row).daxpy(1., sat_row2);
      row++;
    },
    {
      /* akk part of 2h1p loop */
      calc_c22_1_cols(inp,scf,symtab,k,l,a,row,c_col0,c_col1);
      Rectangular_proxymatrix sat_row(c_col0,1,row);
      mat(dim_1h+row,1,dim_1h,row).daxpy(1., sat_row);
      row++;
    });
  
  delete [] c_col0;
  delete [] c_col1;

  // Build the diagonal elements
  Blas_matrix c22_diag(dim_2h1p);
  c22_diag = 0.;

  calc_k2(inp,scf,symtab, &c22_diag(0,0));
  //c22_diag = -1.;
  calc_c22_1_diag(inp,scf,symtab, &c22_diag(0,0));

  mat(dim_1h,dim_2h1p,dim_1h,dim_2h1p).add_diag(c22_diag);
}


int ND_ADC3_matrix::operator()(Blas_matrix** block_out, Blas_matrix** block_in, int count) 
{
  // Emulate Lanczos subspace iterations, see adc_analyzer.cpp
 if (count_mult < 2*dim_1h) {
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



void ND_ADC3_matrix::multiply_main_block(Blas_matrix** block_out,  Blas_matrix** block_in, int count) 
{
  double *c11=0,*f_matrix=0;
  Blas_matrix main_block(dim_1h, dim_1h);
  main_block =  0.;
  c11 = new double[dim_1h*(dim_1h+1)/2];
  f_matrix = new double [dim_1h*dim_1h];
  for (int i = 0; i < dim_1h*(dim_1h+1)/2; i++) c11[i] = 0.;
  for (int i = 0; i < dim_1h*dim_1h; i++) f_matrix[i] = 0.;
//  make_matrix(dim_1h,dim_1h,'p','z',&c11);
//  make_matrix(dim_1h,dim_1h,'g','z',&f_matrix);  
  calc_k1(inp,scf,symtab,c11,f_matrix);
  calc_c11_2(inp,scf,symtab,c11,f_matrix);
  if (third_order) {
//#ifdef THIRD
      calc_c11_3(inp,scf,symtab,c11,f_matrix);
//#endif
  }
  for(unsigned col = 0; col < dim_1h; col++) {
    Rectangular_proxymatrix main_row(c11+col*(col+1)/2, 1, col+1);
    main_block(col,1,0,col+1).daxpy(1., main_row);
    
    Transposed_proxymatrix main_col(c11+col*(col+1)/2, 1, col);
    main_block(0,col,col,1).daxpy(1., main_col);
    
  }
//  free(c11); free(f_matrix);
  delete [] c11; delete [] f_matrix;
  main_block.daxpy(-1., *sigma_);
  
  Blas_matrix** vec_in = block_in; 
  Blas_matrix** vec_out = block_out;
  unsigned block_rows = main_block.rows();
  unsigned block_cols = main_block.cols();
  for (int i = 0; i < count; i++) {
    (**vec_out)(0, block_rows, 0, 1)
      .dgemv('N', main_block, (**vec_in)(0, block_cols, 0, 1));
    vec_in++; vec_out++;
  }
}

void ND_ADC3_matrix::multiply_main_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, 
					     int count) 
{

  Blas_matrix c12(dim_1h, dim_2h1p);
  c12 = 0.;

  calc_c12_1(inp,scf,symtab,&c12(0,0));
  if (third_order) { 
//#ifdef THIRD
      calc_c12_2(inp,scf,symtab,&c12(0,0));
//#endif
  }
  // a sign error in Joerg's code
  c12 *= -1.;
  Blas_matrix** vec_in = block_in; 
  Blas_matrix** vec_out = block_out;
  unsigned block_rows = c12.rows();
  unsigned block_cols = c12.cols();
  
  for (int i = 0; i < count; i++) {
    (**vec_out)(dim_1h, dim_2h1p, 0, 1)
      .dgemv('T', c12, (**vec_in)(0, dim_1h, 0, 1));
    (**vec_out)(0, dim_1h, 0, 1)
      .dgemv('N', c12, (**vec_in)(dim_1h, dim_2h1p, 0, 1));
    vec_out++; vec_in++;
  }
  
}

void ND_ADC3_matrix::multiply_sat_block(Blas_matrix** block_out,  Blas_matrix** block_in, 
					int count) 
{
  
  // Multiply with the off diagonal elements first
  Blas_matrix c_col0(dim_2h1p+1);
  Blas_matrix c_col1(dim_2h1p+1);

  int k_sym,l_sym,a_sym,kI_sym;
  int *k,*l,*a;
  
  int row = 0;
  // Note here that the diagonal elements are not included
  FOR_ALL_2H1P_AKL({
      /* akl part of 2h1p loop (k!=l) */
      
      calc_c22_1_cols(inp,scf,symtab,k,l,a,row,&c_col0(0,0),&c_col1(0,0));
      
      multiply_sat_row(block_out,block_in,count, c_col0(0,row,0,1), row++);
      
      multiply_sat_row(block_out,block_in,count, c_col1(0,row,0,1), row++);

    },
    {
      /* akk part of 2h1p loop */
      calc_c22_1_cols(inp,scf,symtab,k,l,a,row,&c_col0(0,0),&c_col1(0,0));
      
      multiply_sat_row(block_out,block_in,count, c_col0(0,row,0,1), row++);

     });
  

  // Multiply with the  diagonal elements
  Blas_matrix c22_diag(dim_2h1p);
  c22_diag = 0.;
  calc_k2(inp,scf,symtab,&c22_diag(0,0));
  calc_c22_1_diag(inp,scf,symtab,&c22_diag(0,0));
  
  for(int col = dim_1h; col < dim_1h+dim_2h1p; col++) {
    Blas_matrix** vec_in = block_in, **vec_out = block_out;
    
    for(int vec = 0; vec < count; vec++) {
      (**vec_out)(col) += c22_diag(col-dim_1h) * (**vec_in)(col);
      vec_out++, vec_in++;
    }
  }
}

void ND_ADC3_matrix::multiply_sat_row(Blas_matrix** block_out,  Blas_matrix** block_in,  int count,
				      Submatrix row, unsigned block_rows) 
{
  Blas_matrix** vec_in = block_in, **vec_out = block_out;
  for(int vec = 0; vec < count; vec++) {
    (**vec_out)(dim_1h, block_rows, 0 , 1)
      .daxpy((**vec_in)(dim_1h+block_rows), row);
    
    (**vec_out)(dim_1h+block_rows) += 
      row
      * (**vec_in)(dim_1h, block_rows, 0, 1);
    vec_out++; vec_in++;
  }
}
