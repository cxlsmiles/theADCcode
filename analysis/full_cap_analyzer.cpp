#include "full_cap_analyzer.hpp"
#include "adc2_prop/adc2_cap_matrix.hpp"
#include "ndadc3_prop/nd_adc3_cap_matrix.hpp"
#include "ndadc3_ip/nd_adc3_matrix.hpp"
#include "adc2_dip/adc2_matrix.hpp"
#include "../blas_matrix.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "../input_struct.hpp"
#include "../matrices.hpp"


#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include <cmath>

//23 Feb 2012

struct Complex {
   double a;
   double b;
   Complex() : a(0.0), b(0.0) {}
   Complex(double real, double imag) : a(real), b(imag) {}
   double real() const { return a; }
   bool operator() (const Complex &c1, const Complex &c2) { return (c1.real() < c2.real()); }
};

// END

extern "C" {
#if defined INTEL
#include "mkl_lapack.h"
#elif defined PGI
#include "acml.h"
#elif defined GCC
struct _dcomplex { double real, imag; };
typedef struct _dcomplex dcomplex;
void zgeev_( char* jobvl, char* jobvr, int* n, dcomplex* a,
                int* lda, dcomplex* w, dcomplex* vl, int* ldvl, dcomplex* vr, int* ldvr,
                dcomplex* work, int* lwork, double* rwork, int* info );
#endif
}

using namespace std;

#define LANCZOS_LIMIT 1
#if defined INTEL
#define COMPLEX MKL_Complex16
#elif defined PGI
#define COMPLEX doublecomplex
#elif defined GCC
#define COMPLEX dcomplex
#endif


const static double au2eV = 27.211396;


// This method builds full-CAP trajectories of a given propagator. 
// It obtains the MO CAP matrix from phis, then builds the FULL propagator and property matrices.
// Depending on the CAP settings, it builds a complex matrix out of the two full matrices and diagonalizes
// it at each CAP step.

int Full_CAP_analyzer::analyze( ADC_matrix& adc_mat)
{ 

  // Get the cap matrix
  int dim = phis_.number_orbitals();
  Triangular_matrix<double>* mat;
  mat = new Triangular_matrix<double>(dim);
  Blas_matrix capmo(dim, dim);
  phis_.get_mo_cap(&cap_.boxx, &cap_.boxy, &cap_.boxz, &dim, &capmo(0,0));
  if (dim < 0)
    throw string (" Full_CAP_analyzer:: error reading CAP matrix.");
  for(int i = 0; i < dim; i++)
    for(int j = 0;  j < dim; j++)
      (*mat)(i,j) = capmo(i,j);

//cout << "CAP MO\n";
//capmo.print();

  ADC_analyzer::analyze(adc_mat);
  Blas_matrix real_part, imag_part; 
  
  // TODO: find a cleaner way to do this:
  // Determine what type of ISR representation we need (N-1) or (N-2):
  Adc2_matrix* is_adc2 = dynamic_cast<Adc2_matrix*>(&adc_mat);
  ND_ADC3_matrix* is_ndadc3 = dynamic_cast<ND_ADC3_matrix*>(&adc_mat);
  ADC_matrix* cap;
  if (is_adc2){
    cap = new Adc2_CAP_matrix(phis_, tab_, mat, adc_mat.symmetry(), adc_mat.spin());
  } else if (is_ndadc3){
    cap = new ND_ADC3_CAP_matrix(phis_, mat, adc_mat.symmetry(), 0);
  } else {
    cout << "Full_CAP_analyzer:Requested an ISR representation of CAP that is not yet implemented.";
    return 1;
  }
  

  // Build the lower triangles of the matrices
  cap->build_matrix(imag_part);
  delete mat;
  delete cap;
  adc_mat.build_matrix(real_part);
   cout << "Print imag part\n";
   imag_part.print();
   cout << "Print real part\n";
   real_part.print();
  
  // copy the lower triangle to the upper
  for(unsigned row = 0; row < real_part.rows(); row++)
    for(unsigned col = 0; col < row; col++) {
      real_part(col,row) = real_part(row,col);
      imag_part(col,row) = imag_part(row,col);
    }
  
  Rectangular_matrix<COMPLEX> CAP_mat(real_part.rows(), real_part.cols());
  
  Rectangular_matrix<COMPLEX> w(real_part.rows(), 1), vl(1,1), vr(1,1);
  Rectangular_matrix<COMPLEX> wp(real_part.rows(), 1);
  Rectangular_matrix<COMPLEX> wm(real_part.rows(), 1);
  int info = 0;
  int rows = CAP_mat.rows();
  int lwork = 2*rows;
  Rectangular_matrix<COMPLEX> work(lwork,1);
  Rectangular_matrix<double>  rwork(lwork,1);
  int one = 1;
  
  // THE CAP loop

  int iter = cap_.nicap;
  int ioffset = cap_.ioffset;
  double deltaeta = cap_.deltaeta;

  ofstream outfile(cap_.output);
  cout << " Building CAP trajectories.\n";

  outfile << "#Real en.  Imag. en.  Real+ en.  Imag.+ en.  Real- en.  Imag.- en. Eta" << endl;
  for(unsigned row = w.rows(); row > 0; row--) {
    wp(row-1,0).real=0.0;
    wp(row-1,0).imag=0.0;
    wm(row-1,0).real=0.0;
    wm(row-1,0).imag=0.0;
  }

  for (int count = ioffset; count < iter + ioffset; count++) {
    time_t t1,t2;
    time(&t1);
    double eta_value;
    double eta = cap_.sicap;
    // Build the imag part each iteration
    for(unsigned row = 0; row < CAP_mat.rows(); row++ ) {
      for(unsigned col = 0; col < CAP_mat.cols(); col++){
        CAP_mat(row,col).real = real_part(row,col) + cap_.srcap*imag_part(row,col);

        switch (cap_.increase) {
        case 1:                   // add real part if requested     //exponential
          CAP_mat(row,col).imag = -(eta*count*(count+1)) * imag_part(row,col);
          eta_value = eta*count*(count + 1); break;
        case 2:                                                     // linear
          CAP_mat(row,col).imag = -(eta*count) * imag_part(row,col);
          eta_value = eta*count; break;
        case 3:                                                     // power
          CAP_mat(row,col).imag = -eta*(pow(1.1, count)-1.0)/0.1 * imag_part(row,col);
          eta_value = eta*(pow(1.1, count) - 1.0)/0.1; break;
        default:;
        }
      }
    }

#if defined INTEL
    zgeev_("N", "N", &rows, &CAP_mat(0,0), &rows, &w(0,0), &vl(0,0), &one,
         &vr(0,0), &one, &work(0,0), &lwork, &rwork(0,0), &info);
#elif defined PGI
    zgeev('N', 'N', rows, &CAP_mat(0,0), rows, &w(0,0), &vl(0,0), one,
	      &vr(0,0), one, &info);
#elif defined GCC
    zgeev_("N", "N", &rows, &CAP_mat(0,0), &rows, &w(0,0), &vl(0,0), &one,
         &vr(0,0), &one, &work(0,0), &lwork, &rwork(0,0), &info);
#endif
    time(&t2);
    cout << "step: "<< count + 1<< ", time: " << t2 - t1 << " s" <<endl;

    if (deltaeta){
	    time(&t1);
	    // Build the imag part each iteration
	    for(unsigned row = 0; row < CAP_mat.rows(); row++) {
	      for(unsigned col = 0; col < CAP_mat.cols(); col++){
		CAP_mat(row,col).real = real_part(row,col) + cap_.srcap*imag_part(row,col);

		switch (cap_.increase) {
		case 1:                   // add real part if requested     //exponential
		  CAP_mat(row,col).imag = -(eta*count*(count+1) + deltaeta) * imag_part(row,col); break;
		case 2:                                                     // linear
		  CAP_mat(row,col).imag = -(eta*count + deltaeta) * imag_part(row,col); break;
		case 3:                                                     // power
		  CAP_mat(row,col).imag = -(eta*(pow(1.1, count)-1.0)/0.1 + deltaeta) * imag_part(row,col); break;
		default:;
	      }
            }
           }
            
	#if defined INTEL
	    zgeev_("N", "N", &rows, &CAP_mat(0,0), &rows, &wp(0,0), &vl(0,0), &one,
		 &vr(0,0), &one, &work(0,0), &lwork, &rwork(0,0), &info);
	#elif defined PGI
	    zgeev('N', 'N', rows, &CAP_mat(0,0), rows, &wp(0,0), &vl(0,0), one,
		      &vr(0,0), one, &info);
	#elif defined GCC
	    zgeev_("N", "N", &rows, &CAP_mat(0,0), &rows, &wp(0,0), &vl(0,0), &one,
		 &vr(0,0), &one, &work(0,0), &lwork, &rwork(0,0), &info);
	#endif

	    time(&t2);
	    cout << "step+: "<< count + 1<< ", time: " << t2 - t1 << " s" <<endl;

	    time(&t1);
	    // Build the imag part each iteration
	    for(unsigned row = 0; row < CAP_mat.rows(); row++) {
	      for(unsigned col = 0; col < CAP_mat.cols(); col++){
		CAP_mat(row,col).real = real_part(row,col) + cap_.srcap*imag_part(row,col);

		switch (cap_.increase) {
		case 1:                   // add real part if requested     //exponential
		  CAP_mat(row,col).imag = -(eta*count*(count+1) - deltaeta) * imag_part(row,col); break;
		case 2:                                                     // linear
		  CAP_mat(row,col).imag = -(eta*count - deltaeta) * imag_part(row,col); break;
		case 3:                                                     // power
		  CAP_mat(row,col).imag = -(eta*(pow(1.1, count)-1.0)/0.1 - deltaeta) * imag_part(row,col); break;
		default:;
	      }
	  }
         }


	#if defined INTEL
	    zgeev_("N", "N", &rows, &CAP_mat(0,0), &rows, &wm(0,0), &vl(0,0), &one,
		 &vr(0,0), &one, &work(0,0), &lwork, &rwork(0,0), &info);
	#elif defined PGI
	    zgeev('N', 'N', rows, &CAP_mat(0,0), rows, &wm(0,0), &vl(0,0), one,
		      &vr(0,0), one, &info);
	#elif defined GCC
	    zgeev_("N", "N", &rows, &CAP_mat(0,0), &rows, &wm(0,0), &vl(0,0), &one,
		 &vr(0,0), &one, &work(0,0), &lwork, &rwork(0,0), &info);
	#endif
    
	    time(&t2);
	    cout << "step-: "<< count + 1<< ", time: " << t2 - t1 << " s" <<endl;
	}
    //Sort the elements of w, wp, wm by ascending order in order to be able to obtain the derivative FEB 23, 2012
    vector<Complex> w_values, wm_values, wp_values;
    struct Complex condition ;
    
    for (int r = w.rows(); r > 0; --r)
    {
        w_values.push_back(Complex(w(r-1,0).real, w(r-1,0).imag));
        wm_values.push_back(Complex(wm(r-1,0).real, wm(r-1,0).imag));
        wp_values.push_back(Complex(wp(r-1,0).real, wp(r-1,0).imag));
    };
    sort(w_values.begin(), w_values.end(), condition);
    sort(wm_values.begin(), wm_values.end(), condition);
    sort(wp_values.begin(), wp_values.end(), condition);

    // Print the results into a file
// Feb 24 2012
    for(int i = 0; i < w_values.size(); ++i)
    {
      outfile << left << setprecision(15) << setw(20) << w_values[i].a * au2eV << setw(23) << w_values[i].b * au2eV << ' ';
      outfile << left << setprecision(15) << setw(20) << wm_values[i].a * au2eV << setw(23) << wm_values[i].b * au2eV << ' ';
      outfile << left << setprecision(15) << setw(20) << wp_values[i].a * au2eV << setw(23) << wp_values[i].b * au2eV << ' ';
      outfile << left << setprecision(10) << setw(20) << eta_value << endl;
    }
// END

  }


  outfile.close();
  return 0;
}

Full_CAP_analyzer::Full_CAP_analyzer(SCF_data_reader& phis, Integral_table& tab, struct Diag_info& diag, struct Eigen_info& eig, struct Cap_info& cap)
  : phis_(phis), tab_(tab), ADC_analyzer(diag, eig), cap_(cap) 
{

//  if (!cap_.deltaeta) throw string("Full_CAP_analyzer:: delta eta for the derivative not selected.");
  if (!cap_.sicap) throw string("Full_CAP_analyzer:: step size of CAP not selected.");
  if (!cap_.nicap) throw string("Full_CAP_analyzer:: iteration number of CAP not selected.");
  if (cap_.ioffset < 0) throw string("Full_CAP_analyzer:: offset of the CAP iterations not selected.");
  if (!cap_.boxx || !cap_.boxy || !cap_.boxz) throw string("Full_CAP_analyzer:: box size of CAP not selected.");

}
