#include "isr_dipole_analyzer.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "../blas_matrix.hpp"
#include "ndadc3_prop/nd_adc3_cap_matrix.hpp"
#include "adc2_prop/adc2_cap_matrix.hpp"
#include "adc2_dip/adc2_matrix.hpp"
#include "ndadc3_ip/nd_adc3_matrix.hpp"
#include "../matrices.hpp"
#include "ndadc3_ip/self_energy.hpp"
#include "../input_struct.hpp"

#include <iostream>

using namespace std;



int ISR_dipole_analyzer::analyze(ADC_matrix& adc_mat)
{

  cout << "Dipole analyzer\n";

  int n =  0;
  phis_.get_dip(0,0,0,&n);

  n *= -1;
  Blas_matrix x(n,n), y(n,n),z(n,n) ;

  phis_.get_dip(&x(0,0), &y(0,0), &z(0,0), &n);

  Triangular_matrix<double>* mat = new Triangular_matrix<double>(n);

  for(int i = 0; i < mat->rows(); i++) {
    for(int j = 0; j <= i; j++)
      (*mat)(j,i)  =  z(j,i);
  }

  cout <<" Z mat:\n";z.print();


  // TODO: find a cleaner way to do this:
  // Determine what type of ISR representation we need (N-1) or (N-2):
  Adc2_matrix* is_adc2 = dynamic_cast<Adc2_matrix*>(&adc_mat);
  ND_ADC3_matrix* is_ndadc3 = dynamic_cast<ND_ADC3_matrix*>(&adc_mat);
   
  ADC_matrix* cap;
  if (is_adc2){
    cap = new Adc2_CAP_matrix(phis_, tab_, mat, adc_mat.symmetry(), adc_mat.spin());
  } else if (is_ndadc3){
    cap = new ND_ADC3_CAP_matrix(phis_, mat, adc_mat.symmetry(), adc_mat.spin());
  } else {
    cout << "Full_CAP_analyzer:Requested an ISR representation of CAP that is not yet implemented.";
    return 1;
  }


  Blas_matrix isr_mat;
//   adc_mat.build_matrix(isr_mat);
  
//    cout << "ADC mat:\n"; 
//    isr_mat.print();



  
  Blas_matrix eig_vecs, eig_vals, norms;

  diagonalize(adc_mat,eig_vecs, eig_vals, norms);
  print_eigenvals(adc_mat, eig_vecs, eig_vals, norms);

  cap->build_matrix(isr_mat);
//   for(unsigned row = 0; row < isr_mat.rows(); row++)
//     for(unsigned col = 0; col < row; col++) 
//       isr_mat(col,row) = isr_mat(row,col);

cout << "ISR mat:\n";
isr_mat.print();


  int num_eigvals = eig_vals.size();
  for(int i = 0; i < num_eigvals; i++) {

    Blas_matrix vec(num_eigvals, 1, &eig_vecs(0,i));
    
    Blas_matrix tempvec(num_eigvals);

    tempvec = 0.;
    tempvec.dsymv(isr_mat,vec);
    

    double dip_mom =  vec * tempvec;
    double d_nuc =dip_.nucl ;
    cout << i+1 << "eval(au):" << eig_vals(i) << "\teval(eV): " << eig_vals(i)*27.2114<< " " 
	 << "\tDel(au):"<< dip_mom << "\tDtot(Debye):"<< (-dip_mom+d_nuc) * 2.54158025293807<< endl;
    
    //break;
  }

  delete cap;
  delete mat;
  return 0;

}
