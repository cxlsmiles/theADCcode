#include "adc_analyzer.hpp"
#include "../adc_matrix.hpp"
#include "../blas_matrix.hpp"
#include "../input_struct.hpp"

#include "./Lanczos/lanczos.h"

#include <ctime>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;
using namespace LanczosProjection;


#define LANCZOS_LIMIT 1
const static double au2eV = 27.211396;




// Prints the eigenvalues of the ADC_matrix. It also prints the sorted
// contributions of the main space configuration.
void ADC_analyzer::print_eigenvals(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals,  Blas_matrix& resnorms)
{  

  cout << "\n Eigenvalue (eV), ps (%), residue\n"
       << " ----------------------------------";
  
  bool has_printed = false;  
  for (unsigned i = 0; i < eig_vals.size(); i++) { // For each eigenvalue
 
    std::ostringstream pop_info;

    double last_pop_max = 1.1;

//        TSVETA May 30, 2012 --> the if statement is added
    if (adc_mat.main_block_size())
    {
        double ps = 0.;
        for(unsigned conf = 0; conf < adc_mat.main_block_size(); conf++) // Find the pole strength (the sum of the main space contributions)
          ps += eig_vecs(conf,i)*eig_vecs(conf,i);
    
        ps *= 100;
        if (ps < eig_.ps) {eig_vals(i) = 0.; continue;} // Check if the pole strength is above the selected threshold
    
        cout << endl << endl;
        cout << ' ' << i+1 << ": " << setiosflags(ios::fixed) << setprecision(6)
	     << eig_vals(i) * au2eV << ", " << setprecision(2) << ps << ", "
	     << setprecision(6) << resnorms(i);
    
    
        has_printed = true;

//        std::ostringstream pop_info;
    
//        double last_pop_max = 1.1;
    
        for(unsigned conf_count = 0; conf_count < adc_mat.main_block_size(); conf_count++) { // Sort the contributions
      
          double pop_max = 0.;
          unsigned conf = 0;
          //Sorts the population coefficients accordingly.
          for(unsigned conf_max = 0; conf_max < adc_mat.main_block_size(); conf_max++) 
	
	
	    if ((abs(pop_max) < abs(eig_vecs(conf_max, i))) && (last_pop_max > abs(eig_vecs(conf_max, i)))) {
	      pop_max = eig_vecs(conf_max, i);
	      conf = conf_max;
	    }
          if (abs(pop_max) < eig_.thresh) break;
          last_pop_max = abs(pop_max);
      
      
          if (!(conf_count % 6))
	    pop_info << endl;
          // Collect the info that is to be printed
          pop_info << setw(7) << adc_mat.get_conf(conf) << ':'  
	           << setw(9) << setiosflags(ios::fixed) << setprecision(6) << eig_vecs(conf,i); 
      
        }
    
        if (!pop_info.str().empty()){
          cout << "\n Overlaps with main-space configurations:";
          cout << pop_info.str();
        }
 
    } else {
        cout << "\n " << i+1 << ": " << setiosflags(ios::fixed) << setprecision(6)
	    << eig_vals(i) * au2eV;    
    }

    // Add the satellite states
    last_pop_max = 1.1;
    pop_info.str("");
    pop_info.clear();
    for(unsigned conf_count = adc_mat.main_block_size(); conf_count < adc_mat.size();  conf_count++) { // Sort the contributions
        // Diagonalized with Lanczos:
        if (eig_vecs.rows() < adc_mat.size()) break;

        double pop_max = 0.;
        unsigned conf = 0;
        //Sorts the population coefficients.
        for(unsigned conf_max = adc_mat.main_block_size(); conf_max < adc_mat.size(); conf_max++) 
	
	
    if ((abs(pop_max) < abs(eig_vecs(conf_max, i))) && (last_pop_max > abs(eig_vecs(conf_max, i)))) {
	    pop_max = eig_vecs(conf_max, i);
	    conf = conf_max;
    }
        if (abs(pop_max) < eig_.thresh) break;
        last_pop_max = abs(pop_max);
      
      
        if (!((conf_count - adc_mat.main_block_size()) % 4))
    pop_info << endl;
        // Collect the info that is to be printed
        pop_info << setw(17) << adc_mat.get_conf(conf) << ':'  
	        << setw(9) << setiosflags(ios::fixed) << setprecision(6) << eig_vecs(conf,i); 
      
    }
    
    if (!pop_info.str().empty()){
        cout << "\n Overlaps with satellite-space configurations:";
        cout << pop_info.str();
    }
  }
  
  cout << endl;
  
//  if (!has_printed) 
//    if (!eig_vecs.size())
//      cout << " No configurations in the main space.\n";
//    else
//      cout << " The pole strength of all roots is below " << eig_.ps << "%.\n";

  if (!has_printed && eig_vecs.size())
    cout << " The pole strength of all roots is below " << eig_.ps << "%.\n";
}



int ADC_analyzer::analyze( ADC_matrix& adc_mat)
{

  if (!adc_mat.size()) return 1;
  
  cout << "\n Computing spectrum for ";
  cout << "symmetry " << adc_mat.symmetry()+1 << ", spin " << adc_mat.spin()+1 << endl;
  cout << " Number of ISR configurations: " << adc_mat.size() << endl;  

// May 7, 2013
//Blas_matrix ADC_mat_print;
//adc_mat.build_matrix(ADC_mat_print);
//std::fstream ADC_mat_file;

//ADC_mat_file.open("ADC_Matrix.dat");
//for (int row = 0; row < ADC_mat_print.rows(); row++)
//{
//   for (int col = 0; col <  ADC_mat_print.cols(); col++)
//         ADC_mat_file <<  ADC_mat_print(col,row) << "  ";
//   ADC_mat_file << endl;

//}
//ADC_mat_file.close();

//ADC_mat_print.print();


// Tsveta

  Blas_matrix eig_vecs, eig_vals, resnorms;
  diagonalize(adc_mat, eig_vecs,  eig_vals,  resnorms);
  print_eigenvals(adc_mat, eig_vecs,  eig_vals,  resnorms);
  
  return 0;
  
}


void ADC_analyzer::diagonalize(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals,  Blas_matrix& resnorms)
{
  if (diag_.type == 1) {

    lapack_diagonalize(adc_mat, eig_vecs, eig_vals);
    resnorms.allocate(eig_vecs.rows(), 1);
    resnorms = 0.;
    
  } else {

      lanczos_diagonalize(adc_mat, eig_vecs, eig_vals, resnorms);
  }  

}



// Uses LAPACK to diagonalize the matrix. This involves the construction of the whole matrix
void ADC_analyzer::lapack_diagonalize(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals)
{

  cout << " Full diagonalizing:\n";

//       TSVETA May, 30 --> This if statement is commented in order for the matrix
//                          to be built and diagonalized even though the main block is zero
//
//  if (!adc_mat.main_block_size()) return;
  
  adc_mat.build_matrix(eig_vecs);
  
  if (eig_vecs.size())
    eig_vecs.dsyev('V', eig_vals);

}


// Uses Lanczos to diagonalize the matrix. 
// Subspace iteration is used including only the satellite part  (see Tarantelli's Chem. Phys. 329, 11).
// If not explicitly requested, only the mains-space vector components and the residuals are produced (requiring only partial diagonalizing
// of the band matrix).
void ADC_analyzer::lanczos_diagonalize(ADC_matrix& adc_mat, Blas_matrix& eig_vecs,  Blas_matrix& eig_vals,  Blas_matrix& resnorms)
{
  
  bool full_vecs;
  if (diag_.vecs.empty()) full_vecs = false;
  else full_vecs = true;
  
  if (full_vecs)
    cout << " Lanczos diagonalizing: ";
  else 
    cout << " Lanczos with fast band-matrix diagonalizing: ";
  
  
  //////////// Initialize the Lanczos diagonalizer
  
  int m_nroots = adc_mat.main_block_size();//adc_mat.main_block_size(); 
  if (!m_nroots) return;
  cout << "block size " << m_nroots << endl;
 
  int m_maxiter = diag_.iter;  
  
  time_t t1,t2;

  Blas_matrix **guess = new Blas_matrix*[m_nroots];
  for(int j = 0; j< m_nroots;j++) { // Build the initial guess (A set of Cartesian vectors)
    guess[j] = new Blas_matrix(adc_mat.size(), 1);
    *guess[j] = 0.;
    (*guess[j])(j) = 1.;
  }

  // Construct the Lanczos diagonalizer
  Lanczos <double,Blas_matrix,ADC_matrix> projection (adc_mat, CRIT_EIGENVALUE, 1e-6, false, !full_vecs);
  
  
  // Copy the guess vector pointers to the Lanczos object.  
  // These vectors will be overwritten during Lanczos iterations.
  projection.AddStartBlock(guess,m_nroots);

  /////////////// Do the Lanczos procedure (Build the band-diagonal matrix)
  int mat_dim = 0;
  // Do the Lanczos iterations!
  int iterI = 0;    
  do {      
    
    //boost::timer t;
    time(&t1);
    cout << " Iteration:" << iterI;
    flush(cout);
    ++iterI;

    mat_dim = projection.Iterate();

    time(&t2);
    cout << " time: " << t2 - t1 << " s" <<endl;
    flush(cout);
  } while(iterI <= m_maxiter);		      


  int lanc_space = mat_dim - m_nroots;
  cout << " Size of Lanczos space: " << lanc_space << endl;
  

  time(&t1);
  cout << " Diagonalizing band matrix, ";  flush(cout);
  // Diagonalize the band matrix
  projection.Diagonalize();
  time(&t2);
  cout << " time: " << t2 - t1 << " s" <<endl;
  
  ///////////////////// Get the data from the Lanczos diagonalizer
  
  time(&t1);
  // Get the eigenvectors
  Blas_matrix energies(lanc_space, 1, projection.GetEigenvals());
  // Get the residuals
  Blas_matrix resid(lanc_space, 1, projection.GetResidualNorms());

  int vec_size;
  if (full_vecs) vec_size = lanc_space;
  else vec_size = 2*m_nroots;
  // Get the eigenvectors of the band-diagonal matrix (These are equivalent with the Lanczos vectors in the short version).
  Blas_matrix vecs(vec_size, lanc_space, projection.GetEigenvecs());
  vector<unsigned> non_spurious;
  // Find the spurious vectors using eigenvectors of the band-diagonal matrix
  find_spurious(m_nroots, vecs, non_spurious);  
  if (!non_spurious.empty())
    cout << " Warning: " << lanc_space - non_spurious.size() << " spurious Lanczos roots will be dropped.\n";
  
  if (!full_vecs) {
    
    eig_vals.allocate(non_spurious.size());
    resnorms.allocate(non_spurious.size());
    eig_vecs.allocate(adc_mat.main_block_size(), non_spurious.size());
    
    eig_vecs = 0.;
    
    for (int i = 0; i<non_spurious.size(); i++) {
      eig_vals(i) = energies(non_spurious[i]);
      resnorms(i) = resid(non_spurious[i]);
      Submatrix submat = vecs(0,adc_mat.main_block_size(), non_spurious[i], 1);
      eig_vecs(0,adc_mat.main_block_size(),i,1).daxpy(1.,submat);
    }
    
    
  } else{ // Additional steps are needed for obtaining the full Lanczos vectors    
    
    // Remove the number of requested vectors that are outside the reasonable range
    vector<unsigned> out_list;
    set<unsigned int>::iterator it;
    for (it=diag_.vecs.begin(); it!=diag_.vecs.end(); it++) 
      if (*it >= non_spurious.size()) 
	out_list.push_back(*it);
    
    if (out_list.size()) {
      cout << " Warning: The following numbers of requested vectors are outside the range: ";
	for (int i = 0; i < out_list.size();i++) {
	  diag_.vecs.erase ( out_list[i] );
	  cout << out_list[i]+1 << ' ';
	}
      cout << endl;
    }
    
    // Leave only those non-spurious vectors that have been requested
    out_list.clear();
    for (int i = 0; i < non_spurious.size(); i++) {
      it = diag_.vecs.find(i);
      if (it != diag_.vecs.end()) {
	out_list.push_back(non_spurious[i]);
      }
    }
    non_spurious = out_list;

    // Get the full eigenvectors (do all needed to call BuildVectors)
    double** vec_list = new double*[non_spurious.size()];
    Blas_matrix** vec_addr = new Blas_matrix*[non_spurious.size()];
    for(int i = 0; i<non_spurious.size(); i++) {
      vec_list[i] = &vecs(0,non_spurious[i]);
      vec_addr[i] = new Blas_matrix(adc_mat.size(),1);
    }
    
    //Rebuild the initial guess before resetting
    for(int j = 0; j< m_nroots;j++) { 
      *guess[j] = 0.;
      (*guess[j])(j) = 1.;
    }

    projection.Reset(guess);
    adc_mat.reset();//reset the internal counting of the matrix-vector multiplies 
    projection.BuildVectors(vec_addr, vec_list, diag_.vecs.size(), lanc_space, m_nroots);

    
    
    // Get only the requested eigenvalues, residuals and eigenvectors
    eig_vals.allocate(non_spurious.size());
    resnorms.allocate(non_spurious.size());
    eig_vecs.allocate(adc_mat.size(), non_spurious.size());
    eig_vecs = 0.;
    for (int i = 0; i<non_spurious.size(); i++) {
      eig_vals(i) = energies(non_spurious[i]);
      resnorms(i) = resid(non_spurious[i]);
      for(int j = 0; j < eig_vecs.rows();j++)
	eig_vecs(j,i) = (*vec_addr[i])(j);
    }
    
      
    // Delete the auxiliary data
    delete [] vec_list;
    for(int i = 0; i<non_spurious.size(); i++) {
      delete vec_addr[i];
    }
    delete [] vec_addr;
    
  }
    
  projection.CleanUp();   
  

  time(&t2);
  cout << " Eigenvalues and eigenvectors fetched for ";  
  cout << " time: " << t2 - t1 << " s" <<endl;  flush(cout);
  
  
  for(int j = 0; j< m_nroots;j++) {
    delete guess[j]; 
  }
  delete [] guess;
  
}



// This method looks for non-spurious vectors. The vectors are obtained with subspace iterations
// and all spurious vectors have zero components in the mains space. (see Tarantelli's Chem. Phys. 329, 11)
void ADC_analyzer::find_spurious(int main_space, Blas_matrix& vecs, vector<unsigned>& nspur)
{
  static const double spur_thresh = 10e-10;
  
  // for each vector
  for(int i = 0; i < vecs.cols(); i++ ) {
    bool is_spurious = true;

    // for each component in the main space
    for (int j=0; j < main_space; j++) 
      if (abs(vecs(j,i)) > spur_thresh )
	is_spurious = false;

    if (!is_spurious) {nspur.push_back(i);}

  }
}
