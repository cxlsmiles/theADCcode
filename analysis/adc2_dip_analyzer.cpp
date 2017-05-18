#include "adc2_dip_analyzer.hpp"
#include "matrices.hpp"
#include "adc2_dip/adc2_matrix.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "blas_matrix.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

static const double SQRT_2 =  1.41421356237310;
static const double SQRT_1_2 = 0.707106781186547;
const static double au2eV = 27.211396;


ADC2_DIP_analyzer::ADC2_DIP_analyzer(SCF_data_reader& phis,  OrbSet groups, struct Diag_info& diag, struct Eigen_info& eig) 
  : _phis(phis), _groups(groups), ADC_analyzer(diag, eig)  
{


  ostringstream oss;
  set<unsigned>::iterator it;
  if (!groups.empty()) {
    
    pair<set<unsigned>::iterator,bool> inserted;
    set<unsigned> orbs;
    
    list<pair<string, set<unsigned> > >::iterator lt;
    // test if the define sets overlap
    for (lt=groups.begin(); lt!=groups.end(); lt++) {
      
      for (it=lt->second.begin(); it!=lt->second.end(); it++) {
	
	inserted = orbs.insert(*it);
	if (!inserted.second) {
	  
	  oss.str("");
	  oss << "Orbital " << *it+1 << " from set " << lt->first 
	      << " is already assinged to another  group" << endl;
	  
	  throw oss.str();
	}
	
      }
    }
    
    // call phis_get_scfvec to find the number of atomic basis functions
    // So use an inderect way to do that.
    int n = 0, len = 0;
    _phis.get_scfvec(0, &n, &len);
    len = abs(len);
 

    for (int i = 0; i < len; i++) {
      
      inserted = orbs.insert(i);
      
      if (inserted.second) {
	oss.str("");
	oss << "Orbital " << i+1 << " not assigned to any set" << endl;
	
	throw oss.str();
      }
      
      
    }
    
    if (orbs.size() > len) {
      oss.str("");
	oss << "More orbitals requested than the basis size " 
	    << len << "." << endl;
	
	throw oss.str();
    }
    
  }
  
}


int ADC2_DIP_analyzer::analyze( ADC_matrix& adc_mat)
{


  Adc2_matrix* is_adc2 = dynamic_cast<Adc2_matrix*>(&adc_mat);

  if (!is_adc2)
    throw string("ADC2_DIP_analyzer:Requested two-hole population analysis of a non-DIP propagator.");
  
  Blas_matrix eig_vals, eig_vecs, norms;

  cout << " Computing the spectrum for ";
  cout << "symmetry " << adc_mat.symmetry()+1 << ", spin " << adc_mat.spin()+1 << endl;

  diagonalize(adc_mat,eig_vecs, eig_vals, norms);
  print_eigenvals(adc_mat, eig_vecs, eig_vals, norms);


  if (!_groups.empty())
    twohole_mulliken(adc_mat, eig_vecs, eig_vals);


  return 0;
}

//See Tarantelli JCP, 94 (1991) 523
//eq. 3a - 8
void ADC2_DIP_analyzer::twohole_mulliken(ADC_matrix& adc2_mat, Blas_matrix& eig_vecs, Blas_matrix& eig_vals)
{
  
  cout << "\n ADC two-hole population analysis\n";

  int nSym, nBas = 0, nCent;
  
  _phis.get_info(&nSym, &nBas, &nCent);
  
  int n = nBas, len = 0;
  // Get the number of basis functions "len"
  _phis.get_scfvec(0, &n, &len);
  len = abs(len);
  n = abs(n);

  Blas_matrix* C = new Blas_matrix(len, n);
  _phis.get_scfvec(&(*C)(0,0), &n, &len);

  // This is an array of indexes that serves as a map
  // between the pair of indexes of atomic orbitals
  // and the physical number of the rows and columns
  Triangular_matrix<int> index(len);
  // indexes in the TRIPLET case for TESTING
  int count = 0;
  // The indexes obey the restrictions specified in eq. 3a-3d.
  // Those are different for triplet and singlet.
  for(int row = 0; row < index.rows(); row++) 
    for(int col = 0; col <= row; col++) {
      index(row,col) = count++;
      //cout << row << ',' << col << ':' << index(col,row) << endl;
    }
  

  // Do the same here but for the pairs of indexes of molecular orbitals
  Triangular_matrix<int> MO_index(_phis.number_occupied());
  count = 0;
  // Never forget symmetry
  
  
  for(int row = 0; row < _phis.number_occupied(); row++) {
    if (adc2_mat.symmetry() != _phis.irrep_product(_phis.irrep(row),_phis.irrep(row)) )
      continue;
    if (adc2_mat.spin() == 2)
      continue;
    
    MO_index(row,row) = count++;
    
    
    //cout << MO_index(row,row) <<':' << adc2_mat.get_conf(MO_index(row,row)) <<endl;
  }
  
  for(int row = 0; row < _phis.number_occupied(); row++) 
    for(int col = 0; col < row; col++) {
      if (adc2_mat.symmetry() != _phis.irrep_product(_phis.irrep(col),_phis.irrep(row)) )
	continue;
      if ((col == row) && (adc2_mat.spin() == 2))
	continue;
      MO_index(row,col) = count++;
//              cout << row << ',' << col << ':' << MO_index(col,row) << endl;
//              cout << "conf:" << adc2_mat.get_conf(count - 1) << endl;
      
//             cout << MO_index(row,col) <<':' << adc2_mat.get_conf(MO_index(row,col)) <<endl;
      
      
    }
  
  //cout << endl << "count:" << count << " blcksize:" <<  adc2_mat.main_block_size() << endl;
  
  int AO_indexes = index.size();
  //int MO_indexes = count;
  int MO_indexes = adc2_mat.main_block_size();
  
  
  Blas_matrix U(AO_indexes, MO_indexes);
  U = 0.;
  //eq. (3a)
  //loop over the cols of U in the triplet case

  int fact1, fact2;
  if (adc2_mat.spin() == 0) 
    fact1 = fact2 = 1;
  else {
    fact1 = -1; fact2 = 0;
  }

  cout << endl;

  //#pragma omp parallel for
  for(int i = 0; i < _phis.number_occupied(); i++)
    for(int j = 0; j <= i; j++) {
      if (adc2_mat.symmetry() != _phis.irrep_product(_phis.irrep(i),_phis.irrep(j)) )
	continue;
      if ((i == j) && (adc2_mat.spin() == 2))
	continue;
      
      // loop over the rows of U
      for(int p = 0; p < len; p++ )
	for(int q = 0; q <= p; q++) {
	  
	    if ((i != j) && (p != q)){
	      U(index(p,q),MO_index(i,j)) 
		= (*C)(p,i)*(*C)(q,j) + fact1 * (*C)(q,i)*(*C)(p,j);
	    } else if ((i != j) && (p == q)) {
  
	      U(index(p,p),MO_index(i,j))
		= (*C)(p,i)*(*C)(p,j) * fact2;
	    
	    } else if ((i == j) && (p != q)){

	      U(index(p,q),MO_index(i,i))
		= SQRT_2 * (*C)(p,i)*(*C)(q,i) * fact2;
	    } else {

	      U(index(p,p),MO_index(i,i))
		= SQRT_1_2 * (*C)(p,i)*(*C)(p,i) * fact2;
            }
	  }
      }
  
  //deallocate C
  delete C;
  
  Blas_matrix* S = new Blas_matrix(len, len);
  _phis.get_overlap(&(*S)(0,0),&len);

   ////Test eq. C^T . S . C = 1, don't deallocate C above
   //Blas_matrix T((*C).cols(),(*S).cols());
   //T = 0.;
   //T.dgemm('T','N',(*C),(*S));
   //cout << "n:" << (*C).cols(); 
   //Blas_matrix N((*C).cols(),(*C).cols());
   //N = 0.;
   //N.dgemm('N','N',T,(*C));
   //N.print();



  //eq. (4)
  Blas_matrix O(AO_indexes, AO_indexes);

#pragma omp parallel for
  for(int p = 0; p < len; p++)
    for(int q = 0; q <= p; q++)
      // loop over the cols of O
      for(int r = 0; r < len; r++ )
	for(int s = 0; s <= r; s++) 
	  O(index(p,q),index(r,s)) 
	    = (*S)(p,r)*(*S)(q,s) + fact1 * (*S)(p,s)*(*S)(q,r);
  delete S;
  

   //Test eq. 5
   //Blas_matrix T(U.cols(),O.cols());
   //cout << "Ucols:" << U.cols() << " O.cols:" << O.cols();
   //T = 0.;
   //T.dgemm('T','N',U,O);
   //Blas_matrix N(U.cols(),U.cols());
   //N = 0.;
   //N.dgemm('N','N',T,U);
   //N.print();

  cout << "    Eigenvalue                           Orbital contributions\n"
       << "    -------------------------------------------------------------------------------\n";




  OrbSet::iterator gr, gr1;  
  int pop_columns = _groups.size()*(_groups.size()+1)/2;
  int printed = 0;
  const int cols_per_row = 8; // prints up to 8 columns per row
  do { 

    int limit = (pop_columns - printed < cols_per_row) ? pop_columns - printed : cols_per_row;
    int count = 0;
    bool stop_printing = false;

    cout << setw(11) << ' ';
    for(gr = _groups.begin(); gr != _groups.end(); gr++) 
      for(gr1 = gr; gr1 != _groups.end(); gr1++) {
	
	if (count < printed) {count++; continue;}
	if (stop_printing) break;
	
	string label(gr->first);
	if (gr!=gr1)
	  label.append( "/" + gr1->first);
	cout << setw(9)<< right <<label;  

	count++;
	if (!((count-printed)%limit)) stop_printing = true;

      }
    cout << endl;
    

    // loop over all eigenvectors
    for (int i = 0; i < eig_vecs.cols(); i++) {
      
      if (!eig_vals(i)) continue;
      
      
      // use submatrix because this interface is not defined for blas_matrix
      //eq.6
      Blas_matrix Y(AO_indexes, 1);
      Y = 0;
      
      Y(0, AO_indexes, 0 , 1)
	.dgemv('N', U, eig_vecs(0, adc2_mat.main_block_size(), i, 1));
      
      cout << setw(11) << setiosflags(ios::fixed) << setprecision(4) << right  << eig_vals(i) * au2eV;
      
      double Qpqn;
      set<unsigned>::iterator p, q;
      
      stop_printing = false;
      count = 0;
      for(gr = _groups.begin(); gr != _groups.end(); gr++) {
	for(gr1 = gr; gr1 != _groups.end(); gr1++) {
	  
	  if (count < printed) {count++; continue;}
	  if (stop_printing) break;
	  
	  
	  
	  Qpqn = 0.;
	 	  
#pragma omp parallel for reduction(+:Qpqn)
	  //	  for(p = gr->second.begin(); p != gr->second.end(); p++) { //can't be split by omp
	  for(unsigned  elem = 0; elem < gr->second.size(); elem++) {
            set<unsigned>::iterator pp = gr->second.begin();
	    for(int i = 0; i < elem; i++) pp++;
	    set<unsigned>::iterator end;
	    if (gr == gr1) 
	      {end = pp; end++;}
	    else 
	      end = gr1->second.end();
	    
	    for(set<unsigned>::iterator qq = gr1->second.begin(); qq != end ; qq++) {
	      double temp = Y * O( 0, O.cols(),index(*pp,*qq), 1);	
	      temp *= Y(index(*pp,*qq));
	      Qpqn += temp;
	      
	    }
	    

	    
	  }

	  
	  cout << setw(9)<< right << Qpqn;

	  count++;
	  if (! ((count-printed) % limit)) stop_printing = true;	  
	}
	


      }
      
      cout << endl;
      
    }
    cout << endl;
    printed += limit;
  } while (printed < pop_columns);
}
