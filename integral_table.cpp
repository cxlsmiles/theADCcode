#include "integral_table.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "integral_blocks.hpp"
#include <iostream>
#include <string>

using namespace std;

const unsigned Integral_table::pairs2table[][3] =
  {
    {OOOO, OVOO, VVOO},
    {OVOO, OVOV, VVOV},
    {VVOO, VVOV, VVVV}
  };


//The function vpqrs_ and and the global pointer (allocated by adc_selector)
// allow to the Fortran and C code to use the integral table.

extern Integral_table* integral_table;

extern "C" double vpqrs_(int* a, int* b, int* c, int* d)
{
  // we have different numbering of the orbitals and symmetries
  // Joerg's starts from 1
  return integral_table->integral(*a-1,*b-1,*c-1,*d-1);
}


  //Returns a pointer to the table that holds the integral for the requested orbital pairs
Matrix<double>* Integral_table::get_table(struct Pair_info& pair1, struct Pair_info& pair2) const 
{
  if (pair1.sym != pair2.sym) throw 0; // The requested integral is zero due to symmetry
  Matrix<double>* table = int_tables_[pairs2table[pair1.type][pair2.type]][pair2.sym];
  if(!table) throw 1;                  // Checks whether the table has been loaded
  if (pair1.type < pair2.type) swap_pos(pair1.pos, pair2.pos); // The "higher" type of pair corresponds to the row "index"
  return table;                                                // of the tables (see above). Here it's sufficient to swap the
}                                                              // the positions instead of the pairs themselves.



// Builds the integral tables using the SCF data provided by the SCF_data_reader.
// The user may specify which tables need be build; that information is stored in 'include':
// 00111111
//   ||||||
//   |||||OOOO to be loaded
//   ||||OVOO to be loaded
//   |||OVOV to be loaded
//   ||VVOO to be loaded
//   |VVOV to be loaded
//   VVVV to be loaded
Integral_table::Integral_table(SCF_data_reader& phis, char include) : scf_(phis)
{
  // Initialize the auxiliary data
  pair_groups_.push_back(vector<unsigned>(scf_.number_irreps(),0));
  pair_groups_.push_back(vector<unsigned>(scf_.number_irreps(),0));
  pair_groups_.push_back(vector<unsigned>(scf_.number_irreps(),0));

  lowest_vir_in_group_.assign(scf_.number_irreps(), 0);
  vir_group_sizes_.assign(scf_.number_irreps(), 0);
  for(unsigned i = scf_.number_occupied(); i < scf_.number_orbitals(); i++){
    if (!vir_group_sizes_[scf_.irrep(i)])
      lowest_vir_in_group_[scf_.irrep(i)] = i;
    vir_group_sizes_[scf_.irrep(i)]++;
  }

  //Initialize the tables
  create_lookup();          
  allocate_tables(include); 
  load_integrals();         

  // Allocate and initialize the pointers to the integral blocks
  rkij_blocks_ = new Triangular_container<Daxpy_argument>(scf_.number_occupied(), scf_.number_occupied());
  rsij_blocks_ = new Triangular_container<Daxpy_argument>(scf_.number_irreps(), scf_.number_occupied());
  risj_blocks_ = new Rectangular_container<Daxpy_argument>(scf_.number_irreps(), scf_.number_occupied());

  // The (rk|ij) and (rs|ij) blocks
  for(unsigned int i = 0; i < scf_.number_occupied(); i++)
    for(unsigned int j = 0; j <=i; j++) {
      
      for(unsigned int k = 0; k < scf_.number_occupied(); k++) 
	(*rkij_blocks_)(i,j,k) = load_rkij(i,j,k);      

      for(unsigned int sym = 0; sym < scf_.number_irreps(); sym++) 
	(*rsij_blocks_)(i,j,sym) = load_rsij(i,j,sym);
    }

  // The (ri|sj) blocks
  for(unsigned int i = 0; i < scf_.number_occupied(); i++)
    for(unsigned int j = 0; j < scf_.number_occupied(); j++) 
      for(unsigned int sym = 0; sym < scf_.number_irreps(); sym++) 
	(*risj_blocks_)(i,j,sym) = load_risj(i,j,sym);

  // Clear the auxiliary data
  vir_group_sizes_.clear();
  lowest_vir_in_group_.clear();
  pair_groups_.clear();

}

// Initializes the lookup table of orbital pairs
inline void Integral_table::create_lookup()
{
  
  lookup_ = new Triangular_matrix<Pair_info>(scf_.number_orbitals());
  // Initialize the information for all but the VV pairs
  for(unsigned int col_orb = 0; col_orb < scf_.number_orbitals(); col_orb++) 
    for(unsigned int row_orb = col_orb; row_orb < scf_.number_orbitals(); row_orb++) {
      
      unsigned pair_sym = scf_.irrep_product(scf_.irrep(col_orb), scf_.irrep(row_orb));
      
      // Find the type of the pair
      Pair_type pair_type;
      if (row_orb >= scf_.number_occupied()) {
	if (col_orb >= scf_.number_occupied()) {
	  pair_type = VV;
	  continue; // skip the VV pairs
	}else 
	  pair_type = OV;
	
      } else 
	pair_type = OO;
      
      (*lookup_)(row_orb, col_orb).type = pair_type;
      (*lookup_)(row_orb, col_orb).sym  = pair_sym;
      (*lookup_)(row_orb, col_orb).pos  = pair_groups_[pair_type][pair_sym]++;
    }
  // In order to ensure the contiguity of the (rs|ij) blocks, 
  // the VV pairs must be order by symmetry
  // For example, say, there are 2 occupied and 4 virtual orbitals:
  // SYM1: O1, V1, V3; SYM2: O2, V2, V4
  // If there were no symmetry separation of the tables the (rs|O1O2) and (rs|O1O1) block
  // would look like that:
  // SYM2 B(O1,O2):V1           V3          V2          V4      SYM1 B(O1,O1):    V1          V3          V2          V4
  //       V1      0            0      (V1V2|O1O2) (V1V4|O1O2)            V1 (V1V1|O1O1) (V1V3|O1O1)      0           0
  //       V3      0            0      (V3V2|O1O2) (V3V4|O1O2)            V3 (V3V1|O1O1) (V3V3|O1O1)      0           0
  //       V2 (V2V1|O1O2)  (V2V3|O1O2)      0           0                 V2      0           0      (V2V2|O1O1) (V2V4|O1O1)
  //       V4 (V4V1|O1O2)  (V4V3|O1O2)      0           0                 V4      0           0      (V2V4|O1O1) (V4V4|O1O1)
  // The "subblocks" are defined by the symmetry of the virtual orbital index.

  // Thus the goal is to achieve the following ordering of the VVOO symmetry tables:
  // SYM1             SYM2
  // OOVV: O1O1 O2O2  OOVV: O1O2 
  // V1V1   A    C    V1V2   E   
  // V1V3   A    C    V1V4   E   
  // V3V3   A    C    V2V3   E   
  // V2V2   B    D    V3V4   E
  // V2V4   B    D       
  // V4V4   B    D       

  // Note that when the symmetry of the column virtual orbital of B(O1,O2) changes from SYM1 -> SYM2
  // the block gets transposed
  //        (V2V1|O1O2)  (V2V3|O1O2)    ->    (V1V2|O1O2)  (V1V4|O1O2)  
  //        (V4V1|O1O2)  (V4V3|O1O2)          (V3V2|O1O2)  (V3V4|O1O2)
  // and both blocks are "contained" in the second table OOVV(SYM2)

  // Initialize the totally symmetric VV pairs
  for(int col_sym = 0; col_sym < scf_.number_irreps(); col_sym++)
    for(int row_sym = 0; row_sym < scf_.number_irreps(); row_sym++)
      
      for(int col_orb = 0; col_orb < scf_.number_orbitals(); col_orb++) {
	
	if (scf_.irrep(col_orb) != col_sym) continue;
	
	for(int row_orb = col_orb; row_orb < scf_.number_orbitals(); row_orb++) {
      
	  if (scf_.irrep(row_orb) != row_sym) continue;
	  
	  int pair_sym = scf_.irrep_product(scf_.irrep(col_orb), scf_.irrep(row_orb));
	  
	  if (pair_sym != 0)
	    continue; // Skip non-totally symmetric
	  
	  // Finds the type of the pair
	  Pair_type pair_type;
	  if (row_orb >= scf_.number_occupied()) {
	    if (col_orb >= scf_.number_occupied())
	      pair_type = VV;
	    else {
	      pair_type = OV;
	      continue; // Skip the OV pairs
	    }
	  } else {
	    pair_type = OO;
	    continue; // Skip the OO pairs
	  }
	  
  	  (*lookup_)(row_orb, col_orb).type = pair_type;
	  (*lookup_)(row_orb, col_orb).sym  = pair_sym;
	  (*lookup_)(row_orb, col_orb).pos  = pair_groups_[pair_type][pair_sym]++;
	}
      }
  
  // Initialize the rest of the VV pairs
  // Here the order of the symmetry loop is noteworthy:
  // The row virtual orbital index has lower symmetry than the column virtual orbital index.
  // That is important to know when "transposing" the memory block (see above)
  for(int row_sym = 0; row_sym < scf_.number_irreps(); row_sym++)
    for(int col_sym = row_sym+1; col_sym < scf_.number_irreps(); col_sym++)
      
      for(int col_orb = 0; col_orb < scf_.number_orbitals(); col_orb++) {
	
	if (scf_.irrep(col_orb) != col_sym) continue;
	
	for(int row_orb = 0; row_orb < scf_.number_orbitals(); row_orb++) {
	  
	  if (scf_.irrep(row_orb) != row_sym) continue;
	  
	  int pair_sym = scf_.irrep_product(scf_.irrep(col_orb), scf_.irrep(row_orb));
	  
	  if (pair_sym == 0)
	    continue; // Skip the totally symmetric
	  
	  Pair_type pair_type;
	  if (row_orb >= scf_.number_occupied()) {
	    if (col_orb >= scf_.number_occupied())
	      pair_type = VV;
	    else {
	      pair_type = OV;
	      continue; //Skip the OV pairs
	    }
	  } else {
	    pair_type = OO;
	    continue; // Skip the OO pairs
	  }
	  
	  (*lookup_)(row_orb, col_orb).type = pair_type;
	  (*lookup_)(row_orb, col_orb).sym = pair_sym;
	  (*lookup_)(row_orb, col_orb).pos = pair_groups_[pair_type][pair_sym]++;
	}
      }
}

// Computes the powers of 2
unsigned char Integral_table::power2(int power) const
{ 
  unsigned char res = 1; 
  for(int i = 1; i <= power; i++) res *= 2; 
  return res; 
}

// Allocates the tables requested by the user
// and initializes them to zero. 
inline void Integral_table::allocate_tables(char include)
{
  
  // Set all table pointers to zero.
  for (unsigned type = OOOO; type < TABLE_TYPES; type++) 
    for (unsigned sym = 0; sym < MAX_SYMS; sym++)
      int_tables_[type][sym] = 0;

  for (int type = OOOO; type < TABLE_TYPES; type++) {    

    if (power2(type) & include)  // if the user has specified the type
      
      for (unsigned sym = 0; sym < scf_.number_irreps(); sym++) 
	switch (type) {
	case OOOO: 
	  int_tables_[OOOO][sym]
	    = new Triangular_matrix<double>(pair_groups_[OO][sym]);
	  break;
	case OVOO:
	  int_tables_[OVOO][sym]
	    = new Rectangular_matrix<double>(pair_groups_[OV][sym],
					     pair_groups_[OO][sym]);
	  break;
	case OVOV:
	  int_tables_[OVOV][ sym]
	    = new Triangular_matrix<double>(pair_groups_[OV][sym]);
	  break;
	case VVOO:
	  int_tables_[VVOO][ sym]
	    = new Rectangular_matrix<double>(pair_groups_[VV][sym],
					     pair_groups_[OO][sym]);
	  break;
	case VVOV:
	  int_tables_[VVOV][ sym]
	    = new Rectangular_matrix<double>(pair_groups_[VV][sym],
					     pair_groups_[OV][sym]);
	  break;
	case VVVV:
	  int_tables_[VVVV][ sym]
	    = new Triangular_matrix<double>(pair_groups_[VV][sym]);
	}
  }
  
  // Initialize the allocated tables to zero.
  for(unsigned col = 0; col < MAX_SYMS; col++)
    for(unsigned row = 0; row < TABLE_TYPES; row++) 
      if (int_tables_[row][col]) 
	for(unsigned int ii = 0; ii < int_tables_[row][col]->rows(); ii++)
	  for(unsigned int jj = 0; jj < int_tables_[row][col]->cols(); jj++) 
	    (*int_tables_[row][ col])(ii,jj) = 0.;
}


// Distribute the integrals delivered by the SCF_data_reader to the proper tables
inline void Integral_table::load_integrals()
{
  // Allocate an array to store the number of integrals read
  Rectangular_matrix<unsigned int> integral_count(TABLE_TYPES, MAX_SYMS);
  for(unsigned int i = 0; i < integral_count.rows(); i++ )
    for(unsigned int j = 0; j < integral_count.cols(); j++ )
      integral_count(i,j) = 0;
  
  // Loop until the SCF reader sets l to -1, see phis.hpp
  Matrix<double>* table;
  double integral;
  int i,j,k,l;
  Pair_info first_pair, second_pair;
  do {

    integral = scf_.get_integral(i, j, k, l);
    if (l == -1) break;

    first_pair  = lookup_->operator()(i,j);
    second_pair = lookup_->operator()(k,l);
    
    try {
      
      table = get_table(first_pair, second_pair);
      
    } catch (int a) {
      
      // If the integral index is not totally symmetric 
      // (though that will mean that the SCF reader is not working properly)
      // or the table for that type has not been created
      // just continue to the next integral

      continue;
    }

    table->operator()(first_pair.pos, second_pair.pos) = integral;
    integral_count(pairs2table[first_pair.type][second_pair.type], second_pair.sym) ++; // Keep track of what is being read
    
  } while (1);
  
  // Print the gathered information
  cout << " Integral table loaded:\n";
  for(unsigned type = 0; type < TABLE_TYPES; type++) {
    unsigned int sum_loaded = 0, sum_all = 0;
    for(unsigned sym = 0; sym < MAX_SYMS; sym++)
      if (int_tables_[type][sym]) {
	sum_loaded += integral_count(type,sym);
	sum_all += int_tables_[type][sym]->size();
      }

    string int_type;
    switch (type) {
    case OOOO: int_type = "(OO|OO)"; break;
    case OVOO: int_type = "(OV|OO)"; break;
    case OVOV: int_type = "(OV|OV)"; break;
    case VVOO: int_type = "(VV|OO)"; break;
    case VVOV: int_type = "(VV|OV)"; break;
    case VVVV: int_type = "(VV|VV)"; break;
    }

    if (sum_all)
      cout << " type " << int_type << " integrals: " 
	   << sum_loaded << " of " << sum_all << " found.\n";
  }
  cout << endl;
}

 
Integral_table::~Integral_table()
{
  delete lookup_;

  for(int sym = 0; sym < MAX_SYMS; sym++)
    for(int type = 0; type < TABLE_TYPES; type++)
      delete int_tables_[type][sym];

  delete risj_blocks_; delete rkij_blocks_; delete rsij_blocks_;

}


// Loads a Vij,k block (a vector) of integrals
Daxpy_argument* Integral_table::load_rkij(unsigned i,unsigned j,unsigned k)
{
 
  unsigned sym = scf_.irrep_product(scf_.irrep(i),scf_.irrep_product(scf_.irrep(j),scf_.irrep(k)));
  
  unsigned first_virtual = lowest_vir_in_group_[sym];
  if (!first_virtual) return 0; // There is no virtual orbital for the given symmetry
  
  unsigned num_virtuals = vir_group_sizes_[sym];

  Pair_info first_pair = (*lookup_)(i,j),
    second_pair = (*lookup_)(k,first_virtual);

  Matrix<double>& mat = *int_tables_[OVOO][first_pair.sym];
  
  unsigned col_start = first_pair.pos;
  unsigned row_start = second_pair.pos;
  unsigned col_slice = 1;
  unsigned row_slice = num_virtuals;
  
  return new Rectangular_submatrix(mat, row_start, row_slice, col_start, col_slice);
}

// Loads a B(i,j) block of integrals for a specified symmetry of the virtual orbital
// column index (see above)
Daxpy_argument* Integral_table::load_rsij(unsigned i,unsigned j,unsigned sym)
{
  unsigned  first_virtual_col =  lowest_vir_in_group_[sym];
  if (!first_virtual_col) return 0;

  unsigned num_virtual_col =  vir_group_sizes_[sym];

  unsigned sym_row = scf_.irrep_product(scf_.irrep_product(scf_.irrep(i),scf_.irrep(j)),sym);
  unsigned int first_virtual_row = lowest_vir_in_group_[sym_row];
  if (!first_virtual_row) return 0;

  unsigned num_virtual_row = vir_group_sizes_[sym_row];

  Pair_info first_pair = (*lookup_)(i,j),
    second_pair = (*lookup_)(first_virtual_row,first_virtual_col);
  
  unsigned row_start, col_start, col_slice, row_slice;
  if (sym > sym_row) {// The choice here (and below) depends on the ordering of the VV pairs
    col_slice = num_virtual_col;
    row_slice = num_virtual_row;
  } else {
    col_slice = num_virtual_row;
    row_slice = num_virtual_col;
  }
  col_start = first_pair.pos;
  row_start = second_pair.pos;

  // The cast is OK, since the VVOO tables are always rectangular
  Rectangular_matrix<double>& mat = dynamic_cast<Rectangular_matrix<double>&>(*int_tables_[VVOO][first_pair.sym]);  
  if (first_pair.sym == 0) 
    return new Lower_triangular_argument(mat, row_start, col_start, row_slice);
  if (sym > sym_row)
    return new Rectangular_proxymatrix(mat,row_start,row_slice,col_start,col_slice);
  else 
    return new Transposed_proxymatrix(mat,row_start,row_slice,col_start,col_slice);
}


// Loads a A(i,j) block of integrals for a specified symmetry of the virtual orbital
// column index 
Daxpy_argument* Integral_table::load_risj(unsigned i,unsigned j,unsigned sym)
{
  
  unsigned  first_virtual_col =  lowest_vir_in_group_[sym];
  if (!first_virtual_col) return 0;
  unsigned num_virtual_col =  vir_group_sizes_[sym];

  int sym_row = scf_.irrep_product(scf_.irrep_product(scf_.irrep(i),scf_.irrep(j)),sym);
  unsigned int first_virtual_row = lowest_vir_in_group_[sym_row];
  if (!first_virtual_row) return 0;
  unsigned num_virtual_row = vir_group_sizes_[sym_row];

  Pair_info first_pair = (*lookup_)(i,first_virtual_row),
    second_pair = (*lookup_)(j,first_virtual_col);
 
  unsigned col_start,row_start,col_slice,row_slice;
  if (i > j) { // Here (and below) the choice depends on the ordering of the OO pairs
    col_start = first_pair.pos;
    row_start = second_pair.pos;
    col_slice = num_virtual_row;
    row_slice = num_virtual_col;
  } else {
    col_start = second_pair.pos;
    row_start = first_pair.pos;
    col_slice = num_virtual_col;
    row_slice = num_virtual_row;
  }

  // The OVOV tables are triangular
  Triangular_matrix<double>& mat = dynamic_cast<Triangular_matrix<double>&>(*int_tables_[OVOV][first_pair.sym]);
  if (i == j) 
    return  new Upper_triangular_argument(mat, row_start, col_start, col_slice);
  if (i < j) 
    return  new Rectangular_submatrix(mat,row_start,row_slice,col_start,col_slice);
  else 
    return  new Transposed_submatrix(mat,row_start,row_slice,col_start,col_slice);
}




