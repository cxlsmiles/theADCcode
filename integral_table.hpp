#ifndef __INTEGRAL_TABLE_HPP__
#define __INTEGRAL_TABLE_HPP__


// This file contains the declaration of the Integral_table class.
// It uses the SCF_data_reader to fetch the two electron integrals from the files produced by the SCF calculation
// and sorts them out. It furthermore allows for access to blocks of integrals,
// see also integral_blocks.hpp

// All pairs of orbitals are mapped by a triangular matrix lookup_ which
// gives information about the symmetry of the pair, its type and the relevant row/column
// in the tables of two-electron integrals. This information is contained in the struct Pair_info.

// The types of orbital pairs are: OO - two occupied, OV - a virtual and an occupied, VV two virtuals

// The rows and columns of the two-electron tables correspond to pairs of orbitals.
// Thus, the two-electron tables are separated according to the type of pairs: OOOO, OVOO, OVOV, VVOO, VVOV, VVVV;
// The tables are further separated by symmetry.

// For example, for only one symmetry the tables have the following form (3 occupied and 2 virtual orbitals):

// OOOO:  O1,O1 O1,O2 O1,O3 O2,O2 O2,O3 O3,O3    OVOO:  O1,O1 O1,O2 O1,O3 O2,O2 O2,O3 O3,O3
// O1,O1    V                                    O1,V1    A     D     G     J     M     P
// O1,O2    V     V                              O1,V2    A     D     G     J     M     P
// O1,O3    V     V     V                        O2,V1    B     E     H     K     N     Q
// O2,O2    V     V     V     V                  O2,V2    B     E     H     K     N     Q
// O2,O3    V     V     V     V     V            O3,V1    C     F     I     L     O     R
// O3,O3    V     V     V     V     V     V      O3,V2    C     F     I     L     O     R

// The ordering of the pairs is essential when a block of integrals need be fetched.
// Note that integrals with the same index triplet of occupied orbitals (only the virtual orbital index varies),
// e.g. (O2V1|O1O1), (O2V2|O1O1), form contiguous chunks of memory that can be fed to a BLAS routine.

// VVOO:  O1,O1 O1,O2 O1,O3 O2,O2 O2,O3 O3,O3    OVOV:  O1,V1 O1,V2 O2,V1 O2,V2 O3,V1 O3,V2
// V1,V1    A     B     C     D     E     F      O1,V1    A         
// V1,V2    A     B     C     D     E     F      O1,V2    A     A     
// V2,V2    A     B     C     D     E     F      O2,V1    B     B     D 
//                                               O2,V2    B     B     D     D
//                                               O3,V1    C     C     E     E     F
//                                               O3,V2    C     C     E     E     F     F

// Integrals with the same index doublet of occupied orbitals (only the two virtual orbital indexes vary) 
// form matrices whose columns are contiguous chunks of memory, e.g.:
//   (V1V1|O1O2)                   and          (O2V1|O1V1)  (O2V1|O1V2)    
//   (V1V2|O1O2)  (V2V2|O1O2)                   (O2V2|O1V1)  (O2V2|O1V2)

// The last block has to be transposed when the index doublet is inverted (O2O1 -> O1O2), i.e. the block is then
// contiguous along its rows.

// In order to form similar contiguous blocks when several symmetries are present, for the VV pairs it is required
// to use somewhat more complex ordering scheme dependent on the symmetry (see the create_lookup method).

// The VVOV  and VVVV tables have not been tested for bugs!!!

// integral_blocks.hpp declares several classes that point to the areas defined above. For a definition
// of the integral blocks with the same index doublets and triplets see Chem. Phis. 329, p. 13.
// Pointers to these classes are initialized and stored in three members: rkij_blocks_, risj_blocks_, rsij_blocks_. When a
// particular block is requested, its pointer is returned.


#include "matrices.hpp"
#include "block_container.hpp"
#include <vector>


#define TABLE_TYPES 6
#define MAX_SYMS    8

class SCF_data_reader;
class Daxpy_argument;


class Integral_table {
  SCF_data_reader& scf_;                     // Provides the SCF data (including the two-electron integrals)

  enum Pair_type  {OO = 0, OV, VV};          // Types of orbital pairs
  enum Table_type {OOOO = 0, OVOO, OVOV, VVOO, VVOV, VVVV}; // Types of integral tables
  
  struct Pair_info {                         // Contains information for each pair of orbitals:
    unsigned short int pos;                  // row/column index in the integral tables;
    unsigned char type, sym;                 // type of the pair and its symmetry.
  };
  
  Triangular_matrix<Pair_info> *lookup_;                // Maps the orbitals to the pair info and the pairs to the tables
  static const unsigned pairs2table[3][3];              // Maps the orbital types to the table types
  Matrix<double>* int_tables_[TABLE_TYPES][MAX_SYMS];   // Contains all tables: maximum 6 types and 8 symmetries


                                                            // See Chem. Phis. 329, p.13:
  Triangular_container<Daxpy_argument>*  rkij_blocks_;       // Stores pointers to the (rk|ij) blocks
  Rectangular_container<Daxpy_argument>* risj_blocks_;       // Stores pointers to the (ri|sj) blocks
  Triangular_container<Daxpy_argument>*  rsij_blocks_;       // Stores pointers to the (rs|ij) blocks
  
  // Auxiliary variables used in the initialization of the tables
  std::vector<unsigned int> vir_group_sizes_;     // Stores the number of virtual orbitals in a given symmetry.
  std::vector<unsigned int> lowest_vir_in_group_; // Stores the number of the first virtual orbital in a symmetry group
  std::vector<std::vector<unsigned int> > pair_groups_; // Stores the ordering of the orbital pairs grouped by type and symmetry
  // Initializing methods used by the constructor
  void create_lookup();                                        // Initializes the orbital-to-pair mapping matrix
  void allocate_tables(char include);                          // Allocates the requested integral tables
  void load_integrals();                                       // Loads the integrals passed by the SCF_reader
  Daxpy_argument* load_rkij(unsigned i, unsigned j, unsigned k);// Loads the proxy classes of the Vij,k blocks
  Daxpy_argument* load_rsij(unsigned i,unsigned j,unsigned sym);// Loads the proxy classes of the Bij blocks
  Daxpy_argument* load_risj(unsigned i,unsigned j,unsigned sym);// Loads the proxy classes of the Aij blocks
  
  unsigned char power2(int a) const;              // An auxiliary function used to check what type of tables are requested

  inline void swap_pos(unsigned short& pos1, unsigned short& pos2) const   // Swaps two short ints
  {unsigned short temp = pos1; pos1 = pos2; pos2 = temp;}

  Matrix<double>* get_table(struct Pair_info& pair1, struct Pair_info& pair2) const;



public:
  
  Integral_table(SCF_data_reader& phis, char include);
  ~Integral_table();
  
//   inline double integral(unsigned i, unsigned j, unsigned k, unsigned l) const
//   { 
//     Pair_info first_pair = lookup_->operator()(i,j), 
//       second_pair = lookup_->operator()(k,l);
//     if (first_pair.sym != second_pair.sym) return 0.;
//     Matrix<double>* table = int_tables_[pairs2table[first_pair.type][second_pair.type]][second_pair.sym];
//     if(!table) return 0.;
//     if (first_pair.type < second_pair.type) swap_pos(first_pair.pos, second_pair.pos);

//     return table->operator()(first_pair.pos, second_pair.pos);
//   }


  inline double integral(unsigned i, unsigned j, unsigned k, unsigned l) const
  {
    Pair_info first_pair = lookup_->operator()(i,j), 
      second_pair = lookup_->operator()(k,l);
    if (first_pair.sym != second_pair.sym) return 0.; 
    Matrix<double>* table = int_tables_[pairs2table[first_pair.type][second_pair.type]][second_pair.sym];
    if(!table) return 0.;
    if (first_pair.type < second_pair.type)     return table->operator()(second_pair.pos,first_pair.pos);
    return table->operator()(first_pair.pos, second_pair.pos);
  }


  inline Daxpy_argument* get_rkij_block(unsigned int i, unsigned int j, unsigned int k) 
  { return rkij_blocks_->operator()(i,j,k);   } // Get a (rk|ij) block 
  inline Daxpy_argument* get_rsij_block(unsigned int i, unsigned int j, unsigned int sym) 
  { return rsij_blocks_->operator()(i,j,sym); } // Get a (rs|ij) block for a specified symmetry of the column virtual orbital 's'
  inline Daxpy_argument* get_risj_block(unsigned int i, unsigned int j, unsigned int sym) 
  { return risj_blocks_->operator()(i,j,sym); } // Get a (ri|sj) block for a specified symmetry of the column virtual orbital 's'

};


#endif //#ifndef __INTEGRAL_TABLE_HPP__
