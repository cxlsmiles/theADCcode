#include "singlet.hpp"
#include "config.hpp"

static const double SQRT_2   = 1.41421356237310;
static const double SQRT_3   = 1.73205080756888;
static const double SQRT_1_2 = 0.707106781186548;
static const double SQRT_3_2 = 1.22474487139159;
static const double SQRT_3_4 = 0.866025403784439;
static const double HALF     = 0.5;
static const double ONE      = 1.;
static const double TWO      = 2.;
static const double THREEHALVES = 1.5;

// Chem. Phis. 329, p. 19, eq. A.3
double Singlet::w_term(unsigned i, unsigned k, unsigned s, unsigned r) 
{ 
  
  double result = 0.;
  double term;
 
  for(unsigned int m = 0; m < number_occupied(); m++) {
    
    if (sym(m) != sym_product(sym(i), sym(s), sym(r)))
      continue;
    
    term = (energy(r) + energy(s)) 
      - energy(m) 
      - HALF * (energy(i) + energy(k));
    
    term /= ((energy(r) + energy(s))
	     - (energy(i) + energy(m)))
      * ((energy(r) + energy(s))
	 - (energy(k) + energy(m)));

    term *= V1122(r, i, s, m) * V1122(r, k, s, m)
      + V1122(r, m, s, i) * V1122(r, m, s, k)
      + V1122_MINUS(r, i, s, m) * V1122_MINUS(r, k, s, m);

    result += (r == s) ? (HALF * term) : term; 
  } 
  
  return result;
}

// Chem. Phis. 329, p. 19, eq. A.2
double Singlet::u_term(unsigned i, unsigned j, unsigned k, unsigned l, unsigned s, unsigned r)
{
  double term;
  
  term = (energy(s) + energy(r))
    - HALF * (energy(i) + energy(j) 
	    + energy(k) + energy(l));

  term /= ((energy(r) + energy(s))
	   - (energy(i) + energy(j)))
    * ((energy(r) + energy(s))
       - (energy(k) + energy(l)));
  //Differes from triplet
  term *= V1122_PLUS(r, i, s, j)
    * V1122_PLUS(r, k, s, l);
  
  return (r == s) ? (HALF * term) : term;
}


//2h closed shell, 2h closed shell block
// Chem. Phis. 329, p. 19, eq. A.5
bool Singlet::block_ii_jj(const Config &row, const Config &col, double &element)
{

  double W = 0., U = 0.;
  
  unsigned i = row.occ[0];
  unsigned j = col.occ[0];

  bool delta_ij = i == j;
  

  for (unsigned int r = number_occupied(); r < number_orbitals(); r++)
    for (unsigned int s = number_occupied(); s <= r; s++) {
      
      if (delta_ij) 
	W += w_term(i, i, s, r);
      
      if (sym(r) == sym(s)) {
	U += u_term(i, i, j, j, s, r);
      }
    }


  element = V1122(i, j, i, j) - HALF * U;

  
  if(delta_ij)
    element += TWO * (W - energy(i));

  return true;
  
}

//2h open shell, 2h closed shell block
// Chem. Phis. 329, p. 19, eq. A.4
bool Singlet::block_ij_kk(const Config &row, const Config &col, double &element)
{

  double W = 0., U = 0.;
  
  unsigned i = row.occ[0];
  unsigned j = row.occ[1];  
  unsigned k = col.occ[0];

  bool delta_ik = i == k;
  bool delta_jk = j == k;
  

  for (unsigned int r = number_occupied(); r < number_orbitals(); r++)
    for (unsigned int s = number_occupied(); s <= r; s++) {
      
      if( delta_ik || delta_jk ) 
	W += w_term(i, j, s, r);
      
      if (sym(r) == sym(s))
	U += u_term(i, j, k, k, s, r);
      
    }

  element = V1122(i, k, j, k) - HALF * U + W;
  element *= SQRT_2;

  return true;
  
}

//2h-2h block
// Chem. Phis. 329, p. 19, eq. A.1
bool Singlet::block_ij_kl(const Config &row, const Config &col, double &element)
{

  double W = 0., U = 0.;
  
  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  
  unsigned k = col.occ[0];
  unsigned l = col.occ[1];
  
  bool delta_ik = i == k;
  bool delta_jk = j == k;
  bool delta_jl = j == l;
 

  for (unsigned int r = number_occupied(); r < number_orbitals(); r++)
    for (unsigned int s = number_occupied(); s <= r; s++) {
	  
      if( delta_ik && (sym(j) == sym(l)) )
	W += w_term(j, l, s, r);
      
      //Differs from triplet
      if( delta_jk && (sym(i) == sym(l)) )
	W += w_term(i, l, s, r);
      
      if( delta_jl && (sym(i) == sym(k)) )
	W += w_term(i, k, s, r);
      
      if (sym_product(sym(i), sym(j))
	  == sym_product(sym(r), sym(s)))
	U += u_term(i, j, k, l, s, r);
      
    }
  
  //Differs from triplet
  element = V1122_PLUS(i, k, j, l);
  element += W - U;
  
  if(delta_ik && delta_jl)
    element += -(energy(i) + energy(j));
  
  return true;
  
}

// Chem. Phis. 329, p. 17, Table A.2 
bool Singlet::block_lkk_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned int l = row.occ[0];
  unsigned int k = row.occ[1];
  unsigned int i = col.occ[0];
  
  bool delta_ik = i == k;
  bool delta_il = i == l;
  
  if (!delta_ik && !delta_il )
    return false;
  

  unsigned int block_size = size_vir_group(row.vir);
  block.allocate(block_size);

  block = 0.;
  
  if (delta_ik){
    block.daxpy( TWO * SQRT_2, V(i,i,l));
    block.daxpy(-SQRT_2, V(i,l,i));
  }

  if (delta_il) {
    block.daxpy(-SQRT_2, V(i,k,k));
  }

  return true;
}
// Chem. Phis. 329, p. 17, Table A.2 
bool Singlet::block_lkk_ij(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned l = row.occ[0];
  unsigned k = row.occ[1];
  
  unsigned i = col.occ[0];
  unsigned j = col.occ[1];
  
  bool delta_ik = i == k;
  bool delta_il = i == l;
  bool delta_jk = j == k;
  bool delta_jl = j == l;
  
  if (!delta_ik && !delta_il && !delta_jk && !delta_jl)
    return false;

  unsigned int block_size = size_vir_group(row.vir);
  block.allocate(block_size);

  block = 0.;
  
  if (delta_ik) {
    block.daxpy( TWO, V(i,j,l));
    block.daxpy(-ONE, V(j,l,i));
  }
  if (delta_il) {
    block.daxpy(-ONE, V(j,k,k));
  }
  if (delta_jk) {
    block.daxpy( TWO, V(i,j,l));
    block.daxpy(-ONE, V(i,l,j));
  }
  if (delta_jl) {
    block.daxpy(-ONE, V(i,k,k));
  }
  return true;
}

// Chem. Phis. 329, p. 17, Table A.2 
bool Singlet::block_klm_ii(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned k = row.occ[0];
  unsigned l = row.occ[1];
  unsigned m = row.occ[2];
  
  unsigned i = col.occ[0];
  
  bool delta_ik = i == k;
  bool delta_il = i == l;
  bool delta_im = i == m;
  
  if (!delta_ik && !delta_il && !delta_im )
    return false;
  
  unsigned int block_size = size_vir_group(row.vir);
  block.allocate(2 * block_size);
  block = 0.;
  
  Submatrix part[] = {block(0, block_size, 0, 1),
		      block(block_size, block_size, 0, 1)};
  
  if (delta_ik) {
    part[0].daxpy( TWO, V(k,m,l));
    part[0].daxpy(-ONE, V(k,l,m));
    part[1].daxpy( SQRT_3, V(k,l,m));
  }

  if (delta_il) {
    part[0].daxpy(-ONE, V(k,l,m));
    part[0].daxpy(-ONE, V(l,m,k));
    part[1].daxpy( SQRT_3, V(k,l,m));
    part[1].daxpy(-SQRT_3, V(l,m,k));
  }
  
  if (delta_im) {
    part[0].daxpy( TWO, V(k,m,l));
    part[0].daxpy(-ONE, V(l,m,k));
    part[1].daxpy(-SQRT_3, V(l,m,k));
  }

  return true;
}

// Chem. Phis. 329, p. 17, Table A.2 
bool Singlet::block_klm_ij(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned k = row.occ[0];
  unsigned l = row.occ[1];
  unsigned m = row.occ[2];
  
  unsigned i = col.occ[0];
  unsigned j = col.occ[1];
  
  bool delta_ik = i == k;
  bool delta_il = i == l;
  bool delta_im = i == m;
  bool delta_jk = j == k;
  bool delta_jl = j == l;
  bool delta_jm = j == m;
  
  if (!delta_ik && !delta_il && !delta_im
      && !delta_jk && !delta_jl && !delta_jm)
    return false;
  
  unsigned int block_size = size_vir_group(row.vir);
  block.allocate(2 * block_size);
  block = 0.;

  Submatrix part[] = {block(0, block_size, 0, 1),
		      block(block_size, block_size, 0, 1)};
  if (delta_ik) {
    part[0].daxpy( SQRT_2, V(j,m,l));
    part[0].daxpy(-SQRT_1_2, V(j,l,m));
    part[1].daxpy( SQRT_3_2, V(j,l,m));
  }

  if (delta_il) {
    part[0].daxpy(-SQRT_1_2, V(j,k,m));
    part[0].daxpy(-SQRT_1_2, V(j,m,k));
    part[1].daxpy( SQRT_3_2, V(j,k,m));
    part[1].daxpy(-SQRT_3_2, V(j,m,k));
  }
  
  if (delta_im) {
    part[0].daxpy( SQRT_2, V(j,k,l));
    part[0].daxpy(-SQRT_1_2, V(j,l,k));
    part[1].daxpy(-SQRT_3_2, V(j,l,k));
  }

  if (delta_jk) {
    part[0].daxpy( SQRT_2, V(i,m,l));
    part[0].daxpy(-SQRT_1_2, V(i,l,m));
    part[1].daxpy( SQRT_3_2, V(i,l,m));

  }

  if (delta_jl) {
    part[0].daxpy(-SQRT_1_2, V(i,k,m));
    part[0].daxpy(-SQRT_1_2, V(i,m,k));
    part[1].daxpy( SQRT_3_2, V(i,k,m));
    part[1].daxpy(-SQRT_3_2, V(i,m,k));
  }

  if (delta_jm) {
    part[0].daxpy( SQRT_2, V(i,k,l));
    part[0].daxpy(-SQRT_1_2, V(i,l,k));
    part[1].daxpy(-SQRT_3_2, V(i,l,k));
  }

  return true;
}

// Chem. Phis. 329, p. 18, Table A.3 
bool Singlet::block_jii_lkk(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned j = row.occ[0];
  unsigned i = row.occ[1];
  
  unsigned l = col.occ[0];
  unsigned k = col.occ[1];
  
  bool delta_il = i == l;
  bool delta_jl = j == l;
  bool delta_ik = i == k;
  bool delta_jk = j == k;
  bool delta_sym = sym(row.vir) == sym(col.vir);
  
  if (
      !( 
	delta_ik 
	|| (delta_ik && delta_jl) 
	|| (delta_il && delta_jk) ||
	( delta_sym
	  && (delta_ik || delta_il || delta_jk || delta_jl)	  
	  )
	)	 
      )
    return false;

  
  unsigned int row_block_size = size_vir_group(row.vir);
  unsigned int col_block_size = size_vir_group(col.vir);
  block.allocate(row_block_size, col_block_size);
  block = 0.;
 
  unsigned int col_sym = sym(col.vir);

  if (delta_ik) {
    block.daxpy( TWO, A(j,l,col_sym));
    block.daxpy(-ONE, B(j,l,col_sym));
  }
  
  if (delta_ik && delta_jl) {
    block.daxpy( ONE, A(i,i,col_sym));
    block.daxpy(-TWO, B(i,i,col_sym));
  }
  
  if (delta_il && delta_jk) {
    block.daxpy( ONE, A(l,j,col_sym));
    block.daxpy( ONE, B(j,l,col_sym));    

  }
  
  
  if (delta_sym) {
    double diag_term = 0.;
    
    if (delta_ik)
      diag_term += TWO * V1122(j,l,i,i) - V1122(i,j,i,l);
    if (delta_il)
      diag_term -= V1122(j,k,l,k);
    if (delta_jk)
      diag_term -= V1122(i,j,i,l);
    if (delta_jl)
      diag_term += V1122(i,k,i,k);
    
    if (delta_jl && delta_ik) {
      diag_term -= energy(j) + energy(i) + energy(i);
    }
    
    block.add_diag(diag_term);
    
    if (delta_jl && delta_ik) {
      block.add_diag(diag_energies(col_sym));
    }
  }
  
  return true;
}

// Chem. Phis. 329, p. 18, Table A.4
bool Singlet::block_ijk_mll(const Config &row, const Config &col, Blas_matrix& block)
{
  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  unsigned k = row.occ[2];
  
  unsigned m = col.occ[0];
  unsigned l = col.occ[1];
  
  bool delta_im = i == m;
  bool delta_jm = j == m;
  bool delta_km = k == m;
  bool delta_il = i == l;
  bool delta_jl = j == l;
  bool delta_kl = k == l;
  bool delta_sym = sym(row.vir) == sym(col.vir);
  
  if (
      !( 
	(delta_im && delta_jl) ||
	(delta_im && delta_kl) ||
	(delta_jm && delta_il) ||
	(delta_jm && delta_kl) ||
	(delta_km && delta_il) ||
	(delta_km && delta_jl) ||
	/* Check if any of the diagonal terms exists.*/
	(
	 delta_sym 
	 && (delta_im || delta_il || delta_jm || 
	     delta_jl || delta_km || delta_kl)
	 )
	)
      )
    return false;
  

  unsigned int row_block_size = size_vir_group(row.vir);
  unsigned int col_block_size = size_vir_group(col.vir);
  block.allocate(2 * row_block_size, col_block_size);
  block = 0.;
  
  Submatrix part[] = 
    {block(0, row_block_size, 0, col_block_size),
     block(row_block_size, row_block_size, 0, col_block_size)};

  unsigned int col_sym = sym(col.vir);
  
  if (delta_im && delta_jl) {
    part[0].daxpy( SQRT_1_2, A(k,j,col_sym));
    part[0].daxpy( SQRT_1_2, B(j,k,col_sym));
    part[1].daxpy(-SQRT_3_2, A(k,j,col_sym));
    part[1].daxpy( SQRT_3_2, B(j,k,col_sym));
  }
  
  if (delta_im && delta_kl) {
    part[0].daxpy(-SQRT_2, A(j,k,col_sym));
    part[0].daxpy( SQRT_1_2, B(j,k,col_sym));
    part[1].daxpy( SQRT_3_2, B(j,k,col_sym));
  }
  
  if (delta_jm && delta_il) {
    part[0].daxpy( SQRT_1_2, A(k,i,col_sym));
    part[0].daxpy(-SQRT_2, B(i,k,col_sym));
    part[1].daxpy(-SQRT_3_2, A(k,i,col_sym));
  }
  
  if (delta_jm && delta_kl) {
    part[0].daxpy( SQRT_1_2, A(i,k,col_sym));
    part[0].daxpy(-SQRT_2, B(i,k,col_sym));
    part[1].daxpy( SQRT_3_2, A(i,k,col_sym));
  }
  
  if (delta_km && delta_il) {
    part[0].daxpy(-SQRT_2, A(j,i,col_sym));
    part[0].daxpy( SQRT_1_2, B(i,j,col_sym));
    part[1].daxpy(-SQRT_3_2, B(i,j,col_sym));
  }
  
  if (delta_km && delta_jl) {
    part[0].daxpy( SQRT_1_2, A(i,j,col_sym));
    part[0].daxpy( SQRT_1_2, B(i,j,col_sym));
    part[1].daxpy( SQRT_3_2, A(i,j,col_sym));
    part[1].daxpy(-SQRT_3_2, B(i,j,col_sym));
  }
  
  /* The diagonal terms, Table A.4. */
  if (delta_sym) {
    double diag_term[] = {0., 0.};
    if (delta_im) {
      diag_term[0] -= V1122(j,l,k,l);
      diag_term[1] -= V1122(j,l,k,l);
    }
    if (delta_il) {
      diag_term[0] += TWO * V1122(i,k,j,m) - V1122(i,j,k,m);
      diag_term[1] += V1122(i,j,k,m);
    }
    if (delta_jm) {
      diag_term[0] +=  TWO * V1122(i,l,k,l);
      //noop diag1
    }
    if (delta_jl) {
      diag_term[0] -= V1122_PLUS(i,j,k,m);
      diag_term[1] += V1122_MINUS(i,j,k,m);
    }
    if (delta_km) {
      diag_term[0] -= V1122(i,l,j,l);
      diag_term[1] += V1122(i,l,j,l);
    }
    if (delta_kl) {
      diag_term[0] += TWO * V1122(i,k,j,m) - V1122(i,m,j,k);
      diag_term[1] -= V1122(i,m,j,k);
    }
    
    part[0].add_diag(diag_term[0] * SQRT_1_2);
    part[1].add_diag(diag_term[1] * SQRT_3_2);
  }
  
  return true;
}

// Chem. Phis. 329, p. 19, Table A.6
bool Singlet::block_ijk_lmn(const Config &row, const Config &col, Blas_matrix &block)
{
  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  unsigned k = row.occ[2];
  
  unsigned l = col.occ[0];
  unsigned m = col.occ[1];
  unsigned n = col.occ[2];
  
  bool delta_il = i == l;
  bool delta_jl = j == l;
  bool delta_kl = k == l;
  bool delta_jm = j == m;
  bool delta_km = k == m;
  bool delta_jn = j == n;
  bool delta_kn = k == n;
  bool delta_sym = sym(row.vir) == sym(col.vir);
  
  if (
      !(
	(delta_il && delta_jm) ||
	(delta_il && delta_km) ||
	(delta_il && delta_kn) ||
	(delta_jl && delta_km) ||
	(delta_jl && delta_kn) ||
	(delta_jm && delta_kn) ||
	/* Check if any of the diagonal terms exists.*/
	(
	 delta_sym
	 && (delta_il || delta_jl || delta_jm ||
	     delta_jn || delta_kl || delta_km || delta_kn)
	 )
	)
      )
    return false;
  
  unsigned int row_block_size = size_vir_group(row.vir);
  unsigned int col_block_size = size_vir_group(col.vir);
  block.allocate(2 * row_block_size, 2 * col_block_size);
  block = 0.;
  
  Submatrix part[][2] = {
    {block(0, row_block_size, 0, col_block_size),
     block(0, row_block_size, col_block_size, col_block_size)},
    
    {block(row_block_size, row_block_size, 0, col_block_size),
     block(row_block_size, row_block_size, col_block_size, col_block_size)}};

  unsigned int col_sym = sym(col.vir);
    
  if (delta_il && delta_jm) {
    part[0][0].daxpy( HALF, A(k,n,col_sym));
    part[0][0].daxpy(-ONE, B(k,n,col_sym));
    

    part[0][1].daxpy(-SQRT_3_4, A(k,n,col_sym));

    part[1][0].daxpy(-SQRT_3_4, A(k,n,col_sym));

    part[1][1].daxpy(THREEHALVES, A(k,n,col_sym));
    part[1][1].daxpy(-ONE, B(k,n,col_sym));
  }
  
  if (delta_il && delta_km) {
    part[0][0].daxpy(-ONE, A(j,n,col_sym));
    part[0][0].daxpy( HALF, B(j,n,col_sym));
    
    part[0][1].daxpy( SQRT_3, A(j,n,col_sym));
    part[0][1].daxpy(-SQRT_3_4, B(j,n,col_sym));
	
    part[1][0].daxpy(-SQRT_3_4, B(j,n,col_sym));

    part[1][1].daxpy(-HALF, B(j,n,col_sym));
  }
			  
  if (delta_il && delta_kn) {
    part[0][0].daxpy( TWO, A(j,m,col_sym));
    part[0][0].daxpy(-ONE, B(j,m,col_sym));

    // noop part[0][1] 
    // noop part[1][0] 

    part[1][1].daxpy(-ONE, B(j,m,col_sym));
  }
  
  if (delta_jl && delta_km) {
    part[0][0].daxpy( HALF, A(i,n,col_sym));
    part[0][0].daxpy( HALF, B(i,n,col_sym));
    
    part[0][1].daxpy(-SQRT_3_4, A(i,n,col_sym));
    part[0][1].daxpy(SQRT_3_4, B(i,n,col_sym));
    
    part[1][0].daxpy(SQRT_3_4, A(i,n,col_sym));
    part[1][0].daxpy(-SQRT_3_4, B(i,n,col_sym));
    
    part[1][1].daxpy(-THREEHALVES, A(i,n,col_sym));
    part[1][1].daxpy( HALF, B(i,n,col_sym));
  }
  
  if (delta_jl && delta_kn) {
    part[0][0].daxpy(-ONE, A(i,m,col_sym));
    part[0][0].daxpy( HALF, B(i,m,col_sym));

    part[0][1].daxpy(SQRT_3_4, B(i,m,col_sym));
    
    part[1][0].daxpy(-SQRT_3, A(i,m,col_sym));
    part[1][0].daxpy(SQRT_3_4, B(i,m,col_sym));

    part[1][1].daxpy(-HALF, B(i,m,col_sym));
  }
  
  if (delta_jm && delta_kn) {
    part[0][0].daxpy( HALF, A(i,l,col_sym));
    part[0][0].daxpy(-ONE, B(i,l,col_sym));
    
    part[0][1].daxpy(SQRT_3_4, A(i,l,col_sym));

    part[1][0].daxpy(SQRT_3_4, A(i,l,col_sym));
    
    part[1][1].daxpy(THREEHALVES, A(i,l,col_sym));
    part[1][1].daxpy(-ONE, B(i,l,col_sym));
  }
  
  if (delta_sym) {
    
    double diag_term[][2] = {{0., 0.},
			     {0., 0.}};
    
    if (delta_il) {
      diag_term[0][0] += V1122(j,m,k,n) - HALF * V1122(j,n,k,m);
      diag_term[0][1] += V1122(j,n,k,m);
      diag_term[1][0] += V1122(j,n,k,m);
      diag_term[1][1] += V1122(j,m,k,n) + HALF * V1122(j,n,k,m);
    }
    
    if (delta_jl) {
      diag_term[0][0] -= HALF * V1122_PLUS(i,m,k,n);
      diag_term[0][1] -= V1122_PLUS(i,m,k,n);
      diag_term[1][0] += V1122_MINUS(i,n,k,m); 
      diag_term[1][1] += HALF * V1122_MINUS(i,m,k,n);
    }
    
    if (delta_jm) {
      diag_term[0][0] += V1122_PLUS(i,l,k,n);
      //noop diag_term[0][1]
      //noop diag_term[1][0] 
      diag_term[1][1] += V1122_MINUS(i,l,k,n);
    }
    
    if (delta_jn) {
      diag_term[0][0] -= HALF * V1122_PLUS(i,l,k,m);
      diag_term[0][1] += V1122_PLUS(i,l,k,m);
      diag_term[1][0] += V1122_MINUS(i,l,k,m);
      diag_term[1][1] += HALF * V1122_MINUS(i,l,k,m);
    }
    
    if (delta_kl) {
      diag_term[0][0] += V1122(i,n,j,m) - HALF * V1122(i,m,j,n);
      diag_term[0][1] += V1122(i,m,j,n);
      diag_term[1][0] -= V1122(i,m,j,n);
      diag_term[1][1] -= V1122(i,n,j,m) + HALF * V1122(i,m,j,n);
    }
    
    if (delta_km) {
      diag_term[0][0] -= HALF * V1122_PLUS(i,l,j,n);
      diag_term[0][1] += V1122_MINUS(i,l,j,n);
      diag_term[1][0] += V1122_PLUS(i,l,j,n);
      diag_term[1][1] += HALF * V1122_MINUS(i,l,j,n);
    }
    
    if (delta_kn) {
      diag_term[0][0] += V1122(i,l,j,m) - HALF * V1122(i,m,j,l);
      diag_term[0][1] -= V1122(i,m,j,l);
      diag_term[1][0] -= V1122(i,m,j,l);
      diag_term[1][1] += V1122(i,l,j,m) + HALF * V1122(i,m,j,l);
   }

    part[0][0].add_diag(diag_term[0][0]);
    part[0][1].add_diag(diag_term[0][1] * SQRT_3_4);
    part[1][0].add_diag(diag_term[1][0] * SQRT_3_4);
    part[1][1].add_diag(diag_term[1][1]);

    if (delta_il && delta_jm && delta_kn){
      
      double epsi_ijk = -(energy(i) + energy(j) 
			  + energy(k));
      
      part[0][0].add_diag(diag_energies(col_sym));    
      part[0][0].add_diag(epsi_ijk);
      
      
      part[1][1].add_diag(diag_energies(col_sym));    
      part[1][1].add_diag(epsi_ijk);
      
      
    }
  }
  
  return true;  
}
