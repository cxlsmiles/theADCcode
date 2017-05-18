module adc2_pol_adapter

  use parameters
  use constants
 
  implicit none  



contains
  
  subroutine init_data(nb, no, orbs, syms, en, ni, sym, hcs, hc)
    integer, intent(in) :: nb, no, ni, sym, hcs;
    integer, dimension(nb), intent(in) :: orbs, syms;
    real*8,  dimension(nb), intent(in) :: en;
    integer, dimension(hcs), intent(in) :: hc;


  integer, dimension(2) :: shp

!!$ Setting up the multiplication table

  shp(:)=(/ 8,8 /)
  MT=reshape(mtrow, shp)


  nBas = nb;
  nIrr = ni;
  nirrep = sym;
  nOcc = no;
  nVirt = nBas - nOcc;
  
  allocate(roccnum(nBas), orbSym(nBas), e(nBas));

  roccnum(:) = orbs(:);
  e(:) = en(:);
  orbSym(:) = syms(:);
  hcentre(0) = hcs;
  hcentre(1:hcs)  = hc(:);

  hinit=hcentre(1)


  end subroutine init_data
  
  
  
  
end module adc2_pol_adapter
          
          
       
       
       
       
    
