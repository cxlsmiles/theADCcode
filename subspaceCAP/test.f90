module test
implicit none

contains

  subroutine test_me(nmo,dens)
    
    integer, intent(in) :: nmo
    double precision, intent(in):: dens(nmo,nmo)
    
    integer :: i,j

    double precision :: vpqrs
    external vpqrs

    write(*,*) 'this is a test'

    write(*,*) 'print one integral:', vpqrs(1,1,1,1)

    write(*,*) 'print the density matrix'


    do i=1,nmo
       do j=1,nmo
          write(*,*) i, j, dens(i,j)
       enddo
    enddo
    

  end subroutine test_me
  
end module test
