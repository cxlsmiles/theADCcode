module test

  
  use parameters
  use iso_c_binding
 
  implicit none  

  type(c_ptr) :: al;

contains
  
  subroutine test_me()

    
    write(*,*) "test me"
    
    write(*,*) roccnum(10)
   


  end subroutine test_me



  subroutine init_ptr()

    write(*,*) "initptr me"
 

    al = c_loc(roccnum);

  end subroutine init_ptr

  subroutine copy_ptr(cp, sz)

    type(c_ptr) :: cp
    integer, dimension(:), pointer :: loc
    integer :: sz

    write(*,*) "copy ptr"
    
    call c_f_pointer(cp, roccnum, (/sz/))


    
    
    write(*,*) roccnum


  end subroutine copy_ptr


    
end module test
          
          
       
       
       
       
    
    
    
       
       
       
       
       
       
    
    
    
    

    
    
    
    
