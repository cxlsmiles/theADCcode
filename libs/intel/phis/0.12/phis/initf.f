  function phis_init(arg,project)
!#######Interface to C-function's
        interface 
        function phis_init_old(arg)
          integer,intent(in)            :: arg
          integer                       :: phis_init_old
        end function phis_init_old
        
        function check(path)
          logical                       :: check
          character(len=*)              :: path
        end function check
        end interface        
!#######Global variables   
        integer,intent(in)                              :: arg           !! Should be always present
        character(len=*),intent(in),optional            :: project       !! New feature, can absent in calling 
!#######Local variables
        character(len=4096)                             :: buffer
!#######Main code
        if (present(project).and.(len_trim(project)/=0)) then 
           buffer=project
        else
           call getcwd(buffer)
        endif

        if (check(buffer)) then
            call chdir(buffer)
            phis_init=phis_init_old(arg)                     
        endif

        buffer=""
  end function phis_init

