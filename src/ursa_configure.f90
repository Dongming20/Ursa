program ursa_configure

    use tools

    implicit none


    integer, allocatable :: iarg,clen,stat
    character(len=:), allocatable :: file_id,file_name,charI,name,dgr_fem,line,typeofmeanfield,gw_sigma

    integer, allocatable :: degree

    integer :: i,p,nbat0,Nstates0
    logical :: file_exist
    character(len=2) :: atname


    allocate(iarg,clen,stat)
    allocate(character(len=10) :: file_id)
    allocate(character(len=50) :: file_name)
    allocate(character(len=10) :: charI)
    allocate(character(len=50) :: name)
    allocate(character(len=2) :: dgr_fem)
    allocate(character(len=50) :: line)
    allocate(character(len=3) :: typeofmeanfield)
    allocate(character(len=6) :: gw_sigma)
    allocate(degree)



        !!!!! Can read up to 1 argument command line
    iarg = command_argument_count()
    if (iarg/=4) then
        print *,'only 4 argument possible'
        stop
    end if
    call get_command_argument(1,name,clen,stat) ! name of the file
    call get_command_argument(2,typeofmeanfield,clen,stat) ! DFT or Hartree Fock
    call get_command_argument(3,gw_sigma,clen,stat) ! self-energy calculation methods
    call get_command_argument(4,dgr_fem,clen,stat) ! degree of polynomial of FEM







    INQUIRE(FILE=trim(name)//'.xyz', EXIST=file_exist)
    if (.not.(file_exist)) then
    nbat0=0
    print *,""//trim(name)//".xyz file not found"
    stop
    endif

!!!! Read atom from xyz file 
    open(10,file=trim(name)//'.xyz',status='old')
    read(10,*) nbat0
    if (nbat0>0) then !! no atom case also possible
    allocate(at(nbat0))
    read(10,*)
    end if
    do i=1,nbat0
    read(10,*) atname!,x,y,z

    p=1
    do while ((trim(PTABLE(p)%name)/=trim(atname)).and.(p<=LENGTH_PTABLE)) 
        p=p+1
    end do
    if (p>LENGTH_PTABLE) then
        print *,'Atom ',trim(atname),' unknown in database !!!'
        stop
    end if

    if (.not.(PTABLE(p)%def)) then
        print *,'Atom ',trim(atname),' not defined in database !!!'
        stop
    end if
    at(i)%core=p
    end do
    close(10)


    Nstates0=0
    do i=1,nbat0
        Nstates0=Nstates0+at(i)%core
    enddo
    Nstates0=Nstates0/2





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!!!!!!!!!!!!! Create the name.gw file    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *,'**Create .gw file'

    open(10,file=trim(adjustl(name))//'.gw',status='replace')


    !!!!! SIMULATION
    write(10,*) '&MODEL'
    write(10,*) "meanfield='"//trim(adjustl(typeofmeanfield))//"'                    !Options: 'dft', 'hf'(hartree-fock)"
    write(10,*) "gw_sigma='"//trim(adjustl(gw_sigma))//"'                      !Options: 'cd', 'casida', 'none'(no GW)"
    write(10,*) "dgr_fem='"//trim(adjustl(dgr_fem))//"'                       !Options: 'p2', 'p3'"
    write(10,*) "offset=6.0                           !Vaccum space of isolated system"
    write(10,*) '/'
    write(10,*) 


    !!!!! DFT
    write(10,*) '&DFT'
    write(10,*) "xc='pbe'                             !Options: 'lda','pbe'"
    write(10,*) "ks_M0="//trim(str1(2*Nstates0+10))//"                              !Estimated number of Eigenvalues"
    write(10,*) "ks_Emin=-50.0d0                      !Minimum of Eigenvalue interval in a.u."
    write(10,*) "ks_Emax=0.05d0                       !Maximum of Eigenvalue interval in a.u."
    write(10,*) "ks_alpha=0.5d0             "
    write(10,*) "ks_eps=1.0d-8             "
    write(10,*) "ks_poisson_bc=0                      !Option: 1 for integral, 0 for neutral "
    write(10,*) '/'
    write(10,*) 


    !!!!! Hartree-Fock
    write(10,*) '&HartreeFock'
    write(10,*) "hf_poisson_bc=0                      !Option: 1 for integral, 0 for neutral "
    write(10,*) '/'
    write(10,*) 


    !!!!! Contour-Deformation approach
    write(10,*) '&ContourDeformation'
    write(10,*) "cd_method='gs'                       !Options: 'gs','spectrum','nlfeast'"
    write(10,*) "cd_Hadamard='no'                     !Options: 'no','yes'"
    write(10,*) "cd_orbital=1                         !The orbital wants to solve"
    write(10,*) "cd_Emin_grid="//trim(str2(-24.5d0))//"                  !Energy grid minimum in (eV)"
    write(10,*) "cd_Emax_grid="//trim(str2(2.5d0))//"                    !Energy grid maximum in (eV)"
    write(10,*) "cd_N_grid=20                              "
    write(10,*) "cd_Emid=(-24.5,0.0)                  !Midpoint of the contour in FEAST in (eV)"
    write(10,*) "cd_r="//trim(str2(-24.5d0+2.5d0))//"                          !Radius of the contour in FEAST in (eV)"
    write(10,*) "cd_M0=3                              !Estimated number of Eigenvalues"
    write(10,*) "cd_igorbital=1                    "
    write(10,*) "cd_fpm1=1                     "
    write(10,*) "cd_fpm3=10                    "
    write(10,*) "cd_fpm8=8                     "
    write(10,*) "cd_fpm45=2                    "
    write(10,*) '/'
    write(10,*) 


    !!!!! Fully-analytic approach
    write(10,*) '&FullyAnalytic'
    write(10,*) "casida_method='gs'                   !Options: 'gs','nlfeast'"
    write(10,*) "casida_N_empty=1500                  !Default is 1500 "
    write(10,*) "casida_orbital1=1                     !The orbital wants to solve"
    write(10,*) "casida_Emin_grid1="//trim(str2(-24.5d0))//"              !Energy grid minimum in (eV)"
    write(10,*) "casida_Emax_grid1="//trim(str2(-22.5d0))//"                !Energy grid maximum in (eV)"
    write(10,*) "casida_orbital2=2                     !The orbital wants to solve"
    write(10,*) "casida_Emin_grid2="//trim(str2(0.0d0))//"              !Energy grid minimum in (eV)"
    write(10,*) "casida_Emax_grid2="//trim(str2(2.5d0))//"                !Energy grid maximum in (eV)"
    write(10,*) "casida_N_grid=20                              "
    write(10,*) "casida_Emid=(-24.5,0.0)              !Midpoint of the contour in FEAST in (eV)"
    write(10,*) "casida_r="//trim(str2(-24.5d0+2.5d0))//"                      !Radius of the contour in FEAST in (eV)"
    write(10,*) "casida_M0=3                          !Estimated number of Eigenvalues"
    write(10,*) "casida_igorbital=1                    "
    write(10,*) "casida_fpm1=1                     "
    write(10,*) "casida_fpm3=10                    "
    write(10,*) "casida_fpm8=8                     "
    write(10,*) "casida_fpm45=2                    "
    write(10,*) '/'
    write(10,*) 


    close(10)   






end program
