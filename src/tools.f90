module tools

    implicit none

    !!!! constant 2018
    !!!! https://physics.nist.gov/cuu/Constants/index.html
    double precision,parameter :: pi=4.0d0*atan(1.0d0)
    !double precision,parameter :: echarge=1.602176634d-19 ! 1.602176565d-19 
    !double precision,parameter :: epsilon0=8.8541878128d-12 !8.854187817620d-12
    !double precision,parameter :: ao=5.29177210903d-11
    !double precision,parameter :: Rydberg=13.605693122994d0
    !double precision,parameter :: hplanck=6.62607015d-34
    !double precision,parameter :: hbar=hplanck/(2.0d0*pi)
    !double precision,parameter :: me=9.1093837015d-31 !9.10938291d-31 !9.1093897D-31

    double precision,parameter :: bohr=0.529177249d0 !angstrom for atomic units
    double precision,parameter :: hartree=27.2113961d0 !eV for atomic units
    double precision,parameter :: radius=50.0d0/bohr !atoms' radius limit 

    double precision,parameter :: a0=1.0d-30

    complex(kind=(kind(1.0d0))),parameter :: ZZERO=(0.0d0,0.0d0),ZONE=(1.0d0,0.0d0)

    integer,parameter :: LENGTH_PTABLE=86 ! incomplete list (118 ! number of known atoms)
    ! integer,dimension(20),parameter :: FILLED_NLSTATE=[2,2,6,2,6,2,10,6,2,10,6,2,14,10,6,2,14,10,6,2]
    ! !1s2,2s2,2p6,3s2,3p6,4s2,3d10,4p6,5s2,4d10,5p6,6s2,4f14,5d10,6p6,7s2,5f14,6d10,7p6,8s2

    double precision,dimension(:,:),allocatable :: r_init,n_init
    integer,dimension(:),allocatable :: Nn_init
    integer,allocatable :: io,lb,ub

    double precision,allocatable :: r0,t


    !!! atomic number -name- weigth - radius calculated - radius empirical - filled-core states (considered)
    !!!  # https://www.ptable.com/
    !!! Rq: when empirical radius is unknown, we make it equal to calculated (3 digits with 0)
    type atomic_element
        character(len=2) :: name
        logical :: def=.false. ! defined (in database) yes/no
        double precision :: w=0.0d0,rcalc=0.0d0,remp=0.0d0
        integer :: filled_core=0
    end type atomic_element

    type(atomic_element),dimension(LENGTH_PTABLE),parameter :: PTABLE = &
        [(atomic_element("H",.true.,1.008d0,0.53d0,0.25d0,0)),&     !1
        (atomic_element("He",.true.,4.0026d0,0.31d0,0.310d0,0)),&   !2
        (atomic_element("Li",.true., 6.94d0,1.67d0, 1.45d0,2)),& !3
        (atomic_element("Be",.true., 9.0122d0,1.12d0, 1.05d0,2)),& !4
        (atomic_element("B",.true., 10.81d0, 0.87d0, 0.85d0,2)),&  !5
        (atomic_element("C",.true., 12.011d0, 0.67d0, 0.70d0,2)),&  !6
        (atomic_element("N",.true., 14.007d0, 0.56d0, 0.65d0,2)),&  !7
        (atomic_element("O",.true., 15.999d0, 0.48d0, 0.60d0,2)),&  !8
        (atomic_element("F",.true., 18.998d0, 0.42d0, 0.50d0,2)),&  !9
        (atomic_element("Ne",.true., 20.180d0, 0.38d0, 0.380d0,2)),& !10
        (atomic_element("Na",.true., 22.990d0, 1.9d0, 1.80d0,10)),&   !11
        (atomic_element("Mg",.true., 24.305d0, 1.45d0, 1.50d0,10)),&  !12
        (atomic_element("Al",.true., 26.982d0,1.18d0, 1.25d0,10)),&  !13
        (atomic_element("Si",.true., 28.085d0,1.11d0, 1.10d0,10)),&  !14
        (atomic_element("P",.true., 30.974d0,0.98d0, 1.00d0,10)),&   !15
        (atomic_element("S",.true., 32.06d0,0.88d0, 1.00d0,10)),&   !16
        (atomic_element("Cl",.true., 35.45d0,0.79d0, 1.00d0,10)),&  !17
        (atomic_element("Ar",.true., 39.948d0,0.71d0, 0.71d0,10)),&  !18
        (atomic_element("K",.true., 39.098d0,2.43d0, 2.20d0,18)),&   !19
        (atomic_element("Ca",.true., 40.078d0,1.94d0, 1.80d0,18)),&  !20
        (atomic_element("Sc")),&  !21
        (atomic_element("Ti")),&  !22
        (atomic_element("V")),&  !23
        (atomic_element("Cr")),&  !24
        (atomic_element("Mn")),&  !25
        (atomic_element("Fe")),&  !26
        (atomic_element("Co")),&  !27
        (atomic_element("Ni")),&  !28
        (atomic_element("Cu")),& !29
        (atomic_element("Zn")),& !30
        (atomic_element("Ga",.true.,69.723d0, 1.36d0, 1.30d0,18)),&  !31
        (atomic_element("Ge",.true., 72.630,1.25d0, 1.25d0,18)),&  !32
        (atomic_element("As",.true., 74.922d0,1.14d0, 1.15d0,18)),&  !33
        (atomic_element("Se",.true., 78.971d0,1.03d0, 1.15d0,18)),&  !34
        (atomic_element("Br",.true., 79.904d0,0.94d0, 1.15d0,18)),&  !35
        (atomic_element("Kr",.true.,83.798d0,0.88d0, 0.880d0,18)),&  !36
        (atomic_element("Rb",.true., 85.468d0,2.65d0, 2.35d0,28)),&  !37
        (atomic_element("Sr",.true., 87.62d0,2.19d0, 2.00d0,28)),&  !38
        (atomic_element("Y")),&  !39
        (atomic_element("Zr")),&  !40
        (atomic_element("Nb")),&   !41
        (atomic_element("Mo")),&  !42
        (atomic_element("Tc")),&  !43
        (atomic_element("Ru")),&  !44
        (atomic_element("Rh")),&  !45
        (atomic_element("Pd")),&  !46
        (atomic_element("Ag")),& !47
        (atomic_element("Cd")),& !48
        (atomic_element("In",.true., 114.82d0,1.56d0, 1.55d0,28)),&  !49
        (atomic_element("Sn",.true., 118.71d0,1.45d0, 1.45d0,28)),&  !50
        (atomic_element("Sb",.true., 121.76d0,1.33d0, 1.45d0,28)),&  !51
        (atomic_element("Te",.true., 127.60d0,1.23d0, 1.40d0,28)),&  !52
        (atomic_element("I",.true., 126.90d0,1.15d0, 1.40d0,28)),&   !53
        (atomic_element("Xe",.true., 131.29d0,1.08d0, 1.080d0,28)),& !54
        (atomic_element("Cs",.true., 132.91d0,2.98d0, 2.60d0,36)),&  !55
        (atomic_element("Ba",.true., 137.33d0,2.53d0, 2.15d0,36)),&  !56
        (atomic_element("La")),&  !57
        (atomic_element("Ce")),&  !58
        (atomic_element("Pr")),&   !59
        (atomic_element("Nd")),&  !60
        (atomic_element("Pm")),&  !61
        (atomic_element("Sm")),&  !62
        (atomic_element("Eu")),&  !63
        (atomic_element("Gd")),&  !64
        (atomic_element("Tb")),& !65
        (atomic_element("Dy")),& !66
        (atomic_element("Ho")),&  !67
        (atomic_element("Er")),&  !68
        (atomic_element("Tm")),&   !69
        (atomic_element("Yb")),&  !70
        (atomic_element("Lu")),&  !71
        (atomic_element("Hf")),&  !72
        (atomic_element("Ta")),&  !73
        (atomic_element("W")),&  !74
        (atomic_element("Re")),& !75
        (atomic_element("Os")),& !76
        (atomic_element("Ir")),&  !77
        (atomic_element("Pt")),&  !78
        (atomic_element("Au")),&  !79
        (atomic_element("Hg")),& !80
        (atomic_element("Tl",.true., 204.38d0,1.56d0, 1.90d0,46)),&  !81
        (atomic_element("Pb",.true., 207.2d0,1.54d0, 1.80d0,46)),&  !82
        (atomic_element("Bi",.true., 208.98d0,1.43d0, 1.60d0,46)),&  !83
        (atomic_element("Po",.true., 209d0, 1.35d0,1.90d0,46)),&  !84
        (atomic_element("At",.true., 210d0, 1.27d0, 1.270d0,46)),& !85
        (atomic_element("Rn",.true., 222d0, 1.20d0,1.200d0,46))]  !86



    type atom
        double precision,dimension(3) :: c
        integer :: core
    end type atom

    type(atom),dimension(:),allocatable :: at !! number of atoms
    integer :: Nbat
    double precision,dimension(:,:),allocatable :: point_muffin,rbox,dedge
    integer,dimension(:,:),allocatable :: face_muffin
    double precision :: offset ! 6.0d0
    double precision :: r,ro
    character(len=10) ::qtetgen="1.5" !"1.5" !(mesh refinement default is 2 in tetgen see Fig 15 tetgen documentation)
    !!! Personal remark: look at fil vpot..*pdf, gaussian 3d for c atom mesh
    !!! qtetgen=1.5 P2 mesh is ~1500 ->not so good result for gaussian reproduction
    !!!qtetgen=1 P2 is ~7500 great but too big
    !!! qtetgen=1.25 P2 is ~3000 great result




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!  .gw file !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    double precision :: offset0
    character(len=3) :: xc0
    integer :: ks_M00, ks_poisson_bc0
    double precision :: ks_Emin0,ks_Emax0,ks_alpha0,ks_eps0
    integer :: hf_poisson_bc0

    character(len=8) :: cd_method0,cd_Hadamard0
    integer :: cd_N_grid0,cd_M00,cd_fpm10,cd_fpm30,cd_fpm80,cd_fpm450,cd_igorbital0,cd_orbital0
    double precision :: cd_Emin_grid0,cd_Emax_grid0,cd_r0
    complex(kind=(kind(1.0d0))) :: cd_Emid0

    character(len=7) :: casida_method0
    integer :: casida_N_grid0,casida_M00,casida_fpm10,casida_fpm30,casida_fpm80,casida_fpm450,casida_igorbital0,casida_N_empty0
    integer :: casida_orbital10,casida_orbital20
    double precision :: casida_Emin_grid10,casida_Emax_grid10,casida_Emin_grid20,casida_Emax_grid20,casida_r0
    complex(kind=(kind(1.0d0))) :: casida_Emid0
    
    

    contains

    subroutine load_xyz(name,scaling)
        character(len=*),intent(in) ::name
        double precision,intent(in) :: scaling
        ! integer,intent(out) :: nbat
        ! type(atom),dimension(:),allocatable,intent(out) :: at
        integer :: i,p
        character(len=2) :: atname
        double precision :: x,y,z
        logical :: file_exist

        !print *,trim(name)//'.xyz'
        !!! check if file xyz exists    
            INQUIRE(FILE=trim(name)//'.xyz', EXIST=file_exist)
            if (.not.(file_exist)) then
            nbat=0
            return
            endif

        !!!! Read atom from xyz file 
            open(10,file=trim(name)//'.xyz',status='old')
            read(10,*) nbat
            if (nbat>0) then !! no atom case also possible
            allocate(at(nbat))
            read(10,*)
            end if
            do i=1,nbat
            read(10,*) atname,x,y,z
            at(i)%c(1)=x*scaling
            at(i)%c(2)=y*scaling
            at(i)%c(3)=z*scaling
            !find atom ID in database
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
    end subroutine load_xyz



    


    subroutine load_gw(name,dgr_fem0,meanfield0,gw_sigma0)

      !!! PARAMETER MODEL !!!
      character(len=*), intent(in) :: name
      character(len=*), intent(inout) :: dgr_fem0,meanfield0,gw_sigma0
      integer :: rc
      character(len=25) :: dgr_fem,meanfield,gw_sigma
      double precision :: offset
      namelist /MODEL/ dgr_fem,meanfield,gw_sigma,offset


      !!! PARAMETER DFT !!!
      character(len=25) :: xc                              !Options: 'lda,pbe'
      integer :: ks_M0, ks_poisson_bc
      double precision :: ks_Emin,ks_Emax,ks_alpha,ks_eps
      namelist /DFT/ xc,ks_M0,ks_poisson_bc,ks_Emin,ks_Emax,ks_alpha,ks_eps


      !!! PARAMETER Hartree-Fock !!!
      integer :: hf_poisson_bc                        !Option: 1 for integral, 0 for neutral
      namelist /HartreeFock/ hf_poisson_bc


      !!! PARAMETER Contour-Deformation approach !!!
      character(len=25) :: cd_method,cd_Hadamard                         !Options: 'gs,spectrum,nlfeast'
      integer :: cd_N_grid,cd_M0,cd_fpm1,cd_fpm3,cd_fpm8,cd_fpm45,cd_igorbital,cd_orbital
      double precision :: cd_Emin_grid,cd_Emax_grid,cd_r
      complex(kind=(kind(1.0d0))) :: cd_Emid
      namelist /ContourDeformation/ cd_method,cd_Hadamard,cd_N_grid,cd_M0,cd_fpm1,cd_fpm3,cd_fpm8,cd_fpm45,cd_igorbital&
      ,cd_Emin_grid,cd_Emax_grid,cd_r,cd_Emid,cd_orbital


      !!! PARAMETER Fully-Analytic approach !!!
      character(len=25) :: casida_method                       !Options: 'gs,spectrum,nlfeast'
      integer :: casida_N_grid,casida_M0,casida_fpm1,casida_fpm3,casida_fpm8,casida_fpm45,casida_igorbital,casida_N_empty
      integer :: casida_orbital1,casida_orbital2
      double precision :: casida_Emin_grid1,casida_Emax_grid1,casida_Emin_grid2,casida_Emax_grid2,casida_r
      complex(kind=(kind(1.0d0))) :: casida_Emid
      namelist /FullyAnalytic/ casida_method,casida_N_empty,casida_N_grid,casida_M0,casida_fpm1,casida_fpm3,casida_fpm8&
      ,casida_fpm45,casida_igorbital,casida_Emin_grid1,casida_Emax_grid1,casida_r,casida_Emid,casida_orbital1&
      ,casida_orbital2,casida_Emin_grid2,casida_Emax_grid2



      open(10,file=trim(adjustl(name))//'.gw',status='old')
        read (unit=10, nml=MODEL, iostat=rc)
        if (rc /= 0) print *,"Error: invalid Namelist format for MODEL"
        read (unit=10, nml=DFT, iostat=rc)
        if (rc /= 0) print *,"Error: invalid Namelist format for DFT"
        read (unit=10, nml=HartreeFock, iostat=rc)
        if (rc /= 0) print *,"Error: invalid Namelist format for HartreeFock"
        read (unit=10, nml=ContourDeformation, iostat=rc)
        if (rc /= 0) print *,"Error: invalid Namelist format for ContourDeformation"
        read (unit=10, nml=FullyAnalytic, iostat=rc)
        if (rc /= 0) print *,"Error: invalid Namelist format for FullyAnalytic"
      close(10)


          !!! MODEL !!!
          dgr_fem0=trim(dgr_fem)
          meanfield0=trim(meanfield)
          gw_sigma0=trim(gw_sigma)
          offset0=offset



          !!! DFT !!!
          xc0=trim(xc)
          ks_M00=ks_M0
          ks_poisson_bc0=ks_poisson_bc
          ks_Emin0=ks_Emin
          ks_Emax0=ks_Emax
          ks_alpha0=ks_alpha
          ks_eps0=ks_eps

          !!! HF !!!
          hf_poisson_bc0=hf_poisson_bc

          !!! CD !!!
          cd_method0=trim(cd_method)
          cd_Hadamard0=trim(cd_Hadamard)
          cd_N_grid0=cd_N_grid
          cd_orbital0=cd_orbital
          cd_Emin_grid0=cd_Emin_grid
          cd_Emax_grid0=cd_Emax_grid
          cd_M00=cd_M0
          cd_Emid0=cd_Emid
          cd_r0=cd_r
          cd_fpm10=cd_fpm1
          cd_fpm30=cd_fpm3
          cd_fpm80=cd_fpm8
          cd_fpm450=cd_fpm45
          cd_igorbital0=cd_igorbital
          

          !!! Casida !!!
          casida_method0=trim(casida_method)
          casida_N_empty0=casida_N_empty
          casida_orbital10=casida_orbital1
          casida_Emin_grid10=casida_Emin_grid1
          casida_Emax_grid10=casida_Emax_grid1
          casida_orbital20=casida_orbital2
          casida_Emin_grid20=casida_Emin_grid2
          casida_Emax_grid20=casida_Emax_grid2
          casida_N_grid0=casida_N_grid
          casida_M00=casida_M0
          casida_Emid0=casida_Emid
          casida_r0=casida_r
          casida_igorbital0=casida_igorbital
          casida_fpm10=casida_fpm1
          casida_fpm30=casida_fpm3
          casida_fpm80=casida_fpm8
          casida_fpm450=casida_fpm45


    end subroutine load_gw





    function comp(a,b) result(c)
        integer,dimension(:),intent(in) :: a,b
        integer :: m,n,c
        c=0
        do m=1,2
            do n=1,2
            if (b(n)==a(m)) then
                c=c+1
            endif
            enddo
        enddo
    end function comp


    recursive subroutine quicksort_eigenpairs(a, b, c, n, left, right)
    complex(kind=(kind(1.0d0))), dimension(:), intent(inout) :: a
    complex(kind=(kind(1.0d0))), dimension(:,:), intent(inout) :: b
    double precision, dimension(:), intent(inout) :: c
    integer, intent(in) :: n, left, right
    integer, allocatable :: i, j
    double precision, allocatable :: pivot
    complex(kind=(kind(1.0d0))),allocatable :: temp
    complex(kind=(kind(1.0d0))), dimension(:), allocatable :: temp2
    double precision, allocatable :: temp3

    allocate(i,j)
    allocate(pivot)
    allocate(temp)
    allocate(temp2(n))
    allocate(temp3)

    if (left < right) then
       pivot = dble(a((left + right) / 2))
       i = left
       j = right
       do
          do while (dble(a(i)) < pivot)
             i = i + 1
          end do
          do while (dble(a(j)) > pivot)
             j = j - 1
          end do
          if (i <= j) then
             temp = a(i)
             a(i) = a(j)
             a(j) = temp
             temp2 = b(:,i)
             b(:,i) = b(:,j)
             b(:,j) = temp2
             temp3 = c(i)
             c(i) = c(j)
             c(j) = temp3
             i = i + 1
             j = j - 1
          end if
          if (i > j) exit
       end do
       call quicksort_eigenpairs(a, b, c, n, left, j)
       call quicksort_eigenpairs(a, b, c, n, i, right)
    end if

    ! deallocate(temp2)

    deallocate(i,j)
    deallocate(pivot)
    deallocate(temp)
    deallocate(temp2)
    deallocate(temp3)

  end subroutine quicksort_eigenpairs


  function outer_product(u, v) result(res)
    double precision, intent(in) :: u(:), v(:)
    double precision :: res(size(u),size(v))
    integer :: col
    do col = 1, size(v)
      res(:,col) = v(col) * u
    end do
end function outer_product


function outer_product_complex(u, v) result(res)
    complex(kind=(kind(1.0d0))), intent(in) :: u(:),v(:)
    complex(kind=(kind(1.0d0))) :: res(size(u),size(v))
    integer :: col
    do col = 1, size(v)
      res(:,col) = v(col) * u
    end do
end function outer_product_complex





subroutine iterative_solver(B0, rhs, x, n, m0, tol, max_iter)
    implicit none
    double precision, dimension(:,:), intent(in) :: B0  ! Matrix B
    double precision, dimension(:,:), intent(in) :: rhs    ! Right-hand side vector
    double precision, dimension(:,:), intent(inout) :: x  ! Solution vector
    integer, intent(in) :: n,m0            ! Size of matrix and vectors
    double precision, intent(in) :: tol                ! Tolerance for convergence
    integer, intent(in) :: max_iter           ! Maximum number of iterations
    double precision, dimension(n,m0) :: x_old,x_temp ! Store previous iteration solution
    double precision, dimension(n,m0) :: residual        ! Residual error
    integer :: iter                           ! Iteration counter
    integer :: i, j                          ! Loop indices
    
    ! Initialize x_old
    ! allocate(x_old(n))
    ! allocate(x_temp(n))
    x_old = 0.0
    
    ! Perform iterative method
    do iter = 1, max_iter
        ! Update x using the equation x = B*x + b
        ! do i = 1, n
        !     x(i) = b(i)
        !     do j = 1, n
        !         ! if (j /= i) then
        !             x(i) = x(i) + B0(i,j) * x_old(j)
        !         ! end if
        !     end do
        ! end do

        call DSYMM('L','L',n,1,1.0d0,B0,n,x_old,n,0.0d0,x_temp,n)

        x=x_temp+rhs
        
        ! Check convergence
        residual = abs(x - x_old)
        if (maxval(residual) < tol) exit
        
        ! Update x_old for next iteration
        x_old = x
    end do
    
    ! Deallocate temporary array
    ! deallocate(x_old)
    
    ! Output convergence information
    if (iter <= max_iter) then
        print *, "Converged in ", iter, " iterations."
    else
        print *, "Did not converge within maximum iterations."
    end if
    
end subroutine iterative_solver


subroutine find_interval(arr, n, value, lower_bound, upper_bound)
    ! Input variables
    double precision, intent(in) :: arr(:)  ! A sorted array of double precision numbers
    integer, intent(in) :: n       ! The size of the array
    double precision, intent(in) :: value   ! The number to be located within the array
    
    ! Output variables
    integer, intent(out) :: lower_bound, upper_bound  ! The bounds between which 'value' lies

    ! Local variables
    integer :: i

    ! Check if the value is out of the array bounds
    if (value < arr(1)) then
      lower_bound = -1.0
      upper_bound = 1
      return
    elseif (value > arr(n)) then
      lower_bound = n
      upper_bound = -1.0
      return
    end if

    ! Initialize the bounds to default values
    lower_bound = -1.0
    upper_bound = -1.0

    ! Loop through the array to find the interval
    do i = 1, n-1
      if (value >= arr(i) .and. value <= arr(i+1)) then
        ! lower_bound = arr(i)
        ! upper_bound = arr(i+1)
        lower_bound = i
        upper_bound = i+1
        return
      end if
    end do

  end subroutine find_interval


  function inv(A) result(Ainv)
    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(size(A,1),size(A,2)) :: Ainv
  
    double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
  
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI
  
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
  
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)
  
    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if
  
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
  
    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv


  recursive subroutine quicksort_integer(arr, left, right)
    integer, intent(inout) :: arr(:)
    integer, intent(in) :: left, right
    integer :: pivot_index

    if (left < right) then
      ! Partition the array and get the pivot index
      pivot_index = partition(arr, left, right)

      ! Recursively apply quicksort to the two sub-arrays
      call quicksort_integer(arr, left, pivot_index - 1)
      call quicksort_integer(arr, pivot_index + 1, right)
    end if
  end subroutine quicksort_integer

  function partition(arr, left, right) result(pivot_index)
    integer, intent(inout) :: arr(:)
    integer, intent(in) :: left, right
    integer :: pivot_index, pivot_value, i, j, temp

    pivot_value = arr(right)
    i = left - 1

    ! Loop through the array to partition it
    do j = left, right - 1
      if (arr(j) <= pivot_value) then
        i = i + 1
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
      end if
    end do

    ! Place the pivot element in its correct position
    temp = arr(i + 1)
    arr(i + 1) = arr(right)
    arr(right) = temp

    pivot_index = i + 1
  end function partition


  character(len=20) function str1(k)
  !   "Convert an integer to string."
      integer, intent(in) :: k
      write (str1, *) k
      str1 = adjustl(str1)
  end function str1

  character(len=20) function str2(k)
  !   "Convert an integer to string."
      double precision, intent(in) :: k
      write (str2, '(f10.2)') k
      str2 = adjustl(str2)
  end function str2


end module tools
