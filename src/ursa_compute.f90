program ursa_compute

    use class_linkedlist
    use tools
    use basisfunctions
    use potentials
    use potentials_gw
    

    implicit none

    ! include 'mpif.h'

    integer, allocatable :: dummy,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,Ne,Nn,Nlocal,Nquadrature,Nbc
    integer, allocatable :: i,k,l,ii,ll,m,n,jj,kk
    integer, allocatable :: Nquadrature1D,Nquadrature1D_range
    integer, allocatable :: Nmuffin,Nface,Nmuffin_total,Nface_total
    logical, allocatable :: mtest
    logical :: module_exists,found

    double precision,allocatable :: dummy_real,value0
    


    !!!!!!!! DFT initial guess relates variables 
    double precision,dimension(:), allocatable :: psi_0,psi_100,psi_200,psi_210,psi_21x,psi_21y,psi_21z
    double precision,dimension(:), allocatable :: psi_300,psi_31x,psi_31y,psi_31z
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: psi_211_minus,psi_211_plus
    double precision,dimension(:,:), allocatable :: psi_00,psi_01
    


    double precision, allocatable :: Z0, rxyz

    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: xy_AB,xy_V!,ab



    !!!!!!! reading file name related variables
    ! character(len=10) :: file_id
    ! character(len=50) :: file_name
    ! character(len=10) :: charI

    integer, allocatable :: iarg,clen,stat
    ! character(len=50) ::name,line

    character(len=:), allocatable :: file_id,file_name,charI,name,dgr_fem,line,meanfield,gw_sigma

    integer, allocatable :: degree
    

    
    





    !! input parameters for FEAST

    ! character(len=1) :: UPLO='F' 
    character(len=:), allocatable :: UPLO

    integer, dimension(:), allocatable :: fpm
    integer, allocatable :: M0,M0_GW ! search subspace dimension
    double precision, allocatable :: Emin,Emax,Emin_GW,Emax_GW ! search interval
    !! output variables for FEAST
    double precision, dimension(:), allocatable :: E,res,E_1,res_1,E_2,res_2,res_GW
    double precision, allocatable :: epsout
    integer , allocatable:: loop, info, M00


    !!!!! complex symmetric FEAST !!!!
    ! double precision :: r_zs
    complex(kind=(kind(1.0d0))), allocatable :: Emid
    complex(kind=(kind(1.0d0))), dimension(:), allocatable :: E_complex,res_complex

    !!!!!!!!!!!!!!!!! Others
    ! integer,dimension(64) :: fpm 
    ! double precision :: epsout
    ! integer :: loop
    ! integer :: i,k,j
    integer, allocatable :: M0_nl,M_nl!,info
    double precision, allocatable :: r_nl
    complex(kind=(kind(1.0d0))), allocatable :: Emid_nl
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: E_nl,E_nl_temp ! eigenvalues
    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: X_nl,X_nl_temp! eigenvectors
    double precision,dimension(:),allocatable :: res_nl,res_nl_temp ! eigenvalue
    !!! RCI
     integer, allocatable :: id,ijob,infoloc2,infoloc1,dmax,nfact,nnza
     complex(kind=(kind(1.0d0))) ,dimension(:,:),allocatable :: Aztemp_nl
     complex(kind=(kind(1.0d0))) ,dimension(:,:,:),allocatable :: Az_nl
     complex(kind=(kind(1.0d0))) ,dimension(:),allocatable :: sAztemp
     complex(kind=(kind(1.0d0))) ,dimension(:,:),allocatable :: saz,psaz
     complex ,dimension(:,:,:),allocatable :: Ac_nl
     complex(kind=(kind(1.0d0))), allocatable :: Ze_nl
     complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork_nl,work_nl,zBq_nl,zAq_nl,zaux
     double precision, dimension(:),allocatable :: dwork_nl
     character(len=:), allocatable :: UPLO2
     integer ,dimension(:,:),allocatable :: ipivloc_nl
     double precision,external :: zlange
     
    !  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne
    complex, dimension(:,:),allocatable ::cwork_nl
    real ,dimension(:,:,:),allocatable :: sA_nl



    !!!!!!! iterative solver
    !!!! for bicgstab
    logical, allocatable :: comb,com
    double precision, dimension(:),allocatable::nres,norm_nl
    integer, allocatable :: linloops
    double precision, allocatable :: lintargeterror
    integer(8), allocatable  :: fout
    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: bwork1,bwork2
    integer, allocatable :: ijob2,jj1,jj2,infoloc3
    




    !!!!!!!!!! DSYEV 
    double precision,dimension(:),allocatable :: DSYEV_work
    integer, allocatable :: DSYEV_lwork
    integer :: DSYEV_info
    double precision,dimension(:),allocatable :: E_hf




    double precision, dimension(:), allocatable :: V,V_i, H

    integer, allocatable :: neigh_size,Nstates0,Nstates,Nempty
    integer,dimension(:),allocatable :: Nempty_list


    !!!!!!! SCF loop !!!!!!!
    integer, allocatable :: itmax,M9,it,itm,nrestart,info_scf,ite
    double precision, allocatable :: alpha,alpha0,eps,norm,normf
    double precision,dimension(:),allocatable :: ni_old,ri,ri_old,res_scf,f,work,nq0,temp1,ri_g,ri_old_g,ni_old_g
    double precision,dimension(:,:),allocatable :: B00,dN,dR,dN_trans,dN_g,dR_g
    integer,dimension(:),allocatable:: ipiv
    double precision,dimension(:,:),allocatable :: nq_temp,psi_temp
    double precision,dimension(:),allocatable :: temp2


    double precision, allocatable :: Eh,Ex,Ec,Ehfx
    integer, allocatable :: HOMO_state,LUMO_state,orbital
    double precision, allocatable :: E_HOMO_DFT



    

    !!!! convert to p2 basis
    integer, allocatable :: Nn0,Ne0,Ng0,Nf0,g01,g02,g03,g04,edge1,edge2,face1,face2,face3,dummy0,Nt0
    integer,dimension(:,:),allocatable :: ele0,edge0,face0,node0
    double precision, dimension(:,:),allocatable :: point0,pointe0
    integer,dimension(:),allocatable :: color0, colore0, colore_face0
    double precision, allocatable :: x_min,x_max,y_min,y_max,z_min,z_max
    
    CHARACTER(len=:), allocatable :: rootdir

    



    !!!! MPI parameters 
     
    !  integer :: code, nb_procs, rank,NEW_COMM_WORLD
    integer, allocatable :: code,rank,nb_procs,lrank,lnb_procs,color_mpi,key,NEW_COMM_WORLD
    complex(kind=(kind(1.0d0))), allocatable :: sum1,sum2
    character(len=:), allocatable :: cnL3
    integer, allocatable :: nL3
    double precision, allocatable :: start_time,finish_time,start_time0,finish_time0,start_time1,finish_time1&
    ,start_time000,finish_time000

    ! integer, dimension(MPI_STATUS_SIZE) :: mpi_recv_status
    integer, allocatable :: M000
    integer, parameter :: tag=100


    double precision,dimension(:),allocatable :: psi_ks00,psi_gw00,psi_ks01,psi_gw01,psi_ks02,psi_gw02
    integer :: dft00,gw00

    

    !!!!!!!!!! dynamic allocation for scalar 
    allocate(dummy,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,Ne,Nn,Nlocal,Nquadrature,Nbc)
    allocate(i,k,l,ii,ll,m,n,jj,kk)
    allocate(Nquadrature1D,Nquadrature1D_range)
    allocate(Nmuffin,Nface,Nmuffin_total,Nface_total)
    allocate(mtest)
    allocate(dummy_real)


    allocate(d,psi_square_100,psi_square_200,psi_square_21x,psi_square_21y,psi_square_21z)
    allocate(psi_square_300,psi_square_31x,psi_square_31y,psi_square_31z)
    allocate(Z0, rxyz)


    allocate(character(len=10) :: file_id)
    allocate(character(len=50) :: file_name)
    allocate(iarg,clen,stat)
    allocate(character(len=10) :: charI)
    allocate(character(len=50) :: name)
    allocate(character(len=50) :: line)
    allocate(character(len=2) :: dgr_fem)
    allocate(character(len=3) :: meanfield)
    allocate(character(len=6) :: gw_sigma)
    allocate(degree)


    allocate(character(len=1) :: UPLO)
    UPLO='F'
    allocate(fpm(64))

    allocate(M0,M0_GW)
    allocate(Emin,Emax,Emin_GW,Emax_GW)
    allocate(epsout)
    allocate(loop, info, M00)

    allocate(Emid)

    allocate(M0_nl,M_nl)

    allocate(r_nl)
    allocate(Emid_nl)
    allocate(id,ijob,infoloc2,infoloc1,dmax,nfact,nnza)

    allocate(Ze_nl)
    allocate(character(len=1) :: UPLO2)
    UPLO2='U'

    allocate(comb,com)

    allocate(linloops)
    allocate(lintargeterror)
    allocate(fout)
    allocate(ijob2,jj1,jj2,infoloc3)

    allocate(DSYEV_lwork)

    allocate(neigh_size,Nstates0,Nstates,Nempty)

    allocate(itmax,M9,it,itm,nrestart,info_scf,ite)
    allocate(alpha,alpha0,eps,norm,normf)

    allocate(Eh,Ex,Ec,Ehfx)
    allocate(HOMO_state,LUMO_state,orbital)
    allocate(E_HOMO_DFT)

    allocate(Nn0,Ne0,Ng0,Nf0,g01,g02,g03,g04,edge1,edge2,face1,face2,face3,dummy0,Nt0)
    allocate(x_min,x_max,y_min,y_max,z_min,z_max)

    allocate(character(len=255) :: rootdir)

    allocate(code,rank,nb_procs,lrank,lnb_procs,color_mpi,key,NEW_COMM_WORLD)
    allocate(sum1,sum2)

    allocate(character(len=3) :: cnL3)
    allocate(nL3)

    allocate(start_time,finish_time,start_time0,finish_time0,start_time1,finish_time1,start_time000,finish_time000)
    allocate(M000)


    ! !!!!!! tools module scalar dynamic allocation 
    ! allocate(Nn_init)


    !!!!!! basisfunction module scalar dynamic allocation 
    allocate(Jdet,Jdet1,Jdet2,Jdet20)

    allocate(phi_p2(10))
    allocate(phi_p2_del(10,3))
    allocate(phi_p3(20))
    allocate(phi_p3_del(20,3))

    allocate(J0(3,3),Jit(3,3),J1(3,3),Jit1(3,3),J2(3,3),Jit2(3,3),J20(3,3),Jit20(3,3))
    allocate(p1_matrix(4,4))
    allocate(p2_matrix(10,10))
    allocate(p3_matrix(20,20))
    allocate(gpoint(11,3))
    allocate(gweight(11))
    allocate(localpointp2(10,3))

    allocate(character(len=1) :: matdescra(6))
    allocate(character(len=1) :: matdescrb(6))

    matdescra=(/'G','X','X','F','X','X'/)
    matdescrb=(/'S','L','N','F','X','X'/)
    


    !!!!!! potentials module scalar dynamic allocation 
    allocate(character(len=5) :: dummyc)

    !!!!!! potentials_gw_casida module scalar dynamic allocation 
    allocate(ii10,ii11,ii30)
    allocate(MAXFCT,MNUM,MTYPE,MSGLVL,PHASE,idum,pardiso_info,MTYPE_real)

    allocate(j_imag,j_real,j_zero)
    j_imag=(0.0d0,1.0d0)
    j_real=(1.0d0,0.0d0)
    j_zero=(0.0d0,0.0d0)

    allocate(E_guess_min,E_guess_max)
    allocate(Nguess_complex,N_energy)
    allocate(alpha_gw)
    allocate(eta_omega1,eta_omega2,eta_time)



    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! call MPI_INIT(code)

    ! call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code )
    ! call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code )
    
    !!!!!!!!!!!!!!!!!!!!!!! READ INPUT # of procs in L3
    ! call GET_COMMAND_ARGUMENT(1,cnL3) !! number of L3 processors
    ! read(cnL3,'(I3)') nL3

    ! nL3=1

    

    

    

    



    

    !!!!! Can read up to 1 argument command line
    iarg = command_argument_count()
    ! if (iarg/=4) then
    !     print *,'only 4 argument possible'
    !     stop
    ! end if
    ! call get_command_argument(1,name,clen,stat) ! name of the file
    ! call get_command_argument(2,meanfield,clen,stat) ! DFT or Hartree Fock -- Options: dft, hf
    ! call get_command_argument(3,gw_sigma,clen,stat) ! gw self-energy frequency integral  -- Options: cd, casida
    ! call get_command_argument(4,dgr_fem,clen,stat) ! degree of polynomial of FEM  -- Options: p2, p3
    if (iarg/=1) then
        print *,'only 1 argument possible'
        stop
    end if
    call get_command_argument(1,name,clen,stat) ! name of the file



    call load_gw(name,dgr_fem,meanfield,gw_sigma)




    !!!! Read atoms from xyz file
    ! call load_xyz(name,1.0d0,nbat,at)
    call load_xyz(name,1.0d0)


    if (dgr_fem=='p2') then
        Nlocal = 10
        degree = 2
    else if (dgr_fem=='p3') then
        Nlocal = 20
        degree = 3
    end if
    Nquadrature=11
    



    


    ! then

        print *," -------------------------------------------"
        print *,"|                 set Mesh                 |" !!!!! node number is ZERO-base for TETGEN, but ONE-base for DFT code
        print *," -------------------------------------------"
	
	! call get_environment_variable('DFTROOT',rootdir) 
	
	! print *,trim(rootdir)



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! CREATE .poly MESH (easy way-- Could be periodic) !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! legacy option
    allocate(dedge(2,3)) !distance to edges - dedge(1,:) is left(negative); dedge(2,:) is right(positive)
    offset=offset0
    dedge=offset

    ! 8 points of the real box (computational domain)
    allocate(rbox(8,3))

    rbox(1,1)=minval(at(:)%c(1))-dedge(1,1)
    rbox(1,2)=minval(at(:)%c(2))-dedge(1,2)
    rbox(1,3)=minval(at(:)%c(3))-dedge(1,3)


    rbox(2,1)=minval(at(:)%c(1))-dedge(1,1)
    rbox(2,2)=maxval(at(:)%c(2))+dedge(2,2)
    rbox(2,3)=minval(at(:)%c(3))-dedge(1,3)

    rbox(3,1)=minval(at(:)%c(1))-dedge(1,1)
    rbox(3,2)=maxval(at(:)%c(2))+dedge(2,2)
    rbox(3,3)=maxval(at(:)%c(3))+dedge(2,3)

    rbox(4,1)=minval(at(:)%c(1))-dedge(1,1)
    rbox(4,2)=minval(at(:)%c(2))-dedge(1,2)
    rbox(4,3)=maxval(at(:)%c(3))+dedge(2,3)

    rbox(5,1)=maxval(at(:)%c(1))+dedge(2,1)
    rbox(5,2)=minval(at(:)%c(2))-dedge(1,2)
    rbox(5,3)=minval(at(:)%c(3))-dedge(1,3)

    rbox(6,1)=maxval(at(:)%c(1))+dedge(2,1)
    rbox(6,2)=maxval(at(:)%c(2))+dedge(2,2)
    rbox(6,3)=minval(at(:)%c(3))-dedge(1,3)

    rbox(7,1)=maxval(at(:)%c(1))+dedge(2,1)
    rbox(7,2)=maxval(at(:)%c(2))+dedge(2,2)
    rbox(7,3)=maxval(at(:)%c(3))+dedge(2,3)

    rbox(8,1)=maxval(at(:)%c(1))+dedge(2,1)
    rbox(8,2)=minval(at(:)%c(2))-dedge(1,2)
    rbox(8,3)=maxval(at(:)%c(3))+dedge(2,3)



    Nmuffin_total=0
    do k=1,nbat
       i=at(k)%core
       write(charI,"(I5)"), i
       open(10,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3n',status='old')
       read(10,*) Nmuffin
       close(10)
       Nmuffin_total=Nmuffin_total+Nmuffin
    enddo

    Nface_total=0
    do k=1,nbat
        i=at(k)%core
        write(charI,"(I5)"), i
        open(10,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3f',status='old')
        read(10,*) Nface
        close(10)
        Nface_total=Nface_total+Nface
    enddo

    ! Nmuffin_total=0
    ! do k=1,nbat
    !    i=at(k)%core
    !    write(charI,"(I5)"), i
    !    open(10,file='../../database_muffins/He.1.node',status='old')
    !    read(10,*) Nmuffin
    !    close(10)
    !    Nmuffin_total=Nmuffin_total+Nmuffin
    ! enddo

    ! Nface_total=0
    ! do k=1,nbat
    !     i=at(k)%core
    !     write(charI,"(I5)"), i
    !     open(10,file='../../database_muffins/He.1.face',status='old')
    !     read(10,*) Nface
    !     close(10)
    !     Nface_total=Nface_total+Nface
    ! enddo



    open(10,file=trim(name)//'.poly',status='replace')
    !!!!!!!!!!!!!!!!!!!!
    !!!!!!!! NODES
    !!!!!!!!!!!!!!!!!!!!
     write(10,'(A)') "# Part 1 - node list"
     write(10,'(A)') "# node count,[dim], [attribute], [boundary marker]"

     

    !  write(10,*) 8+Nmuffin*nbat,3,0,1
     write(10,*) 8+Nmuffin_total,3,0,1
     !! box 3d
     20 FORMAT(I0,X,3(XF10.6),2X,I0)
     do i=1,8
        write(10,20) i-1,dble(rbox(i,:)),1
     enddo

     Nmuffin_total=0
     !! atoms muffin surface
    do k=1,nbat
        ! radius is average between calculated- empirical atomic radius
        ! radius muffin is 3/4 of that
        r=(PTABLE(at(k)%core)%rcalc+PTABLE(at(k)%core)%remp)/2.0d0 
        ! if (p==1) r=PTABLE(p)%remp !! special cae H atom (too close from other atom.. ex SiH4 overlap)
        ro=r*(3.0d0/4.0d0)

        ! print *,ro

        if (k>1) then
            i=at(k-1)%core
            write(charI,"(I5)"), i
            open(12,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3n',status='old')
            read(12,*) Nmuffin  
            close(12)
            Nmuffin_total=Nmuffin_total+Nmuffin
        endif

        ! print *,Nmuffin

        i=at(k)%core
        write(charI,"(I5)"), i
        open(12,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3n',status='old')
        read(12,*) Nmuffin  
        allocate(point_muffin(1:Nmuffin,1:3))
        do ii=1,Nmuffin
        read(12,*) dummyc,point_muffin(ii,1),point_muffin(ii,2),point_muffin(ii,3)
        enddo
        close(12)

        ! if (k>1) then
        !     i=at(k-1)%core
        !     write(charI,"(I5)"), i
        !     open(12,file='../../database_muffins/He.1.node',status='old')
        !     read(12,*) Nmuffin  
        !     close(12)
        !     Nmuffin_total=Nmuffin_total+Nmuffin
        ! endif

        ! ! print *,Nmuffin

        ! i=at(k)%core
        ! write(charI,"(I5)"), i
        ! open(12,file='../../database_muffins/He.1.node',status='old')
        ! read(12,*) Nmuffin  
        ! allocate(point_muffin(1:Nmuffin,1:3))
        ! do ii=1,Nmuffin
        ! read(12,*) dummyc,point_muffin(ii,1),point_muffin(ii,2),point_muffin(ii,3)
        ! enddo
        ! close(12)

        ! print *,Nmuffin

        ro=1.0d0!0.85d0!1.0d0

        do i=1,Nmuffin
            !!! according to atom number place different radius of muffin
            write(10,20) 8+Nmuffin_total+i-1,point_muffin(i,:)*ro+at(k)%c(:),0 
            ! write(10,20) 8+Nmuffin_total+i-1,point_muffin(i,:)*1.0d0+at(k)%c(:),0 
            ! write(10,20) 8+Nmuffin_total+i-1,point_muffin(i,:)*0.6d0+at(k)%c(:),0 
        end do

        deallocate(point_muffin)

    end do


    
    

    ! stop



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !!!!!!!!!!!!!! FACES 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(10,'(A)') ""
    write(10,'(A)') "# Part 2 - facet list"
    write(10,'(A)') "# facet count, [boundary marker]"

    

    ! write(10,*) 6+Nface*Nbat,1
     write(10,*) 6+Nface_total,1
    
    write(10,'(A)') "#list all facets one by one<<<<<<<<"
    write(10,'(A)') "#First line: <# of polygons> [# of holes] [boundary marker]"
    write(10,'(A)') "#Following lines list # of polygons:"
    write(10,'(A)') "#<# of corners> <corner 1> <corner 2> ... <corner #>"
    write(10,'(A)') "#..."
    write(10,'(A)') "#Following lines list # of holes (if any):"
    write(10,'(A)') "#<hole #> <x> <y> <z>"

    !! Domain box 6 faces
    write(10,*) 1,0,1
    write(10,*) 4,0,1,5,4 !1,2,6,5 
    write(10,*) 1,0,1
    write(10,*) 4,0,1,2,3 !1,2,3,4 
    write(10,*) 1,0,1
    write(10,*) 4,4,5,6,7 !5,6,7,8 
    write(10,*) 1,0,1
    write(10,*) 4,2,3,7,6 !3,4,8,7 
    write(10,*) 1,0,1
    write(10,*) 4,0,3,7,4 !1,4,8,5 
    write(10,*) 1,0,1
    write(10,*) 4,1,2,6,5 !2,3,7,6 

    ! do i=1,Nbat
    !     do ii=1,Nface
    !         write(10,*) 1, 0, 0
    !         write(10,*) 3,face_muffin(ii,1)+8+(i-1)*Nmuffin,face_muffin(ii,2)+8+(i-1)*Nmuffin,face_muffin(ii,3)+8+(i-1)*Nmuffin
    !     enddo
    ! enddo

    Nmuffin_total=0
    do i=1,Nbat

        if (i>1) then
            l=at(i-1)%core
            write(charI,"(I5)"), l
            open(12,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3n',status='old')
            read(12,*) Nmuffin
            close(12)
            Nmuffin_total=Nmuffin_total+Nmuffin
        endif

        ! print *,Nmuffin_total

        k=at(i)%core
        write(charI,"(I5)"), k
        open(12,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3f',status='old')
        read(12,*) Nface
        allocate(face_muffin(1:Nface,1:3))
        do ii=1,Nface
        read(12,*) dummyc,face_muffin(ii,1),face_muffin(ii,2),face_muffin(ii,3)

        ! if (i>1) then
        !     l=at(i-1)%core
        !     write(charI,"(I5)"), l
        !     open(12,file='../../database_muffins/He.1.node',status='old')
        !     read(12,*) Nmuffin
        !     close(12)
        !     Nmuffin_total=Nmuffin_total+Nmuffin
        ! endif

        ! ! print *,Nmuffin_total

        ! k=at(i)%core
        ! write(charI,"(I5)"), k
        ! open(12,file='../../database_muffins/He.1.face',status='old')
        ! read(12,*) Nface
        ! allocate(face_muffin(1:Nface,1:3))
        ! do ii=1,Nface
        ! read(12,*) dummyc,face_muffin(ii,1),face_muffin(ii,2),face_muffin(ii,3)

        write(10,*) 1, 0, 0
        write(10,*) 3,face_muffin(ii,1)+8+Nmuffin_total,face_muffin(ii,2)+8+Nmuffin_total,face_muffin(ii,3)+8+Nmuffin_total

        ! if (ii==160) print *,face_muffin(ii,1)+8+Nmuffin_total,face_muffin(ii,2)+8+Nmuffin_total,face_muffin(ii,3)+8+Nmuffin_total
        
        enddo
        close(12)

        deallocate(face_muffin)

    enddo


    ! stop



    !!!!!!!!!!!!!!!!!!!!
    !!!! Holes !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!
    write(10,'(A)') 
    write(10,'(A)') "# Part 3 - onle hole count"
    write(10,'(A)') "#Following lines list # of holes:"
    write(10,'(A)') "#<hole #> <x> <y> <z>"
    ! write(10,*) nbat
    ! do i=1,nbat
    ! write(10,'(I0,X,3(XF10.6))') i, at(i)%c(:)
    ! enddo
    write(10,*) 0
    

    !!!!!!!!!!!!!!!!!!!!
    !!!! Region !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!
    write(10,'(A)') 
    write(10,'(A)') "# Part 4 - one line region count"
    write(10,'(A)') "#Following lines list # of region attributes:"
    write(10,'(A)') "#<region #> <x> <y> <z> <region number> <region attribute>"
    write(10,*) 0

    close(10)






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!  tetgen  !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! q refine  - h for help
    !! V for statistics !! k is vtk 
    line='tetgen -pq'//trim(qtetgen)//'efkV '//trim(name)//'.poly'
    ! line='tetgen -pq'//trim(qtetgen)//'ef '//trim(name)//'.poly'  
    print *,line
    call execute_command_line(line)

    ! stop

    ! call execute_command_line('mv '//trim(name)//".1.node"//" "//trim(name)//".node")
    ! call execute_command_line('mv '//trim(name)//".1.ele" //" "//trim(name)//".ele")
    ! call execute_command_line('mv '//trim(name)//".1.face"//" "//trim(name)//".face")
    ! call execute_command_line('mv '//trim(name)//".1.edge"//" "//trim(name)//".edge")
    ! call execute_command_line('mv '//trim(name)//".1.vtk"//" "//trim(name)//".vtk")

    !!!!!!!!!!!!! copy xyz file to keep track of the input file for the generated mesh
    ! call execute_command_line('cp '//trim(name)//".xyz"//" "//trim(name)//"_itmesh.xyz")


    ! end if ! rank0

    ! call MPI_BARRIER(MPI_COMM_WORLD ,code)


    open(10,file=trim(name)//'.1.edge',status='old')
    read (10,*) Ng0
    allocate(edge0(1:Ng0,1:2))
    do i=1,Ng0
        read(10,*) dummy0,edge1,edge2
        edge0(i,1)=edge1+1
        edge0(i,2)=edge2+1
    enddo
    close(10)


    open(10,file=trim(name)//'.1.node',status='old')
    read(10,*) Nn0
    allocate(point0(1:Nn0,1:3))
    ! allocate(color(1:Nn0))
    do i=1,Nn0
        read(10,*) dummy0, point0(i,1), point0(i,2), point0(i,3)!, color(i)
    enddo
    close(10)


    open(10,file=trim(name)//'.1.ele',status='old')
    read(10,*) Ne
    allocate(ele(1:Ne,1:Nlocal))
    do i=1,Ne
    read(10,*) dummy0,g01,g02,g03,g04
    ele(i,1)=g01+1
    ele(i,2)=g02+1
    ele(i,3)=g03+1
    ele(i,4)=g04+1
    end do
    close(10)


    open(10,file=trim(name)//'.1.face',status='old')
    read (10,*) Nf0
    allocate(face0(1:Nf0,1:3))
    do i=1,Nf0
        read(10,*) dummy0,face1,face2,face3
        face0(i,1)=face1+1
        face0(i,2)=face2+1
        face0(i,3)=face3+1
    enddo
    close(10)



    x_min = minval(point0(:,1))
    x_max = maxval(point0(:,1))
    y_min = minval(point0(:,2))
    y_max = maxval(point0(:,2))
    z_min = minval(point0(:,3))
    z_max = maxval(point0(:,3))

    allocate(color0(1:Nn0))
    color0=0


    allocate(BC(0))
    Nbc=0
    do i=1,Nn0
        if ((abs(point0(i,1)-x_min)<1.0d-8).or.(abs(point0(i,1)-x_max)<1.0d-8).or.(abs(point0(i,2)-y_min)<1.0d-8).or.&
            (abs(point0(i,2)-y_max)<1.0d-8).or.(abs(point0(i,3)-z_min)<1.0d-8).or.(abs(point0(i,3)-z_max)<1.0d-8)) then
            color0(i)=111
            BC=[BC,(/i/)]
            Nbc=Nbc+1
        end if
    enddo

    ! Nbc=0
    ! do i=1,Nn0
    !     if (color(i)==111) then
    !         color0(i)=111
    !         BC=[BC,(/i/)]
    !         Nbc=Nbc+1
    !     end if
    ! enddo


    ! do jj=1,2

    !     if (jj==1) then

    !         Nlocal=10
    !         degree=2
    !         dgr_fem='p2'
    !     else
    !         Nlocal=20
    !         degree=3
    !         dgr_fem='p3'
    !     end if


    if (dgr_fem=='p2') then

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!! p2 !!!!!!!!!!!!!!!!!!!! !!! Could be accelerated by MPI or CUDA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    allocate(pointe0(1:Ng0,1:3))
    allocate(colore0(1:Ng0))
    ! allocate(colore_face(1:Ng))


    colore0(:)=0
    ! colore_face0(:)=0

    do i=1,Ng0
        pointe0(i,1)=(point0(edge0(i,1),1)+point0(edge0(i,2),1))/2.0d0
        pointe0(i,2)=(point0(edge0(i,1),2)+point0(edge0(i,2),2))/2.0d0
        pointe0(i,3)=(point0(edge0(i,1),3)+point0(edge0(i,2),3))/2.0d0

        if ((color0(edge0(i,1))==111).and.(color0(edge0(i,2))==111)) then
            if ((abs(pointe0(i,1)-x_min)<1.0d-8).or.(abs(pointe0(i,1)-x_max)<1.0d-8).or.(abs(pointe0(i,2)-y_min)<1.0d-8).or.&
                (abs(pointe0(i,2)-y_max)<1.0d-8).or.(abs(pointe0(i,3)-z_min)<1.0d-8).or.(abs(pointe0(i,3)-z_max)<1.0d-8)) then
                    colore0(i)=111
                    BC=[BC,(/i+Nn0/)]
                    Nbc=Nbc+1
            end if
        end if
    
    enddo


    Nn=Nn0+Ng0
    allocate(point(1:Nn,1:3))
    allocate(color(1:Nn))
    ! allocate(color0_face(1:Nn))
    do i=1,Nn0
        point(i,1)=point0(i,1)
        point(i,2)=point0(i,2)
        point(i,3)=point0(i,3)
        color(i)=color0(i)
        !color0_face(i)=color_face(i)
    enddo
    do i=1,Ng0
        point(i+Nn0,1)=pointe0(i,1)
        point(i+Nn0,2)=pointe0(i,2)
        point(i+Nn0,3)=pointe0(i,3)
        color(i+Nn0)=colore0(i)
        !color0_face(i+Nn)=colore_face(i)
    enddo


    do i=1,Ne
        do k=1,Ng0
            if (comp((/ele(i,1),ele(i,2)/),edge0(k,:))==2) then
                ele(i,5)=Nn0+k
            else if (comp((/ele(i,2),ele(i,3)/),edge0(k,:))==2) then
                ele(i,6)=Nn0+k
            else if (comp((/ele(i,1),ele(i,3)/),edge0(k,:))==2) then
                ele(i,7)=Nn0+k
            else if (comp((/ele(i,1),ele(i,4)/),edge0(k,:))==2) then
                ele(i,8)=Nn0+k
            else if (comp((/ele(i,2),ele(i,4)/),edge0(k,:))==2) then
                ele(i,9)=Nn0+k
            else if (comp((/ele(i,3),ele(i,4)/),edge0(k,:))==2) then
                ele(i,10)=Nn0+k
            endif
        enddo
    enddo


    else if (dgr_fem=='p3') then

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!! p3 !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocate(pointe(1:(Ng*2+Nf),1:3))
    allocate(pointe0(1:(Ng0*2+Nf0),1:3))
    allocate(colore0(1:(Ng0*2+Nf0)))
    ! allocate(colore_face(1:Ng))

    do i=0,Ng0*2-1,2
        pointe0(i+1,1)=point0(edge0(i/2+1,1),1)+1.0d0/3.0d0*(point0(edge0(i/2+1,2),1)-point0(edge0(i/2+1,1),1))
        pointe0(i+1,2)=point0(edge0(i/2+1,1),2)+1.0d0/3.0d0*(point0(edge0(i/2+1,2),2)-point0(edge0(i/2+1,1),2))
        pointe0(i+1,3)=point0(edge0(i/2+1,1),3)+1.0d0/3.0d0*(point0(edge0(i/2+1,2),3)-point0(edge0(i/2+1,1),3))
        pointe0(i+2,1)=point0(edge0(i/2+1,1),1)+2.0d0/3.0d0*(point0(edge0(i/2+1,2),1)-point0(edge0(i/2+1,1),1))
        pointe0(i+2,2)=point0(edge0(i/2+1,1),2)+2.0d0/3.0d0*(point0(edge0(i/2+1,2),2)-point0(edge0(i/2+1,1),2))
        pointe0(i+2,3)=point0(edge0(i/2+1,1),3)+2.0d0/3.0d0*(point0(edge0(i/2+1,2),3)-point0(edge0(i/2+1,1),3))

        if ((color0(edge0(i/2+1,1))==111).and.(color0(edge0(i/2+1,2))==111)) then
            if ((abs(pointe0(i+1,1)-x_min)<1.0d-8).or.(abs(pointe0(i+1,1)-x_max)<1.0d-8).or.(abs(pointe0(i+1,2)-y_min)<1.0d-8).or.&
                (abs(pointe0(i+1,2)-y_max)<1.0d-8).or.(abs(pointe0(i+1,3)-z_min)<1.0d-8).or.(abs(pointe0(i+1,3)-z_max)<1.0d-8)) then
                    colore0(i+1)=111
                    BC=[BC,(/i+1+Nn0/)]
                    Nbc=Nbc+1
            end if
            if ((abs(pointe0(i+2,1)-x_min)<1.0d-8).or.(abs(pointe0(i+2,1)-x_max)<1.0d-8).or.(abs(pointe0(i+2,2)-y_min)<1.0d-8).or.&
                (abs(pointe0(i+2,2)-y_max)<1.0d-8).or.(abs(pointe0(i+2,3)-z_min)<1.0d-8).or.(abs(pointe0(i+2,3)-z_max)<1.0d-8)) then
                    colore0(i+2)=111
                    BC=[BC,(/i+2+Nn0/)]
                    Nbc=Nbc+1
            end if
        end if
    enddo

    do i=1,Nf0
        pointe0(2*Ng0+i,1)=1.0d0/3.0d0*(point0(face0(i,1),1)+point0(face0(i,2),1)+point0(face0(i,3),1))
        pointe0(2*Ng0+i,2)=1.0d0/3.0d0*(point0(face0(i,1),2)+point0(face0(i,2),2)+point0(face0(i,3),2))
        pointe0(2*Ng0+i,3)=1.0d0/3.0d0*(point0(face0(i,1),3)+point0(face0(i,2),3)+point0(face0(i,3),3))

        if ((color0(face0(i,1))==111).and.(color0(face0(i,2))==111).and.(color0(face0(i,3))==111)) then
            if ((abs(pointe0(i,1)-x_min)<1.0d-8).or.(abs(pointe0(i,1)-x_max)<1.0d-8).or.(abs(pointe0(i,2)-y_min)<1.0d-8).or.&
                (abs(pointe0(i,2)-y_max)<1.0d-8).or.(abs(pointe0(i,3)-z_min)<1.0d-8).or.(abs(pointe0(i,3)-z_max)<1.0d-8)) then

                colore0(i+Ng0*2)=111
                BC=[BC,(/i+Nn0+Ng0*2/)]
                Nbc=Nbc+1

            end if
        end if

    enddo

    Nn=Nn0+Ng0*2+Nf0
    allocate(point(1:Nn,1:3))
    allocate(color(1:Nn))
    ! allocate(color0_face(1:Nn))
    do i=1,Nn0
        point(i,1)=point0(i,1)
        point(i,2)=point0(i,2)
        point(i,3)=point0(i,3)
        color(i)=color0(i)
        !color0_face(i)=color_face(i)
    enddo
    do i=1,Ng0*2+Nf0
        point(i+Nn0,1)=pointe0(i,1)
        point(i+Nn0,2)=pointe0(i,2)
        point(i+Nn0,3)=pointe0(i,3)
        color(i+Nn0)=colore0(i)
        !color0_face(i+Nn)=colore_face(i)
    enddo


    allocate(pointex(1:Ne,1:12))
    allocate(pointey(1:Ne,1:12))
    allocate(pointez(1:Ne,1:12))

    pointex(:,1)=point0(ele(:,1),1)+1.0d0/3.0d0*(point0(ele(:,2),1)-point0(ele(:,1),1))
    pointey(:,1)=point0(ele(:,1),2)+1.0d0/3.0d0*(point0(ele(:,2),2)-point0(ele(:,1),2))
    pointez(:,1)=point0(ele(:,1),3)+1.0d0/3.0d0*(point0(ele(:,2),3)-point0(ele(:,1),3))
    pointex(:,2)=point0(ele(:,1),1)+2.0d0/3.0d0*(point0(ele(:,2),1)-point0(ele(:,1),1))
    pointey(:,2)=point0(ele(:,1),2)+2.0d0/3.0d0*(point0(ele(:,2),2)-point0(ele(:,1),2))
    pointez(:,2)=point0(ele(:,1),3)+2.0d0/3.0d0*(point0(ele(:,2),3)-point0(ele(:,1),3))
    pointex(:,3)=point0(ele(:,2),1)+1.0d0/3.0d0*(point0(ele(:,3),1)-point0(ele(:,2),1))
    pointey(:,3)=point0(ele(:,2),2)+1.0d0/3.0d0*(point0(ele(:,3),2)-point0(ele(:,2),2))
    pointez(:,3)=point0(ele(:,2),3)+1.0d0/3.0d0*(point0(ele(:,3),3)-point0(ele(:,2),3))
    pointex(:,4)=point0(ele(:,2),1)+2.0d0/3.0d0*(point0(ele(:,3),1)-point0(ele(:,2),1))
    pointey(:,4)=point0(ele(:,2),2)+2.0d0/3.0d0*(point0(ele(:,3),2)-point0(ele(:,2),2))
    pointez(:,4)=point0(ele(:,2),3)+2.0d0/3.0d0*(point0(ele(:,3),3)-point0(ele(:,2),3))
    pointex(:,5)=point0(ele(:,3),1)+1.0d0/3.0d0*(point0(ele(:,1),1)-point0(ele(:,3),1))
    pointey(:,5)=point0(ele(:,3),2)+1.0d0/3.0d0*(point0(ele(:,1),2)-point0(ele(:,3),2))
    pointez(:,5)=point0(ele(:,3),3)+1.0d0/3.0d0*(point0(ele(:,1),3)-point0(ele(:,3),3))
    pointex(:,6)=point0(ele(:,3),1)+2.0d0/3.0d0*(point0(ele(:,1),1)-point0(ele(:,3),1))
    pointey(:,6)=point0(ele(:,3),2)+2.0d0/3.0d0*(point0(ele(:,1),2)-point0(ele(:,3),2))
    pointez(:,6)=point0(ele(:,3),3)+2.0d0/3.0d0*(point0(ele(:,1),3)-point0(ele(:,3),3))
    pointex(:,7)=point0(ele(:,1),1)+1.0d0/3.0d0*(point0(ele(:,4),1)-point0(ele(:,1),1))
    pointey(:,7)=point0(ele(:,1),2)+1.0d0/3.0d0*(point0(ele(:,4),2)-point0(ele(:,1),2))
    pointez(:,7)=point0(ele(:,1),3)+1.0d0/3.0d0*(point0(ele(:,4),3)-point0(ele(:,1),3))
    pointex(:,8)=point0(ele(:,1),1)+2.0d0/3.0d0*(point0(ele(:,4),1)-point0(ele(:,1),1))
    pointey(:,8)=point0(ele(:,1),2)+2.0d0/3.0d0*(point0(ele(:,4),2)-point0(ele(:,1),2))
    pointez(:,8)=point0(ele(:,1),3)+2.0d0/3.0d0*(point0(ele(:,4),3)-point0(ele(:,1),3))
    pointex(:,9)=point0(ele(:,2),1)+1.0d0/3.0d0*(point0(ele(:,4),1)-point0(ele(:,2),1))
    pointey(:,9)=point0(ele(:,2),2)+1.0d0/3.0d0*(point0(ele(:,4),2)-point0(ele(:,2),2))
    pointez(:,9)=point0(ele(:,2),3)+1.0d0/3.0d0*(point0(ele(:,4),3)-point0(ele(:,2),3))

    pointex(:,10)=point0(ele(:,2),1)+2.0d0/3.0d0*(point0(ele(:,4),1)-point0(ele(:,2),1))
    pointey(:,10)=point0(ele(:,2),2)+2.0d0/3.0d0*(point0(ele(:,4),2)-point0(ele(:,2),2))
    pointez(:,10)=point0(ele(:,2),3)+2.0d0/3.0d0*(point0(ele(:,4),3)-point0(ele(:,2),3))
    pointex(:,11)=point0(ele(:,3),1)+1.0d0/3.0d0*(point0(ele(:,4),1)-point0(ele(:,3),1))
    pointey(:,11)=point0(ele(:,3),2)+1.0d0/3.0d0*(point0(ele(:,4),2)-point0(ele(:,3),2))
    pointez(:,11)=point0(ele(:,3),3)+1.0d0/3.0d0*(point0(ele(:,4),3)-point0(ele(:,3),3))
    pointex(:,12)=point0(ele(:,3),1)+2.0d0/3.0d0*(point0(ele(:,4),1)-point0(ele(:,3),1))
    pointey(:,12)=point0(ele(:,3),2)+2.0d0/3.0d0*(point0(ele(:,4),2)-point0(ele(:,3),2))
    pointez(:,12)=point0(ele(:,3),3)+2.0d0/3.0d0*(point0(ele(:,4),3)-point0(ele(:,3),3))



    do i=1,Ne
        do l=1,12
            do k=1,2*Ng0
                if ((abs(pointex(i,l)-point(k+Nn0,1))<1.0d-5).and.(abs(pointey(i,l)-point(k+Nn0,2))<1.0d-5)&
                                        .and.(abs(pointez(i,l)-point(k+Nn0,3))<1.0d-5)) then
                    ele(i,l+4)=k+Nn0
                endif
            enddo
        enddo
    enddo





    allocate(pointfx(1:Ne,1:4))
    allocate(pointfy(1:Ne,1:4))
    allocate(pointfz(1:Ne,1:4))


    pointfx(:,1)=1.0d0/3.0d0*(point0(ele(:,2),1)+point0(ele(:,3),1)+point0(ele(:,4),1))
    pointfy(:,1)=1.0d0/3.0d0*(point0(ele(:,2),2)+point0(ele(:,3),2)+point0(ele(:,4),2))
    pointfz(:,1)=1.0d0/3.0d0*(point0(ele(:,2),3)+point0(ele(:,3),3)+point0(ele(:,4),3))
    pointfx(:,2)=1.0d0/3.0d0*(point0(ele(:,1),1)+point0(ele(:,3),1)+point0(ele(:,4),1))
    pointfy(:,2)=1.0d0/3.0d0*(point0(ele(:,1),2)+point0(ele(:,3),2)+point0(ele(:,4),2))
    pointfz(:,2)=1.0d0/3.0d0*(point0(ele(:,1),3)+point0(ele(:,3),3)+point0(ele(:,4),3))
    pointfx(:,3)=1.0d0/3.0d0*(point0(ele(:,1),1)+point0(ele(:,2),1)+point0(ele(:,4),1))
    pointfy(:,3)=1.0d0/3.0d0*(point0(ele(:,1),2)+point0(ele(:,2),2)+point0(ele(:,4),2))
    pointfz(:,3)=1.0d0/3.0d0*(point0(ele(:,1),3)+point0(ele(:,2),3)+point0(ele(:,4),3))
    pointfx(:,4)=1.0d0/3.0d0*(point0(ele(:,1),1)+point0(ele(:,2),1)+point0(ele(:,3),1))
    pointfy(:,4)=1.0d0/3.0d0*(point0(ele(:,1),2)+point0(ele(:,2),2)+point0(ele(:,3),2))
    pointfz(:,4)=1.0d0/3.0d0*(point0(ele(:,1),3)+point0(ele(:,2),3)+point0(ele(:,3),3))


    ! allocate(node(1:Ne,1:4))

    do i=1,Ne
        do l=1,4
            do k=1,Nf0
                if ((abs(pointfx(i,l)-point(k+Nn0+Ng0*2,1))<1.0d-5).and.(abs(pointfy(i,l)-point(k+Nn0+Ng0*2,2))<1.0d-5)&
                                        .and.(abs(pointfz(i,l)-point(k+Nn0+Ng0*2,3))<1.0d-5)) then
                    ele(i,l+16)=k+Nn0+Ng0*2
                endif
            enddo
        enddo
    enddo

    deallocate(pointex)
    deallocate(pointey)
    deallocate(pointez)
    deallocate(pointfx)
    deallocate(pointfy)
    deallocate(pointfz)

    
    end if !!! end p2/p3 selection



    ! print *,Nn,Nn0,Ng0,Nf0
    ! print *,Nbc


    point=point/bohr

    do i=1,Nbat
        at(i)%c(:)=at(i)%c(:)/bohr
    enddo

    ! print *,point(9,:)
    ! print *,point(10,:)

    print *,"Nn -", Nn,"Ne -",Ne


    Nstates=0
    do i=1,Nbat
        Nstates=Nstates+at(i)%core
    enddo
    Nstates=Nstates/2

    print *,'# of occupied orbitals',Nstates


    deallocate(dedge)
    deallocate(rbox)
    deallocate(edge0)
    deallocate(point0)
    deallocate(face0)
    deallocate(pointe0)
    deallocate(colore0)
    

    ! allocate(point(1:Nn,1:3))
    ! allocate(ele(1:Ne,1:(4+6)))
    ! allocate(color(1:Nn))



    

    ! stop







    ! open(10,file=trim(name)//'_p2.node',status='replace')
    ! write(10,fmt='(4i7)') Nn,3,0,1
    ! do i=1,Nn
    !     write(10,fmt='(i7,3f12.5,i7)') i,real(point(i,:)*bohr),color(i)
    ! enddo
    ! close(10)

    ! 11 format (11i6)

    ! open(10,file=trim(name)//'_p2.ele',status='replace')
    ! write(10,fmt='(1i7)') Ne
    ! do i=1,Ne
    !     write(10,11) i,ele(i,1),ele(i,2),ele(i,3),ele(i,4),ele(i,5),ele(i,6),ele(i,7),ele(i,8),ele(i,9),ele(i,10)!&
    !     !,ele(i,11),ele(i,12),ele(i,13),ele(i,14),ele(i,15),ele(i,16),ele(i,17),ele(i,18),ele(i,19),ele(i,20)
    ! enddo
    ! close(10)

    
    ! stop


    ! stop


    ! ! open(10,file='Ne.1_p2.pxyz',status='replace')
    ! open(10,file='He_p2_extendedBC.pxyz',status='replace')

    ! do i=1,Nn
    !     write(10,*) point(i,:)*bohr
    ! enddo
    ! close(10)


    ! stop





















    ! ! open(10,file='CO_atommesh.1_p2.node',status='old')
    ! ! open(10,file='Ne.1_p2.node',status='old')
    ! ! open(10,file='He_p2.node',status='old')
    ! ! open(10,file='He.1_p2.node',status='old')
    ! open(10,file='He.1_p2 (copy).node',status='old')
    ! ! open(10,file='at_6.1_p2.node',status='old')
    ! read(10,*) Nn
    
    ! allocate(point(1:Nn,1:3))
    ! allocate(color(1:Nn))

    ! allocate(BC(0))
    ! Nbc=0

    ! do i=1,Nn
    !     read(10,*) dummy,point(i,1),point(i,2),point(i,3),color(i)!,color_new(i),color_new1(i)
    !     if (color(i)==111) then
    !         BC=[BC,(/i/)]
    !         Nbc=Nbc+1
    !     end if
    ! enddo
    ! close(10)

    

    ! point=12.0d0*point ! scaling
    ! ! point=point/0.529177249d0!*1.5446d0 ! convert to atomic units


    ! ! do i=1,Nn
    ! !     if (color(i)==1) color(i)=111
    ! ! enddo





    ! ! open(10,file='CO_atommesh.1_p2.ele',status='old')
    ! ! open(10,file='Ne.1_p2.ele',status='old')
    ! ! open(10,file='He_p2.ele',status='old')
    ! ! open(10,file='He.1_p2.ele',status='old')
    ! open(10,file='He.1_p2 (copy).ele',status='old')
    ! ! open(10,file='at_6.1_p2.ele',status='old')
    
    ! read(10,*) Ne
    ! allocate(ele(1:Ne,1:10))

    ! do i=1,Ne
    ! read(10,*) dummy,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10
    ! ele(i,1)=g1+1
    ! ele(i,2)=g2+1
    ! ele(i,3)=g3+1
    ! ele(i,4)=g4+1
    ! ele(i,5)=g5+1
    ! ele(i,6)=g6+1
    ! ele(i,7)=g7+1
    ! ele(i,8)=g8+1
    ! ele(i,9)=g9+1
    ! ele(i,10)=g10+1
    ! end do

    ! close(10)


    ! do i=1,Nbat
    !     at(i)%c(:)=at(i)%c(:)/bohr
    ! enddo


    ! print *,"Nn -", Nn,"Ne -",Ne


    ! Nstates=0
    ! do i=1,Nbat
    !     Nstates=Nstates+at(i)%core
    ! enddo
    ! Nstates=Nstates/2

    ! print *,'# of occupied orbitals',Nstates

    ! print *,Nbc





    

    

    

    if (meanfield=='dft') then

        print *," -------------------------------------------"
        print *,"| compute initial guess of electron density |"
        print *,"| & diagonalize Kohn-Sham equation          |"
        print *," -------------------------------------------"
    
    else if (meanfield=='hf') then

        print *," -------------------------------------------"
        print *,"| compute initial guess of electron density |"
        print *,"| & diagonalize Hartree-Fock equation       |"
        print *," -------------------------------------------"

    end if

    ! if (jj==1) then

    allocate(psi_0(1:Nn))
    ! allocate(psi_100(1:Nn))
    ! allocate(psi_200(1:Nn))
    ! allocate(psi_21x(1:Nn))
    ! allocate(psi_21y(1:Nn))
    ! allocate(psi_21z(1:Nn))
    ! allocate(psi_300(1:Nn))
    ! allocate(psi_31x(1:Nn))
    ! allocate(psi_31y(1:Nn))
    ! allocate(psi_31z(1:Nn))

    psi_0=0.0d0

    psi_square_100=0.0d0
    psi_square_200=0.0d0
    psi_square_21x=0.0d0
    psi_square_21y=0.0d0
    psi_square_21z=0.0d0
    psi_square_300=0.0d0
    psi_square_31x=0.0d0
    psi_square_31y=0.0d0
    psi_square_31z=0.0d0



    allocate(r_init(10000,Nbat),n_init(10000,Nbat))
    allocate(Nn_init(Nbat))
    allocate(io)
    allocate(lb,ub)

    
    r_init=0.0d0
    n_init=0.0d0



    do i=1,Nbat

        Nn_init(i)=0

    k=at(i)%core
    write(charI,"(I5)"), k
    open(12,file='../../database_muffins/at_'//trim(adjustl(charI))//'.lda_p2_rnv',status='old')
        
        ! do
        !     read(12,*,end=10) r_init(Nn_init+1),dummy_real,n_init(Nn_init+1)
        !     Nn_init=Nn_init+1
        ! enddo
        do 
            READ(12,*,iostat=io) r_init(Nn_init(i)+1,i),dummy_real,n_init(Nn_init(i)+1,i)
            IF (io/=0) EXIT
            Nn_init(i)=Nn_init(i)+1
        enddo

    close(12)

    enddo

    r_init=r_init/1.0d-10!/bohr
    n_init=n_init*1.0d-30!/bohr**3
        
    do i=1,Nn
        do ll=1,Nbat

            rxyz=distance(point(i,:),at(ll)%c(:))

            call find_interval(r_init(:,ll), Nn_init(ll), rxyz*bohr, lb, ub)

            ! if (i==9) then

            ! print *,lb,ub,rxyz*bohr,rxyz
            ! print *,point(i,:)
            ! print *,at(ll)%c(:)
            ! print *,(n_init(lb,ll)+n_init(ub,ll))/2.0d0

            ! stop

            ! end if

            

            psi_0(i)=psi_0(i)+((n_init(ub,ll)-n_init(lb,ll))/(r_init(ub,ll)-r_init(lb,ll))*(rxyz*bohr-r_init(lb,ll))&
            +n_init(lb,ll))*bohr**3

        enddo
    enddo





!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !!!!!!!!!  compute intitial guess of electron density  !!!!!!!!!!!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !!!! Could be optimized by producing database in advance and read from it. 

!     do i=1,Nn
!         do ll=1,Nbat
            
!             rxyz=distance(point(i,:),at(ll)%c(:))

!             Z0=dble(at(ll)%core)
        
!             if (rxyz>radius) then

!                 psi_square_100=0.0d0 
!                 psi_square_200=0.0d0 
!                 psi_square_21x=0.0d0 
!                 psi_square_21y=0.0d0 
!                 psi_square_21z=0.0d0
!                 psi_square_300=0.0d0 
!                 psi_square_31x=0.0d0 
!                 psi_square_31y=0.0d0 
!                 psi_square_31z=0.0d0 
!             else
                
!                 call psi_1s(at(ll)%c(:),point(i,:),Z0)
!                 call psi_2s(at(ll)%c(:),point(i,:),Z0) 
!                 call psi_2x(at(ll)%c(:),point(i,:),Z0) 
!                 call psi_2y(at(ll)%c(:),point(i,:),Z0) 
!                 call psi_2z(at(ll)%c(:),point(i,:),Z0) 
!                 call psi_3s(at(ll)%c(:),point(i,:),Z0)
!                 call psi_3x(at(ll)%c(:),point(i,:),Z0)
!                 call psi_3y(at(ll)%c(:),point(i,:),Z0)
!                 call psi_3z(at(ll)%c(:),point(i,:),Z0)    



!                 if (PTABLE(at(ll)%core)%filled_core==0) then

!                     ! psi_100(i)=psi_100(i)+dble(at(ll)%core)*psi_100_temp**2

!                     ! psi_0(i)=psi_100(i)

!                     psi_0(i)=psi_0(i)+dble(at(ll)%core)*psi_square_100**2
                
!                 else if (PTABLE(at(ll)%core)%filled_core==2) then
    
!                     ! psi_100(i)=psi_100(i)+2.0d0*psi_100_temp**2
!                     ! psi_200(i)=psi_200(i)+(dble(at(ll)%core)-2.0d0)/4.0d0*psi_200_temp**2
!                     ! psi_21x(i)=psi_21x(i)+(dble(at(ll)%core)-2.0d0)/4.0d0*psi_21x_temp**2
!                     ! psi_21y(i)=psi_21y(i)+(dble(at(ll)%core)-2.0d0)/4.0d0*psi_21y_temp**2
!                     ! psi_21Z(i)=psi_21Z(i)+(dble(at(ll)%core)-2.0d0)/4.0d0*psi_21z_temp**2
!                     ! ! psi_200(i)=psi_200(i)+2.0d0*psi_200_temp**2
!                     ! ! psi_21x(i)=psi_21x(i)+(dble(at(ll)%core)-4.0d0)/3.0d0*psi_21x_temp**2
!                     ! ! psi_21y(i)=psi_21y(i)+(dble(at(ll)%core)-4.0d0)/3.0d0*psi_21y_temp**2
!                     ! ! psi_21Z(i)=psi_21Z(i)+(dble(at(ll)%core)-4.0d0)/3.0d0*psi_21z_temp**2

!                     ! psi_0(i)=psi_100(i)+psi_200(i)+psi_21x(i)+psi_21y(i)+psi_21z(i)

!                     ! if (at(ll)%core<=4) then
!                     !     psi_0(i)=psi_0(i)+2.0d0*psi_square_100**2
!                     ! else
!                         psi_0(i)=psi_0(i)+2.0d0*psi_square_100**2&
!                         +(dble(at(ll)%core)-2.0d0)*psi_square_200**2
!                         ! +(dble(at(ll)%core)-2.0d0)/4.0d0*psi_square_200**2&
!                         ! +(dble(at(ll)%core)-2.0d0)/4.0d0*psi_square_21x**2&
!                         ! +(dble(at(ll)%core)-2.0d0)/4.0d0*psi_square_21y**2&
!                         ! +(dble(at(ll)%core)-2.0d0)/4.0d0*psi_square_21z**2
!                     ! end if
    
!                 else if (PTABLE(at(k)%core)%filled_core==10) then
    
!                     ! psi_100(i)=psi_100(i)+2.0d0*psi_100_temp**2
!                     ! psi_200(i)=psi_200(i)+2.0d0*psi_200_temp**2
!                     ! psi_21x(i)=psi_21x(i)+2.0d0*psi_21x_temp**2
!                     ! psi_21y(i)=psi_21y(i)+2.0d0*psi_21y_temp**2
!                     ! psi_21Z(i)=psi_21Z(i)+2.0d0*psi_21z_temp**2
!                     ! psi_300(i)=psi_300(i)+(dble(at(ll)%core)-10.0d0)/4.0d0*psi_300_temp**2
!                     ! psi_31x(i)=psi_31x(i)+(dble(at(ll)%core)-10.0d0)/4.0d0*psi_31x_temp**2
!                     ! psi_31y(i)=psi_31y(i)+(dble(at(ll)%core)-10.0d0)/4.0d0*psi_31y_temp**2
!                     ! psi_31Z(i)=psi_31Z(i)+(dble(at(ll)%core)-10.0d0)/4.0d0*psi_31z_temp**2

!                     ! psi_0(i)=psi_100(i)+psi_200(i)+psi_21x(i)+psi_21y(i)+psi_21z(i)+psi_300(i)+psi_31x(i)+psi_31y(i)+psi_31z(i)
                    
                    
!                     psi_0(i)=psi_0(i)+2.0d0*psi_square_100**2&
!                     +2.0d0*psi_square_200**2&
!                     +2.0d0*psi_square_21x**2&
!                     +2.0d0*psi_square_21y**2&
!                     +2.0d0*psi_square_21z**2&
!                     +(dble(at(ll)%core)-10.0d0)/4.0d0*psi_square_300**2&
!                     +(dble(at(ll)%core)-10.0d0)/4.0d0*psi_square_31x**2&
!                     +(dble(at(ll)%core)-10.0d0)/4.0d0*psi_square_31y**2&
!                     +(dble(at(ll)%core)-10.0d0)/4.0d0*psi_square_31z**2

    
!                 end if

!             end if

!         enddo

!         ! psi_0(i)=psi_100(i)+psi_200(i)+psi_21x(i)+psi_21y(i)+psi_21z(i)+psi_300(i)+psi_31x(i)+psi_31y(i)+psi_31z(i) 
!         psi_0(i)=sqrt(psi_0(i))

!     enddo

! ! else if (jj==2) then
! !     psi_0(:)=sqrt(nq(:))

! ! end if


!     ! allocate(psi_0(1:Nn))
!     ! deallocate(psi_100)
!     ! deallocate(psi_200)
!     ! deallocate(psi_21x)
!     ! deallocate(psi_21y)
!     ! deallocate(psi_21z)
!     ! deallocate(psi_300)
!     ! deallocate(psi_31x)
!     ! deallocate(psi_31y)
!     ! deallocate(psi_31z)







    ! allocate(psi_i(1:Nn,1:Nstates+1))
    ! do k=1,Nstates+1
    ! write (file_id, '(I0)') k
    ! ! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
    ! open (unit=20,file='Ne.1_p2_nessie.psi'//trim(file_id),status='old')
    ! do i=1,Nn
    !     read(20,*) psi_i(i,k)
    ! enddo
    ! close (20)
    ! enddo
    ! psi_i=psi_i*sqrt(a0*bohr**3)
    
    
    

    





    ! stop







    if (meanfield=="hf") then


        allocate(psi_i(1:Nn,1:Nstates+1))
        do k=1,Nstates+1
        write (file_id, '(I0)') k
        ! open (unit=20,file='CO.1_p2_nessie.psi'//trim(file_id),status='old')
        ! open (unit=20,file='CO_p2_dft.psi'//trim(file_id),status='old')
        ! open (unit=20,file='CO_p2_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='CO_0.6415.1_p2_nessie.psi'//trim(file_id),status='old')
        ! open (unit=20,file='He.1_p2_nessie.psi'//trim(file_id),status='old')
        ! open (unit=20,file='He_p2_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='He_p3_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='Ne.1_p2_nessie.psi'//trim(file_id),status='old')
        ! open (unit=20,file='Ne_p2_dft.psi'//trim(file_id),status='old')
        ! open (unit=20,file='LiH.1_p2_nessie.psi'//trim(file_id),status='old')
        open (unit=20,file='LiH_p2_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='Li2_p2_dft.psi'//trim(file_id),status='old')
        ! open (unit=20,file='H2.1_p2_nessie.psi'//trim(file_id),status='old')
        ! open (unit=20,file='H2_p2_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='F2_p2_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='P2_p2_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='NP_p2_dft_1.psi'//trim(file_id),status='old')
        ! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
        do i=1,Nn
            read(20,*) psi_i(i,k)
        enddo
        close (20)
        enddo
        ! psi_i=psi_i*sqrt(a0*bohr**3)

        ! do k=1,Nstates
        ! write (file_id, '(I0)') k
        ! ! open (unit=20,file='CO.1_p2_nessie.psi'//trim(file_id),status='old')
        ! ! open (unit=20,file='CO_0.6415.1_p2_nessie.psi'//trim(file_id),status='old')
        ! open (unit=20,file='He_p2_dft.psi'//trim(file_id),status='old')
        ! ! open (unit=20,file='LiH.1_p2_nessie.psi'//trim(file_id),status='old')
        ! ! open (unit=20,file='H2.1_p2_nessie.psi'//trim(file_id),status='old')
        ! ! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
        ! do i=1,Nn
        !     read(20,*) psi_i(i,k)
        ! enddo
        ! close (20)
        ! enddo


        ! psi_i = 0.0d0

        ! do i=1,Nn
        !     do ll=1,Nbat

        !         rxyz=distance(point(i,:),at(ll)%c(:))

        !         psi_i(i,1) = psi_i(i,1) + 1.0d0/sqrt(pi)*at(ll)%core**(3.0d0/2.0d0)*exp(-at(ll)%core*rxyz)

        !     enddo
        ! enddo



        ! do i=1,Nn
        !     psi_i(i,1)=psi_1s(at(1)%c(:),point(i,:),dble(at(1)%core))
        ! enddo




        print *,"+++++++++"

    end if





































!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!     construct Kinetic energy and mass matrix      !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(neigh_AB(Nn))
allocate(neigh_V(Nn))

call set_linkedlistAB(Ne,Nlocal,Nquadrature)

allocate(volumegweight(1:Ne*Nquadrature))
call construct_volumegweight(Ne,Nquadrature)

neigh_size = 0
do i=1,Nn
    neigh_size=neigh_size+neigh_AB(i)%size
enddo

allocate(A(neigh_size))
allocate(B(neigh_size))
! allocate(A00(neigh_size))
allocate(IA(Nn+1))
allocate(JA(neigh_size))
! allocate(IA0(Nn+1))
! allocate(JA0(neigh_size))
allocate(V(neigh_size))
allocate(V_i(neigh_size))


allocate(xy_AB(1:2))
allocate(xy_V(1:2))
ii=0
do i=1,Nn
    IA(i)=ii+1
    do l=1,neigh_AB(i)%size
        ii=ii+1
        JA(ii)=neigh_AB(i)%get(l)
        xy_AB=neigh_AB(i)%getab(l)
        A(ii)=dble(xy_AB(1))
        B(ii)=dble(xy_AB(2))
    enddo
enddo
IA(Nn+1)=neigh_size+1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!    quaternion     !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate(quaternion(4*Nn))
! allocate(quaternion_A(16*neigh_size))
! allocate(quaternion_IA(4*Nn+1))
! allocate(quaternion_JA(16*neigh_size))

! call build_quaternionmatrix(Nn)

! ii=0
! do i=1,4*Nn
!     quaternion_IA(i)=ii+1
!     do l=1,quaternion(i)%size
!         ii=ii+1
!         quaternion_JA(ii)=quaternion(i)%get(l)
!     enddo
! enddo

! ! print *,ii,neigh_size

! ! stop


! quaternion_IA(4*Nn+1)=16*neigh_size+1

! allocate(IdenMat_quaternion(4*Nn,Nn))

! do i=1,Nn
!     IdenMat_quaternion((i-1)*4+1,i)=1.0d0
!     ! IdenMat_quaternion((i-1)*4+3,i)=1.0d0!sqrt(1.0d0/3.0d0)
!     ! IdenMat_quaternion((i-1)*4+4,i)=sqrt(1.0d0/3.0d0)
! enddo

! allocate(G_quaternion(4*Nn,Nn))

! ! stop
    



allocate(IA_pntrb(Nn+1))
allocate(IA_pntre(Nn+1))
IA_pntrb=IA(1:Nn)
IA_pntre=IA(2:Nn+1)


allocate(H_dense(Nn,Nn))
allocate(S_dense(Nn,Nn))
allocate(H_dense_psi(Nn,Nn))
allocate(B_dense(Nn,Nn))
H_dense=0.0d0
S_dense=0.0d0

do i=1,Nn
    do l=1,neigh_AB(i)%size
        xy_AB=neigh_AB(i)%getab(l)
        ii=neigh_AB(i)%get(l)
        H_dense(i,ii)=H_dense(i,ii)+dble(xy_AB(1))
        H_dense(i,ii)=H_dense(i,ii)*0.5d0
        S_dense(i,ii)=S_dense(i,ii)+dble(xy_AB(2))
    enddo
enddo

! open(10,file='S_dense.txt',status='replace')
! do i=1,Nn
!     write(10,*) S_dense(i,:)
! enddo

! stop

! print *,B(1)

! stop

allocate(neigh_NnToNg(1:Ne*Nquadrature))
call construct_NnToNg(Ne,Nquadrature,Nlocal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!     FEAST eigensolver parameters setting      !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

M0=ks_M00!40 ! 2*Nstates+10!+40 ! 40 is an arbitrary number 
                ! it could be even optimized by solving each eigenvalue interval provided by FEAST, which is much faster than this.
Emin=ks_Emin0
Emax=ks_Emax0


!!! Allocate memory for eigenvalues. eigenvectors, residual 
allocate(E(M0), psi(Nn,M0), res(M0))

allocate(E_hf(Nn))











! allocate(snq1(Nn))
! allocate(psi_ks00(Nn))
! allocate(psi_gw00(Nn))


! ! do ii=1,Nstates+1

! write (file_id, '(I0)') Nstates
! ! open (unit=20,file=trim(name)//'_p2_dft_1.psi'//trim(file_id),status='old')
! open (unit=20,file=trim(name)//'_p2_hf_1.psi'//trim(file_id),status='old')
! ! open (unit=20,file='CO_p2_dft_2.psi'//trim(file_id),status='old')
! ! open (unit=20,file='CO_p2_hf_1.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Li2_p2_hf_2.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Ne_p2_hf_2.psi'//trim(file_id),status='old')
! ! open (unit=20,file='H2_p2_dft_1.psi'//trim(file_id),status='old')

! do i=1,Nn
!     read(20,*) psi_ks00(i)
! enddo
! close (20)

! write (file_id, '(I0)') Nstates+1
! ! open (unit=20,file=trim(name)//'_p2_dft_1.psi'//trim(file_id),status='old')
! open (unit=20,file=trim(name)//'_p2_hf_1.psi'//trim(file_id),status='old')
! ! open (unit=20,file='CO_p2_dft_2.psi'//trim(file_id),status='old')
! ! open (unit=20,file='CO_p2_hf_1.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Li2_p2_hf_2.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Ne_p2_hf_2.psi'//trim(file_id),status='old')
! ! open (unit=20,file='H2_p2_dft_1.psi'//trim(file_id),status='old')
! allocate(psi_ks01(Nn))
! do i=1,Nn
!     read(20,*) psi_ks01(i)
! enddo
! close (20)

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_ks00,snq1)
! ! ! print *,dot_product(psi_ks00,snq1)

! ! psi_ks00=psi_ks00/sqrt(dot_product(psi_ks00,snq1))

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_ks00,snq1)
! ! ! print *,dot_product(psi_ks00,snq1)


! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_ks01,snq1)
! ! print *,dot_product(psi_ks01,snq1)

! ! psi_ks01=psi_ks01/sqrt(dot_product(psi_ks01,snq1))

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_ks01,snq1)
! ! print *,dot_product(psi_ks01,snq1)

! ! write (file_id, '(I0)') 9
! ! open (unit=20,file='CO_p2_dft_2.psi'//trim(file_id),status='old')
! ! ! open (unit=20,file='Ne_p2_dft_1.psi'//trim(file_id),status='old')
! ! ! open (unit=20,file='H2_p2_dft_1.psi'//trim(file_id),status='old')
! ! allocate(psi_ks02(Nn))
! ! do i=1,Nn
! !     read(20,*) psi_ks02(i)
! ! enddo
! ! close (20)



! write (file_id, '(I0)') 3
! ! open (unit=20,file=trim(name)//'_p2_hf_1.psi'//trim(file_id),status='old')
! open (unit=20,file=trim(name)//'_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='CO_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Li2_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Ne_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='H2_p2_gw.psi'//trim(file_id),status='old')

! do i=1,Nn
!     read(20,*) psi_gw00(i)
! enddo
! close (20)

! write (file_id, '(I0)') 4
! ! open (unit=20,file=trim(name)//'_p2_hf_1.psi'//trim(file_id),status='old')
! open (unit=20,file=trim(name)//'_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='CO_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Li2_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='Ne_p2_gw.psi'//trim(file_id),status='old')
! ! open (unit=20,file='H2_p2_gw.psi'//trim(file_id),status='old')
! allocate(psi_gw01(Nn))
! do i=1,Nn
!     read(20,*) psi_gw01(i)
! enddo
! close (20)

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_gw00,snq1)
! ! ! print *,dot_product(psi_gw00,snq1)

! ! psi_gw00=psi_gw00/sqrt(dot_product(psi_gw00,snq1))

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_gw00,snq1)
! ! ! print *,dot_product(psi_gw00,snq1)


! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_gw01,snq1)
! ! print *,dot_product(psi_gw01,snq1)

! ! psi_gw01=psi_gw01/sqrt(dot_product(psi_gw01,snq1))

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_gw01,snq1)
! ! print *,dot_product(psi_gw01,snq1)

! ! write (file_id, '(I0)') 5
! ! open (unit=20,file='Ne_p2_gw.psi'//trim(file_id),status='old')
! ! ! open (unit=20,file='H2_p2_gw.psi'//trim(file_id),status='old')
! ! allocate(psi_gw02(Nn))
! ! do i=1,Nn
! !     read(20,*) psi_gw02(i)
! ! enddo
! ! close (20)




! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi_ks00-psi_gw00),snq1)
! print *,dot_product((psi_ks00-psi_gw00),snq1)

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi_ks00+psi_gw00),snq1)
! print *,dot_product((psi_ks00+psi_gw00),snq1)

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi_ks01-psi_gw01),snq1)
! print *,dot_product((psi_ks01-psi_gw01),snq1)

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi_ks01+psi_gw01),snq1)
! print *,dot_product((psi_ks01+psi_gw01),snq1)

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi_ks02-psi_gw01),snq1)
! ! print *,dot_product((psi_ks02-psi_gw01),snq1)

! ! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi_ks02+psi_gw01),snq1)
! ! print *,dot_product((psi_ks02+psi_gw01),snq1)


! ! enddo


! deallocate(snq1)
! deallocate(psi_ks00)
! deallocate(psi_gw00)


! stop




! stop


! allocate(Etotal(1:Nstates,1:2))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pre-compute nuclei potential and store it for later use !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(Vi(Nn))

do i=1,Nn
    Vi(i) = nucleipotential((/point(i,1),point(i,2),point(i,3)/),radius)   !!! for isolated system
enddo

! print *,Vi(1),Vi(2)

! stop

call set_linkedlistV(Ne,Nlocal,Nquadrature,1)



ii=0
do i=1,Nn
    do l=1,neigh_V(i)%size
        ii=ii+1
        xy_V=neigh_V(i)%getab(l)
        V_i(ii)=dble(xy_V(1))   
    enddo
enddo

do i=1,Nn
    do l=1,neigh_V(i)%size
        xy_V=neigh_V(i)%getab(l)
        ii=neigh_AB(i)%get(l)
        H_dense(i,ii)=H_dense(i,ii)+dble(xy_V(1))
        ! H_dense(i,ii)=H_dense(i,ii)*0.5d0
        ! S_dense(i,ii)=S_dense(i,ii)+dble(xy_AB(2))
    enddo
enddo

do i=1,Nn
    do ii=1,neigh_V(i)%size
        call neigh_V(i)%deleteFirst()
    enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!             DIIS SCF loop               !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


itmax = 50
M9 = 50
alpha = ks_alpha0!0.5d0!0.45d0
alpha0 = 1.0d0
eps = ks_eps0!1.0d-8!1.0d-8!1E-8
it=0
itm = 0
nrestart = 0
ite=0


!!!! allocate related memory for DIIS-Pulay iteration
allocate(nq(Nn))
allocate(ni(Nn))
allocate(ni_old(Nn))
allocate(ni_previous(Nn,itmax))
allocate(ni_xc(Nn))
allocate(ri(Nn))
allocate(ri_old(Nn))
allocate(res_scf(itmax))

allocate(ri_g(Ne*Nquadrature))
allocate(ri_old_g(Ne*Nquadrature))
allocate(ni_g(Ne*Nquadrature))
allocate(ni_old_g(Ne*Nquadrature))




allocate(B00(M9,M9))
allocate(f(M9))

allocate(dN(M9,Nn))
allocate(dR(M9,Nn))

allocate(dN_g(M9,Ne*Nquadrature))
allocate(dR_g(M9,Ne*Nquadrature))

allocate(dN_trans(Nn,M9))

allocate(work(M9))
allocate(ipiv(M9))

allocate(temp1(Nn))

allocate(Vh(Nn))
allocate(Vx(Nn))
allocate(Vc(Nn))
allocate(V_total(Nn))

allocate(Vh_g(Ne*Nquadrature))
allocate(Vx_g(Ne*Nquadrature))
allocate(Vc_g(Ne*Nquadrature))
allocate(V_total_g(Ne*Nquadrature))

allocate(vh_temp(1:Nn))
allocate(H(neigh_size))
allocate(H_dft(neigh_size))


allocate(point_g(1:Ne*Nquadrature,1:3))


allocate(psi_gauss0(1:Ne*Nquadrature,1:Nstates))


allocate(nq_g(1:Ne*Nquadrature))
allocate(nq_g_gradient(1:Ne*Nquadrature,1:3))

allocate(ex_g_pbe(1:Ne*Nquadrature))
allocate(ec_g_pbe(1:Ne*Nquadrature))
allocate(vx_g_pbe_n(1:Ne*Nquadrature))
allocate(vc_g_pbe_n(1:Ne*Nquadrature))
allocate(vx_g_pbe_g(1:Ne*Nquadrature,1:3))
allocate(vc_g_pbe_g(1:Ne*Nquadrature,1:3))

allocate(nq_gradient(1:Nn,1:3))
allocate(psi_gradient(1:Nn,1:3,1:Nstates))

allocate(ex_pbe(1:Nn))
allocate(ec_pbe(1:Nn))
allocate(vx_pbe_n(1:Nn))
allocate(vc_pbe_n(1:Nn))
allocate(vx_pbe_g(1:Nn,1:3))
allocate(vc_pbe_g(1:Nn,1:3))


if (gw_sigma=='cd') then

    allocate(psi_point_g(1:Ne*Nquadrature,1:Nstates))

else if (gw_sigma=='none') then

    allocate(psi_point_g(1:Ne*Nquadrature,1:Nstates))

else if (gw_sigma=='casida') then

    allocate(psi_point_g(1:Ne*Nquadrature,1:Nn))

end if



!!!!!!!!!!!!!! Hartree Fock variables !!!!!!!!!!!!!!!

if (meanfield=="hf") then

allocate(hf_Vhfx(Nn,Ne*Nquadrature))
allocate(Vhfx_final(Nn,Nn))
allocate(NNmatrix00(Nn,Nn))
allocate(psii(1:Nn,1:Nstates))
allocate(psig(1:Nn,1:Nstates))
allocate(psi_previous(1:Nn,1:Nstates,1:M9))
allocate(error_previous(1:Nn,1:Nstates,1:M9))
do i=1,Nstates
    psii(:,i)=psi_i(:,i)
enddo
call psi_interpolation_HF(Ne,Nquadrature,Nstates,1)
call construct_hf_v_kernel(Nn,Ne,Nquadrature)

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(matBC(Nn,Nbc))

call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,0.0d0,0,0,0)
nq(1:Nn)=psi_0(1:Nn)!**2

!!!!!!!!!!!!!! Hartree Fock variables !!!!!!!!!!!!!!!
if (meanfield=="hf") then
nq(1:Nn)=0.0d0
do i=1,Nstates
    nq(:)=nq(:)+2.0d0*psi_i(:,i)**2
enddo
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call mkl_dcsrmm('N',Ne*Nquadrature,1,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!     ,psi_0,Nn,0.0d0,nq_g,Ne*Nquadrature)

! do i=1,Ne*Nquadrature
!     if (nq_g(i)<0.0d0) nq_g(i) = 0.0d0
! enddo

deallocate(psi_0)

! allocate(psii(1:Nn,1:Nstates))
! allocate(psig(1:Nn,1:Nstates))
! allocate(psi_previous(1:Nn,1:Nstates,1:M9))
! allocate(error_previous(1:Nn,1:Nstates,1:M9))

! psi_previous=0.0d0
! error_previous=0.0d0

! allocate(error_previous0(1:Nn,1:M9))
allocate(temp2(Nn))
! do i=1,neigh_size
!     if (isnan(A00(i))) stop 'NaN number happened'
!     if (A00(i)-1.0d0==A00(i)) stop 'infinity happened'
! enddo

! stop
! nq=0.0d0
! do i=1,Nstates
!     nq(:)=nq(:)+2.0d0*psi_i(:,i)**2
! enddo


! allocate(apt_v(64))
! allocate(aiparm_v(64))

if (meanfield=="dft") then



print *, '***'
print *," ------------------------------------------"
print *,"|           SCF-iteration starts           |"
print *," ------------------------------------------"
print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'



! call mkl_dcsrgemv('N',Nn,B,IA,JA,psi_0(:),Vh)
! print *,dot_product(psi_0(:),Vh)

! stop



do while (.true.)
    
    print *, 'SCF-Cycle ----',it

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!       initial setting for electron density        !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (it>0) ni_old=ni
    
    ni=nq

    if (it>0) call psi_interpolation(Ne,Nquadrature,Nstates,1,degree)

    if (xc0=='lda') then

        if (it>0) call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,2.0d0,ks_poisson_bc0,0,0)

        ! call hartreepotential(Nn,Ne,Nquadrature,Nbc,0)
        call hartreepotential(Nn,Ne,Nquadrature,Nstates,Nbc,2.0d0,it,0)

        ! call mkl_dcsrmm('N',Ne*Nquadrature,1,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        ! ,Vh,Nn,0.0d0,Vh_g,Ne*Nquadrature)

        
        call exchangepotential(Nn)
        call correlationpotential(Nn)


        V_total(1:Nn) = Vh(1:Nn)+Vx(1:Nn)+Vc(1:Nn)!+Vi(1:Nn)
        ! V_total_g(1:Ne*Nquadrature) = Vh_g(1:Ne*Nquadrature)+Vx_g(1:Ne*Nquadrature)+Vc_g(1:Ne*Nquadrature)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!  construct Hamiltonian  !!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (it>0) then
            do i=1,Nn
                ! do ii=1,neigh_V(i)%size
                !     call neigh_V(i)%deleteFirst()
                ! enddo
                call neigh_V(i)%initializeAB()
            enddo
        end if
        
        call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

        ii=0
        do i=1,Nn
            do l=1,neigh_V(i)%size
                ii=ii+1
                xy_V=neigh_V(i)%getab(l)
                V(ii)=dble(xy_V(1))
            enddo
        enddo


    else if (xc0=='pbe') then

        print *,'PBE functional is not supported yet, LIBXC needs to be installed.'
        print *,'Please contact the developer for more information.'

        stop

            ! if (it==0) then

            ! !!!!!!!!!!!!!!!!!!
            ! !!!!!! LDA !!!!!!!
            ! !!!!!!!!!!!!!!!!!!

            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! !!!!!!  compute_Hartree/Exchange/Correlation Potentials  !!!!!!
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! if (it>0) call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,2.0d0,ks_poisson_bc0,0,0)

            ! ! call hartreepotential(Nn,Ne,Nquadrature,Nbc,0)
            ! call hartreepotential(Nn,Ne,Nquadrature,Nstates,Nbc,2.0d0,it,0)

            ! ! call mkl_dcsrmm('N',Ne*Nquadrature,1,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
            ! ! ,Vh,Nn,0.0d0,Vh_g,Ne*Nquadrature)

            
            ! call exchangepotential(Nn)
            ! call correlationpotential(Nn)


            ! V_total(1:Nn) = Vh(1:Nn)+Vx(1:Nn)+Vc(1:Nn)!+Vi(1:Nn)
            ! ! V_total_g(1:Ne*Nquadrature) = Vh_g(1:Ne*Nquadrature)+Vx_g(1:Ne*Nquadrature)+Vc_g(1:Ne*Nquadrature)

            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! !!!!!!!!!!!!!!!!  construct Hamiltonian  !!!!!!!!!!!!!!!!!!!!!!
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! if (it>0) then
            !     do i=1,Nn
            !         ! do ii=1,neigh_V(i)%size
            !         !     call neigh_V(i)%deleteFirst()
            !         ! enddo
            !         call neigh_V(i)%initializeAB()
            !     enddo
            ! end if
            
            ! call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

            ! ii=0
            ! do i=1,Nn
            !     do l=1,neigh_V(i)%size
            !         ii=ii+1
            !         xy_V=neigh_V(i)%getab(l)
            !         V(ii)=dble(xy_V(1))
            !     enddo
            ! enddo


            ! else

            ! !!!!!!!!!!!!!!!!!!
            ! !!!!!! PBE !!!!!!!
            ! !!!!!!!!!!!!!!!!!!

            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! !!!!!!  compute_Hartree/Exchange/Correlation Potentials  !!!!!!
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! if (it>0) call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,2.0d0,ks_poisson_bc0,0,0)

            ! ! call hartreepotential(Nn,Ne,Nquadrature,Nbc,0)
            ! call hartreepotential(Nn,Ne,Nquadrature,Nstates,Nbc,2.0d0,it,0)

            ! V_total(1:Nn) = Vh(1:Nn)!+Vx(1:Nn)+Vc(1:Nn)!+Vi(1:Nn)

            ! ! inquire(file='xc_f90_lib_m.mod', exist=module_exists)

            ! ! found = check_file_in_path("xc_f90_lib_m.mod")

            ! ! if (found) then
            !     ! call compute_nq_g_gradient(Nn,Ne,Nlocal,Nquadrature,Nstates)
            !     call compute_nq_gradient(Nn,Ne,Nlocal,Nstates)

            !     ! print *, nq_g(1)
            !     ! print *, nq_g_gradient(1,:)
            !     ! print *, psi_point_g(1,1)

            !     ! call compute_nq_gradient(Nn,Ne,Nlocal)

            !     ! call libxc_exchange_g_pbe(Ne,Nquadrature)
            !     ! call libxc_correlation_g_pbe(Ne,Nquadrature)
            !     call libxc_exchange_pbe(Nn)
            !     call libxc_correlation_pbe(Nn)
            !     ! call libxc_exchange_pbe(Nn)
            !     ! call libxc_correlation_pbe(Nn)
            !     ! call libxc_exchange_lda(Ne,Nquadrature)
            !     ! call libxc_correlation_lda(Ne,Nquadrature)
            ! ! else
            ! !     print *, "LIBXC .mod file does not exist."
            ! !     stop
            ! ! end if

            


            ! V_total(1:Nn)=V_total(1:Nn)+vx_pbe_n(1:Nn)+vc_pbe_n(1:Nn)

            ! ! print *,Vh(1),Vx(1),Vc(1)

            ! ! exit


            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! !!!!!!!!!!!!!!!!  construct Hamiltonian  !!!!!!!!!!!!!!!!!!!!!!
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! if (it>0) then
            !     do i=1,Nn
            !         ! do ii=1,neigh_V(i)%size
            !         !     call neigh_V(i)%deleteFirst()
            !         ! enddo
            !         call neigh_V(i)%initializeAB()
            !     enddo
            ! end if
            
            ! call set_linkedlistV(Ne,Nlocal,Nquadrature,2)

            ! ii=0
            ! do i=1,Nn
            !     do l=1,neigh_V(i)%size
            !         ii=ii+1
            !         xy_V=neigh_V(i)%getab(l)
            !         V(ii)=dble(xy_V(1))
            !     enddo
            ! enddo

            ! end if !!! it>0

    end if !!! LDA/PBE 'if' ends


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! build total Hamiltonian !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! V=V+V_i

    H=0.5d0*A+V_i+V

    ! print *,'toto'

    ! call MPI_BARRIER(MPI_COMM_WORLD ,code)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!    solve Kohn-Sham equation using FEAST    !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! M0=2*Nstates+15!+20!+10!+40
    call feastinit(fpm)
    ! fpm(1)=1
    ! call pfeastinit(fpm,MPI_COMM_WORLD,nL3) !! for pfeast which provides solving each interval in parallel

    call dfeast_scsrgv(UPLO, Nn, H, IA, JA, B, IA, JA, fpm, epsout, loop, Emin, Emax, M0, E, psi, M00, res, info)

    ! if (info/=0) then
    !     print *,'FEAST_error info --------',info
    !     stop
    ! end if

        print *, '--- Energy states up to occupied + LUMO (eV) ---'
        do i=1,Nstates+1
            print *, i,E(i)*hartree
        enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!    post-process wavefunctions and density    !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    nq=0.0d0
    ! nq_g=0.0d0
    do i=1,Nstates
        nq(:)=nq(:)+2.0d0*psi(:,i)**2
        ! nq_g(:)=nq_g(:)+2.0d0*psi_gauss0(:,i)**2
        ! psig(:,i)=psi(:,i)
    enddo





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!  using Pulay mixing to accelarate SCF convergence   !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (it>0) ri_old=ri
    ri=nq-ni

    norm=norm2(ri)        
    normf=norm2(ni)       
    res_scf(it+1)=norm/normf

    print *, 'SCF-loop residual ------- ',res_scf(it+1)!,norm,normf


    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    if (res_scf(it+1)<eps) then
        print *, 'Ground State Convergence Reached!!!'

        H_dft = H



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!     compute Eh,Ex,Ec        !!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(snq2(1:Nn))

        print *,"=============================================="

        if (xc0=='lda') then

            print *,"-----  Orbital,  DFT_Ekinetic (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,0.5d0*A,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            print *,"-----  Orbital,  DFT_Eext (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V_i,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            !!!!!!!!!! Eh !!!!!!!!!!

            V_total(1:Nn)=Vh(1:Nn)!+vx_pbe_n(1:Nn)+vc_pbe_n(1:Nn)

            if (it>0) then
                do i=1,Nn
                    ! do ii=1,neigh_V(i)%size
                    !     call neigh_V(i)%deleteFirst()
                    ! enddo
                    call neigh_V(i)%initializeAB()
                enddo
            end if
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo

            ! print *,"DFT_Eh ---"
            print *,"-----  Orbital,  DFT_Eh (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo


            !!!!!!!!!! Ex !!!!!!!!!!

            V_total(1:Nn)=Vx(1:Nn)!+vx_pbe_n(1:Nn)+vc_pbe_n(1:Nn)

            if (it>0) then
                do i=1,Nn
                    ! do ii=1,neigh_V(i)%size
                    !     call neigh_V(i)%deleteFirst()
                    ! enddo
                    call neigh_V(i)%initializeAB()
                enddo
            end if
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo

            ! print *,"DFT_Eh ---"
            print *,"-----  Orbital,  DFT_Ex (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            !!!!!!!!!! Ec !!!!!!!!!!

            V_total(1:Nn)=Vc(1:Nn)!+vx_pbe_n(1:Nn)+vc_pbe_n(1:Nn)

            if (it>0) then
                do i=1,Nn
                    ! do ii=1,neigh_V(i)%size
                    !     call neigh_V(i)%deleteFirst()
                    ! enddo
                    call neigh_V(i)%initializeAB()
                enddo
            end if
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo

            ! print *,"DFT_Eh ---"
            print *,"-----  Orbital,  DFT_Ec (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

        else if (xc0=='pbe') then

            print *,"-----  Orbital,  DFT_Ekinetic (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,0.5d0*A,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            print *,"-----  Orbital,  DFT_Eext (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V_i,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            !!!!!!!!!! Eh !!!!!!!!!!

            V_total(1:Nn)=Vh(1:Nn)!+vx_pbe_n(1:Nn)+vc_pbe_n(1:Nn)

            if (it>0) then
                do i=1,Nn
                    ! do ii=1,neigh_V(i)%size
                    !     call neigh_V(i)%deleteFirst()
                    ! enddo
                    call neigh_V(i)%initializeAB()
                enddo
            end if
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo

            ! print *,"DFT_Eh ---"
            print *,"-----  Orbital,  DFT_Eh (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            !!!!!!!!!! Ex !!!!!!!!!!

            V_total(1:Nn)=vx_pbe_n(1:Nn)!+vc_pbe_n(1:Nn)

            if (it>0) then
                do i=1,Nn
                    ! do ii=1,neigh_V(i)%size
                    !     call neigh_V(i)%deleteFirst()
                    ! enddo
                    call neigh_V(i)%initializeAB()
                enddo
            end if
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo

            ! print *,"DFT_Ex ---"
            print *,"-----  Orbital,  DFT_Ex (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            !!!!!!!!!! Ec !!!!!!!!!!

            V_total(1:Nn)=vc_pbe_n(1:Nn)

            if (it>0) then
                do i=1,Nn
                    ! do ii=1,neigh_V(i)%size
                    !     call neigh_V(i)%deleteFirst()
                    ! enddo
                    call neigh_V(i)%initializeAB()
                enddo
            end if
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,0)

            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo

            ! print *,"DFT_Ec ---"
            print *,"-----  Orbital,  DFT_Ec (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo


            !!!!!!!!!! Exc_pbe_g !!!!!!!!!!

            do i=1,Nn
                ! do ii=1,neigh_V(i)%size
                !     call neigh_V(i)%deleteFirst()
                ! enddo
                call neigh_V(i)%initializeAB()
            enddo
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,3)
        
            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo

            ! print *,"DFT_Ex+Ec pbe_g ---"
            print *,"-----  Orbital,  DFT_Ex+Ec pbe_g (eV) -----"

            do i=1,Nstates+1

            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

        end if !!! lda/pbe 
        
        print *,"=============================================="
        



        if (gw_sigma=='cd') then 

            print *,'--- solve more unoccupied states ---'

            ! Nempty=1500!1375!Nn*1/2!-1!2000!Nn*1/2!2500!3279!1000!1496!1000!1496!250!1496!1000!1496
            M0=2*Nstates+25!Nstates+Nempty+50!Nn!*2/3!3280!1497 ! 4  !10 !300
            Emin=-40.0d0 ! -30.0d0 
            Emax=0.25d0!324.0d0!2.0d7!2.0d4!2.0d7!26000.0d0 ! -0.05d0  !-0.05d0
            deallocate(E,psi,res)
            allocate(E(M0), psi(Nn,M0), res(M0))


            ! ! Nempty=1500!1375!Nn*1/2!-1!2000!Nn*1/2!2500!3279!1000!1496!1000!1496!250!1496!1000!1496
            ! M0=Nn!Nstates+Nempty+50!Nn!*2/3!3280!1497 ! 4  !10 !300
            ! Emin=-50.0d0 ! -30.0d0 
            ! Emax=2.0d7!324.0d0!2.0d7!2.0d4!2.0d7!26000.0d0 ! -0.05d0  !-0.05d0
            ! deallocate(E,psi,res)
            ! allocate(E(M0), psi(Nn,M0), res(M0))

            call feastinit(fpm)
            call dfeast_scsrgv(UPLO, Nn, H, IA, JA, B, IA, JA, fpm, epsout, loop, Emin, Emax, M0, E, psi, M00, res, info)

            ! if (info/=0) then
            !     print *,'FEAST_error info --------',info
            !     stop
            ! end if
    
    
            print *, '--- Energy states up to all solved unoccupied states (eV) ---'
            do i=1,M00
                print *, i,E(i)*hartree
            enddo
            
            ! print *,M00,E(M00)!,E(1000)!,E_GW(2000)

        end if

        if (gw_sigma=='casida') then 

            print *,'--- solve all unoccupied states ---'

            ! Nempty=1500!1375!Nn*1/2!-1!2000!Nn*1/2!2500!3279!1000!1496!1000!1496!250!1496!1000!1496
            M0=Nn!Nstates+Nempty+50!Nn!*2/3!3280!1497 ! 4  !10 !300
            Emin=-50.0d0 ! -30.0d0 
            Emax=2.0d7!324.0d0!2.0d7!2.0d4!2.0d7!26000.0d0 ! -0.05d0  !-0.05d0
            deallocate(E,psi,res)
            allocate(E(M0), psi(Nn,M0), res(M0))

            call feastinit(fpm)
            call dfeast_scsrgv(UPLO, Nn, H, IA, JA, B, IA, JA, fpm, epsout, loop, Emin, Emax, M0, E, psi, M00, res, info)

            if (info/=0) then
                print *,'FEAST_error info --------',info
                stop
            end if
    
    
            ! print *, '--- Energy states up to occupied (eV) ---'
            ! do i=1,Nstates+1
            !     print *, i,E_GW(i)*hartree
            ! enddo
            
            print *,M00,E(M00)!,E(1000)!,E_GW(2000)

        end if


        V_total(1:Nn) = Vh(1:Nn)
        do i=1,Nn
            ! do ii=1,neigh_V(i)%size
            !     call neigh_V(i)%deleteFirst()
            ! enddo
            call neigh_V(i)%initializeAB()
        enddo
        call set_linkedlistV(Ne,Nlocal,Nquadrature,0)
        do i=1,Nn
            do l=1,neigh_V(i)%size
                xy_V=neigh_V(i)%getab(l)
                ii=neigh_V(i)%get(l)
                H_dense(i,ii)=H_dense(i,ii)+dble(xy_V(1))
            enddo
        enddo

        ! H_dense=H_dense+Vhfx_final_GW

        ii=0
        do i=1,Nn
            do l=1,neigh_V(i)%size
                ii=ii+1
                xy_V=neigh_V(i)%getab(l)
                V(ii)=dble(xy_V(1))
            enddo
        enddo

        V=V+V_i

        H=0.5d0*A+V

        ! print *,H(1:10)

        ! stop


        

        exit
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!     DIIS - Pulay mixing     !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (it==0) then
        nq=nq-alpha*ri
    else
        !! build subspace
        dN(itm,:)=ni-ni_old
        dR(itm,:)=ri-ri_old

        do i=1,itm
            do ii=1,itm
               B00(i,ii)=dot_product(dR(i,:),dR(ii,:)) 
            end do
            f(i)=dot_product(dR(i,:),ri)
        enddo
        !! solve linear system Bx=f   
        call  DSYSV('L',itm,1, B00, M9, IPIV, f, M9, WORK,M9, info_scf )
        nq=0.0d0
        do i=1,itm
            nq=nq+f(i)*(dN(i,:)+alpha*dR(i,:))
        enddo
        nq=ni+alpha*ri-alpha0*nq
        do i=1,Nn
            if (nq(i)<0.0d0) then
                nq(i) = 0.0d0
            end if
        enddo
    end if

    it=it+1
    itm=itm+1



enddo  !!!! dft SCF while loop



else if (meanfield=="hf") then


print *, '***'
print *," ------------------------------------------"
print *,"|           SCF-iteration starts           |"
print *," ------------------------------------------"
print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'



    DSYEV_lwork = 3*Nn-1
    allocate(DSYEV_work(DSYEV_lwork))

    do while (.true.)
    
        print *, 'SCF-Cycle ----',it
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!       initial setting for electron density        !!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        if (it>0) ni_old=ni
        
        
        ni=nq
        ! psii=psig
    
        if (it>0) then
            do i=1,Nstates
                psii(:,i)=psig(:,i)
            enddo
        end if
    
    
        do i=1,Nstates
            psi_previous(:,i,itm+1)=psii(:,i)
        enddo
    
        if (it>0) call psi_interpolation_HF(Ne,Nquadrature,Nstates,1)
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!  compute_Hartree/Exchange/Correlation Potentials  !!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        ! call set_hf_BC(Nn,Nbc,Ne,Nquadrature,Nstates)
        ! call hf_hartreepotential(Nn,Nbc)

        call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,2.0d0,hf_poisson_bc0,0,0)
        ! call hartreepotential(Nn,Ne,Nquadrature,Nbc,0)
        call hartreepotential(Nn,Ne,Nquadrature,Nstates,Nbc,2.0d0,it,0)
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! compute exchange kernel and diagonalize exchange potential  !!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        call construct_hf_exchangekernel(Nstates)
    
        do i=1,Nn
            hf_Vhfx(i,:)=hf_Vhfx(i,:)*volumegweight
        enddo
    
        call mkl_dcsrmm('T',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,transpose(hf_Vhfx),Ne*Nquadrature,0.0d0,NNmatrix00,Nn)

        call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,transpose(NNmatrix00),Nn,0.0d0,Vhfx_final,Nn)
    

        NNmatrix00=transpose(Vhfx_final)

        Vhfx_final = (Vhfx_final + NNmatrix00)/2.0d0
    
        V_total(1:Nn) = Vh(1:Nn)!/2.0d0!+Vx(1:Nn)+Vc(1:Nn)!+Vi(1:Nn)
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!  construct Hamiltonian  !!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        if (it>0) then
            do i=1,Nn
                do ii=1,neigh_V(i)%size
                    call neigh_V(i)%deleteFirst()
                enddo
            enddo
        end if
        
        call set_linkedlistV(Ne,Nlocal,Nquadrature,0)
    
        ii=0
        do i=1,Nn
            do l=1,neigh_V(i)%size
                ii=ii+1
                xy_V=neigh_V(i)%getab(l)
                V(ii)=dble(xy_V(1))
            enddo
        enddo
        V=V+V_i
    
        H=0.5d0*A+V
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!  Hartree Fock !!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        H_dense_psi=H_dense
    
        do i=1,Nn
            do l=1,neigh_V(i)%size
                xy_V=neigh_V(i)%getab(l)
                ii=neigh_AB(i)%get(l)
                H_dense_psi(i,ii)=H_dense_psi(i,ii)+dble(xy_V(1))
            enddo
        enddo
    
    
        H_dense_psi=H_dense_psi+Vhfx_final
        B_dense = S_dense
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!    solve Kohn-Sham equation using FEAST    !!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! allocate(DSYEV_work(2*Nn))  !!!!! slightly overestimate the work space

        ! !!!!! Query the optimal workspace
        ! DSYEV_lwork = -1
        ! call dsygv(1,'V','U',Nn,H_dense_psi,Nn,S_dense,Nn,E,DSYEV_work,DSYEV_lwork,info)
        ! DSYEV_lwork = int(DSYEV_work(1))


        
        !!!!! solve eigenvalue problem
        call dsygv(1,'V','U',Nn,H_dense_psi,Nn,B_dense,Nn,E_hf,DSYEV_work,DSYEV_lwork,info)

        
    
        ! M0=2*Nstates+10
        ! call feastinit(fpm)
        ! ! fpm(1)=1
        ! ! call pfeastinit(fpm,MPI_COMM_WORLD,nL3) !! for pfeast which provides solving each interval in parallel
    
        ! ! call dfeast_scsrgv(UPLO, Nn, H, IA, JA, B, IA, JA, fpm, epsout, loop, Emin, Emax, M0, E, psi, M00, res, info)
    
        ! call dfeast_sygv('F', Nn, H_dense, Nn, S_dense, Nn, fpm, epsout, loop, Emin, Emax, M0, E, psi, M00, res, info)
    
        ! ! if ((rank==0).and.(info/=0)) then
        ! !     print *,'FEAST_error info --------',info
        ! !     stop
        ! ! end if
    
        print *, '--- Energy states up to occupied + LUMO (eV) ---'
        do i=1,Nstates+1
            print *, i,E_hf(i)*hartree
        enddo
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!    post-process wavefunctions and density    !!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        nq=0.0d0
        psig(:,:)=0.0d0
        do i=1,Nstates
            nq(:)=nq(:)+2.0d0*H_dense_psi(:,i)**2
            psig(:,i)=H_dense_psi(:,i)
        enddo
    
        ! H_dense=H_dense-Vhfx_final
    
        ! do i=1,Nn
        !     do l=1,neigh_V(i)%size
        !         xy_V=neigh_V(i)%getab(l)
        !         ii=neigh_AB(i)%get(l)
        !         H_dense(i,ii)=H_dense(i,ii)-dble(xy_V(1))
        !     enddo
        ! enddo
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!  using Pulay mixing to accelarate SCF convergence   !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (it>0) ri_old=ri
        ri=nq-ni
    
        norm=norm2(ri)        
        normf=norm2(ni)       
        res_scf(it+1)=norm/normf
    
        print *, 'SCF-loop residual ------- ',res_scf(it+1)!,norm,normf
    
    
        print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    
        if (res_scf(it+1)<eps) then
        ! if (res_scf(it+1)<1.0d-3) then
            print *, 'Ground State Convergence Reached!!!'

            do i=1,Nstates+40
                print *, i,E_hf(i)*hartree
            enddo

            print *,'----'
            print *,E_hf(100) *hartree
            print *,E_hf(200) *hartree
            print *,E_hf(300) *hartree
            print *,E_hf(400) *hartree
            print *,E_hf(500) *hartree
            print *,E_hf(600) *hartree
            print *,E_hf(700) *hartree
            print *,E_hf(800) *hartree
            print *,E_hf(900) *hartree
            print *,E_hf(1000)*hartree
            print *,E_hf(1100)*hartree
            print *,E_hf(1200)*hartree
            print *,E_hf(1300)*hartree
            print *,E_hf(1400)*hartree
            print *,E_hf(1500)*hartree
            print *,'----'

            deallocate(E,psi)
            allocate(E(Nn),psi(Nn,Nn))

            E(1:Nn)=E_hf(1:Nn)
            psi(:,1:Nn)=H_dense_psi(:,1:Nn)

            call psi_interpolation_HF(Ne,Nquadrature,Nstates,1)
            call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,2.0d0,hf_poisson_bc0,0,0)
            ! call hartreepotential(Nn,Ne,Nquadrature,Nbc,0)
            call hartreepotential(Nn,Ne,Nquadrature,Nstates,Nbc,2.0d0,it,0)
            V_total(1:Nn) = Vh(1:Nn)!/2.0d0!+Vx(1:Nn)+Vc(1:Nn)!+Vi(1:Nn)
            if (it>0) then
                do i=1,Nn
                    do ii=1,neigh_V(i)%size
                        call neigh_V(i)%deleteFirst()
                    enddo
                enddo
            end if
            
            call set_linkedlistV(Ne,Nlocal,Nquadrature,0)
        
            ii=0
            do i=1,Nn
                do l=1,neigh_V(i)%size
                    ii=ii+1
                    xy_V=neigh_V(i)%getab(l)
                    V(ii)=dble(xy_V(1))
                enddo
            enddo
            V=V+V_i
        
            H=0.5d0*A+V

            call construct_hf_exchangekernel(Nstates)
    
            do i=1,Nn
                hf_Vhfx(i,:)=hf_Vhfx(i,:)*volumegweight
            enddo
        
            call mkl_dcsrmm('T',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                        ,transpose(hf_Vhfx),Ne*Nquadrature,0.0d0,NNmatrix00,Nn)

            call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,transpose(NNmatrix00),Nn,0.0d0,Vhfx_final,Nn)

            NNmatrix00=transpose(Vhfx_final)

            Vhfx_final = (Vhfx_final + NNmatrix00)/2.0d0

            H_dense_psi=H_dense

            do i=1,Nn
                do l=1,neigh_V(i)%size
                    xy_V=neigh_V(i)%getab(l)
                    ii=neigh_AB(i)%get(l)
                    H_dense(i,ii)=H_dense(i,ii)+dble(xy_V(1))
                    H_dense_psi(i,ii)=H_dense_psi(i,ii)+dble(xy_V(1))
                enddo
            enddo
        
            allocate(Hhf_dense(Nn,Nn))
        
            Hhf_dense = H_dense_psi+Vhfx_final

            exit

        end if
    
    
        do i=1,Nstates
            if (psig(1,i)*psii(1,i)<0) psig(:,i)=-psig(:,i)
            error_previous(:,i,itm+1)=psig(:,i)-psii(:,i)
        enddo
    
    
        if (it>2) then
    
            print *,"use DIIS(Pulay mixing) to extrapolate the solution"
            print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    
            B00(:,:)=0.0d0
    
            value0=0.0d0
            !! solve linear system Bx=f
            psig(:,:)=0.0d0
            
            do i=itm+1-2,itm+1 !1,itm+1
                do ii=itm+1-2,itm+1 !1,itm+1
                    do k=1,Nstates
                    ! B000(i,ii,k)=dot_product(error_previous(:,k,i),error_previous(:,k,ii)) 
                    call mkl_dcsrgemv('N', Nn, B, IA, JA, error_previous(:,k,ii) , temp1)
                    B00(i-itm+2,ii-itm+2)=B00(i-itm+2,ii-itm+2)+dot_product(error_previous(:,k,i),temp1)
                    enddo
                end do
            enddo
    
    
            
            
            value0=0.0d0
            !! solve linear system Ax=f
            psig(:,:)=0.0d0
            f(:)=0.0d0
            B00(:,4)=1.0d0
            B00(4,:)=1.0d0
            B00(4,4)=0.0d0
            f(4)=1.0d0
            call  DSYSV('L',4,1, B00(1:4,1:4), 4, IPIV, f(1:4), 4, WORK,4, info_scf )

            do k=1,Nstates
            do i=itm+1-2,itm+1
                temp2(:)=psi_previous(:,k,i)+error_previous(:,k,i)
                psig(:,k)=psig(:,k)+f(i-itm+2)*temp2(:)
            enddo
            enddo
            
        end if
    
        nq=0.0d0
        do k=1,Nstates
            nq=nq+2.0d0*psig(:,k)**2
        enddo
    
        it=it+1
        itm=itm+1

        M0=2*Nstates+10
        M00=2*Nstates+10
    
    enddo  !!!!! end scf while loop

    ! E(1:M0)=E_hf(1:M0)
    ! psi(:,1:M0)=H_dense_psi(:,1:M0)

    ! if (gw_sigma=='casida') then

        

    ! end if

    

    deallocate(DSYEV_work)

end if  !!!!! end DFT/HF selection

if (meanfield=='cd') then

    Nempty = Nn-Nstates

else if (meanfield=='casida') then

    Nempty = casida_N_empty0!Nn-Nstates
    
end if

! open(12,file=trim(name)//'_p2_density.node',status='replace')
! do i=1,Nn
! write(12,*) H_dense_psi(i,1)
! end do
! close(12)


! open(10,file=trim(name)//'_p2_dft.psi1',status='replace')
! do i=1,Nn
!     write(10,*) psi(i,1)
! enddo
! close(10)

do k=1,Nstates+7!1
write (file_id, '(I0)') k
! open (unit=20,file='CO.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='CO_0.6415.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='He.1_p2_nessie.psi'//trim(file_id),status='old')
! open(unit=20,file=trim(name)//'_p2_dft.psi'//trim(file_id),status='replace')
open(unit=20,file=trim(name)//'_'//trim(dgr_fem)//'_'//trim(meanfield)//'_1.psi'//trim(file_id),status='replace')
! open (unit=20,file='LiH.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
do i=1,Nn
    write(20,*) psi(i,k)
enddo
close (20)
enddo






do i=1,Nn
    call neigh_AB(i)%deallocate_linked_list()
    call neigh_V(i)%deallocate_linked_list()
    call neigh_NnToNg(i)%deallocate_linked_list()
enddo







!!!!!!!!!!!! deallocate scf-DIIS related memory !!!!!!!!!!!!!

deallocate(nq)
deallocate(ni)
deallocate(ni_old)
deallocate(ni_previous)
deallocate(ni_xc)
deallocate(ri)
deallocate(ri_old)
deallocate(res_scf)
deallocate(B00)
deallocate(f)
deallocate(dN)
deallocate(dR)
deallocate(dN_trans)
deallocate(work)
deallocate(ipiv)

deallocate(nq_g)
deallocate(nq_g_gradient)

deallocate(nq_gradient)
deallocate(psi_gradient)



! allocate(temp1(Nn))

!!!!!!!!!!!! deallocate DFT potential related memory !!!!!!!!!!!!!

deallocate(Vh)
deallocate(Vx)
deallocate(Vc)
deallocate(V_total)
deallocate(vh_temp)

! allocate(H(neigh_size))
! allocate(H_dft(neigh_size))


! allocate(point_g(1:Ne*Nquadrature,1:3))
! deallocate(psi_point_g)



deallocate(ex_g_pbe)
deallocate(ec_g_pbe)
deallocate(vx_g_pbe_n)
deallocate(vc_g_pbe_n)
deallocate(vx_g_pbe_g)
deallocate(vc_g_pbe_g)

deallocate(ex_pbe)
deallocate(ec_pbe)
deallocate(vx_pbe_n)
deallocate(vc_pbe_n)
deallocate(vx_pbe_g)
deallocate(vc_pbe_g)

deallocate(matBC)





if (gw_sigma=='none') then

    stop

end if





































!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    GW starts                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


print *,' '
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'!!!!!!!!!!!!!  GW starts  !!!!!!!!!!!!!!!'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!     compute GW <Vc> in dft      !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    

print *, "*** allocate memory for GW calculations ***"


if (gw_sigma=='cd') then

    ! allocate(psi_point_g(1:Ne*Nquadrature,1:Nstates))
    call psi_interpolation(Ne,Nquadrature,Nstates,1,degree)

else if (gw_sigma=='casida') then

    ! allocate(psi_point_g(1:Ne*Nquadrature,1:Nn))
    call psi_interpolation(Ne,Nquadrature,Nn,1,degree)
    ! call psi_interpolation(Ne,Nquadrature,Nstates,1,degree)

end if




print *, "*** construct bare Coulomb potential -- 1/|r-r'| ***"

! allocate(coulombmatrix00(Nn,Ne*Nquadrature))
! call compute_coulombmatrix00(Nn,Ne,Nquadrature)
allocate(v_kernel(Ne*Nquadrature,Nn))
call construct_v_kernel(Nn,Ne,Nquadrature)

! do i=1,Nn
!     coulombmatrix00(i,:)=coulombmatrix00(i,:)*volumegweight(:)
! enddo

do i=1,Nn
    v_kernel(:,i)=v_kernel(:,i)*volumegweight(:)
enddo

! allocate(v_kernel_transpose(Ne*Nquadrature,Nn))

allocate(v00(Nn,Nn))
call mkl_dcsrmm('T',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                            ,v_kernel,Ne*Nquadrature,0.0d0,v00,Nn)


                                            

! print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
! print *,'!!!!! construct GW Exchange Self-energy -- Sigma_x operator !!!!!'
! print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

print *, "*** construct GW Exchange Self-energy -- Sigma_x ***"

allocate(Vhfx_final_GW(Nn,Nn))

if (meanfield=="hf") then

    Vhfx_final_GW = Vhfx_final

else if (meanfield=="dft") then

    call construct_exchangematrix(Nn,Ne,Nquadrature,Nstates)

end if




allocate(snq1(Nn))

print *,"=============================================="
! print *,'Exchange orbital energy -- <psi_i|Sigma_x|psi_i>'
! print *,' '


print *,"----  Orbital,  GW <psi|Sigma_x|psi> (eV) ----"

do i=1,Nstates+1

    orbital=i

! CALL DGEMM('N','N',Nn,1,Nn,1.0d0,Vhfx_final_GW,Nn,psi(:,HOMO_state),Nn,0.0d0,snq1,Nn)
call DSYMM('L','L',Nn,1,1.0d0,Vhfx_final_GW,Nn,psi(:,orbital),Nn,0.0d0,snq1,Nn)

print *,i,dot_product(psi(:,orbital),snq1)*hartree

enddo

print *,"=============================================="

! stop





! 5000 continue

! deallocate(v_kernel)
! deallocate(v00)
! deallocate(Vhfx_GW)

! deallocate(NNmatrix_temp02)
! deallocate(psi_point_g)

! deallocate(snq1)

deallocate(volume)

100 format (I4,I4,ES25.16,ES25.16,ES25.16)


! stop






































if (gw_sigma=='cd') then


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!      Set energy guess of GW Sigma     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

print *,'*** compute the energy grid points ***'

N_energy = cd_N_grid0

allocate(E_guess(N_energy))

E_guess_min = cd_Emin_grid0/hartree
E_guess_max = cd_Emax_grid0/hartree

do i=1,N_energy
    E_guess(i)=E_guess_min+(i-1)*(E_guess_max-E_guess_min)/(N_energy-1)
enddo

! ! E_guess(1)=-22.9d0/hartree
! E_guess(1)=-21.9d0/hartree
! E_guess(2)=-21.7d0/hartree
! E_guess(3)=-21.5d0/hartree
! E_guess(4)=  2.0d0/hartree
! E_guess(5)=  2.2d0/hartree
! E_guess(6)=  2.35d0/hartree





allocate(E_Sigma(1:size(E_guess)*2))




! allocate(E_guess(1))
! ! E_guess=(/-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.95d0,-0.9d0,-0.85d0,-0.8d0,-0.75d0,-0.7d0,-0.65d0,-0.6d0,-0.5d0,-0.4d0,-0.853d0/)
! ! E_guess=(/-1.89d0,-1.72d0,-1.20d0/)
! ! E_guess=(/-1.2d0,-0.878308482d0,-0.575938296056237d0-0.01d0,-0.575938296056237d0+0.01d0/)
! ! E_guess=(/-0.878308482d0/)
! ! E_guess=(/-38.5d0/hartree,-33.0d0/hartree,-28.0d0/hartree,-23.0d0/hartree,&
! !                 -18.0d0/hartree,-14.8d0/hartree,-13.8d0/hartree,-13.6d0/hartree/)
! ! E_guess=(/-23.84d0/hartree/)
! ! E_guess=(/-0.545727237d0,-0.536539909d0/)
! ! E_guess=(/(E(1)*hartree-0.1)/hartree/)
! E_guess=(/E(1)-0.003d0/)!(/-23.84/hartree/)!(/-16.05/hartree/)
! allocate(E_Sigma(1:size(E_guess)))




print *,E(1:Nstates)
print *,'toto'
print *,E_guess(:)


Nguess_complex=N_energy!40
allocate(E_guess_complex(1:Nguess_complex))

print *,pi

do ll=1,Nguess_complex
    E_guess_complex(ll)%re=-23.5d0/hartree-28.0d0/hartree*dcos(dble((ll-1))/Nguess_complex*2*pi)
    E_guess_complex(ll)%im=28.0d0/hartree*dsin(dble((ll-1))/Nguess_complex*2*pi)
    ! E_guess_complex(ll)%re=E_guess(ll)
    ! E_guess_complex(ll)%im=0.2
enddo

do ll=1,Nguess_complex
    E_guess_complex(ll)%re=-23.5d0/hartree-28.0d0/hartree*dcos(dble((ll-1))/Nguess_complex*2*pi)
    E_guess_complex(ll)%im=28.0d0/hartree*dsin(dble((ll-1))/Nguess_complex*2*pi)
enddo

allocate(E_Sigma_complex(1:size(E_guess_complex)))


print *,'================'
do i=1,size(E_guess_complex)
    print *,E_guess_complex(i)
enddo
print *,'================'


! stop




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  compute energy integral gauus nodes  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


print *,'*** compute imaginary axis energy integral gauus nodes ***'

eta_omega1=0.001d0!0.000000000000001d0
eta_omega2=0.001d0!0.000001d0
eta_time =0.01d0


call gaussian_integral_1D



Nquadrature1D=9!19!15!9
Nquadrature1D_range=1

allocate(gpoint_1D(1:Nquadrature1D))
allocate(gweight_1D(1:Nquadrature1D))

gpoint_1D(:)=gpoint9_1D(:)
gweight_1D(:)=gweight9_1D(:)

allocate(omega_range(1:Nquadrature1D,1:Nquadrature1D_range))
allocate(omega_min(1:Nquadrature1D_range))
allocate(omega_max(1:Nquadrature1D_range))

omega_min(1)=0.0d0
omega_max(1)=1.0d3!5.0d4

allocate(Xi_range(1:Nquadrature1D,1:Nquadrature1D_range))
allocate(Xi_min(1:Nquadrature1D_range))
allocate(Xi_max(1:Nquadrature1D_range))

Xi_min(1)=0.0d0
Xi_max(1)=1.0d0

alpha_gw=0.2d0!0.4d0!0.2d0

do ii=1,Nquadrature1D_range
    do i=1,Nquadrature1D
        ! omega_range(i,ii)=(omega_max(ii)-omega_min(ii))/2.0d0*gpoint_1D(i)+(omega_max(ii)+omega_min(ii))/2.0d0
        ! print *,omega_range(i,ii)
        Xi_range(i,ii)=(Xi_max(ii)-Xi_min(ii))/2.0d0*gpoint_1D(i)+(Xi_max(ii)+Xi_min(ii))/2.0d0
        ! omega_range(i,ii)=Xi_range(i,ii)/(1.0d0-Xi_range(i,ii))
        ! omega_range(i,ii)=tan(pi/2.0d0*Xi_range(i,ii))*gamma
        omega_range(i,ii)=exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))-1.0d0
        print *,Xi_range(i,ii),omega_range(i,ii)
    enddo
enddo


print *,"+++++++++++++++++"





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




allocate(IdentityMatrix(1:Nn,1:Nn))
IdentityMatrix(:,:)=cmplx(0.0d0,0.0d0,8)
do i=1,Nn
    IdentityMatrix(i,i)=cmplx(1.0d0,0.0d0,8)
enddo



allocate(IdentityMatrix_real(1:Nn,1:Nn))
IdentityMatrix_real(:,:)=0.0d0
do i=1,Nn
    IdentityMatrix_real(i,i)=1.0d0
enddo




if ((cd_method0=='gs').or.(cd_method0=='spectrum')) then

    print *,'--------------------'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnzA=neigh_size/2+Nn ! slightly overestimated
    allocate(uisa(Nn+1))
    allocate(ujsa(nnzA))
    allocate(usa(nnzA))
    allocate(usb(nnzA))
    call dcsr_convert_upper(Nn,UPLO,H,IA,JA,usa,uisa,ujsa)
    call dcsr_convert_upper(Nn,UPLO,B,IA,JA,usb,uisa,ujsa)
    nnza=uisa(Nn+1)-1 ! by definition


    allocate(usa_dft(nnzA))
    call dcsr_convert_upper(Nn,UPLO,H_dft,IA,JA,usa_dft,uisa,ujsa)
    allocate(psaz_dft(nnza))



    allocate(bwork1_vectemp1(1:Nn,1:1))
    allocate(bwork1_vectemp2(1:Ne*Nquadrature,1:1))
    allocate(bwork1_vectemp3(1:Nn,1:1))


    allocate(gw_W_complex(1:Nn,1:Nn,Nquadrature1D*Nquadrature1D_range))

    

    


    ! allocate(SigmaY(1:Nn,1:1),Y_primed(1:Nn,1:1))


    print *,'*** compute dynamically screened Coulomb interaction - W ***'


    call compute_gw_W_CD(Nn,Nstates,Nempty,E(1:Nstates),omega_range(:,:),Nquadrature1D_range,Nquadrature1D,nnza,Ne,Nquadrature&
    ,E_guess_complex(1),meanfield,cd_Hadamard0)


 
    ! allocate(quaternion(4*Nn))
    ! allocate(quaternion_A(16*neigh_size))
    ! allocate(quaternion_IA(4*Nn+1))
    ! allocate(quaternion_JA(16*neigh_size))

    ! call build_quaternionmatrix(Nn)

    ! ii=0
    ! do i=1,4*Nn
    !     quaternion_IA(i)=ii+1
    !     do l=1,quaternion(i)%size
    !         ii=ii+1
    !         quaternion_JA(ii)=quaternion(i)%get(l)
    !     enddo
    ! enddo

    ! quaternion_IA(4*Nn+1)=16*neigh_size+1


    allocate(IdenMat_quaternion(4*Nn,Nn))

    IdenMat_quaternion=0.0d0
    do i=1,Nn
        IdenMat_quaternion((i-1)*4+1,i)=1.0d0
        ! IdenMat_quaternion((i-1)*4+3,i)=1.0d0!sqrt(1.0d0/3.0d0)
        ! IdenMat_quaternion((i-1)*4+4,i)=sqrt(1.0d0/3.0d0)
    enddo

    allocate(G_quaternion(4*Nn,Nn))

    ! allocate(NNmatrix_temp02(Nn,Nn))
    ! allocate(NnNgtemp(Nn,Ne*Nquadrature))



    if (cd_method0=='spectrum') then


    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! spectrum method !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (cd_Hadamard0=='no') then

    allocate(NNmatrix_temp005(Nn,Nn))
    allocate(NNmatrix_temp006(Nn,Nn))

    allocate(gw_Sigma_complex(1:Nn,1:Nn))
    allocate(gw_Sigma_complex_poles(1:Nn,1:Nn))


    do m=1,Nguess_complex
    
call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,E_guess(m)*j_real,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
        ,1,psi(:,Nstates)*j_real,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)
    
call compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,E_guess(m)*j_real,E,M00,1,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)



call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,dble(gw_Sigma_complex+gw_Sigma_complex_poles)*j_real,Nn&
            ,ZZERO,NNmatrix_temp005,Nn)

call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,transpose(NNmatrix_temp005),Nn&
            ,ZZERO,NNmatrix_temp006,Nn)

        NNmatrix_temp005=E_guess(m)*j_real*S_dense-H_dense-Vhfx_final_GW-NNmatrix_temp006

        NNmatrix_temp006=dense_matrix_inverse(Nn,NNmatrix_temp005)

        ! do k=1,5
        !     print *,dble(NNmatrix_temp006(k,k))
        ! enddo

        vtemp=0.0d0
        do k=1,Nn
            vtemp=vtemp+abs(dble(NNmatrix_temp006(k,k)))
            ! vtemp=vtemp+dble(NNmatrix_temp006(k,k))
        enddo

        print *,E_guess(m)*hartree,vtemp

    enddo


    deallocate(NNmatrix_temp005)
    deallocate(NNmatrix_temp006)


        else if (cd_Hadamard0=='yes') then


    allocate(NNmatrix_temp02(Nn,Nn))
    allocate(NNmatrix_temp005(Nn,Nn))
    allocate(NNmatrix_temp006(Nn,Nn))

    allocate(NnNgtemp(Nn,Ne*Nquadrature))

    allocate(gw_Sigma_g(1:Ne*Nquadrature,1:Ne*Nquadrature))
    allocate(gw_Sigma_g_poles(1:Ne*Nquadrature,1:Ne*Nquadrature))


    do m=1,Nguess_complex
    
call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,E_guess(m)*j_real,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
        ,1,psi(:,Nstates)*j_real,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)
    
call compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,E_guess(m)*j_real,E,M00,1,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)


call outer_product('ge',volumegweight(:),volumegweight(:),NgNgtemp)
    gw_Sigma_g=(gw_Sigma_g+gw_Sigma_g_poles)*NgNgtemp!outer_product(volumegweight(:),volumegweight(:))

    call mkl_dcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ,dble(gw_Sigma_g),Ne*Nquadrature,0.0d0,NnNgtemp,Nn)


    call mkl_dcsrmm('T',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ,transpose(NnNgtemp),Ne*Nquadrature,0.0d0,NNmatrix_temp02,Nn)


        NNmatrix_temp005=E_guess(m)*j_real*S_dense-H_dense-Vhfx_final_GW-NNmatrix_temp02

        NNmatrix_temp006=dense_matrix_inverse(Nn,NNmatrix_temp005)

        ! do k=1,5
        !     print *,dble(NNmatrix_temp006(k,k))
        ! enddo

        vtemp=0.0d0
        do k=1,Nn
            vtemp=vtemp+abs(dble(NNmatrix_temp006(k,k)))
        enddo

        print *,E_guess(m)*hartree,vtemp


    enddo

    deallocate(NNmatrix_temp02)
    deallocate(NNmatrix_temp005)
    deallocate(NNmatrix_temp006)

    deallocate(NnNgtemp)

        end if  !!! hadamard



    else if (cd_method0=='gs') then

        if (cd_Hadamard0=='no') then

    allocate(gw_Sigma_complex(1:Nn,1:Nn))
    allocate(gw_Sigma_complex_poles(1:Nn,1:Nn))

    do m=1,Nguess_complex

call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,E_guess(m)*j_real,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
        ,1,psi(:,Nstates)*j_real,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)
    
call compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,E_guess(m)*j_real,E,M00,1,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! approximated Hadamard product !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,psi(:,cd_orbital0)*j_real,Nn&
    ,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)
    call ZSYMM('L','L',Nn,1,ZONE,dble(gw_Sigma_complex+gw_Sigma_complex_poles)*j_real&
    ,Nn,bwork1_vectemp1(:,1),Nn,ZZERO,bwork1_vectemp3(:,1),Nn)
    call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,bwork1_vectemp3(:,1),Nn&
    ,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

    print *,dot_product(psi(:,cd_orbital0)*j_real,bwork1_vectemp1(:,1))

    enddo

        else if (cd_Hadamard0=='yes') then

    allocate(gw_Sigma_g(1:Ne*Nquadrature,1:Ne*Nquadrature))
    allocate(gw_Sigma_g_poles(1:Ne*Nquadrature,1:Ne*Nquadrature))


    do m=1,Nguess_complex

call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,E_guess(m)*j_real,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
        ,1,psi(:,Nstates)*j_real,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)
    
call compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,E_guess(m)*j_real,E,M00,1,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Hadamard product !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call ZSYMM('L','L',Ne*Nquadrature,1,ZONE,dble(gw_Sigma_g+gw_Sigma_g_poles)*j_real&
    ,Ne*Nquadrature,psi_point_g(:,cd_orbital0)*volumegweight(:)*j_real,Ne*Nquadrature,ZZERO,bwork1_vectemp2(:,1),Ne*Nquadrature)



    print *,dot_product(psi_point_g(:,cd_orbital0)*volumegweight(:)*j_real,bwork1_vectemp2(:,1))

    enddo

        end if !!! hadamard


    

    end if









































else if (cd_method0=='nlfeast') then





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! CALL NFEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!! easier case
Emid_nl=cmplx(E(1),0.0d0,8)
! r_nl=abs(E(1))+3.0d0/hartree!abs(E_GW(1))+3.0d0/hartree!19.0d0/hartree!51.5d0/hartree!9.0d0/hartree!19.5d0/hartree!503.0d0/hartree!28.0d0/hartree
r_nl=abs(E(1))+5.0d0/hartree!+3.0d0/hartree
M0_nl=15!M0!15!15!20


!!!!!!!!!!!!  FEAST
! fpm(30)=121118 ! code name for direct call to this routine   
  call feastinit(fpm)
  fpm(1)=1 ! change from default (printing info on screen)
  fpm(8)=8 ! number of contour points
  fpm(15)=1 !! 1 needed for synev
  fpm(16)=1 !! Trapezoida
  fpm(42)=0 ! double (0) or single (1) precision

  fpm(45)=2 ! precision bicgstab
  fpm(46)=1000 ! number ofiterations

  fpm(5)=1

  fpm(18)=50

  fpm(3)=8

 ! fpm(4)=0
  
!! feast output  
  allocate(E_nl(1:M0_nl))     ! Eigenvalue
  allocate(X_nl(1:Nn,1:M0_nl)) ! Eigenvectors
  allocate(res_nl(1:M0_nl))   ! Residual

!   X_nl(1:Nn,1:M0_nl)=psi(1:Nn,1:M0_nl)
!   E_nl(1)=cmplx(E(1),0.0d0,8)!cmplx(-23.87d0/hartree,0.0d0,8)
!   do i=2,M0_nl
!     E_nl(i)=cmplx(E(i),0.0d0,8)!(0.0d0,0.0d0)!cmplx(2.6d0/hartree+(i-1)*2.0d-2,0.0d0,8)
!   enddo


  X_nl=ZZERO
  E_nl=ZZERO

    if (M0_nl>M0+5) then

  X_nl(1:Nn,1:M0)=psi(1:Nn,1:M0)
  E_nl(1)=cmplx(E(1),0.0d0,8)!cmplx(-23.87d0/hartree,0.0d0,8)
  do i=2,M0
    E_nl(i)=cmplx(E(i),0.0d0,8)!(0.0d0,0.0d0)!cmplx(2.6d0/hartree+(i-1)*2.0d-2,0.0d0,8)
  enddo

    else

  X_nl(1:Nn,1:M0_nl)=psi(1:Nn,1:M0_nl)
  E_nl(1)=cmplx(E(1),0.0d0,8)!cmplx(-23.87d0/hartree,0.0d0,8)
  do i=2,M0_nl
    E_nl(i)=cmplx(E(i),0.0d0,8)!(0.0d0,0.0d0)!cmplx(2.6d0/hartree+(i-1)*2.0d-2,0.0d0,8)
  enddo

    end if
  
 
!!!  store factorization
 nfact=1
 if (fpm(10)==1) nfact=fpm(8)

!  print *,nfact

 allocate(ipivloc_nl(Nn,nfact))

!  allocate(Az_nl(Nn,Nn,nfact))
 allocate(Aztemp_nl(Nn,Nn))

!! RCI work-array
allocate(zAq_nl(M0_nl,M0_nl))
allocate(zBq_nl(M0_nl,M0_nl))
allocate(work_nl(Nn,M0_nl))  
allocate(zwork_nl(Nn,2*M0_nl))



nnzA=neigh_size/2+Nn ! slightly overestimated
allocate(uisa(Nn+1))
allocate(ujsa(nnzA))
allocate(usa(nnzA))
allocate(usb(nnzA))
call dcsr_convert_upper(Nn,'F',H,IA,JA,usa,uisa,ujsa)
call dcsr_convert_upper(Nn,'F',B,IA,JA,usb,uisa,ujsa)
nnza=uisa(Nn+1)-1 ! by definition





if (meanfield=='dft') then
allocate(usa_dft(nnzA))
call dcsr_convert_upper(Nn,'F',H_dft,IA,JA,usa_dft,uisa,ujsa)
allocate(psaz_dft(nnza))
end if


!!!! mixed precision set-up (only double for now)
allocate(saz(nnza,nfact))
allocate(psaz(nnza,nfact))

allocate(saztemp(nnza))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!  use same factorization for (normal+transpose solve)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MAXFCT=1 
MNUM=1
MTYPE=6      ! complex and symmetric
allocate(PT(1:64,1:nfact))
allocate(IPARM(1:64,1:nfact))
do i=1,nfact !! multiple factorization
   pt(:,i)=0
   call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
end do

!IPARM(2,:)=3 ! parallel omp nested dissection (sensible to #threads- no consistency- same runs with different results)
!IPARM(4)=11 !CGS solver
IPARM(25,:)=1 ! parallel omp rhs solve
IPARM(11,:)=0 ! disable scaling (taking care by feast fpm(41))
!!!!!!!!!!!!
if (fpm(64)==1) then
   do i=1,64
      if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
   enddo
endif
!!!!!!!!!!!!
IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used  !!
MSGLVL=0!0 !0- no output, 1- output

! !IPARM(2,:)=3 ! parallel omp nested dissection (sensible to #threads- no consistency- same runs with different results)
! !IPARM(4)=11 !CGS solver
! aIPARM_V(25)=1 ! parallel omp rhs solve
! aIPARM_V(11)=0 ! disable scaling (taking care by feast fpm(41))
! !!!!!!!!!!!!
! if (fpm(64)==1) then
!    do i=1,64
!       if (fpm(64+i)/=-111) aiparm_V(i)=fpm(64+i)
!    enddo
! endif
! !!!!!!!!!!!!
! aIPARM_V(6)=1 ! solution and rhs are input/output, attention zaux is always used  !!
! MSGLVL=0!0 !0- no output, 1- output


allocate(zaux(Nn,M0_nl))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! for bicgstab
  allocate(nres(M0_nl))
  !allocate(norm(M0))



comb=.false.!.true.!.false.
linloops=10000
com=.false. !default
if (fpm(1)/=0) com=.true.
if (fpm(1)<0) then
    fout=abs(fpm(1))+200!fpm(60) !!file id name
else
    fout=6 !screen
endif

!com=.false.

  ! bicgstab-rci     
allocate(bwork1(Nn,7*M0_nl))
allocate(bwork2(M0_nl,4)) 


allocate(bwork1_vectemp1(1:Nn,1:M0_nl))


allocate(bwork1_vectemp2(1:Ne*Nquadrature,1:M0_nl))
allocate(bwork1_vectemp3(1:Nn,1:M0_nl))
allocate(bwork1_vectemp4(1:Ne*Nquadrature,1:M0_nl))
! allocate(bwork1_vectemp1_g(1:Ne*Nquadrature,1:5))
! allocate(bwork1_vectemp2_g(1:Nstates*Nempty,1:M0_nl))
! allocate(bwork1_vectemp3_g(1:Ne*Nquadrature,1:5))





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  





ijob=-1 ! initialization 
do while (ijob/=0)
call zfeast_srcinev(ijob,Nn,Ze_nl,work_nl,zwork_nl,zAq_nl,zBq_nl,fpm,epsout,loop,Emid_nl,r_nl,M0_nl,E_nl,X_nl,M_nl,res_nl,info)
  
   
   select case(ijob)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case(10) !! !! form T(Ze) and factorize preconditioner if any-- dongming does T(Ze)=ZeS-H-Sigma(ze) == Az, or Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

    ! print *,'------------',Ze_nl
    
    ! print *,'case 10'
    
     id=fpm(33) !! id of factorization (for fpm(10) flag) 
     
    !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

        saz(1:nnza,id)=Ze_nl*usb(1:nnza)-usa(1:nnza)
            
     !!!! Factorize (preconditioner) 
            PHASE=12 !include the symbolic factorization
            psaz(1:nnza,id)=saz(1:nnza,id) ! save preconditioner
call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1,bwork1,infoloc3)  
    ! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1,bwork1,infoloc3)  



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(11) !!solve the linear system  T(Ze)x=workc(1:N,1:fpm(23)) result in to workc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! print *,'case 11'

            id=fpm(33)
    
           !!!! direct solver 
            !PHASE=33 ! solve
           !call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
    !!!! iterative bicgstab solver
            lintargeterror=10d0**(-fpm(45))
            linloops=fpm(46)
            nres=1.0d0 ! max rhs norm
            zaux(1:Nn,1:fpm(23))=(0.0d0,0.0d0) !! initial guess
            !call zbicgstab(fpm(44),'U','N',saz(1,id),uisa,ujsa,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! RCI bicgstab using ZeI-A as preconditoner
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      ijob2=-1
      do while (ijob2/=0)
         call zbicgstab_rci(ijob2,Nn,fpm(23),zwork_nl,zaux,bwork1,bwork2,jj1,jj2,nres,linloops,lintargeterror,comb,infoloc2) 
     
         select case(ijob2)
              case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
            !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input
    
          !!!! direct solver 
                 PHASE=33 ! solve
                IPARM(6,:)=0 ! solution and rhs are separated
        call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1(1,jj2)&
        ,bwork1(1,jj1),infoloc3)
        !                    
        !         aIPARM_V(6)=0 ! solution and rhs are separated                                                                             
        ! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1(1,jj2)&
                                                                                                ! ,bwork1(1,jj1),infoloc3)
          
           
      case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)

        call wzcsrmm('U','N',Nn,Nn,fpm(23),(1.0d0,0.0d0),saz(1,id),uisa,ujsa,bwork1(1,jj1),(0.0d0,0.0d0),bwork1(:,jj2))

        ! call mkl_zcsrmm('N',Nn,fpm(23),Nn,(1.0d0,0.0d0),matdescrb,saz(1,id),ujsa,uisa,uisa(2),bwork1(1,jj1),Nn&
        ! ,(0.0d0,0.0d0),bwork1(:,jj2),Nn)

        
        !!!!!!!!!!!!!!! Sigma_x !!!!!!!!!!!!!!!!!
        
        Aztemp_nl=-Vhfx_final_GW*(1.0d0,0.0d0)
        
       call ZSYMM('L','L',Nn,fpm(23),(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)
    !    CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)

       bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)+bwork1_vectemp1

       

           
         end select
      end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
            
            
           call ZCOPY( Nn*fpm(23), zaux, 1, zwork_nl, 1 )
    
      fpm(60)=fpm(60)+linloops
            if (com) then
               write(fout,'(A)',advance='no') '  #it  '
               write(fout,'(I4)',advance='no') linloops
               write(fout,'(A)',advance='no') '; res min='
               write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
               write(fout,'(A)',advance='no') '; res max='
               write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
               write(fout,*)
            endif

! stop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    case(33) !! compute norm at Ze  ||T(Ze)||
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        ! print *,'case 33'


        Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-H_dense-Vhfx_final_GW), Nn, zwork_nl )
        ! Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-Vhfx_final_GW), Nn, zwork_nl )


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(30) !! perform multiplication
        
        !!!!
        do ii=1,fpm(23)
          !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

            ! print 100,ii,fpm(23),E_nl(ii),res(ii)
            print *,ii,fpm(23),E_nl(ii)

            


            saztemp(1:nnza)=E_nl(ii)*usb(1:nnza)-usa(1:nnza)

            call wzcsrmm('U','N',Nn,Nn,1,(-1.0d0,0.0d0),saztemp,uisa,ujsa,X_nl(1,ii),(0.0d0,0.0d0),bwork1_vectemp1(:,1))

            ! call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,saztemp,ujsa,uisa,uisa(2),X_nl(1,ii),Nn&
            ! ,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

            work_nl(:,ii)=bwork1_vectemp1(:,1)

            Aztemp_nl=-Vhfx_final_GW
        
            call ZSYMM('L','L',Nn,1,(-1.0d0,0.0d0),Aztemp_nl(1,1),Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)
            ! CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl,Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

            work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)


            ! Aztemp_nl=E_nl(ii)*S_dense-H_dense-Vhfx_final_GW
        
            ! call ZSYMM('L','L',Nn,1,(-1.0d0,0.0d0),Aztemp_nl(1,1),Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

            ! work_nl(:,ii)=bwork1_vectemp1(:,1)


            

        end do

    end select

end do







! do ii=1,M0_nl
!     E_nl_temp(ii)=E_nl(ii)
!     X_nl_temp(:,ii)=X_nl(:,ii)
!     print *,E_nl(ii)*hartree,res_nl(ii)
! enddo

! print *,'======================'

! do ii=1,M_nl
!     print *,E_nl(ii)*hartree,res_nl(ii)
! enddo

print *,'+++++++++++++++++++++++++++++++++++'


call quicksort_eigenpairs(E_nl,X_nl,res_nl,Nn,1,M0_nl)










allocate(E_nl_temp(1:M0_nl))     ! Eigenvalue
allocate(X_nl_temp(1:Nn,1:M0_nl)) ! Eigenvectors
allocate(res_nl_temp(1:M0_nl))  


print *,'***found eigenvalues with Sigma_x only***'

do ii=1,M0_nl
    E_nl_temp(ii)=E_nl(ii)
    X_nl_temp(:,ii)=X_nl(:,ii)
    res_nl_temp(ii)=res_nl(ii)
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo

print *,'+++++++++++++++++++++++++++++++++++'

! print *,'======================'

! do ii=1,M_nl
!     print *,E_nl(ii)*hartree,res_nl(ii)
! enddo

! stop


do i=1,M0_nl

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    ! print *,dot_product(dble(X_nl(:,i)),snq1)

    vtemp=dot_product(dble(X_nl(:,i)),snq1)

    X_nl(:,i)=X_nl(:,i)/sqrt(vtemp)

    ! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    ! print *,dot_product(dble(X_nl(:,i)),snq1)


    ! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    ! ! print *,dot_product(dble(X_nl_temp(:,i)),snq1)

    ! vtemp=dot_product(dble(X_nl_temp(:,i)),snq1)

    ! X_nl_temp(:,i)=X_nl_temp(:,i)/sqrt(vtemp)

    ! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    ! print *,dot_product(dble(X_nl_temp(:,i)),snq1)

enddo


do k=1,Nstates+3
write (file_id, '(I0)') k
! open (unit=20,file='CO.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='CO_0.6415.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='He.1_p2_nessie.psi'//trim(file_id),status='old')
! open(unit=20,file=trim(name)//'_p2_dft.psi'//trim(file_id),status='replace')
open(unit=20,file=trim(name)//'_'//trim(dgr_fem)//'_'//trim(meanfield)//'_2.psi'//trim(file_id),status='replace')
! open (unit=20,file='LiH.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
do i=1,Nn
    write(20,*) dble(X_nl(i,k))
enddo
close (20)
enddo


! stop

! allocate(ztemp(Nn))

! call mkl_zcsrgemv('N',Nn,B*j_real,IA,JA,X_nl(:,Nstates),ztemp)
! print *,dot_product(X_nl(:,Nstates),ztemp)

! call mkl_zcsrgemv('N',Nn,B*j_real,IA,JA,X_nl(:,Nstates+1),ztemp)
! print *,dot_product(X_nl(:,Nstates+1),ztemp)

! stop








! go to 9191
























































! allocate(NNmatrix_temp(1:Nn,1:Nn))
! allocate(NNmatrix_temp005(1:Nn,1:Nn))
! ! allocate(NNmatrix_temp02(1:Nn,1:Nn))
! allocate(NNmatrix_temp03(1:Nn,1:Nn))

! allocate(gw_W_complex(1:Nn,1:Nn,Nquadrature1D*Nquadrature1D_range))
! allocate(G(1:Nn,1:Nn))
! allocate(gw_Sigma_complex(1:Nn,1:Nn))
! allocate(gw_Sigma_complex_poles(1:Nn,1:Nn))


! allocate(SigmaY(1:Nn,1:1),Y_primed(1:Nn,1:1))


! call compute_gw_W_CD(Nn,Nstates,Nempty,E(1:Nstates),omega_range(:,:),Nquadrature1D_range,Nquadrature1D,nnza,Ne,Nquadrature&
! ,E_guess_complex(1))


    ! allocate(NgNntemp(Ne*Nquadrature,Nn))
    ! allocate(NnNgtemp(Nn,Ne*Nquadrature))
    ! allocate(NnNgtemp2(Nn,Ne*Nquadrature))
    ! allocate(NgNgtemp(Ne*Nquadrature,Ne*Nquadrature))


    ! allocate(NNmatrix_temp(1:Nn,1:Nn))
    ! allocate(NNmatrix_temp005(1:Nn,1:Nn))
    ! allocate(NNmatrix_temp006(1:Nn,1:Nn))
    ! allocate(zNgNntemp(Ne*Nquadrature,Nn))
    ! allocate(zNnNgtemp(Nn,Ne*Nquadrature))
    ! allocate(zNnNgtemp2(Nn,Ne*Nquadrature))

    ! allocate(zNgNgtemp(Ne*Nquadrature,Ne*Nquadrature))
    ! allocate(zNgNgtemp2(Ne*Nquadrature,Ne*Nquadrature))
    
    
    ! allocate(NgNgtemp0(Ne*Nquadrature,Ne*Nquadrature))

    
    

    ! allocate(vectortemp(Ne*Nquadrature))

    ! allocate(chi_matrix(Nn,Nn))
    ! allocate(chi_matrix_real(Nn,Nn))

    allocate(gw_W_complex(1:Nn,1:Nn,Nquadrature1D*Nquadrature1D_range))
    ! allocate(G(1:Nn,1:Nn))
    ! allocate(gw_G_g(1:Ne*Nquadrature,1:Ne*Nquadrature))
    ! allocate(gw_W_g(1:Ne*Nquadrature,1:Ne*Nquadrature))


    if (cd_Hadamard0=='no') then

    allocate(gw_Sigma_complex(1:Nn,1:Nn))
    allocate(gw_Sigma_complex_poles(1:Nn,1:Nn))

    else if (cd_Hadamard0=='yes') then

    allocate(gw_Sigma_g(1:Ne*Nquadrature,1:Ne*Nquadrature))
    allocate(gw_Sigma_g_poles(1:Ne*Nquadrature,1:Ne*Nquadrature))

    ! allocate(NNmatrix_temp02(Nn,Nn))
    allocate(NNmatrix_temp005(Nn,Nn))
    allocate(NNmatrix_temp006(Nn,Nn))

    allocate(NnNgtemp(Nn,Ne*Nquadrature))

    end if


    ! allocate(SigmaY(1:Nn,1:1),Y_primed(1:Nn,1:1))


    call compute_gw_W_CD(Nn,Nstates,Nempty,E(1:Nstates),omega_range(:,:),Nquadrature1D_range,Nquadrature1D,nnza,Ne,Nquadrature&
    ,E_guess_complex(1),meanfield,cd_Hadamard0)

    ! deallocate(NNmatrix_temp02)
    ! deallocate(NNmatrix_temp03)
    ! deallocate(NNmatrix_temp04)
    ! deallocate(NgNntemp)
    ! deallocate(NnNgtemp)
    ! deallocate(NnNgtemp2)
    ! deallocate(NgNgtemp)



allocate(quaternion(4*Nn))
allocate(quaternion_A(16*neigh_size))
allocate(quaternion_IA(4*Nn+1))
allocate(quaternion_JA(16*neigh_size))

call build_quaternionmatrix(Nn)

ii=0
do i=1,4*Nn
    quaternion_IA(i)=ii+1
    do l=1,quaternion(i)%size
        ii=ii+1
        quaternion_JA(ii)=quaternion(i)%get(l)
    enddo
enddo

quaternion_IA(4*Nn+1)=16*neigh_size+1

allocate(IdenMat_quaternion(4*Nn,Nn))

IdenMat_quaternion=0.0d0
do i=1,Nn
    IdenMat_quaternion((i-1)*4+1,i)=1.0d0
    ! IdenMat_quaternion((i-1)*4+3,i)=1.0d0!sqrt(1.0d0/3.0d0)
    ! IdenMat_quaternion((i-1)*4+4,i)=sqrt(1.0d0/3.0d0)
enddo

allocate(G_quaternion(4*Nn,Nn))








!!!!!!!!!!!!  FEAST
! fpm(30)=121118 ! code name for direct call to this routine   
call feastinit(fpm)
fpm(1)=cd_fpm10!1 ! change from default (printing info on screen)
fpm(8)=cd_fpm80!8 ! number of contour points
fpm(15)=1 !! 1 needed for synev
fpm(16)=1 !! Trapezoida
fpm(42)=0 ! double (0) or single (1) precision

fpm(45)=cd_fpm450!2!8!2 ! precision bicgstab
fpm(46)=1000 ! number ofiterations

fpm(5)=1

! fpm(18)=50

! fpm(3)=7

! 3333 continue

! call compute_gw_W_CD(Nn,Nstates,Nempty,E(1:Nstates),omega_range(:,:),Nquadrature1D_range,Nquadrature1D,nnza,Ne,Nquadrature,Ze_nl)

! easier case
Emid_nl=cd_Emid0/hartree!cmplx(dble(E_nl(1)),0.0d0,8)!cmplx(-10.5d0/hartree,0.0d0,8)!cmplx(dble(E_nl(1)),0.0d0,8)
r_nl=cd_r0/hartree!abs(dble(E_nl(1)))+2.9d0/hartree!2.5d0/hartree!13.5d0/hartree!abs(dble(E_nl(1)))+2.5d0/hartree
! r_nl=abs(dble(E_nl(1)))/hartree!4.5d0/hartree!abs(dble(E_nl(1)))+2.0d0/hartree  !!!! Emax is slightly higher than LUMO+1 states 
M0_nl=cd_M00!3!4!5!6!Nstates+3!4!11

! Emid_nl=cmplx(dble(E_nl(2)),0.0d0,8)
! ! r_nl=18.4d0/hartree!51.5d0/hartree!27.0d0/hartree
! r_nl=abs(dble(E_nl(2)))+3.0d0/hartree!E_nl(Nstates+2)!+(E_nl(Nstates+3)-E_nl(Nstates+2))/2.0d0  !!!! Emax is slightly higher than LUMO+1 states 
! M0_nl=5!4!11


fpm(3)=cd_fpm30!12!9!10!8!10!8!7
! fpm(6)=1  ! 0: Using relative error on the trace epsout i.e. epsout / 1: Using relative residual res i.e. maxi res(i) 
! fpm(45)=2 ! precision bicgstab
fpm(16)=1 !!!!! 0: Gauss 1: Trapezoidal; 2: Zolotarev








!! feast output  

deallocate(E_nl,X_nl,res_nl)
allocate(E_nl(1:M0_nl))     ! Eigenvalue
allocate(X_nl(1:Nn,1:M0_nl)) ! Eigenvectors
!   allocate(E_nl_temp(1:M0_nl))     ! Eigenvalue
!   allocate(X_nl_temp(1:Nn,1:M0_nl)) ! Eigenvectors
allocate(res_nl(1:M0_nl))   ! Residual

! print *,'toto7'

X_nl(1:Nn,1:M0_nl)=X_nl_temp(1:Nn,cd_igorbital0:M0_nl+(cd_igorbital0-1))
E_nl(1:M0_nl)=E_nl_temp(cd_igorbital0:M0_nl+(cd_igorbital0-1))
! res_nl(1:M0_nl)=res_nl_temp(1:M0_nl)


! !!!  store factorization
! nfact=1
! if (fpm(10)==1) nfact=fpm(8)

! print *,nfact

! allocate(ipivloc_nl(Nn,nfact))

! allocate(Az_nl(Nn,Nn,nfact))
! allocate(Aztemp_nl(Nn,Nn))

! print *,Nn,M0_nl

!! RCI work-array


! deallocate(zAq_nl,zBq_nl,work_nl,zwork_nl)
! allocate(zAq_nl(M0_nl,M0_nl))
! allocate(zBq_nl(M0_nl,M0_nl))
! allocate(work_nl(Nn,M0_nl)) 
! allocate(zwork_nl(Nn,2*M0_nl))




! nnzA=neigh_size/2+Nn ! slightly overestimated
! allocate(uisa(Nn+1))
! allocate(ujsa(nnzA))
! allocate(usa(nnzA))
! allocate(usb(nnzA))
! call dcsr_convert_upper(Nn,UPLO,H,IA,JA,usa,uisa,ujsa)
! call dcsr_convert_upper(Nn,UPLO,B,IA,JA,usb,uisa,ujsa)
! nnza=uisa(Nn+1)-1 ! by definition


! allocate(usa_dft(nnzA))
! call dcsr_convert_upper(Nn,UPLO,H_dft,IA,JA,usa_dft,uisa,ujsa)
! allocate(psaz_dft(nnza))


!!!! mixed precision set-up (only double for now)
! allocate(saz(nnza,nfact))
! allocate(psaz(nnza,nfact))

! allocate(saztemp(nnza))




! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!  pardiso initialization
! !!!!!!!!  use same factorization for (normal+transpose solve)  
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAXFCT=1 
! MNUM=1
! MTYPE=6      ! complex and symmetric
! allocate(IPARM(1:64,1:nfact), STAT = k)
! print *,k
! allocate(PT(1:64,1:nfact))
! do i=1,nfact !! multiple factorization
!  pt(:,i)=0
!  call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
! end do

! ! allocate(apt_v(64))
! ! allocate(aiparm_v(64))
! ! apt_V=0
! ! call pardisoinit(aPT_V,MTYPE,aIPARM_V)


! !IPARM(2,:)=3 ! parallel omp nested dissection (sensible to #threads- no consistency- same runs with different results)
! !IPARM(4)=11 !CGS solver
! IPARM(25,:)=1 ! parallel omp rhs solve
! IPARM(11,:)=0 ! disable scaling (taking care by feast fpm(41))
! !!!!!!!!!!!!
! if (fpm(64)==1) then
!  do i=1,64
!     if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
!  enddo
! endif
! !!!!!!!!!!!!
! IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used  !!
! MSGLVL=0!0 !0- no output, 1- output

! !IPARM(2,:)=3 ! parallel omp nested dissection (sensible to #threads- no consistency- same runs with different results)
! !IPARM(4)=11 !CGS solver
! aIPARM_V(25)=1 ! parallel omp rhs solve
! aIPARM_V(11)=0 ! disable scaling (taking care by feast fpm(41))
! !!!!!!!!!!!!
! if (fpm(64)==1) then
!    do i=1,64
!       if (fpm(64+i)/=-111) aiparm_V(i)=fpm(64+i)
!    enddo
! endif
! !!!!!!!!!!!!
! aIPARM_V(6)=1 ! solution and rhs are input/output, attention zaux is always used  !!
! MSGLVL=0!0 !0- no output, 1- output



! deallocate(zaux,nres)
! allocate(zaux(Nn,M0_nl))
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! for bicgstab
! allocate(nres(M0_nl))
! !allocate(norm(M0))




! comb=.true.!.false.
! linloops=10000
! com=.false. !default
! if (fpm(1)/=0) com=.true.
! if (fpm(1)<0) then
!   fout=abs(fpm(1))+200!fpm(60) !!file id name
! else
!   fout=6 !screen
! endif

! !com=.false.

! print *,'toto'

! ! bicgstab-rci     
! deallocate(bwork1,bwork2)
! allocate(bwork1(Nn,7*M0_nl))
! allocate(bwork2(M0_nl,4)) 

! deallocate(bwork1_vectemp1)
! allocate(bwork1_vectemp1(1:Nn,1:M0_nl))



allocate(NNmatrix_temp02(Nn,Nn))
allocate(NNmatrix_temp03(Nn,Nn))

allocate(bwork1_vectemp1_g(Ne*Nquadrature,M0_nl))
allocate(bwork1_vectemp3_g(Ne*Nquadrature,M0_nl))

if (cd_Hadamard0=='no') then

allocate(Sigma_c_temp(Nn,Nn,fpm(8)))

else if (cd_Hadamard0=='yes') then

allocate(Sigma_c_temp_g(Ne*Nquadrature,Ne*Nquadrature,fpm(8)))
! allocate(Sigma_c_temp(Ne*Nquadrature,Ne*Nquadrature,fpm(8)))
allocate(NgNgtemp(1:Ne*Nquadrature,1:Ne*Nquadrature))

end if

! allocate(contourpoints(1:fpm(8)))


ii10=0
ii11=0
ii30=0


ijob=-1 ! initialization 
do while (ijob/=0)
call zfeast_srcinev(ijob,Nn,Ze_nl,work_nl,zwork_nl,zAq_nl,zBq_nl,fpm,epsout,loop,Emid_nl,r_nl,M0_nl,E_nl,X_nl,M_nl,res_nl,info)
  
   
   select case(ijob)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case(10) !! !! form T(Ze) and factorize preconditioner if any-- dongming does T(Ze)=ZeS-H-Sigma(ze) == Az, or Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    
    ! print *,'case 10'

    ! print *,Ze_nl
    
    
     id=fpm(33) !! id of factorization (for fpm(10) flag) 
     
    !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

        saz(1:nnza,id)=Ze_nl*usb(1:nnza)-usa(1:nnza)
            
     !!!! Factorize (preconditioner) 
            PHASE=12 !include the symbolic factorization
            psaz(1:nnza,id)=saz(1:nnza,id) ! save preconditioner
call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1,bwork1,infoloc3)   
! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1,bwork1,infoloc3)   



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(11) !!solve the linear system  T(Ze)x=workc(1:N,1:fpm(23)) result in to workc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! print *,'case 11'

            id=fpm(33)




!             if (ii11<fpm(8)/2) then

!             ! call compute_gw_Sigma_CD_times_Y(Nn,Nstates,Nempty,Ze_nl,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
!             ! ,M00,fpm(23),bwork1(:,jj1:jj1+fpm(23)-1),nnza)

!             ! call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,Ze_nl,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
!             ! ,M00,psi(:,HOMO_state-(Nstates-1))*J_real,nnza,Ne,Nquadrature)

!             ! call compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,Ze_nl,E,M00,fpm(23),nnza,Ne,Nquadrature)

! call compute_gw_Sigma_CD_imagintegral_2(Nn,Nstates,Nempty,Ze_nl,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
!             ,M00,psi(:,Nstates)*j_real,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)

! call compute_gw_Sigma_CD_poles_mtx_2(Nn,Nstates,Nempty,Ze_nl,E,M00,fpm(23),nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)


!             if (cd_Hadamard0=='no') then

!                 gw_Sigma_complex=gw_Sigma_complex+gw_Sigma_complex_poles
!                 ! gw_Sigma_g=gw_Sigma_g+gw_Sigma_g_poles

!             ! call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,Ze_nl,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
!             ! ,fpm(23),bwork1(:,jj1:jj1+fpm(23)-1),nnza)

!                 Sigma_c_temp(:,:,ii11+1)=gw_Sigma_complex
!                 ! Sigma_c_temp(:,:,ii11+1)=gw_Sigma_g

!             else if (cd_Hadamard0=='yes') then

!                 gw_Sigma_g=gw_Sigma_g+gw_Sigma_g_poles

!                 Sigma_c_temp_g(:,:,ii11+1)=gw_Sigma_g

!             ! contourpoints(ii11+1)=Ze_nl

!             end if

!             end if

!             if ((ii11>fpm(8)/2-1).and.(ii11<fpm(8))) then

!                 if (cd_Hadamard0=='no') then

!                     Sigma_c_temp(:,:,ii11+1)=Sigma_c_temp(:,:,fpm(8)-ii11)

!                 else if (cd_Hadamard0=='yes') then

!                     Sigma_c_temp_g(:,:,ii11+1)=Sigma_c_temp_g(:,:,fpm(8)-ii11)

!                 end if

!                 ! contourpoints(ii11+1)=Ze_nl

!             end if

            
    
           !!!! direct solver 
            !PHASE=33 ! solve
           !call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
    !!!! iterative bicgstab solver
            lintargeterror=10d0**(-fpm(45))
            linloops=fpm(46)
            nres=1.0d0 ! max rhs norm
            zaux(1:Nn,1:fpm(23))=(0.0d0,0.0d0) !! initial guess
            !call zbicgstab(fpm(44),'U','N',saz(1,id),uisa,ujsa,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! RCI bicgstab using ZeI-A as preconditoner
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      ijob2=-1
      do while (ijob2/=0)
    call zbicgstab_rci(ijob2,Nn,fpm(23),zwork_nl,zaux,bwork1,bwork2,jj1,jj2,nres,linloops,lintargeterror,comb,infoloc2) 
     
         select case(ijob2)
              case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
            !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input
    
          !!!! direct solver 
                 PHASE=33 ! solve
                IPARM(6,:)=0 ! solution and rhs are separated
        call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1(1,jj2)&
                                                                                                        ,bwork1(1,jj1),infoloc3)

        !         aIPARM_V(6)=0 ! solution and rhs are separated  
        ! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1(1,jj2)&
        !                                                                                         ,bwork1(1,jj1),infoloc3)
          
           
      case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)


        

        call wzcsrmm('U','N',Nn,Nn,fpm(23),(1.0d0,0.0d0),saz(1,id),uisa,ujsa,bwork1(1,jj1),(0.0d0,0.0d0),bwork1(:,jj2))

        
        !!!!!!!!!!!!!!! Sigma_x !!!!!!!!!!!!!!!!!
        
        Aztemp_nl=-Vhfx_final_GW*(1.0d0,0.0d0)
        
       call ZSYMM('L','L',Nn,fpm(23),(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)
        ! CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)

       bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)+bwork1_vectemp1(:,1:fpm(23))

       
!         !!!!!!!!!!!!!!! Sigma_c !!!!!!!!!!!!!!!!!

!        if (cd_Hadamard0=='no') then

!        call mkl_zcsrmm('N',Nn,fpm(23),Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,bwork1(1,jj1),Nn,(0.0d0,0.0d0)&
!        ,bwork1_vectemp1,Nn)
! call ZSYMM('L','L',Nn,fpm(23),(1.0d0,0.0d0),dble(Sigma_c_temp(:,:,mod(ii11,fpm(8))+1))*j_real,Nn,bwork1_vectemp1,Nn,(0.0d0,0.0d0)&
!        ,bwork1_vectemp3,Nn)
!        call mkl_zcsrmm('N',Nn,fpm(23),Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,bwork1_vectemp3,Nn,(0.0d0,0.0d0)&
!        ,bwork1_vectemp1,Nn)
   
!            bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)-bwork1_vectemp1


!         else if (cd_Hadamard0=='yes') then


! call mkl_zcsrmm('N',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ,bwork1(1,jj1),Nn,ZZERO,bwork1_vectemp3_g,Ne*Nquadrature)

!        do i=1,fpm(23)
!             bwork1_vectemp3_g(:,i)=bwork1_vectemp3_g(:,i)*volumegweight(:)
!        enddo

! call ZSYMM('L','L',Ne*Nquadrature,fpm(23),ZONE,dble(gw_Sigma_g+gw_Sigma_g_poles)*j_real,Ne*Nquadrature,bwork1_vectemp3_g&
! ,Ne*Nquadrature,ZZERO,bwork1_vectemp1_g,Ne*Nquadrature)

!         do i=1,fpm(23)
!             bwork1_vectemp1_g(:,i)=bwork1_vectemp1_g(:,i)*volumegweight(:)
!         enddo

! call mkl_zcsrmm('T',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!     ,bwork1_vectemp1_g,Ne*Nquadrature,ZZERO,bwork1_vectemp1,Nn)


!             bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)-bwork1_vectemp1


!         end if  !!! hadamard




!     ! call mkl_zcsrmm('N',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!     !                 ,bwork1(1,jj1),Nn,ZZERO,bwork1_vectemp2,Ne*Nquadrature)

!     !                 do i=1,fpm(23)
!     !                     bwork1_vectemp2(:,i)=bwork1_vectemp2(:,i)*volumegweight(:)
!     !                 enddo

!     ! call ZSYMM('L','L',Ne*Nquadrature,fpm(23),ZONE,dble(Sigma_c_temp(:,:,mod(ii11,fpm(8))+1))*j_real,Ne*Nquadrature&
!     !                 ,bwork1_vectemp2,Ne*Nquadrature,ZZERO,bwork1_vectemp4,Ne*Nquadrature)

!     !                 do i=1,fpm(23)
!     !                     bwork1_vectemp4(:,i)=bwork1_vectemp4(:,i)*volumegweight(:)
!     !                 enddo

!     ! call mkl_zcsrmm('T',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!     !                 ,bwork1_vectemp4,Ne*Nquadrature,ZZERO,bwork1_vectemp1,Nn)


!     !                 bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)-bwork1_vectemp1
       





           
         end select
      end do

    !   !    print *,'==================='
    !         print *,bwork1(1,jj2:jj2+fpm(23)-1)
    !   !    print *,bwork1(Nn,jj2:jj2+fpm(23)-1)
    !   !    print *,'==================='
            
            ! do kk=1,fpm(23)
            !     print *,dot_product(bwork1(:,jj1+kk-1),bwork1(:,jj2+kk-1))
            ! enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ii11=ii11+1
            
           call ZCOPY( Nn*fpm(23), zaux, 1, zwork_nl, 1 )
    
      fpm(60)=fpm(60)+linloops
            if (com) then
               write(fout,'(A)',advance='no') '  #it  '
               write(fout,'(I4)',advance='no') linloops
               write(fout,'(A)',advance='no') '; res min='
               write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
               write(fout,'(A)',advance='no') '; res max='
               write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
               write(fout,*)
            endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !     zwork_nl%im=0.0d0

    !     bwork1_vectemp4(:,mod(ii11,fpm(8))*fpm(23)+1:mod(ii11,fpm(8))*fpm(23)+fpm(23))=zwork_nl!bwork1(:,jj2:jj2+fpm(23)-1)


    !     else if ((mod(ii11,fpm(8))>fpm(8)/2-1).and.(mod(ii11,fpm(8))<fpm(8))) then


    ! zwork_nl%re = bwork1_vectemp4(:,((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+1):((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+fpm(23)))%re
    ! zwork_nl%im = 0.0d0   !-bwork1_vectemp4(:,((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+1):((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+fpm(23)))%im


    !     end if 


            

! stop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    case(33) !! compute norm at Ze  ||T(Ze)||
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        ! print *,'case 33'

        call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,Ze_nl,E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
            ,M00,psi(:,Nstates)*j_real,nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)

        call compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,Ze_nl,E,M00,fpm(23),nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)

        

        if (cd_Hadamard0=='no') then

call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,dble(gw_Sigma_complex+gw_Sigma_complex_poles),Nn&
            ,0.0d0,NNmatrix_temp02,Nn)

call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,transpose(NNmatrix_temp02),Nn&
            ,0.0d0,NNmatrix_temp03,Nn)

        Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-H_dense-Vhfx_final_GW-NNmatrix_temp03), Nn, zwork_nl )

        else if (cd_Hadamard0=='yes') then

            call outer_product('ge',volumegweight(:),volumegweight(:),NgNgtemp)
            gw_Sigma_g=(gw_Sigma_g+gw_Sigma_g_poles)*NgNgtemp!outer_product(volumegweight(:),volumegweight(:))

    call mkl_dcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ,dble(gw_Sigma_g),Ne*Nquadrature,0.0d0,NnNgtemp,Nn)


    call mkl_dcsrmm('T',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ,transpose(NnNgtemp),Ne*Nquadrature,0.0d0,NNmatrix_temp02,Nn)

        Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-H_dense-Vhfx_final_GW-NNmatrix_temp02), Nn, zwork_nl )

        end if


! call mkl_dcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                     ,dble(gw_Sigma_g+gw_Sigma_g_poles),Ne*Nquadrature,0.0d0,NnNgtemp,Nn)

! call mkl_dcsrmm('T',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                     ,transpose(NnNgtemp),Ne*Nquadrature,0.0d0,NNmatrix_temp02,Nn)

        ! Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-H_dense-Vhfx_final_GW-NNmatrix_temp03), Nn, zwork_nl )

        
        ! Ze_nl=ZLANGE('F',Nn,Nn,(Ze_nl*S_dense-H_dense-Vhfx_final_GW-dble(gw_Sigma_complex+gw_Sigma_complex_poles)),Nn,zwork_nl)

        ! Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-H_dense-Vhfx_final_GW), Nn, zwork_nl )


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(30) !! perform multiplication
        
        !!form T(lambda(i)) and multiply by x(1:N,i) result in work(1:N,i)

        ! print *,'case 30'
        
        !!!!
        do ii=1,fpm(23)
          !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

            
            ! print 100,ii,fpm(23),E_nl(ii),res(ii)
            print *,ii,fpm(23),E_nl(ii)
            ! print *,X_nl(1:3,ii)

            


            saztemp(1:nnza)=E_nl(ii)*usb(1:nnza)-usa(1:nnza)

            call wzcsrmm('U','N',Nn,Nn,1,(-1.0d0,0.0d0),saztemp,uisa,ujsa,X_nl(1,ii),(0.0d0,0.0d0),bwork1_vectemp1(:,1))

            work_nl(:,ii)=bwork1_vectemp1(:,1)

            Aztemp_nl=-Vhfx_final_GW
        
            call ZSYMM('L','L',Nn,1,(-1.0d0,0.0d0),Aztemp_nl(1,1),Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)
            ! CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl,Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

            work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)



            !!!!!!!!!!!!!!! Sigma_c !!!!!!!!!!!!!!!!!

        call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,E_nl(ii),E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
        ,1,X_nl(1,ii),nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)

        call compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,E_nl(ii),E,M00,fpm(23),nnza,Ne,Nquadrature,meanfield,cd_Hadamard0)

        if (cd_Hadamard0=='no') then

        call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,X_nl(1,ii),Nn,(0.0d0,0.0d0)&
        ,bwork1_vectemp1(:,1),Nn)
call ZSYMM('L','L',Nn,1,(1.0d0,0.0d0),dble(gw_Sigma_complex+gw_Sigma_complex_poles)*j_real,Nn,bwork1_vectemp1(:,1),Nn,(0.0d0,0.0d0)&
        ,bwork1_vectemp3(:,1),Nn)
        call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,bwork1_vectemp3(:,1),Nn,(0.0d0,0.0d0)&
        ,bwork1_vectemp1(:,1),Nn)

        else if (cd_Hadamard0=='yes') then

        call mkl_zcsrmm('N',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ,X_nl(1,ii),Nn,ZZERO,bwork1_vectemp3_g(:,1),Ne*Nquadrature)

       bwork1_vectemp3_g(:,1)=bwork1_vectemp3_g(:,1)*volumegweight(:)

call ZSYMM('L','L',Ne*Nquadrature,1,ZONE,dble(gw_Sigma_g+gw_Sigma_g_poles)*j_real,Ne*Nquadrature,bwork1_vectemp3_g(:,1)&
,Ne*Nquadrature,ZZERO,bwork1_vectemp1_g(:,1),Ne*Nquadrature)

        bwork1_vectemp1_g(:,1)=bwork1_vectemp1_g(:,1)*volumegweight(:)

call mkl_zcsrmm('T',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ,bwork1_vectemp1_g(:,1),Ne*Nquadrature,ZZERO,bwork1_vectemp1(:,1),Nn)



        end if

                work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)



    ! call mkl_zcsrmm('N',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                 ,X_nl(1,ii),Nn,ZZERO,bwork1_vectemp2(:,1),Ne*Nquadrature)

    ! call ZSYMM('L','L',Ne*Nquadrature,1,ZONE,dble(gw_Sigma_g+gw_Sigma_g_poles)*j_real,Ne*Nquadrature&
    !                 ,bwork1_vectemp2(:,1)*volumegweight(:),Ne*Nquadrature,ZZERO,bwork1_vectemp4(:,1),Ne*Nquadrature)

    ! call mkl_zcsrmm('T',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                 ,bwork1_vectemp4(:,1)*volumegweight(:),Ne*Nquadrature,ZZERO,bwork1_vectemp1(:,1),Nn)


    !                 work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)


!             call compute_gw_Sigma_CD_poles(Nn,Nstates,Nempty,E_nl(ii),E,M00,1,X_nl(1,ii),nnza,Ne,Nquadrature)

! call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,dble(bwork1_vectemp3(:,1))*j_real+Y_primed(:,1)&
! ,Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

!         ! print *,bwork1_vectemp1(1:3,1)!+Y_primed(1:3,1)

!                 work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)!Y_primed(:,1)
       



            

        end do

    end select

end do


do ii=1,M0_nl
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo

print *,'======================'

do ii=1,M_nl
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo

print *,'======================'



print *,'+++++++++++++++++++++++++++++++++++'

call quicksort_eigenpairs(E_nl,X_nl,res_nl,Nn,1,M0_nl)

do ii=1,M0_nl
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo

print *,'+++++++++++++++++++++++++++++++++++'


do i=1,M0_nl

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    print *,dot_product(dble(X_nl(:,i)),snq1)

    vtemp=dot_product(dble(X_nl(:,i)),snq1)

    X_nl(:,i)=X_nl(:,i)/sqrt(vtemp)

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    print *,dot_product(dble(X_nl(:,i)),snq1)


    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    ! print *,dot_product(dble(X_nl_temp(:,i)),snq1)

    vtemp=dot_product(dble(X_nl_temp(:,i)),snq1)

    X_nl_temp(:,i)=X_nl_temp(:,i)/sqrt(vtemp)

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    print *,dot_product(dble(X_nl_temp(:,i)),snq1)

enddo


! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,Nstates+1)),snq1)
! print *,dot_product(dble(X_nl(:,Nstates+1)),snq1)

! vtemp=dot_product(dble(X_nl(:,Nstates+1)),snq1)

! X_nl(:,Nstates+1)=X_nl(:,Nstates+1)/sqrt(vtemp)



! call mkl_dcsrgemv('N',Nn,B,IA,JA,(abs(psi(:,Nstates))-abs(dble(X_nl(:,Nstates)))),snq1)
! print *,dot_product((abs(psi(:,Nstates))-abs(dble(X_nl(:,Nstates)))),snq1)

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(abs(psi(:,Nstates+1))-abs(dble(X_nl(:,Nstates+1)))),snq1)
! print *,dot_product((abs(psi(:,Nstates+1))-abs(dble(X_nl(:,Nstates+1)))),snq1)


do i=1,M0_nl

    print *,'========================='

call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,cd_igorbital0+(i-1))+dble(X_nl(:,i))),snq1)
print *,dot_product((psi(:,cd_igorbital0+(i-1))+dble(X_nl(:,i))),snq1)
 
call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,cd_igorbital0+(i-1))-dble(X_nl(:,i))),snq1)
print *,dot_product((psi(:,cd_igorbital0+(i-1))-dble(X_nl(:,i))),snq1)

    print *,'~~~~~~~~~~~~~~~~~~~'

call mkl_dcsrgemv('N',Nn,B,IA,JA,(dble(X_nl_temp(:,cd_igorbital0+(i-1)))+dble(X_nl(:,i))),snq1)
print *,dot_product((dble(X_nl_temp(:,cd_igorbital0+(i-1)))+dble(X_nl(:,i))),snq1)

call mkl_dcsrgemv('N',Nn,B,IA,JA,(dble(X_nl_temp(:,cd_igorbital0+(i-1)))-dble(X_nl(:,i))),snq1)
print *,dot_product((dble(X_nl_temp(:,cd_igorbital0+(i-1)))-dble(X_nl(:,i))),snq1)

    ! print *,'========================='

enddo

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,Nstates+1)+dble(X_nl(:,Nstates+1))),snq1)
! print *,dot_product((psi(:,Nstates+1)+dble(X_nl(:,Nstates+1))),snq1)


! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,Nstates)-dble(X_nl(:,Nstates))),snq1)
! print *,dot_product((psi(:,Nstates)-dble(X_nl(:,Nstates))),snq1)

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,Nstates+1)-dble(X_nl(:,Nstates+1))),snq1)
! print *,dot_product((psi(:,Nstates+1)-dble(X_nl(:,Nstates+1))),snq1)




do k=1,M0_nl!Nstates+1
write (file_id, '(I0)') k
! open (unit=20,file='CO.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='CO_0.6415.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='He.1_p2_nessie.psi'//trim(file_id),status='old')
open(unit=20,file=trim(name)//'_p2_gw.psi'//trim(file_id),status='replace')
! open (unit=20,file='LiH.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
do i=1,Nn
    write(20,*) X_nl(i,k)
enddo
close (20)
enddo



! do ii=1,32
!     print *,contourpoints(ii)
! enddo

! print *,'----------------------'



! do ii=1,32

!     ! ! work_nl(:,:)=(0.0d0,0.0d0)

!     ! saztemp(1:nnza)=contourpoints(ii)*usb(1:nnza)-usa(1:nnza)

!     ! call wzcsrmm('U','N',Nn,Nn,M0_nl,(1.0d0,0.0d0),saztemp,uisa,ujsa,X_nl(1,1),(0.0d0,0.0d0),bwork1_vectemp1(:,1))

!     ! work_nl(:,:)=bwork1_vectemp1(:,:)

!     ! Aztemp_nl=-Vhfx_final_GW

!     ! call ZSYMM('L','L',Nn,M0_nl,(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,X_nl(1,1),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)
!     ! ! CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl,Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

!     ! work_nl(:,:)=work_nl(:,:)+bwork1_vectemp1(:,:)

! !!!!!!!!!!!!!!! Sigma_c !!!!!!!!!!!!!!!!!

! call compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,contourpoints(ii),E,omega_range,Xi_range,Nquadrature1D_range,Nquadrature1D&
! ,1,X_nl(1,ii),nnza,Ne,Nquadrature)

! call mkl_zcsrmm('N',Nn,M0_nl,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,X_nl(1,1),Nn,(0.0d0,0.0d0)&
! ,bwork1_vectemp1(1,1),Nn)
! call ZSYMM('L','L',Nn,M0_nl,(1.0d0,0.0d0),dble(gw_Sigma_complex)*j_real,Nn,bwork1_vectemp1(1,1),Nn,(0.0d0,0.0d0)&
! ,bwork1_vectemp3(1,1),Nn)
! ! call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,bwork1_vectemp3(:,1),Nn,(0.0d0,0.0d0)&
! ! ,bwork1_vectemp1(:,1),Nn)

! !         work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)


!     call compute_gw_Sigma_CD_poles(Nn,Nstates,Nempty,contourpoints(ii),E,M00,M0_nl,X_nl(1,1),nnza,Ne,Nquadrature)

! call mkl_zcsrmm('N',Nn,M0_nl,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
! ,dble(bwork1_vectemp3(:,:))*j_real+Y_primed(:,:)&
! ,Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,:),Nn)

! ! print *,bwork1_vectemp1(1:3,1)!+Y_primed(1:3,1)

!         ! work_nl(:,:)=work_nl(:,:)-bwork1_vectemp1(:,:)!Y_primed(:,1)
!         work_nl(:,:)=bwork1_vectemp1(:,:)!Y_primed(:,1)

!         do i=1,M0_nl

!         print *,dot_product(X_nl(:,i),work_nl(:,i))

!         enddo


! enddo



! 9191 continue



! enddo !! Nemptyloop



end if !!! end cd_method0 selection



else if (gw_sigma=='casida') then




    ! go to 9999


    if (casida_method0=='gs') then

        print *,'*** compute Casida Kx ***'

    Nempty=casida_N_empty0!1000!4000!1500!6500!1500!3000!1500!Nn-Nstates!1000!1500!3000!1500!1000!1500


allocate(psi_ia_g_real(Ne*Nquadrature,Nstates*Nempty))
do k=1,Nstates
    do kk=1,Nempty
        psi_ia_g_real(:,(k-1)*Nempty+kk)=psi_point_g(:,k)*psi_point_g(:,Nstates+kk)!*volumegweight(:)
    enddo
enddo








allocate(NgNntemp(1:Ne*Nquadrature,1:Nstates*Nempty))
! allocate(NnNgtemp(1:Nn,1:Nstates*Nempty))

allocate(casida_Kx(1:Nstates*Nempty,1:Nstates*Nempty))
allocate(casida_R(1:Nstates*Nempty))
allocate(casida_Xs(1:Nstates*Nempty,1:Nstates*Nempty))
allocate(casida_omega(1:Nstates*Nempty))
allocate(casida_R_half(1:Nstates*Nempty))
! allocate(casida_V(1:Nstates*Nempty))
allocate(casida_Kxnj(1:Nempty,1:Nstates*Nempty))
allocate(casida_Kxnj_occupied(1:Nstates,1:Nstates*Nempty))
allocate(casida_cvomega(1:Nstates*Nempty))

allocate(casida_Ys(1:Nstates*Nempty,1:Nstates*Nempty))


allocate(NNmatrix_temp12(1:Nn,1:Nstates*Nempty))
allocate(NNmatrix_temp13(1:Nn,1:Nstates*Nempty))
! allocate(NNmatrix_temp14(1:Nstates*Nempty,1:Nstates*Nempty))
! allocate(vector_temp00(1:Nn))




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   compute casida_Kx   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call compute_gw_Casida_Kx(Nn,Ne,Nquadrature,Nstates,Nempty)

print *,'*** construct Casida matrix  ***'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   solve casida_equation   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call compute_gw_Casida_matrix(Nstates,Nempty,E)

print *,'*** solve Casida eigenvalue equation  ***'

call solve_Casida_matrix(Nstates,Nempty)

! stop


! do i=1,50
!     print *,casida_omega(i)
! enddo

! stop


! do k=1,Nstates
!     do kk=1,Nempty
!         casida_cvomega((k-1)*Nempty+kk)=sqrt(E(Nstates+kk)-E(k))
!     enddo
! enddo






!!!!!!!!! allocate for 'compute_gw_Sigma_Casida'/'compute_gw_Casida_Pi'/'compute_gw_casida_Ec'/'compute_gw_casida_Ec_complex' !!!


Nguess_complex=casida_N_grid0!20!30!20
! allocate(E_guess_complex(1:Nguess_complex))
allocate(E_guess_complex_HOMO(1:Nguess_complex))
allocate(E_guess_complex_LUMO(1:Nguess_complex))
allocate(casida_Ec(1:Nguess_complex))
! allocate(casida_Ec_HOMO(1:Nguess_complex))
! allocate(casida_Ec_LUMO(1:Nguess_complex))

print *,pi

do ll=1,Nguess_complex

    E_guess_complex_HOMO(ll)=casida_Emin_grid10/hartree+(ll-1)*(casida_Emax_grid10/hartree-casida_Emin_grid10/hartree)&
                                /(Nguess_complex-1)

    E_guess_complex_LUMO(ll)=casida_Emin_grid20/hartree+(ll-1)*(casida_Emax_grid20/hartree-casida_Emin_grid20/hartree)&
                                /(Nguess_complex-1)

enddo


print *,'****  1st energy grid  ****'

do m=1,Nguess_complex
    print *,E_guess_complex_HOMO(m)*hartree
enddo

print *,'****  2nd energy grid  ****'

do m=1,Nguess_complex
    print *,E_guess_complex_LUMO(m)*hartree
enddo

















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   compute casida_Kxnj   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(psi_n_g_real(Ne*Nquadrature,Nempty))
allocate(psi_n_g_real_occupied(Ne*Nquadrature,Nstates))


orbital = casida_orbital10!1!Nstates!1!Nstates!HOMO_state!+1

do kk=1,Nempty
    psi_n_g_real(:,kk)=psi_point_g(:,orbital)*psi_point_g(:,Nstates+kk)!*volumegweight(:)
enddo


do kk=1,Nstates
    psi_n_g_real_occupied(:,kk)=psi_point_g(:,orbital)*psi_point_g(:,kk)!*volumegweight(:)
enddo

call compute_gw_Casida_Kxnj(Nn,Ne,Nquadrature,Nstates,Nempty)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!     compute <psi|Sigma_c|psi>    !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! stop

print *,'*** solve',casida_orbital10,'orbital -- Sigma_c  ***'

call compute_gw_casida_Ec_complex(E_guess_complex_HOMO,E,Nguess_complex,Nstates,Nempty)
! call compute_gw_casida_Ec_complex_2(E_guess_complex,E,size(E_guess_complex),Nstates,Nempty)



print *,'=========='
do m=1,Nguess_complex
    print *,casida_Ec(m)*hartree
enddo





orbital = casida_orbital20!1!Nstates!1!Nstates!HOMO_state!+1

do kk=1,Nempty
    psi_n_g_real(:,kk)=psi_point_g(:,orbital)*psi_point_g(:,Nstates+kk)!*volumegweight(:)
enddo


do kk=1,Nstates
    psi_n_g_real_occupied(:,kk)=psi_point_g(:,orbital)*psi_point_g(:,kk)!*volumegweight(:)
enddo

call compute_gw_Casida_Kxnj(Nn,Ne,Nquadrature,Nstates,Nempty)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!     compute <psi|Sigma_c|psi>    !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! stop

print *,'*** solve',casida_orbital20,'orbital -- Sigma_c  ***'

call compute_gw_casida_Ec_complex(E_guess_complex_LUMO,E,Nguess_complex,Nstates,Nempty)
! call compute_gw_casida_Ec_complex_2(E_guess_complex,E,size(E_guess_complex),Nstates,Nempty)



print *,'=========='
do m=1,Nguess_complex
    print *,casida_Ec(m)*hartree
enddo


! stop




! 9999 continue











































else if (casida_method0=='nlfeast') then





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! CALL NFEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!! easier case
Emid_nl=cmplx(E(1),0.0d0,8)!cmplx(E_GW(1),0.0d0,8)!cmplx(-16.0d0/hartree,0.0d0,8)!cmplx(-50.0d0/hartree,0.0d0,8)!cmplx(-7.0d0/hartree,0.0d0,8)!cmplx(-16.0d0/hartree,0.0d0,8)!cmplx(-497.0d0/hartree,0.0d0,8)!cmplx(-23.87d0/hartree,0.0d0,8)
r_nl=abs(E(1))+3.0d0/hartree!abs(E_GW(1))+3.0d0/hartree!19.0d0/hartree!51.5d0/hartree!9.0d0/hartree!19.5d0/hartree!503.0d0/hartree!28.0d0/hartree
M0_nl=20!15!2*Nstates+2!M0!15!15!20


!!!!!!!!!!!!  FEAST
! fpm(30)=121118 ! code name for direct call to this routine   
    call feastinit(fpm)
    fpm(1)=1 ! change from default (printing info on screen)
    fpm(8)=8 ! number of contour points
    fpm(15)=1 !! 1 needed for synev
    fpm(16)=1 !! Trapezoida
    fpm(42)=0 ! double (0) or single (1) precision

    fpm(45)=2 ! precision bicgstab
    fpm(46)=1000 ! number ofiterations

    fpm(5)=1

    fpm(18)=50

    fpm(3)=8


    ! fpm(4)=0
    
!! feast output  
    allocate(E_nl(1:M0_nl))     ! Eigenvalue
    allocate(X_nl(1:Nn,1:M0_nl)) ! Eigenvectors
    allocate(res_nl(1:M0_nl))   ! Residual

    ! X_nl(1:Nn,1:M0_nl)=psi(1:Nn,1:M0_nl)
    ! E_nl(1)=cmplx(E(1),0.0d0,8)!cmplx(-23.87d0/hartree,0.0d0,8)
    ! do i=2,M0_nl
    ! E_nl(i)=cmplx(E(i),0.0d0,8)!(0.0d0,0.0d0)!cmplx(2.6d0/hartree+(i-1)*2.0d-2,0.0d0,8)
    ! enddo

    X_nl=ZZERO
    E_nl=ZZERO

!     if (M0_nl>M0) then

! X_nl(1:Nn,1:M0)=psi(1:Nn,1:M0)
! E_nl(1)=cmplx(E(1),0.0d0,8)!cmplx(-23.87d0/hartree,0.0d0,8)
! do i=2,M0
!     E_nl(i)=cmplx(E(i),0.0d0,8)!(0.0d0,0.0d0)!cmplx(2.6d0/hartree+(i-1)*2.0d-2,0.0d0,8)
! enddo

!     else

X_nl(1:Nn,1:M0_nl)=psi(1:Nn,1:M0_nl)
E_nl(1)=cmplx(E(1),0.0d0,8)!cmplx(-23.87d0/hartree,0.0d0,8)
do i=2,M0_nl
    E_nl(i)=cmplx(E(i),0.0d0,8)!(0.0d0,0.0d0)!cmplx(2.6d0/hartree+(i-1)*2.0d-2,0.0d0,8)
enddo

    ! end if

    ! print *,X_nl(1:3,1:2)

    ! stop

    
    
!!!  store factorization
    nfact=1
    if (fpm(10)==1) nfact=fpm(8)

    ! print *,nfact

    allocate(ipivloc_nl(Nn,nfact))

!  allocate(Az_nl(Nn,Nn,nfact))
    allocate(Aztemp_nl(Nn,Nn))

!! RCI work-array
allocate(zAq_nl(M0_nl,M0_nl))
allocate(zBq_nl(M0_nl,M0_nl))
allocate(work_nl(Nn,M0_nl))  
allocate(zwork_nl(Nn,2*M0_nl))



nnzA=neigh_size/2+Nn ! slightly overestimated
allocate(uisa(Nn+1))
allocate(ujsa(nnzA))
allocate(usa(nnzA))
allocate(usb(nnzA))
call dcsr_convert_upper(Nn,UPLO,H,IA,JA,usa,uisa,ujsa)
call dcsr_convert_upper(Nn,UPLO,B,IA,JA,usb,uisa,ujsa)
nnza=uisa(Nn+1)-1 ! by definition


! allocate(usa_dft(nnzA))
! call dcsr_convert_upper(Nn,UPLO,H_dft,IA,JA,usa_dft,uisa,ujsa)
! allocate(psaz_dft(nnza))


!!!! mixed precision set-up (only double for now)
allocate(saz(nnza,nfact))
allocate(psaz(nnza,nfact))

allocate(saztemp(nnza))




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!  use same factorization for (normal+transpose solve)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MAXFCT=1 
MNUM=1
MTYPE=6      ! complex and symmetric
allocate(PT(1:64,1:nfact))
allocate(IPARM(1:64,1:nfact))
do i=1,nfact !! multiple factorization
    pt(:,i)=0
    call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
end do


!IPARM(2,:)=3 ! parallel omp nested dissection (sensible to #threads- no consistency- same runs with different results)
!IPARM(4)=11 !CGS solver
IPARM(25,:)=1 ! parallel omp rhs solve
IPARM(11,:)=0 ! disable scaling (taking care by feast fpm(41))
!!!!!!!!!!!!
if (fpm(64)==1) then
    do i=1,64
        if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
    enddo
endif
!!!!!!!!!!!!
IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used  !!
MSGLVL=0!0 !0- no output, 1- output

! !IPARM(2,:)=3 ! parallel omp nested dissection (sensible to #threads- no consistency- same runs with different results)
! !IPARM(4)=11 !CGS solver
! aIPARM_V(25)=1 ! parallel omp rhs solve
! aIPARM_V(11)=0 ! disable scaling (taking care by feast fpm(41))
! !!!!!!!!!!!!
! if (fpm(64)==1) then
!    do i=1,64
!       if (fpm(64+i)/=-111) aiparm_V(i)=fpm(64+i)
!    enddo
! endif
! !!!!!!!!!!!!
! aIPARM_V(6)=1 ! solution and rhs are input/output, attention zaux is always used  !!
! MSGLVL=0!0 !0- no output, 1- output


allocate(zaux(Nn,M0_nl))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! for bicgstab
    allocate(nres(M0_nl))
    !allocate(norm(M0))


comb=.false.!.true.!.false.
linloops=10000
com=.false. !default
if (fpm(1)/=0) com=.true.
if (fpm(1)<0) then
    fout=abs(fpm(1))+200!fpm(60) !!file id name
else
    fout=6 !screen
endif

!com=.false.

    ! bicgstab-rci     
allocate(bwork1(Nn,7*M0_nl))
allocate(bwork2(M0_nl,4)) 


allocate(bwork1_vectemp1(1:Nn,1:M0_nl))



! allocate(bwork1_vectemp3(1:Nn,1:M0_nl))
! allocate(bwork1_vectemp4(1:Nn,1:M0_nl*fpm(8)/2))
allocate(bwork1_vectemp1_g(1:Ne*Nquadrature,1:M0_nl))
! allocate(bwork1_vectemp2_g(1:Nstates*Nempty,1:M0_nl))
allocate(bwork1_vectemp3_g(1:Ne*Nquadrature,1:M0_nl))





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  





ijob=-1 ! initialization 
do while (ijob/=0)
call zfeast_srcinev(ijob,Nn,Ze_nl,work_nl,zwork_nl,zAq_nl,zBq_nl,fpm,epsout,loop,Emid_nl,r_nl,M0_nl,E_nl,X_nl,M_nl,res_nl,info)
    
    
    select case(ijob)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(10) !! !! form T(Ze) and factorize preconditioner if any-- dongming does T(Ze)=ZeS-H-Sigma(ze) == Az, or Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    
    ! print *,'case 10'
    
        id=fpm(33) !! id of factorization (for fpm(10) flag) 
        
    !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

        saz(1:nnza,id)=Ze_nl*usb(1:nnza)-usa(1:nnza)
            
        !!!! Factorize (preconditioner) 
            PHASE=12 !include the symbolic factorization
            psaz(1:nnza,id)=saz(1:nnza,id) ! save preconditioner
call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1,bwork1,infoloc3)  
    ! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1,bwork1,infoloc3)  



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(11) !!solve the linear system  T(Ze)x=workc(1:N,1:fpm(23)) result in to workc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! print *,'case 11'

            id=fpm(33)
    
            !!!! direct solver 
            !PHASE=33 ! solve
            !call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
    !!!! iterative bicgstab solver
            lintargeterror=10d0**(-fpm(45))
            linloops=fpm(46)
            nres=1.0d0 ! max rhs norm
            zaux(1:Nn,1:fpm(23))=(0.0d0,0.0d0) !! initial guess
            !call zbicgstab(fpm(44),'U','N',saz(1,id),uisa,ujsa,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! RCI bicgstab using ZeI-A as preconditoner
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        ijob2=-1
        do while (ijob2/=0)
            call zbicgstab_rci(ijob2,Nn,fpm(23),zwork_nl,zaux,bwork1,bwork2,jj1,jj2,nres,linloops,lintargeterror,comb,infoloc2) 
        
            select case(ijob2)
                case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
            !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input
    
            !!!! direct solver 
                    PHASE=33 ! solve
                IPARM(6,:)=0 ! solution and rhs are separated
        call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1(1,jj2)&
        ,bwork1(1,jj1),infoloc3)
        !                    
        !         aIPARM_V(6)=0 ! solution and rhs are separated                                                                             
        ! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1(1,jj2)&
                                                                                                ! ,bwork1(1,jj1),infoloc3)
            
            
        case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)

        call wzcsrmm('U','N',Nn,Nn,fpm(23),(1.0d0,0.0d0),saz(1,id),uisa,ujsa,bwork1(1,jj1),(0.0d0,0.0d0),bwork1(:,jj2))

        ! call mkl_zcsrmm('N',Nn,fpm(23),Nn,(1.0d0,0.0d0),matdescrb,saz(1,id),ujsa,uisa,uisa(2),bwork1(1,jj1),Nn&
        ! ,(0.0d0,0.0d0),bwork1(:,jj2),Nn)

        
        !!!!!!!!!!!!!!! Sigma_x !!!!!!!!!!!!!!!!!
        
        Aztemp_nl=-Vhfx_final_GW*(1.0d0,0.0d0)
        
        call ZSYMM('L','L',Nn,fpm(23),(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)
    !    CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)

        bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)+bwork1_vectemp1

            
            end select
        end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
            
            
            call ZCOPY( Nn*fpm(23), zaux, 1, zwork_nl, 1 )
    
        fpm(60)=fpm(60)+linloops
            if (com) then
                write(fout,'(A)',advance='no') '  #it  '
                write(fout,'(I4)',advance='no') linloops
                write(fout,'(A)',advance='no') '; res min='
                write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
                write(fout,'(A)',advance='no') '; res max='
                write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
                write(fout,*)
            endif

! stop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    case(33) !! compute norm at Ze  ||T(Ze)||
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        ! print *,'case 33'


        Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-H_dense-Vhfx_final_GW), Nn, zwork_nl )
        ! Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-Vhfx_final_GW), Nn, zwork_nl )


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        case(30) !! perform multiplication
        
        !!!!
        do ii=1,fpm(23)
            !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

            ! print 100,ii,fpm(23),E_nl(ii),res(ii)
            print *,ii,fpm(23),E_nl(ii)

            


            saztemp(1:nnza)=E_nl(ii)*usb(1:nnza)-usa(1:nnza)

            call wzcsrmm('U','N',Nn,Nn,1,(-1.0d0,0.0d0),saztemp,uisa,ujsa,X_nl(1,ii),(0.0d0,0.0d0),bwork1_vectemp1(:,1))

            ! call mkl_zcsrmm('N',Nn,1,Nn,(1.0d0,0.0d0),matdescrb,saztemp,ujsa,uisa,uisa(2),X_nl(1,ii),Nn&
            ! ,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

            work_nl(:,ii)=bwork1_vectemp1(:,1)

            Aztemp_nl=-Vhfx_final_GW
        
            call ZSYMM('L','L',Nn,1,(-1.0d0,0.0d0),Aztemp_nl(1,1),Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)
            ! CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl,Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

            work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)


        end do

    end select

end do








! do ii=1,M0_nl
!     E_nl_temp(ii)=E_nl(ii)
!     X_nl_temp(:,ii)=X_nl(:,ii)
!     print *,E_nl(ii)*hartree,res_nl(ii)
! enddo

! print *,'======================'

! do ii=1,M_nl
!     print *,E_nl(ii)*hartree,res_nl(ii)
! enddo

print *,'+++++++++++++++++++++++++++++++++++'


call quicksort_eigenpairs(E_nl,X_nl,res_nl,Nn,1,M0_nl)

allocate(E_nl_temp(1:M0_nl))     ! Eigenvalue
allocate(X_nl_temp(1:Nn,1:M0_nl)) ! Eigenvectors
allocate(res_nl_temp(1:M0_nl))  


print *,'***found eigenvalues with Sigma_x only***'

do ii=1,M0_nl
    E_nl_temp(ii)=E_nl(ii)
    X_nl_temp(:,ii)=X_nl(:,ii)
    res_nl_temp(ii)=res_nl(ii)
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo

print *,'+++++++++++++++++++++++++++++++++++'

! print *,'======================'

! do ii=1,M_nl
!     print *,E_nl(ii)*hartree,res_nl(ii)
! enddo

! stop

do i=1,M0_nl

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    ! print *,dot_product(dble(X_nl(:,i)),snq1)

    vtemp=dot_product(dble(X_nl(:,i)),snq1)

    X_nl(:,i)=X_nl(:,i)/sqrt(vtemp)

    ! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    ! print *,dot_product(dble(X_nl(:,i)),snq1)


    ! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    ! ! print *,dot_product(dble(X_nl_temp(:,i)),snq1)

    ! vtemp=dot_product(dble(X_nl_temp(:,i)),snq1)

    ! X_nl_temp(:,i)=X_nl_temp(:,i)/sqrt(vtemp)

    ! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    ! print *,dot_product(dble(X_nl_temp(:,i)),snq1)

enddo


do k=1,Nstates+3
write (file_id, '(I0)') k
! open (unit=20,file='CO.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='CO_0.6415.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='He.1_p2_nessie.psi'//trim(file_id),status='old')
! open(unit=20,file=trim(name)//'_p2_dft.psi'//trim(file_id),status='replace')
open(unit=20,file=trim(name)//'_'//trim(dgr_fem)//'_'//trim(meanfield)//'_2.psi'//trim(file_id),status='replace')
! open (unit=20,file='LiH.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
do i=1,Nn
    write(20,*) dble(X_nl(i,k))
enddo
close (20)
enddo








! allocate(Nempty_list(1))
! ! Nempty_list = (/1250,1500,1750,2000,2250,2500/)
! Nempty_list = (/1500/)
! ! Nempty_list = (/1496/)

! do l=1,size(Nempty_list) !!! Nemptyloop

    ! Nempty=Nempty_list(l)
Nempty=casida_N_empty0!1000!1500!1000!1500

! allocate(psi_ia(Nn,Nstates*Nempty))
! allocate(psi_ia_real(Nn,Nstates*Nempty))
! do k=1,Nstates
!     do kk=1,Nempty
!         psi_ia(:,(k-1)*Nempty+kk)=sqrt(2.0d0)*psi_GW(:,k)*psi_GW(:,Nstates+kk)*j_real
!         psi_ia_real(:,(k-1)*Nempty+kk)=psi_GW(:,k)*psi_GW(:,Nstates+kk)
!     enddo
! enddo

allocate(psi_ia_g_real(Ne*Nquadrature,Nstates*Nempty))
do k=1,Nstates
    do kk=1,Nempty
        psi_ia_g_real(:,(k-1)*Nempty+kk)=psi_point_g(:,k)*psi_point_g(:,Nstates+kk)!*volumegweight(:)
    enddo
enddo

! allocate(psi_n_g_real(Ne*Nquadrature,Nempty))
! do kk=1,Nempty
!     psi_n_g_real(:,kk)=psi_GW_g(:,orbital)*psi_GW_g(:,Nstates+kk)!*volumegweight(:)
! enddo

! print *,'toto'

! allocate(psi_n_g_real_occupied(Ne*Nquadrature,Nstates))
! do kk=1,Nstates
!     psi_n_g_real_occupied(:,kk)=psi_GW_g(:,orbital)*psi_GW_g(:,kk)!*volumegweight(:)
! enddo

! print *,'toto'


allocate(NgNntemp(1:Ne*Nquadrature,1:Nstates*Nempty))
! allocate(NnNgtemp(1:Nn,1:Nstates*Nempty))

allocate(casida_Kx(1:Nstates*Nempty,1:Nstates*Nempty))
allocate(casida_R(1:Nstates*Nempty))
allocate(casida_Xs(1:Nstates*Nempty,1:Nstates*Nempty))
allocate(casida_omega(1:Nstates*Nempty))
allocate(casida_R_half(1:Nstates*Nempty))
! allocate(casida_V(1:Nstates*Nempty))
! allocate(casida_Kxnj(1:Nempty,1:Nstates*Nempty))
! allocate(casida_Kxnj_occupied(1:Nstates,1:Nstates*Nempty))
! allocate(casida_cvomega(1:Nstates*Nempty))


allocate(NNmatrix_temp12(1:Nn,1:Nstates*Nempty))
allocate(NNmatrix_temp13(1:Nn,1:Nstates*Nempty))
! allocate(NNmatrix_temp14(1:Nstates*Nempty,1:Nstates*Nempty))
! allocate(vector_temp00(1:Nn))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   compute casida_Kx   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call compute_gw_Casida_Kx(Nn,Ne,Nquadrature,Nstates,Nempty)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   solve casida_equation   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call compute_gw_Casida_matrix(Nstates,Nempty,E)

call solve_Casida_matrix(Nstates,Nempty)






print *,"++++++++++++++++++++++++"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     Non-Linear FEAST solver                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! allocate(vector_temp_complex(1:Nn))

do k=1,Nstates*Nempty
casida_Xs(:,k)=sqrt((1.0d0/sqrt(casida_omega(k))))*casida_R_half(:)*casida_Xs(:,k)
enddo

CALL DGEMM('N','N',Nn,Nstates*Nempty,Nstates*Nempty,1.0d0,NNmatrix_temp12,Nn,casida_Xs,Nstates*Nempty,0.0d0,NNmatrix_temp13,Nn)





















!!!!!!!!!!!!  FEAST
! fpm(30)=121118 ! code name for direct call to this routine   
call feastinit(fpm)
fpm(1)=casida_fpm10!1 ! change from default (printing info on screen)
fpm(8)=casida_fpm80!8 ! number of contour points
fpm(15)=1 !! 1 needed for synev
fpm(16)=1 !! Trapezoida
fpm(42)=0 ! double (0) or single (1) precision

fpm(45)=casida_fpm450!2!6!2 ! precision bicgstab
fpm(46)=1000 ! number ofiterations

fpm(5)=1

! fpm(18)=50

! easier case
Emid_nl=casida_Emid0/hartree!cmplx(dble(E_nl(1)),0.0d0,8)!cmplx(-10.5d0/hartree,0.0d0,8)!cmplx(dble(E_nl(1)),0.0d0,8)!cmplx(-16.0d0/hartree,0.0d0,8)!cmplx(-50.0d0/hartree,0.0d0,8)!cmplx(-23.0d0/hartree,0.0d0,8)
r_nl=casida_r0/hartree!13.5d0/hartree!51.5d0/hartree!27.0d0/hartree
! r_nl=abs(dble(E_nl(1)))+2.5d0/hartree!E_nl(Nstates+2)!+(E_nl(Nstates+3)-E_nl(Nstates+2))/2.0d0  !!!! Emax is slightly higher than LUMO+1 states 
M0_nl=casida_M00!3!5!3!Nstates+2!4!11

! Emid_nl=cmplx(dble(E_nl(2)),0.0d0,8)
! ! r_nl=18.4d0/hartree!51.5d0/hartree!27.0d0/hartree
! r_nl=abs(dble(E_nl(2)))+3.0d0/hartree!E_nl(Nstates+2)!+(E_nl(Nstates+3)-E_nl(Nstates+2))/2.0d0  !!!! Emax is slightly higher than LUMO+1 states 
! M0_nl=5!4!11

fpm(3)=casida_fpm30!10!8!7
! fpm(6)=1  ! 0: Using relative error on the trace epsout i.e. epsout / 1: Using relative residual res i.e. maxi res(i) 
! fpm(45)=2 ! precision bicgstab





!! feast output  

deallocate(E_nl,X_nl,res_nl)
allocate(E_nl(1:M0_nl))     ! Eigenvalue
allocate(X_nl(1:Nn,1:M0_nl)) ! Eigenvectors
!   allocate(E_nl_temp(1:M0_nl))     ! Eigenvalue
!   allocate(X_nl_temp(1:Nn,1:M0_nl)) ! Eigenvectors
allocate(res_nl(1:M0_nl))   ! Residual

X_nl(1:Nn,1:M0_nl)=X_nl_temp(1:Nn,casida_igorbital0:M0_nl+(casida_igorbital0-1))
E_nl(1:M0_nl)=E_nl_temp(casida_igorbital0:M0_nl+(casida_igorbital0-1))
! res_nl(1:M0_nl)=res_nl_temp(1:M0_nl)



allocate(bwork1_backup(Nn,M0_nl,fpm(8)/2))
allocate(NNmatrix_temp005(Nn,Nn))
allocate(bwork1_vectemp2(1:Nstates*Nempty,1:M0_nl))
! allocate(bwork1_vectemp2(1:Nstates*Nempty,1:M0_nl))
! ! allocate(bwork1_vectemp3(1:Nn,1:M0_nl))
! ! allocate(bwork1_vectemp4(1:Nn,1:M0_nl*fpm(8)/2))
! allocate(bwork1_vectemp1_g(1:Ne*Nquadrature,1:M0_nl))
! ! allocate(bwork1_vectemp2_g(1:Nstates*Nempty,1:M0_nl))
! allocate(bwork1_vectemp3_g(1:Ne*Nquadrature,1:M0_nl))


ii10=0
ii11=0
ii30=0


ijob=-1 ! initialization 
do while (ijob/=0)
call zfeast_srcinev(ijob,Nn,Ze_nl,work_nl,zwork_nl,zAq_nl,zBq_nl,fpm,epsout,loop,Emid_nl,r_nl,M0_nl,E_nl,X_nl,M_nl,res_nl,info)
  
   
   select case(ijob)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case(10) !! !! form T(Ze) and factorize preconditioner if any-- dongming does T(Ze)=ZeS-H-Sigma(ze) == Az, or Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    
    ! print *,'case 10'
    
     id=fpm(33) !! id of factorization (for fpm(10) flag) 
     
    !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

        saz(1:nnza,id)=Ze_nl*usb(1:nnza)-usa(1:nnza)
            
     !!!! Factorize (preconditioner) 
            PHASE=12 !include the symbolic factorization
            psaz(1:nnza,id)=saz(1:nnza,id) ! save preconditioner
call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1,bwork1,infoloc3)   
! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1,bwork1,infoloc3)   




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(11) !!solve the linear system  T(Ze)x=workc(1:N,1:fpm(23)) result in to workc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! print *,'case 11'

            id=fpm(33)


            ! if (mod(ii11,fpm(8))<fpm(8)/2) then

            
    
           !!!! direct solver 
            !PHASE=33 ! solve
           !call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
    !!!! iterative bicgstab solver
            lintargeterror=10d0**(-fpm(45))
            linloops=fpm(46)
            nres=1.0d0 ! max rhs norm
            zaux(1:Nn,1:fpm(23))=(0.0d0,0.0d0) !! initial guess
            !call zbicgstab(fpm(44),'U','N',saz(1,id),uisa,ujsa,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! RCI bicgstab using ZeI-A as preconditoner
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      ijob2=-1
      do while (ijob2/=0)
    call zbicgstab_rci(ijob2,Nn,fpm(23),zwork_nl,zaux,bwork1,bwork2,jj1,jj2,nres,linloops,lintargeterror,comb,infoloc2) 
        ! call zbicgstab_rci(ijob2,Nn,fpm(23),zwork_nl,zaux,bwork1,bwork2,jj1,jj2,nres,linloops,lintargeterror,.true.,infoloc2) 
     
         select case(ijob2)
              case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
            !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input
    
          !!!! direct solver 
                 PHASE=33 ! solve
                IPARM(6,:)=0 ! solution and rhs are separated
        call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),IPARM(1,id),MSGLVL,bwork1(1,jj2)&
                                                                                                        ,bwork1(1,jj1),infoloc3)

        !         aIPARM_V(6)=0 ! solution and rhs are separated  
        ! call PARDISO(aPT_V,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz(1,id),uisa,ujsa,idum,fpm(23),aIPARM_V,MSGLVL,bwork1(1,jj2)&
        !                                                                                         ,bwork1(1,jj1),infoloc3)
          
            ! print *,bwork1(1:3,jj1)
           
      case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)


        

        call wzcsrmm('U','N',Nn,Nn,fpm(23),(1.0d0,0.0d0),saz(1,id),uisa,ujsa,bwork1(1,jj1),(0.0d0,0.0d0),bwork1(:,jj2))

        ! print *,bwork1(1:3,jj2)

        !!!!!!!!!!!!!!! Sigma_x !!!!!!!!!!!!!!!!!
        
        Aztemp_nl=-Vhfx_final_GW*(1.0d0,0.0d0)
        
       call ZSYMM('L','L',Nn,fpm(23),(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)
        ! CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl(1,1),Nn,bwork1(1,jj1),Nn,(0.0d0,0.0d0),bwork1_vectemp1,Nn)

       bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)+bwork1_vectemp1(:,1:fpm(23))

    !    print *,bwork1(1:3,jj2)
       
!         !!!!!!!!!!!!!!! Sigma_c !!!!!!!!!!!!!!!!!
       
!     !    if (mod(ii11,fpm(8))<fpm(8)/2) then
!     !     bwork1_backup(:,:,mod(ii11,fpm(8))+1)=(0.0d0,0.0d0)

! !    call mkl_zcsrmv('N',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
! !                                                                 ,bwork1(1,jj1),(0.0d0,0.0d0),bwork1_vectemp1_g)

!        call mkl_zcsrmm('N',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ,bwork1(1,jj1),Nn,ZZERO,bwork1_vectemp3_g,Ne*Nquadrature)
!    !    call mkl_zcsrmm('N',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!    !     ,psi_GW(:,1:fpm(23))*j_real,Nn,ZZERO,bwork1_vectemp3_g,Ne*Nquadrature)

!        ! call wzcsrmm('F','N',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1(1,jj1)&
!        !                                                                         ,ZZERO,bwork1_vectemp1_g)

!    !    call wzcsrmm('F','N',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,psi_GW(:,1:fpm(23))&
!    !                                                                             ,ZZERO,bwork1_vectemp1_g)

!        ! bwork1_vectemp1_g=bwork1_vectemp3_g

!        ! print *,norm2(psi_GW_g(:,1)-dble(bwork1_vectemp3_g(:,1)))

!        ! stop

! !    bwork1_vectemp4(:,:)=(0.0d0,0.0d0)
!    ! Ze_nl=(0.0d0,0.0d0)
!        do kk=1,Nstates+Nempty

!            bwork1_vectemp1_g=bwork1_vectemp3_g

!            do i=1,M0_nl
!                ! bwork1_vectemp1_g(:,i)=psi_GW_g(:,kk)*psi_GW_g(:,i)*volumegweight(:)
!                bwork1_vectemp1_g(:,i)=psi_point_g(:,kk)*bwork1_vectemp1_g(:,i)*volumegweight(:)
!            enddo

!            ! if (kk==2) print *,bwork1_vectemp1_g(1,1)


!        ! call mkl_zcsrmv('T',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ! ,bwork1_vectemp1_g,(0.0d0,0.0d0),bwork1_vectemp1)

!        call mkl_zcsrmm('T',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ,bwork1_vectemp1_g,Ne*Nquadrature,ZZERO,bwork1_vectemp1,Nn)

!        ! call wzcsrmm('F','T',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1_vectemp1_g&
!        !                                                                         ,ZZERO,bwork1_vectemp1)

!        ! print *,bwork1_vectemp1(1:3,1)

!        ! stop

!        ! if (kk==2) print *,bwork1_vectemp1(1,1)



!        NNmatrix_temp005=NNmatrix_temp13*j_real

!        ! call zgemv('T',Nn,Nstates*Nempty,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp1,1,(0.0d0,0.0d0),bwork1_vectemp2,1)

!        CALL ZGEMM('T','N',Nstates*Nempty,fpm(23),Nn,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp1,Nn&
!        ,(0.0d0,0.0d0),bwork1_vectemp2,Nstates*Nempty)

!        ! print *,bwork1_vectemp2(1:3,1)

!        ! if (kk==2) print *,bwork1_vectemp2(1,1)

!        ! stop


! ! Ze_nl=Ze_nl+2.0d0*dot_product(bwork1_vectemp2(:,1),bwork1_vectemp2(:,1)*1.0d0/(E_guess(1)-E_GW(kk)+sqrt(casida_omega(:))))

! ! go to 1122

!        if (kk<=Nstates) then
!        do i=1,M0_nl
!        bwork1_vectemp2(:,i)=bwork1_vectemp2(:,i)/(Ze_nl-E(kk)+sqrt(casida_omega(:)))!*2.0d0
!        enddo
!        else
!        do i=1,M0_nl
!        bwork1_vectemp2(:,i)=bwork1_vectemp2(:,i)/(Ze_nl-E(kk)-sqrt(casida_omega(:)))!*2.0d0
!        enddo
!        end if

!        ! call zgemv('N',Nn,Nstates*Nempty,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp2,1,(0.0d0,0.0d0),bwork1_vectemp1,1)
!        CALL ZGEMM('N','N',Nn,fpm(23),Nstates*Nempty,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp2,Nstates*Nempty&
!        ,(0.0d0,0.0d0),bwork1_vectemp1,Nn)

!        ! print *,bwork1_vectemp1(1,1)
!        ! if (kk==2) print *,bwork1_vectemp1(1,1)
!        ! call mkl_zcsrmv('N',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ! ,bwork1_vectemp1,(0.0d0,0.0d0),bwork1_vectemp1_g)
!        call mkl_zcsrmm('N',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ,bwork1_vectemp1,Nn,ZZERO,bwork1_vectemp1_g,Ne*Nquadrature)

!        ! print *,bwork1_vectemp3_g(1,1)
!        ! if (kk==2) print *,bwork1_vectemp1_g(1,1)

!        ! call wzcsrmm('F','N',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1_vectemp1&
!        !                                                                         ,ZZERO,bwork1_vectemp3_g)

!        do i=1,M0_nl
!            bwork1_vectemp1_g(:,i)=psi_point_g(:,kk)*bwork1_vectemp1_g(:,i)*volumegweight(:)
!        enddo

!        ! if (kk==2) print *,bwork1_vectemp1_g(1,1)

!        ! print *,bwork1_vectemp1_g(1,1)

!        ! call mkl_zcsrmv('T',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ! ,bwork1_vectemp1_g,(0.0d0,0.0d0),bwork1_vectemp3)

!        call mkl_zcsrmm('T',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!        ,bwork1_vectemp1_g,Ne*Nquadrature,ZZERO,bwork1_vectemp1,Nn)

!        ! if (kk==2) print *,bwork1_vectemp3(1,1)

!        ! if (kk==2) stop

!        ! print *,bwork1_vectemp3(1,1)

!        ! stop

!        ! call wzcsrmm('F','T',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1_vectemp1_g&
!        !                                                                         ,ZZERO,bwork1_vectemp3)


!        bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)-bwork1_vectemp1(:,1:fpm(23))*2.0d0
!     !    bwork1_backup(:,:,mod(ii11,fpm(8))+1)=bwork1_backup(:,:,mod(ii11,fpm(8))+1)-bwork1_vectemp1(:,1:fpm(23))*2.0d0
!     !    bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)-dble(bwork1_vectemp1(:,1:fpm(23)))*2.0d0
!     !    bwork1_backup(:,:,mod(ii11,fpm(8))+1)=bwork1_backup(:,:,mod(ii11,fpm(8))+1)-dble(bwork1_vectemp1(:,1:fpm(23)))*2.0d0


!        if (mod(kk,500)==0) print *,kk
!        enddo

!     ! !    print *,bwork1(1:3,jj2)

!     ! else if ((mod(ii11,fpm(8))>fpm(8)/2-1).and.(mod(ii11,fpm(8))<fpm(8))) then


!     !     bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)+conjg(bwork1_backup(:,:,fpm(8)-mod(ii11,fpm(8))))!%re
!     !     ! bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)-bwork1_backup(:,:,fpm(8)-mod(ii11,fpm(8)))%im*j_imag
!     !     ! bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)+bwork1_backup(:,:,fpm(8)-mod(ii11,fpm(8)))%re
!     !     ! bwork1(:,jj2:jj2+fpm(23)-1)=bwork1(:,jj2:jj2+fpm(23)-1)-bwork1_backup(:,:,fpm(8)-mod(ii11,fpm(8)))%im*j_imag


!     ! end if  !!! end selection of contour nodes

           
         end select
      end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
           call ZCOPY( Nn*fpm(23), zaux, 1, zwork_nl, 1 )
    
      fpm(60)=fpm(60)+linloops
            if (com) then
            ! if (.true.) then
               write(fout,'(A)',advance='no') '  #it  '
               write(fout,'(I4)',advance='no') linloops
               write(fout,'(A)',advance='no') '; res min='
               write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
               write(fout,'(A)',advance='no') '; res max='
               write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
               write(fout,*)
            endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !     zwork_nl%im=0.0d0

    !     bwork1_vectemp4(:,mod(ii11,fpm(8))*fpm(23)+1:mod(ii11,fpm(8))*fpm(23)+fpm(23))=zwork_nl!bwork1(:,jj2:jj2+fpm(23)-1)


    !     else if ((mod(ii11,fpm(8))>fpm(8)/2-1).and.(mod(ii11,fpm(8))<fpm(8))) then


    ! zwork_nl%re = bwork1_vectemp4(:,((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+1):((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+fpm(23)))%re
    ! zwork_nl%im = 0.0d0   !-bwork1_vectemp4(:,((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+1):((fpm(8)-1-mod(ii11,fpm(8)))*fpm(23)+fpm(23)))%im


    !     end if 


            ii11=ii11+1

! stop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    case(33) !! compute norm at Ze  ||T(Ze)||
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        ! print *,'case 33'

        Ze_nl=ZLANGE('F', Nn, Nn, (Ze_nl*S_dense-H_dense-Vhfx_final_GW), Nn, zwork_nl )


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(30) !! perform multiplication
        
        !!form T(lambda(i)) and multiply by x(1:N,i) result in work(1:N,i)

        ! print *,'case 30'
        
        !!!!
        do ii=1,fpm(23)
          !!! Add the contribution of the energy to the matrix (Ae<=EI-A)

            
            ! print 100,ii,fpm(23),E_nl(ii),res(ii)
            print *,ii,fpm(23),E_nl(ii)
            ! print *,X_nl(1:3,ii)

            


            saztemp(1:nnza)=E_nl(ii)*usb(1:nnza)-usa(1:nnza)

            call wzcsrmm('U','N',Nn,Nn,1,(-1.0d0,0.0d0),saztemp,uisa,ujsa,X_nl(1,ii),(0.0d0,0.0d0),bwork1_vectemp1(:,1))

            work_nl(:,ii)=bwork1_vectemp1(:,1)

            Aztemp_nl=-Vhfx_final_GW
        
            call ZSYMM('L','L',Nn,1,(-1.0d0,0.0d0),Aztemp_nl(1,1),Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)
            ! CALL ZGEMM('N','N',Nn,fpm(23),Nn,(1.0d0,0.0d0),Aztemp_nl,Nn,X_nl(1,ii),Nn,(0.0d0,0.0d0),bwork1_vectemp1(:,1),Nn)

            work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)



            !!!!!!!!!!!!!!! Sigma_c !!!!!!!!!!!!!!!!!
       

!    call mkl_zcsrmv('N',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                                                                 ,bwork1(1,jj1),(0.0d0,0.0d0),bwork1_vectemp1_g)

       call mkl_zcsrmm('N',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ,X_nl(1,ii),Nn,ZZERO,bwork1_vectemp3_g,Ne*Nquadrature)
   !    call mkl_zcsrmm('N',Ne*Nquadrature,fpm(23),Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
   !     ,psi_GW(:,1:fpm(23))*j_real,Nn,ZZERO,bwork1_vectemp3_g,Ne*Nquadrature)

       ! call wzcsrmm('F','N',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1(1,jj1)&
       !                                                                         ,ZZERO,bwork1_vectemp1_g)

   !    call wzcsrmm('F','N',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,psi_GW(:,1:fpm(23))&
   !                                                                             ,ZZERO,bwork1_vectemp1_g)

       ! bwork1_vectemp1_g=bwork1_vectemp3_g

       ! print *,norm2(psi_GW_g(:,1)-dble(bwork1_vectemp3_g(:,1)))

       ! stop

!    bwork1_vectemp4(:,:)=(0.0d0,0.0d0)
   ! Ze_nl=(0.0d0,0.0d0)
       do kk=1,Nstates+Nempty

        ! if (mod(kk,10)==0) print *,kk

           bwork1_vectemp1_g=bwork1_vectemp3_g

            bwork1_vectemp1_g(:,1)=psi_point_g(:,kk)*bwork1_vectemp1_g(:,1)*volumegweight(:)

           ! if (kk==2) print *,bwork1_vectemp1_g(1,1)


       ! call mkl_zcsrmv('T',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ! ,bwork1_vectemp1_g,(0.0d0,0.0d0),bwork1_vectemp1)

       call mkl_zcsrmm('T',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ,bwork1_vectemp1_g,Ne*Nquadrature,ZZERO,bwork1_vectemp1,Nn)

       ! call wzcsrmm('F','T',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1_vectemp1_g&
       !                                                                         ,ZZERO,bwork1_vectemp1)

       ! print *,bwork1_vectemp1(1:3,1)

       ! stop

       ! if (kk==2) print *,bwork1_vectemp1(1,1)



       NNmatrix_temp005=NNmatrix_temp13*j_real

       ! call zgemv('T',Nn,Nstates*Nempty,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp1,1,(0.0d0,0.0d0),bwork1_vectemp2,1)

       CALL ZGEMM('T','N',Nstates*Nempty,1,Nn,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp1,Nn&
       ,(0.0d0,0.0d0),bwork1_vectemp2,Nstates*Nempty)

       ! print *,bwork1_vectemp2(1:3,1)

       ! if (kk==2) print *,bwork1_vectemp2(1,1)

       ! stop


! Ze_nl=Ze_nl+2.0d0*dot_product(bwork1_vectemp2(:,1),bwork1_vectemp2(:,1)*1.0d0/(E_guess(1)-E_GW(kk)+sqrt(casida_omega(:))))

! go to 1122

       if (kk<=Nstates) then

       bwork1_vectemp2(:,1)=bwork1_vectemp2(:,1)/(E_nl(ii)-E(kk)+sqrt(casida_omega(:)))!*2.0d0

       else

       bwork1_vectemp2(:,1)=bwork1_vectemp2(:,1)/(E_nl(ii)-E(kk)-sqrt(casida_omega(:)))!*2.0d0

       end if

       ! call zgemv('N',Nn,Nstates*Nempty,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp2,1,(0.0d0,0.0d0),bwork1_vectemp1,1)
       CALL ZGEMM('N','N',Nn,1,Nstates*Nempty,(1.0d0,0.0d0),NNmatrix_temp005,Nn,bwork1_vectemp2,Nstates*Nempty&
       ,(0.0d0,0.0d0),bwork1_vectemp1,Nn)

       ! print *,bwork1_vectemp1(1,1)
       ! if (kk==2) print *,bwork1_vectemp1(1,1)
       ! call mkl_zcsrmv('N',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ! ,bwork1_vectemp1,(0.0d0,0.0d0),bwork1_vectemp1_g)
       call mkl_zcsrmm('N',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ,bwork1_vectemp1,Nn,ZZERO,bwork1_vectemp1_g,Ne*Nquadrature)

       ! print *,bwork1_vectemp3_g(1,1)
       ! if (kk==2) print *,bwork1_vectemp1_g(1,1)

       ! call wzcsrmm('F','N',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1_vectemp1&
       !                                                                         ,ZZERO,bwork1_vectemp3_g)


        bwork1_vectemp1_g(:,1)=psi_point_g(:,kk)*bwork1_vectemp1_g(:,1)*volumegweight(:)


       ! if (kk==2) print *,bwork1_vectemp1_g(1,1)

       ! print *,bwork1_vectemp1_g(1,1)

       ! call mkl_zcsrmv('T',Ne*Nquadrature,Nn,(1.0d0,0.0d0),matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ! ,bwork1_vectemp1_g,(0.0d0,0.0d0),bwork1_vectemp3)

       call mkl_zcsrmm('T',Ne*Nquadrature,1,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
       ,bwork1_vectemp1_g,Ne*Nquadrature,ZZERO,bwork1_vectemp1,Nn)

       ! if (kk==2) print *,bwork1_vectemp3(1,1)

       ! if (kk==2) stop

       ! print *,bwork1_vectemp3(1,1)

       ! stop

       ! call wzcsrmm('F','T',Ne*Nquadrature,Nn,fpm(23),ZONE,NnToNg*j_real,NnToNg_IA,NnToNg_JA,bwork1_vectemp1_g&
       !                                                                         ,ZZERO,bwork1_vectemp3)


       work_nl(:,ii)=work_nl(:,ii)+bwork1_vectemp1(:,1)*2.0d0


       if (mod(kk,500)==0) print *,kk
       enddo

            

        end do

    end select

end do


do ii=1,M0_nl
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo

print *,'======================'

do ii=1,M_nl
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo




print *,'+++++++++++++++++++++++++++++++++++'

call quicksort_eigenpairs(E_nl,X_nl,res_nl,Nn,1,M0_nl)

do ii=1,M0_nl
    print *,E_nl(ii)*hartree,res_nl(ii)
enddo

print *,'+++++++++++++++++++++++++++++++++++'


do i=1,M0_nl

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    print *,dot_product(dble(X_nl(:,i)),snq1)

    vtemp=dot_product(dble(X_nl(:,i)),snq1)

    X_nl(:,i)=X_nl(:,i)/sqrt(vtemp)

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,i)),snq1)
    print *,dot_product(dble(X_nl(:,i)),snq1)


    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    ! print *,dot_product(dble(X_nl_temp(:,i)),snq1)

    vtemp=dot_product(dble(X_nl_temp(:,i)),snq1)

    X_nl_temp(:,i)=X_nl_temp(:,i)/sqrt(vtemp)

    call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl_temp(:,i)),snq1)
    print *,dot_product(dble(X_nl_temp(:,i)),snq1)

enddo


! call mkl_dcsrgemv('N',Nn,B,IA,JA,dble(X_nl(:,Nstates+1)),snq1)
! print *,dot_product(dble(X_nl(:,Nstates+1)),snq1)

! vtemp=dot_product(dble(X_nl(:,Nstates+1)),snq1)

! X_nl(:,Nstates+1)=X_nl(:,Nstates+1)/sqrt(vtemp)



! call mkl_dcsrgemv('N',Nn,B,IA,JA,(abs(psi(:,Nstates))-abs(dble(X_nl(:,Nstates)))),snq1)
! print *,dot_product((abs(psi(:,Nstates))-abs(dble(X_nl(:,Nstates)))),snq1)

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(abs(psi(:,Nstates+1))-abs(dble(X_nl(:,Nstates+1)))),snq1)
! print *,dot_product((abs(psi(:,Nstates+1))-abs(dble(X_nl(:,Nstates+1)))),snq1)


do i=1,M0_nl

    print *,'========================='

call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,cd_igorbital0+(i-1))+dble(X_nl(:,i))),snq1)
print *,dot_product((psi(:,cd_igorbital0+(i-1))+dble(X_nl(:,i))),snq1)
 
call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,cd_igorbital0+(i-1))-dble(X_nl(:,i))),snq1)
print *,dot_product((psi(:,cd_igorbital0+(i-1))-dble(X_nl(:,i))),snq1)

    print *,'~~~~~~~~~~~~~~~~~~~'

call mkl_dcsrgemv('N',Nn,B,IA,JA,(dble(X_nl_temp(:,cd_igorbital0+(i-1)))+dble(X_nl(:,i))),snq1)
print *,dot_product((dble(X_nl_temp(:,cd_igorbital0+(i-1)))+dble(X_nl(:,i))),snq1)

call mkl_dcsrgemv('N',Nn,B,IA,JA,(dble(X_nl_temp(:,cd_igorbital0+(i-1)))-dble(X_nl(:,i))),snq1)
print *,dot_product((dble(X_nl_temp(:,cd_igorbital0+(i-1)))-dble(X_nl(:,i))),snq1)

    ! print *,'========================='

enddo

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,Nstates+1)+dble(X_nl(:,Nstates+1))),snq1)
! print *,dot_product((psi(:,Nstates+1)+dble(X_nl(:,Nstates+1))),snq1)


! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,Nstates)-dble(X_nl(:,Nstates))),snq1)
! print *,dot_product((psi(:,Nstates)-dble(X_nl(:,Nstates))),snq1)

! call mkl_dcsrgemv('N',Nn,B,IA,JA,(psi(:,Nstates+1)-dble(X_nl(:,Nstates+1))),snq1)
! print *,dot_product((psi(:,Nstates+1)-dble(X_nl(:,Nstates+1))),snq1)




do k=1,M0_nl!Nstates+1
write (file_id, '(I0)') k
! open (unit=20,file='CO.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='CO_0.6415.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='He.1_p2_nessie.psi'//trim(file_id),status='old')
open(unit=20,file=trim(name)//'_p2_gw.psi'//trim(file_id),status='replace')
! open (unit=20,file='LiH.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2.1_p2_nessie.psi'//trim(file_id),status='old')
! open (unit=20,file='H2O.1_p2_nessie.psi'//trim(file_id),status='old')
do i=1,Nn
    write(20,*) X_nl(i,k)
enddo
close (20)
enddo


! stop
















end if !!!! end casida_method selection


end if !!!! end selection CD/Casida


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    GW ends                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
























































    ! print *,"***save node points files in Angstrom"
    ! open(10,file=trim(name)//'.points',status='replace')
    ! do i=1,Nn
    !     write(10,*) point(i,:)*bohr
    ! enddo
    ! close(10)


    
!     ! WRITE(indent, '(A)') ' '
!     ! WRITE(strlev, '(I2)') digits(Nn)
!     ! WRITE(*, '(A)') indent//TRIM(ADJUSTL(strlev))//': this is my level'

!     ! open(10,file=trim(name)//'_p2.node',status='replace')
!     ! write(10,*) Nn,3,0,1
!     ! do i=1,Nn
!     !     write(10,'(A)') indent//TRIM(ADJUSTL(strlev))//i,real(point(i,:)*bohr),color(i)
!     ! enddo
!     ! close(10)


!     ! open(10,file=trim(name)//'_p2.ele',status='replace')
!     ! write(10,*) Ne
!     ! do i=1,Nn
!     !     write(10,fmt='(i4.2)','(A)') indent//TRIM(ADJUSTL(strlev))//i,ele(i,:)
!     ! enddo
!     ! close(10)


!     print *,"***save wavefunction in atomic units"

!     do k=1,Nstates
!         write (file_id, '(I0)') k
!         open(10,file=trim(name)//'.psi'//trim(adjustl(file_id)),status='replace')
!         do i=1,Nn
!             write(10,*) psi(i,k)
!         enddo
!         close(10)
!     enddo

!     print *,"***save electron density file in atomic units"

!     open(10,file=trim(name)//'.n',status='replace')
!     do i=1,Nn
!         do k=1,Nstates
!             write(10,*) 2.0d0*psi(i,k)**2
!         enddo
!     enddo
!     close(10)

!     print *,"***save orbital energy in atomic units"

!     open(10,file=trim(name)//'.orbitalenergy',status='replace')
!     do k=1,Nstates
!         write(10,*) k,E(k)
!     enddo
!     close(10)



!!!! since orbitals and energy states are obtained, then total energy and other quatities can also be calculated.


    print *,">>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<"
    print *," ---------------------------------------- "
    print *,"|  Post-processing is to be continued..  |"
    print *,"|                                        |"
    print *,"|                                        |"
    print *,"|               Stay Tuned!              |"
    print *,"|                                        |"
    print *," ---------------------------------------- "





end program
