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
    if (iarg/=1) then
        print *,'only 1 argument possible'
        stop
    end if
    call get_command_argument(1,name,clen,stat) ! name of the file



    call load_gw(name,dgr_fem,meanfield,gw_sigma)




    !!!! Read atoms from xyz file
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


        i=at(k)%core
        write(charI,"(I5)"), i
        open(12,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3n',status='old')
        read(12,*) Nmuffin  
        allocate(point_muffin(1:Nmuffin,1:3))
        do ii=1,Nmuffin
        read(12,*) dummyc,point_muffin(ii,1),point_muffin(ii,2),point_muffin(ii,3)
        enddo
        close(12)


        ro=1.0d0!0.85d0!1.0d0

        do i=1,Nmuffin
            !!! according to atom number place different radius of muffin
            write(10,20) 8+Nmuffin_total+i-1,point_muffin(i,:)*ro+at(k)%c(:),0 
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

        k=at(i)%core
        write(charI,"(I5)"), k
        open(12,file='../../database_muffins/at_'//trim(adjustl(charI))//'.3f',status='old')
        read(12,*) Nface
        allocate(face_muffin(1:Nface,1:3))
        do ii=1,Nface
        read(12,*) dummyc,face_muffin(ii,1),face_muffin(ii,2),face_muffin(ii,3)


        write(10,*) 1, 0, 0
        write(10,*) 3,face_muffin(ii,1)+8+Nmuffin_total,face_muffin(ii,2)+8+Nmuffin_total,face_muffin(ii,3)+8+Nmuffin_total

        
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
    print *,line
    call execute_command_line(line)



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

    allocate(psi_0(1:Nn))

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

            psi_0(i)=psi_0(i)+((n_init(ub,ll)-n_init(lb,ll))/(r_init(ub,ll)-r_init(lb,ll))*(rxyz*bohr-r_init(lb,ll))&
            +n_init(lb,ll))*bohr**3

        enddo
    enddo



    if (meanfield=="hf") then


        allocate(psi_i(1:Nn,1:Nstates+1))
        do k=1,Nstates+1
        write (file_id, '(I0)') k
        open (unit=20,file=trim(name)//'_p2_dft_1.psi'//trim(file_id),status='old')
        do i=1,Nn
            read(20,*) psi_i(i,k)
        enddo
        close (20)
        enddo
        

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



allocate(neigh_NnToNg(1:Ne*Nquadrature))
call construct_NnToNg(Ne,Nquadrature,Nlocal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!     FEAST eigensolver parameters setting      !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

M0=ks_M00
Emin=ks_Emin0
Emax=ks_Emax0


!!! Allocate memory for eigenvalues. eigenvectors, residual 
allocate(E(M0), psi(Nn,M0), res(M0))

allocate(E_hf(Nn))







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
alpha0 = 1.0d0
if (meanfield=="dft") then
    eps = ks_eps0!1.0d-8!1.0d-8!1E-8
    alpha = ks_alpha0!0.5d0!0.45d0
else if (meanfield=="hf") then
    eps = hf_eps0
    alpha = hf_alpha0
end if
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

    print *,'Contour Deformation method for self-energy calculation will be included in the next release.'

    stop

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



deallocate(psi_0)
allocate(temp2(Nn))


if (meanfield=="dft") then



print *, '***'
print *," ------------------------------------------"
print *,"|           SCF-iteration starts           |"
print *," ------------------------------------------"
print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'




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

        call hartreepotential(Nn,Ne,Nquadrature,Nstates,Nbc,2.0d0,it,0)
        call exchangepotential(Nn)
        call correlationpotential(Nn)


        V_total(1:Nn) = Vh(1:Nn)+Vx(1:Nn)+Vc(1:Nn)!+Vi(1:Nn)

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


    end if !!! LDA/PBE 'if' ends


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! build total Hamiltonian !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    H=0.5d0*A+V_i+V

    ! call MPI_BARRIER(MPI_COMM_WORLD ,code)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!    solve Kohn-Sham equation using FEAST    !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    do i=1,Nstates
        nq(:)=nq(:)+2.0d0*psi(:,i)**2
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

            print *,"-----  Orbital,  E_kinetic (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,0.5d0*A,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            print *,"-----  Orbital,  E_ext (eV) -----"

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
            print *,"-----  Orbital,  E_h (eV) -----"

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
            print *,"-----  Orbital,  E_x (eV) -----"

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


            print *,"-----  Orbital,      E_c (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

        else if (xc0=='pbe') then

            

        end if !!! lda/pbe 
        
        print *,"=============================================="
        



        if (gw_sigma=='cd') then 

            print *,'--- solve more unoccupied states ---'

            M0=2*Nstates+25
            ! Emin=-40.0d0
            Emax=0.25d0
            deallocate(E,psi,res)
            allocate(E(M0), psi(Nn,M0), res(M0))


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

        end if

        if (gw_sigma=='casida') then 

            print *,'--- solve all unoccupied states ---'


            M0=Nn
            ! Emin=-50.0d0
            Emax=2.0d7
            deallocate(E,psi,res)
            allocate(E(M0), psi(Nn,M0), res(M0))

            call feastinit(fpm)
            call dfeast_scsrgv(UPLO, Nn, H, IA, JA, B, IA, JA, fpm, epsout, loop, Emin, Emax, M0, E, psi, M00, res, info)

            if (info/=0) then
                print *,'FEAST_error info --------',info
                stop
            end if
            
            ! print *,M00,E(M00)

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

        call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,2.0d0,hf_poisson_bc0,0,0)
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
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!    solve Hartree-Fock equation using DSYGV    !!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        !!!!! solve eigenvalue problem
        call dsygv(1,'V','U',Nn,H_dense_psi,Nn,B_dense,Nn,E_hf,DSYEV_work,DSYEV_lwork,info)
    
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

            ! do i=1,Nstates+40
            !     print *, i,E_hf(i)*hartree
            ! enddo

            deallocate(E,psi)
            allocate(E(Nn),psi(Nn,Nn))

            E(1:Nn)=E_hf(1:Nn)
            psi(:,1:Nn)=H_dense_psi(:,1:Nn)


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!     compute Eh,Ex,Ec        !!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(snq2(1:Nn))

        print *,"=============================================="


            print *,"-----  Orbital,  E_kinetic (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,0.5d0*A,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo

            print *,"-----  Orbital,  E_ext (eV) -----"

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

            print *,"-----  Orbital,  E_h (eV) -----"

            do i=1,Nstates+1
        
            call mkl_dcsrgemv('N',Nn,V,IA,JA,psi(:,i),snq2)
            Eh=dot_product(snq2,psi(:,i))

            print *,i,Eh*hartree

            enddo


            !!!!!!!!!! Ex !!!!!!!!!!

            print *,"-----  Orbital,  E_x (eV) -----"

            do i=1,Nstates+1
            
                call DSYMM('L','L',Nn,1,1.0d0,Vhfx_final,Nn,psi(:,i),Nn,0.0d0,snq2,Nn)
                
                print *,i,dot_product(psi(:,i),snq2)*hartree
            
            enddo

        
        print *,"=============================================="

























            call psi_interpolation_HF(Ne,Nquadrature,Nstates,1)
            call set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,2.0d0,hf_poisson_bc0,0,0)
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

    deallocate(DSYEV_work)

end if  !!!!! end DFT/HF selection

if (gw_sigma=='cd') then

    Nempty = Nn-Nstates

else if (meanfield=='casida') then

    Nempty = casida_N_empty0!Nn-Nstates
    
end if



do k=1,Nstates+1
write (file_id, '(I0)') k
open(unit=20,file=trim(name)//'_'//trim(dgr_fem)//'_'//trim(meanfield)//'_1.psi'//trim(file_id),status='replace')
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

    call psi_interpolation(Ne,Nquadrature,Nstates,1,degree)

else if (gw_sigma=='casida') then

    call psi_interpolation(Ne,Nquadrature,Nn,1,degree)

end if




print *, "*** construct bare Coulomb potential -- 1/|r-r'| ***"

allocate(v_kernel(Ne*Nquadrature,Nn))
call construct_v_kernel(Nn,Ne,Nquadrature)

do i=1,Nn
    v_kernel(:,i)=v_kernel(:,i)*volumegweight(:)
enddo

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
print *,"----  Orbital,  GW <psi|Sigma_x|psi> (eV) ----"

do i=1,Nstates+1

    orbital=i

! CALL DGEMM('N','N',Nn,1,Nn,1.0d0,Vhfx_final_GW,Nn,psi(:,HOMO_state),Nn,0.0d0,snq1,Nn)
call DSYMM('L','L',Nn,1,1.0d0,Vhfx_final_GW,Nn,psi(:,orbital),Nn,0.0d0,snq1,Nn)

print *,i,dot_product(psi(:,orbital),snq1)*hartree

enddo

print *,"=============================================="

deallocate(volume)

100 format (I4,I4,ES25.16,ES25.16,ES25.16)


! stop






































if (gw_sigma=='cd') then





else if (gw_sigma=='casida') then



    if (casida_method0=='gs') then

        print *,'*** compute Casida Kx ***'

    Nempty=casida_N_empty0


allocate(psi_ia_g_real(Ne*Nquadrature,Nstates*Nempty))
do k=1,Nstates
    do kk=1,Nempty
        psi_ia_g_real(:,(k-1)*Nempty+kk)=psi_point_g(:,k)*psi_point_g(:,Nstates+kk)!*volumegweight(:)
    enddo
enddo








allocate(NgNntemp(1:Ne*Nquadrature,1:Nstates*Nempty))
allocate(casida_Kx(1:Nstates*Nempty,1:Nstates*Nempty))
allocate(casida_R(1:Nstates*Nempty))
allocate(casida_Xs(1:Nstates*Nempty,1:Nstates*Nempty))
allocate(casida_omega(1:Nstates*Nempty))
allocate(casida_R_half(1:Nstates*Nempty))
allocate(casida_Kxnj(1:Nempty,1:Nstates*Nempty))
allocate(casida_Kxnj_occupied(1:Nstates,1:Nstates*Nempty))
allocate(casida_cvomega(1:Nstates*Nempty))

allocate(casida_Ys(1:Nstates*Nempty,1:Nstates*Nempty))


allocate(NNmatrix_temp12(1:Nn,1:Nstates*Nempty))
allocate(NNmatrix_temp13(1:Nn,1:Nstates*Nempty))





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





!!!!!!!!! allocate for 'compute_gw_Sigma_Casida'/'compute_gw_Casida_Pi'/'compute_gw_casida_Ec'/'compute_gw_casida_Ec_complex' !!!


Nguess_complex=casida_N_grid0
allocate(E_guess_complex_HOMO(1:Nguess_complex))
allocate(E_guess_complex_LUMO(1:Nguess_complex))
allocate(casida_Ec(1:Nguess_complex))



do ll=1,Nguess_complex

    E_guess_complex_HOMO(ll)=casida_Emin_grid10/hartree+(ll-1)*(casida_Emax_grid10/hartree-casida_Emin_grid10/hartree)&
                                /(Nguess_complex-1)

    E_guess_complex_LUMO(ll)=casida_Emin_grid20/hartree+(ll-1)*(casida_Emax_grid20/hartree-casida_Emin_grid20/hartree)&
                                /(Nguess_complex-1)

enddo














!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   compute casida_Kxnj   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(psi_n_g_real(Ne*Nquadrature,Nempty))
allocate(psi_n_g_real_occupied(Ne*Nquadrature,Nstates))


orbital = casida_orbital10

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
print *,'=========='
print *,'*** solve orbital',casida_orbital10,' -- Sigma_c  ***'

call compute_gw_casida_Ec_complex(E_guess_complex_HOMO,E,Nguess_complex,Nstates,Nempty)




print *,'  energy grid (eV)','     GW <psi|Sigma_c|psi> (eV)'
do m=1,Nguess_complex
    print *,dble(E_guess_complex_HOMO(m))*hartree,dble(casida_Ec(m))*hartree
enddo





orbital = casida_orbital20

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


print *,'*** solve orbital',casida_orbital20,' -- Sigma_c  ***'

call compute_gw_casida_Ec_complex(E_guess_complex_LUMO,E,Nguess_complex,Nstates,Nempty)



print *,'  energy grid (eV)','     GW <psi|Sigma_c|psi> (eV)'
do m=1,Nguess_complex
    print *,dble(E_guess_complex_LUMO(m))*hartree,dble(casida_Ec(m))*hartree
enddo


! stop




! 9999 continue











































else if (casida_method0=='nlfeast') then


print *,'The upcoming version will include the release of the FEAST nonlinear eigenvalue algorithm.'

stop



















end if !!!! end casida_method selection


end if !!!! end selection CD/Casida


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    GW ends                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





end program
