module basisfunctions

    use class_linkedlist
    use tools

    implicit none

    public

    type :: matrix
        integer,dimension(:),allocatable :: col
        complex(kind=(kind(1.0d0))),dimension(:),allocatable :: A,S
    end type

    type :: psi_gauss
        double precision,dimension(:),allocatable :: state_value,location
    end type

    !!!! Basis !!!!
    double precision, dimension(:,:), allocatable :: point,point_g
    integer, dimension(:), allocatable :: color,color_prime,color_new,color_new1
    integer, dimension(:,:), allocatable :: ele,ele_prime
    double precision, dimension(:), allocatable :: phi_p2,phi_p3
    double precision, dimension(:,:), allocatable :: phi_p2_del,phi_p3_del
    ! double precision, dimension(1:3,1:3) :: J0,Jit,  J1,Jit1,  J2,Jit2, J20,Jit20
    double precision, allocatable :: Jdet,Jdet1,Jdet2,Jdet20
    double precision, dimension(:,:), allocatable :: p1_matrix,p2_matrix,p3_matrix

    double precision, dimension(:), allocatable :: efem_phi
    double precision, dimension(:,:), allocatable :: efem_phi_grad,efem_1s_grad,h_grad
    integer,dimension(:),allocatable :: enriched_node,enriched
    integer,dimension(:,:),allocatable :: enriched_ele

    double precision, dimension(:,:), allocatable :: pointex,pointey,pointez,pointfx,pointfy,pointfz
    ! double precision, dimension(1:4,1:4) :: p1_matrix
    ! double precision, dimension(1:10,1:10) :: p2_matrix
    ! double precision, dimension(1:20,1:20) :: p3_matrix
    ! double precision, dimension(1:11,1:3) :: gpoint
    ! double precision, dimension(1:11) :: gweight

    double precision, allocatable :: psi_square_100,psi_square_200,psi_square_21x,psi_square_21y,psi_square_21z
    double precision, allocatable :: psi_square_300,psi_square_31x,psi_square_31y,psi_square_31z,d

    double precision, dimension(:,:), allocatable :: J0,Jit,J1,Jit1,J2,Jit2,J20,Jit20,localpointp2
    

    double precision, dimension(:,:), allocatable :: gpoint,gpoint5,gpoint10,gpoint14,gpoint15,gpoint24,gpoint31,gpoint45
    double precision, dimension(:), allocatable :: gweight,gweight5,gweight10,gweight14,gweight15,gweight24,gweight31,gweight45
    

    ! double precision, dimension(1:10,1:3) :: localpointp2
    ! double precision, dimension(1:20,1:3) :: localpointp3

    double precision,dimension(:),allocatable :: volume
    type(psi_gauss),dimension(:,:),allocatable :: psi_g,psi_g0,nvi_g,psi_g1
    double precision, dimension(:,:), allocatable :: psi,psi_1,psi_2,psi_i,psi_GW,psi_GW_g,psi_GW0,psi_gauss0
    complex(kind=(kind(1.0d0))), dimension(:,:), allocatable :: psi_complex

    ! double precision,dimension(:),allocatable :: E_GW

    double precision,dimension(:,:),allocatable :: psi_point_g

    ! complex(kind=(kind(1.0d0))),dimension(:,:), allocatable :: psi_GW,psi_GW_g

    double precision, dimension(:), allocatable :: A,B,A00,usa,usb,H_dft,usa_dft,quaternion_A
    complex(kind=(kind(1.0d0))) ,dimension(:),allocatable :: psaz_dft
    double precision, dimension(:,:), allocatable :: Btemp
    integer, dimension(:), allocatable :: IA, JA, IB, JB,IA0, JA0,IA_pntrb,IA_pntre,uisa,ujsa,quaternion_IA,quaternion_JA
    complex(kind=(kind(1.0d0))), dimension(:), allocatable :: B_complex,H_complex
    complex(kind=(kind(1.0d0))), dimension(:,:), allocatable :: C

    double precision, dimension(:,:), allocatable :: H_dense,S_dense,B_dense
    double precision, dimension(:,:), allocatable :: H_dense_psi

    double precision, dimension(:,:), allocatable :: Hhf_dense


    !!!!!!!!!!!!!!!! Hartree Fock variables !!!!!!!!!!!!!!
    double precision, dimension(:,:), allocatable :: psii,psig
    double precision, dimension(:,:,:), allocatable :: psi_previous,error_previous
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type(LinkedList) ,dimension(:),allocatable:: neigh_AB,neigh_V,neigh_AB0,neigh_C,quaternion
    type(matrix), dimension(:), allocatable :: row,row0
    integer,dimension(:),allocatable :: BC

    double precision,dimension(:,:),allocatable :: point_mesh,point_mesh1
    integer, dimension(:,:), allocatable :: ele_mesh

    double precision,dimension(:),allocatable :: nq,ni,nq_g,ni_xc,ni_g
    double precision,dimension(:,:),allocatable :: ni_previous

    double precision,dimension(:,:),allocatable :: nq_gradient,nq_g_gradient
    double precision,dimension(:,:,:),allocatable :: psi_gradient

    integer,dimension(:,:),allocatable :: node_arranged,node_arranged_rest

    double precision,dimension(:,:),allocatable :: cxyz,cxyz_s
    ! integer :: Nxyz,Nxyz_s

    double precision,dimension(:,:),allocatable :: matBC

    CHARACTER(len=:), dimension(:), allocatable :: matdescrb
    CHARACTER(len=:) ,dimension(:), allocatable :: matdescra


    double precision,dimension(:),allocatable :: volumegweight
    type(LinkedList) ,dimension(:),allocatable:: neigh_NnToNg
    double precision, dimension(:), allocatable :: NnToNg,NnToNg_temp
    integer, dimension(:), allocatable :: NnToNg_IA,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre
    


    contains



    subroutine jacobian(e1,e2,e3,e4)

        !! for LAPACK  
        integer,intent(in) :: e1,e2,e3,e4!,lwork
        integer :: info,n
        double precision,dimension(:),allocatable :: work, ipiv
    
        J0(1,1) = point(e2,1)-point(e1,1)
        J0(1,2) = point(e3,1)-point(e1,1)
        J0(1,3) = point(e4,1)-point(e1,1)
        J0(2,1) = point(e2,2)-point(e1,2)
        J0(2,2) = point(e3,2)-point(e1,2)
        J0(2,3) = point(e4,2)-point(e1,2)
        J0(3,1) = point(e2,3)-point(e1,3)
        J0(3,2) = point(e3,3)-point(e1,3)
        J0(3,3) = point(e4,3)-point(e1,3)
    
        Jdet=J0(1,1)*(J0(2,2)*J0(3,3)-J0(2,3)*J0(3,2))&
                 -J0(1,2)*(J0(2,1)*J0(3,3)-J0(2,3)*J0(3,1))&
                 +J0(1,3)*(J0(2,1)*J0(3,2)-J0(2,2)*J0(3,1))
    
        allocate(work(3))
        allocate(ipiv(3))
        Jit = J0
        n = 3
    
        call DGETRF(n, n, Jit, n, ipiv, info)
        call DGETRI(n, Jit, n, ipiv, work, n, info)
        
        Jit = transpose(Jit)

        deallocate(work)
        deallocate(ipiv)
    
    end subroutine jacobian








subroutine jacobian1(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)

    !! for LAPACK  
    integer :: info,n!,lwork
    double precision,intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    double precision,dimension(:),allocatable :: work, ipiv
  
    J1(1,1) = x2-x1 ! point(e2,1)-point(e1,1)
    J1(1,2) = x3-x1 ! point(e3,1)-point(e1,1)
    J1(1,3) = x4-x1 ! point(e4,1)-point(e1,1)
    J1(2,1) = y2-y1 ! point(e2,2)-point(e1,2)
    J1(2,2) = y3-y1 ! point(e3,2)-point(e1,2)
    J1(2,3) = y4-y1 ! point(e4,2)-point(e1,2)
    J1(3,1) = z2-z1 ! point(e2,3)-point(e1,3)
    J1(3,2) = z3-z1 ! point(e3,3)-point(e1,3)
    J1(3,3) = z4-z1 ! point(e4,3)-point(e1,3)
  
    Jdet1=J1(1,1)*(J1(2,2)*J1(3,3)-J1(2,3)*J1(3,2))&
             -J1(1,2)*(J1(2,1)*J1(3,3)-J1(2,3)*J1(3,1))&
             +J1(1,3)*(J1(2,1)*J1(3,2)-J1(2,2)*J1(3,1))
  
    allocate(work(3))
    allocate(ipiv(3))
    Jit1 = J1
    n = 3
  
    call DGETRF(n, n, Jit1, n, ipiv, info)
    call DGETRI(n, Jit1, n, ipiv, work, n, info)
    
    Jit1 = transpose(Jit1)

    deallocate(work)
    deallocate(ipiv)
  
  end subroutine jacobian1

  subroutine jacobian2(e1,e2,e3,e4)

    !! for LAPACK  
    integer,intent(in) :: e1,e2,e3,e4!,lwork
    integer :: info,n
    double precision,dimension(:),allocatable :: work, ipiv
  
    J2(1,1) = point_mesh(e2,1)-point_mesh(e1,1)
    J2(1,2) = point_mesh(e3,1)-point_mesh(e1,1)
    J2(1,3) = point_mesh(e4,1)-point_mesh(e1,1)
    J2(2,1) = point_mesh(e2,2)-point_mesh(e1,2)
    J2(2,2) = point_mesh(e3,2)-point_mesh(e1,2)
    J2(2,3) = point_mesh(e4,2)-point_mesh(e1,2)
    J2(3,1) = point_mesh(e2,3)-point_mesh(e1,3)
    J2(3,2) = point_mesh(e3,3)-point_mesh(e1,3)
    J2(3,3) = point_mesh(e4,3)-point_mesh(e1,3)
  
    Jdet2=abs(J2(1,1)*(J2(2,2)*J2(3,3)-J2(2,3)*J2(3,2))&
             -J2(1,2)*(J2(2,1)*J2(3,3)-J2(2,3)*J2(3,1))&
             +J2(1,3)*(J2(2,1)*J2(3,2)-J2(2,2)*J2(3,1)))
  
    allocate(work(3))
    allocate(ipiv(3))
    Jit2 = J2
    n = 3
  
    call DGETRF(n, n, Jit2, n, ipiv, info)
    call DGETRI(n, Jit2, n, ipiv, work, n, info)
    
    Jit2 = transpose(Jit2)

    deallocate(work)
    deallocate(ipiv)
  
  end subroutine jacobian2

  subroutine jacobian20(e1,e2,e3,e4)

    !! for LAPACK  
    integer,intent(in) :: e1,e2,e3,e4!,lwork
    integer :: info,n
    double precision,dimension(:),allocatable :: work, ipiv
  
    J20(1,1) = point_mesh(e2,1)-point_mesh(e1,1)
    J20(1,2) = point_mesh(e3,1)-point_mesh(e1,1)
    J20(1,3) = point_mesh(e4,1)-point_mesh(e1,1)
    J20(2,1) = point_mesh(e2,2)-point_mesh(e1,2)
    J20(2,2) = point_mesh(e3,2)-point_mesh(e1,2)
    J20(2,3) = point_mesh(e4,2)-point_mesh(e1,2)
    J20(3,1) = point_mesh(e2,3)-point_mesh(e1,3)
    J20(3,2) = point_mesh(e3,3)-point_mesh(e1,3)
    J20(3,3) = point_mesh(e4,3)-point_mesh(e1,3)
  
    Jdet20=abs(J20(1,1)*(J20(2,2)*J20(3,3)-J20(2,3)*J20(3,2))&
             -J20(1,2)*(J20(2,1)*J20(3,3)-J20(2,3)*J20(3,1))&
             +J20(1,3)*(J20(2,1)*J20(3,2)-J20(2,2)*J20(3,1)))
  
    allocate(work(3))
    allocate(ipiv(3))
    Jit20 = J20
    n = 3
  
    call DGETRF(n, n, Jit20, n, ipiv, info)
    call DGETRI(n, Jit20, n, ipiv, work, n, info)
    
    Jit20 = transpose(Jit20)

    deallocate(work)
    deallocate(ipiv)
  
  end subroutine jacobian20
  

  subroutine gaussian_integral()

    !!!! the Keast Rule -- https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html

    gpoint(1,:)   = (/0.2500000000000000d0,  0.2500000000000000d0,  0.2500000000000000d0/)
    gpoint(2,:)   = (/0.7857142857142857d0,  0.0714285714285714d0,  0.0714285714285714d0/)
    gpoint(3,:)   = (/0.0714285714285714d0,  0.0714285714285714d0,  0.0714285714285714d0/)
    gpoint(4,:)   = (/0.0714285714285714d0,  0.0714285714285714d0,  0.7857142857142857d0/)
    gpoint(5,:)   = (/0.0714285714285714d0,  0.7857142857142857d0,  0.0714285714285714d0/)
    gpoint(6,:)   = (/0.1005964238332008d0,  0.3994035761667992d0,  0.3994035761667992d0/)
    gpoint(7,:)   = (/0.3994035761667992d0,  0.1005964238332008d0,  0.3994035761667992d0/)
    gpoint(8,:)   = (/0.3994035761667992d0,  0.3994035761667992d0,  0.1005964238332008d0/)
    gpoint(9,:)   = (/0.3994035761667992d0,  0.1005964238332008d0,  0.1005964238332008d0/)
    gpoint(10,:)  = (/0.1005964238332008d0,  0.3994035761667992d0,  0.1005964238332008d0/)
    gpoint(11,:)  = (/0.1005964238332008d0,  0.1005964238332008d0,  0.3994035761667992d0/)

    gweight(1)   = -0.0789333333333333d0
    gweight(2)   =  0.0457333333333333d0
    gweight(3)   =  0.0457333333333333d0
    gweight(4)   =  0.0457333333333333d0
    gweight(5)   =  0.0457333333333333d0
    gweight(6)   =  0.1493333333333333d0
    gweight(7)   =  0.1493333333333333d0
    gweight(8)   =  0.1493333333333333d0
    gweight(9)   =  0.1493333333333333d0
    gweight(10)  =  0.1493333333333333d0
    gweight(11)  =  0.1493333333333333d0

end subroutine gaussian_integral




subroutine gaussian_integral5()

  !!!! the Keast Rule -- https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html

    gpoint5(1,:) = (/0.2500000000000000d0, 0.2500000000000000d0, 0.2500000000000000d0/)
    gpoint5(2,:) = (/0.5000000000000000d0, 0.1666666666666667d0, 0.1666666666666667d0/)
    gpoint5(3,:) = (/0.1666666666666667d0, 0.1666666666666667d0, 0.1666666666666667d0/)
    gpoint5(4,:) = (/0.1666666666666667d0, 0.1666666666666667d0, 0.5000000000000000d0/)
    gpoint5(5,:) = (/0.1666666666666667d0, 0.5000000000000000d0, 0.1666666666666667d0/)

    gweight5(1) = -0.8000000000000000d0
    gweight5(2) =  0.4500000000000000d0
    gweight5(3) =  0.4500000000000000d0
    gweight5(4) =  0.4500000000000000d0
    gweight5(5) =  0.4500000000000000d0

end subroutine gaussian_integral5

subroutine gaussian_integral10()

    gpoint10(1 ,:) = (/0.5684305841968444d0, 0.1438564719343852d0, 0.1438564719343852d0/)
    gpoint10(2 ,:) = (/0.1438564719343852d0, 0.1438564719343852d0, 0.1438564719343852d0/)
    gpoint10(3 ,:) = (/0.1438564719343852d0, 0.1438564719343852d0, 0.5684305841968444d0/)
    gpoint10(4 ,:) = (/0.1438564719343852d0, 0.5684305841968444d0, 0.1438564719343852d0/)
    gpoint10(5 ,:) = (/0.0000000000000000d0, 0.5000000000000000d0, 0.5000000000000000d0/)
    gpoint10(6 ,:) = (/0.5000000000000000d0, 0.0000000000000000d0, 0.5000000000000000d0/)
    gpoint10(7 ,:) = (/0.5000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0/)
    gpoint10(8 ,:) = (/0.5000000000000000d0, 0.0000000000000000d0, 0.0000000000000000d0/)
    gpoint10(9 ,:) = (/0.0000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0/)
    gpoint10(10,:) = (/0.0000000000000000d0, 0.0000000000000000d0, 0.5000000000000000d0/)

    gweight10(1) = 0.2177650698804054d0
    gweight10(2) = 0.2177650698804054d0
    gweight10(3) = 0.2177650698804054d0
    gweight10(4) = 0.2177650698804054d0
    gweight10(5) = 0.0214899534130631d0
    gweight10(6) = 0.0214899534130631d0
    gweight10(7) = 0.0214899534130631d0
    gweight10(8) = 0.0214899534130631d0
    gweight10(9) = 0.0214899534130631d0
    gweight10(10) = 0.0214899534130631d0

  end subroutine gaussian_integral10

  subroutine gaussian_integral14()

    gpoint14(1 ,:) = (/0.0000000000000000d0, 0.5000000000000000d0, 0.5000000000000000d0/)
    gpoint14(2 ,:) = (/0.5000000000000000d0, 0.0000000000000000d0, 0.5000000000000000d0/)
    gpoint14(3 ,:) = (/0.5000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0/)
    gpoint14(4 ,:) = (/0.5000000000000000d0, 0.0000000000000000d0, 0.0000000000000000d0/)
    gpoint14(5 ,:) = (/0.0000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0/)
    gpoint14(6 ,:) = (/0.0000000000000000d0, 0.0000000000000000d0, 0.5000000000000000d0/)
    gpoint14(7 ,:) = (/0.6984197043243866d0, 0.1005267652252045d0, 0.1005267652252045d0/)
    gpoint14(8 ,:) = (/0.1005267652252045d0, 0.1005267652252045d0, 0.1005267652252045d0/)
    gpoint14(9 ,:) = (/0.1005267652252045d0, 0.1005267652252045d0, 0.6984197043243866d0/)
    gpoint14(10,:) = (/0.1005267652252045d0, 0.6984197043243866d0, 0.1005267652252045d0/)
    gpoint14(11,:) = (/0.0568813795204234d0, 0.3143728734931922d0, 0.3143728734931922d0/)
    gpoint14(12,:) = (/0.3143728734931922d0, 0.3143728734931922d0, 0.3143728734931922d0/)
    gpoint14(13,:) = (/0.3143728734931922d0, 0.3143728734931922d0, 0.0568813795204234d0/)
    gpoint14(14,:) = (/0.3143728734931922d0, 0.0568813795204234d0, 0.3143728734931922d0/)

    gweight14(1 ) = 0.0190476190476190d0
    gweight14(2 ) = 0.0190476190476190d0
    gweight14(3 ) = 0.0190476190476190d0
    gweight14(4 ) = 0.0190476190476190d0
    gweight14(5 ) = 0.0190476190476190d0
    gweight14(6 ) = 0.0190476190476190d0
    gweight14(7 ) = 0.0885898247429807d0
    gweight14(8 ) = 0.0885898247429807d0
    gweight14(9 ) = 0.0885898247429807d0
    gweight14(10) = 0.0885898247429807d0
    gweight14(11) = 0.1328387466855907d0
    gweight14(12) = 0.1328387466855907d0
    gweight14(13) = 0.1328387466855907d0
    gweight14(14) = 0.1328387466855907d0

  end subroutine gaussian_integral14

  subroutine gaussian_integral15()

    gpoint15(1 ,:) = (/0.2500000000000000d0, 0.2500000000000000d0, 0.2500000000000000d0/)
    gpoint15(2 ,:) = (/0.0000000000000000d0, 0.3333333333333333d0, 0.3333333333333333d0/)
    gpoint15(3 ,:) = (/0.3333333333333333d0, 0.3333333333333333d0, 0.3333333333333333d0/)
    gpoint15(4 ,:) = (/0.3333333333333333d0, 0.3333333333333333d0, 0.0000000000000000d0/)
    gpoint15(5 ,:) = (/0.3333333333333333d0, 0.0000000000000000d0, 0.3333333333333333d0/)
    gpoint15(6 ,:) = (/0.7272727272727273d0, 0.0909090909090909d0, 0.0909090909090909d0/)
    gpoint15(7 ,:) = (/0.0909090909090909d0, 0.0909090909090909d0, 0.0909090909090909d0/)
    gpoint15(8 ,:) = (/0.0909090909090909d0, 0.0909090909090909d0, 0.7272727272727273d0/)
    gpoint15(9 ,:) = (/0.0909090909090909d0, 0.7272727272727273d0, 0.0909090909090909d0/)
    gpoint15(10,:) = (/0.4334498464263357d0, 0.0665501535736643d0, 0.0665501535736643d0/)
    gpoint15(11,:) = (/0.0665501535736643d0, 0.4334498464263357d0, 0.0665501535736643d0/)
    gpoint15(12,:) = (/0.0665501535736643d0, 0.0665501535736643d0, 0.4334498464263357d0/)
    gpoint15(13,:) = (/0.0665501535736643d0, 0.4334498464263357d0, 0.4334498464263357d0/)
    gpoint15(14,:) = (/0.4334498464263357d0, 0.0665501535736643d0, 0.4334498464263357d0/)
    gpoint15(15,:) = (/0.4334498464263357d0, 0.4334498464263357d0, 0.0665501535736643d0/)

    gweight15(1 ) = 0.1817020685825351d0
    gweight15(2 ) = 0.0361607142857143d0
    gweight15(3 ) = 0.0361607142857143d0
    gweight15(4 ) = 0.0361607142857143d0
    gweight15(5 ) = 0.0361607142857143d0
    gweight15(6 ) = 0.0698714945161738d0
    gweight15(7 ) = 0.0698714945161738d0
    gweight15(8 ) = 0.0698714945161738d0
    gweight15(9 ) = 0.0698714945161738d0
    gweight15(10) = 0.0656948493683187d0
    gweight15(11) = 0.0656948493683187d0
    gweight15(12) = 0.0656948493683187d0
    gweight15(13) = 0.0656948493683187d0
    gweight15(14) = 0.0656948493683187d0
    gweight15(15) = 0.0656948493683187d0

  end subroutine gaussian_integral15

  subroutine gaussian_integral24()

    gpoint24(1 ,:) = (/0.3561913862225449d0, 0.2146028712591517d0, 0.2146028712591517d0/)
    gpoint24(2 ,:) = (/0.2146028712591517d0, 0.2146028712591517d0, 0.2146028712591517d0/)
    gpoint24(3 ,:) = (/0.2146028712591517d0, 0.2146028712591517d0, 0.3561913862225449d0/)
    gpoint24(4 ,:) = (/0.2146028712591517d0, 0.3561913862225449d0, 0.2146028712591517d0/)
    gpoint24(5 ,:) = (/0.8779781243961660d0, 0.0406739585346113d0, 0.0406739585346113d0/)
    gpoint24(6 ,:) = (/0.0406739585346113d0, 0.0406739585346113d0, 0.0406739585346113d0/)
    gpoint24(7 ,:) = (/0.0406739585346113d0, 0.0406739585346113d0, 0.8779781243961660d0/)
    gpoint24(8 ,:) = (/0.0406739585346113d0, 0.8779781243961660d0, 0.0406739585346113d0/)
    gpoint24(9 ,:) = (/0.0329863295731731d0, 0.3223378901422757d0, 0.3223378901422757d0/)
    gpoint24(10,:) = (/0.3223378901422757d0, 0.3223378901422757d0, 0.3223378901422757d0/)
    gpoint24(11,:) = (/0.3223378901422757d0, 0.3223378901422757d0, 0.0329863295731731d0/)
    gpoint24(12,:) = (/0.3223378901422757d0, 0.0329863295731731d0, 0.3223378901422757d0/)
    gpoint24(13,:) = (/0.2696723314583159d0, 0.0636610018750175d0, 0.0636610018750175d0/)
    gpoint24(14,:) = (/0.0636610018750175d0, 0.2696723314583159d0, 0.0636610018750175d0/)
    gpoint24(15,:) = (/0.0636610018750175d0, 0.0636610018750175d0, 0.2696723314583159d0/)
    gpoint24(16,:) = (/0.6030056647916491d0, 0.0636610018750175d0, 0.0636610018750175d0/)
    gpoint24(17,:) = (/0.0636610018750175d0, 0.6030056647916491d0, 0.0636610018750175d0/)
    gpoint24(18,:) = (/0.0636610018750175d0, 0.0636610018750175d0, 0.6030056647916491d0/)
    gpoint24(19,:) = (/0.0636610018750175d0, 0.2696723314583159d0, 0.6030056647916491d0/)
    gpoint24(20,:) = (/0.2696723314583159d0, 0.6030056647916491d0, 0.0636610018750175d0/)
    gpoint24(21,:) = (/0.6030056647916491d0, 0.0636610018750175d0, 0.2696723314583159d0/)
    gpoint24(22,:) = (/0.0636610018750175d0, 0.6030056647916491d0, 0.2696723314583159d0/)
    gpoint24(23,:) = (/0.2696723314583159d0, 0.0636610018750175d0, 0.6030056647916491d0/)
    gpoint24(24,:) = (/0.6030056647916491d0, 0.2696723314583159d0, 0.0636610018750175d0/)

    gweight24(1 ) = 0.0399227502581679d0
    gweight24(2 ) = 0.0399227502581679d0
    gweight24(3 ) = 0.0399227502581679d0
    gweight24(4 ) = 0.0399227502581679d0
    gweight24(5 ) = 0.0100772110553207d0
    gweight24(6 ) = 0.0100772110553207d0
    gweight24(7 ) = 0.0100772110553207d0
    gweight24(8 ) = 0.0100772110553207d0
    gweight24(9 ) = 0.0553571815436544d0
    gweight24(10) = 0.0553571815436544d0
    gweight24(11) = 0.0553571815436544d0
    gweight24(12) = 0.0553571815436544d0
    gweight24(13) = 0.0482142857142857d0
    gweight24(14) = 0.0482142857142857d0
    gweight24(15) = 0.0482142857142857d0
    gweight24(16) = 0.0482142857142857d0
    gweight24(17) = 0.0482142857142857d0
    gweight24(18) = 0.0482142857142857d0
    gweight24(19) = 0.0482142857142857d0
    gweight24(20) = 0.0482142857142857d0
    gweight24(21) = 0.0482142857142857d0
    gweight24(22) = 0.0482142857142857d0
    gweight24(23) = 0.0482142857142857d0
    gweight24(24) = 0.0482142857142857d0

  end subroutine gaussian_integral24


  subroutine gaussian_integral31()

    gpoint31(1 ,:) = (/0.2500000000000000d0, 0.2500000000000000d0, 0.2500000000000000d0/)
    gpoint31(2 ,:) = (/0.7653604230090441d0, 0.0782131923303186d0, 0.0782131923303186d0/)
    gpoint31(3 ,:) = (/0.0782131923303186d0, 0.0782131923303186d0, 0.0782131923303186d0/)
    gpoint31(4 ,:) = (/0.0782131923303186d0, 0.0782131923303186d0, 0.7653604230090441d0/)
    gpoint31(5 ,:) = (/0.0782131923303186d0, 0.7653604230090441d0, 0.0782131923303186d0/)
    gpoint31(6 ,:) = (/0.6344703500082868d0, 0.1218432166639044d0, 0.1218432166639044d0/)
    gpoint31(7 ,:) = (/0.1218432166639044d0, 0.1218432166639044d0, 0.1218432166639044d0/)
    gpoint31(8 ,:) = (/0.1218432166639044d0, 0.1218432166639044d0, 0.6344703500082868d0/)
    gpoint31(9 ,:) = (/0.1218432166639044d0, 0.6344703500082868d0, 0.1218432166639044d0/)
    gpoint31(10,:) = (/0.0023825066607383d0, 0.3325391644464206d0, 0.3325391644464206d0/)
    gpoint31(11,:) = (/0.3325391644464206d0, 0.3325391644464206d0, 0.3325391644464206d0/)
    gpoint31(12,:) = (/0.3325391644464206d0, 0.3325391644464206d0, 0.0023825066607383d0/)
    gpoint31(13,:) = (/0.3325391644464206d0, 0.0023825066607383d0, 0.3325391644464206d0/)
    gpoint31(14,:) = (/0.0000000000000000d0, 0.5000000000000000d0, 0.5000000000000000d0/)
    gpoint31(15,:) = (/0.5000000000000000d0, 0.0000000000000000d0, 0.5000000000000000d0/)
    gpoint31(16,:) = (/0.5000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0/)
    gpoint31(17,:) = (/0.5000000000000000d0, 0.0000000000000000d0, 0.0000000000000000d0/)
    gpoint31(18,:) = (/0.0000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0/)
    gpoint31(19,:) = (/0.0000000000000000d0, 0.0000000000000000d0, 0.5000000000000000d0/)
    gpoint31(20,:) = (/0.2000000000000000d0, 0.1000000000000000d0, 0.1000000000000000d0/)
    gpoint31(21,:) = (/0.1000000000000000d0, 0.2000000000000000d0, 0.1000000000000000d0/)
    gpoint31(22,:) = (/0.1000000000000000d0, 0.1000000000000000d0, 0.2000000000000000d0/)
    gpoint31(23,:) = (/0.6000000000000000d0, 0.1000000000000000d0, 0.1000000000000000d0/)
    gpoint31(24,:) = (/0.1000000000000000d0, 0.6000000000000000d0, 0.1000000000000000d0/)
    gpoint31(25,:) = (/0.1000000000000000d0, 0.1000000000000000d0, 0.6000000000000000d0/)
    gpoint31(26,:) = (/0.1000000000000000d0, 0.2000000000000000d0, 0.6000000000000000d0/)
    gpoint31(27,:) = (/0.2000000000000000d0, 0.6000000000000000d0, 0.1000000000000000d0/)
    gpoint31(28,:) = (/0.6000000000000000d0, 0.1000000000000000d0, 0.2000000000000000d0/)
    gpoint31(29,:) = (/0.1000000000000000d0, 0.6000000000000000d0, 0.2000000000000000d0/)
    gpoint31(30,:) = (/0.2000000000000000d0, 0.1000000000000000d0, 0.6000000000000000d0/)
    gpoint31(31,:) = (/0.6000000000000000d0, 0.2000000000000000d0, 0.1000000000000000d0/)

    gweight31(1 ) =  0.1095853407966528d0
    gweight31(2 ) =  0.0635996491464850d0
    gweight31(3 ) =  0.0635996491464850d0
    gweight31(4 ) =  0.0635996491464850d0
    gweight31(5 ) =  0.0635996491464850d0
    gweight31(6 ) = -0.3751064406859797d0
    gweight31(7 ) = -0.3751064406859797d0
    gweight31(8 ) = -0.3751064406859797d0
    gweight31(9 ) = -0.3751064406859797d0
    gweight31(10) =  0.0293485515784412d0
    gweight31(11) =  0.0293485515784412d0
    gweight31(12) =  0.0293485515784412d0
    gweight31(13) =  0.0293485515784412d0
    gweight31(14) =  0.0058201058201058d0
    gweight31(15) =  0.0058201058201058d0
    gweight31(16) =  0.0058201058201058d0
    gweight31(17) =  0.0058201058201058d0
    gweight31(18) =  0.0058201058201058d0
    gweight31(19) =  0.0058201058201058d0
    gweight31(20) =  0.1653439153439105d0
    gweight31(21) =  0.1653439153439105d0
    gweight31(22) =  0.1653439153439105d0
    gweight31(23) =  0.1653439153439105d0
    gweight31(24) =  0.1653439153439105d0
    gweight31(25) =  0.1653439153439105d0
    gweight31(26) =  0.1653439153439105d0
    gweight31(27) =  0.1653439153439105d0
    gweight31(28) =  0.1653439153439105d0
    gweight31(29) =  0.1653439153439105d0
    gweight31(30) =  0.1653439153439105d0
    gweight31(31) =  0.1653439153439105d0

  end subroutine gaussian_integral31

  subroutine gaussian_integral45()

    gpoint45(1 ,:) = (/0.2500000000000000d0, 0.2500000000000000d0, 0.2500000000000000d0/)
    gpoint45(2 ,:) = (/0.6175871903000830d0, 0.1274709365666390d0, 0.1274709365666390d0/)
    gpoint45(3 ,:) = (/0.1274709365666390d0, 0.1274709365666390d0, 0.1274709365666390d0/)
    gpoint45(4 ,:) = (/0.1274709365666390d0, 0.1274709365666390d0, 0.6175871903000830d0/)
    gpoint45(5 ,:) = (/0.1274709365666390d0, 0.6175871903000830d0, 0.1274709365666390d0/)
    gpoint45(6 ,:) = (/0.9037635088221031d0, 0.0320788303926323d0, 0.0320788303926323d0/)
    gpoint45(7 ,:) = (/0.0320788303926323d0, 0.0320788303926323d0, 0.0320788303926323d0/)
    gpoint45(8 ,:) = (/0.0320788303926323d0, 0.0320788303926323d0, 0.9037635088221031d0/)
    gpoint45(9 ,:) = (/0.0320788303926323d0, 0.9037635088221031d0, 0.0320788303926323d0/)
    gpoint45(10,:) = (/0.4502229043567190d0, 0.0497770956432810d0, 0.0497770956432810d0/)
    gpoint45(11,:) = (/0.0497770956432810d0, 0.4502229043567190d0, 0.0497770956432810d0/)
    gpoint45(12,:) = (/0.0497770956432810d0, 0.0497770956432810d0, 0.4502229043567190d0/)
    gpoint45(13,:) = (/0.0497770956432810d0, 0.4502229043567190d0, 0.4502229043567190d0/)
    gpoint45(14,:) = (/0.4502229043567190d0, 0.0497770956432810d0, 0.4502229043567190d0/)
    gpoint45(15,:) = (/0.4502229043567190d0, 0.4502229043567190d0, 0.0497770956432810d0/)
    gpoint45(16,:) = (/0.3162695526014501d0, 0.1837304473985499d0, 0.1837304473985499d0/)
    gpoint45(17,:) = (/0.1837304473985499d0, 0.3162695526014501d0, 0.1837304473985499d0/)
    gpoint45(18,:) = (/0.1837304473985499d0, 0.1837304473985499d0, 0.3162695526014501d0/)
    gpoint45(19,:) = (/0.1837304473985499d0, 0.3162695526014501d0, 0.3162695526014501d0/)
    gpoint45(20,:) = (/0.3162695526014501d0, 0.1837304473985499d0, 0.3162695526014501d0/)
    gpoint45(21,:) = (/0.3162695526014501d0, 0.3162695526014501d0, 0.1837304473985499d0/)
    gpoint45(22,:) = (/0.0229177878448171d0, 0.2319010893971509d0, 0.2319010893971509d0/)
    gpoint45(23,:) = (/0.2319010893971509d0, 0.0229177878448171d0, 0.2319010893971509d0/)
    gpoint45(24,:) = (/0.2319010893971509d0, 0.2319010893971509d0, 0.0229177878448171d0/)
    gpoint45(25,:) = (/0.5132800333608811d0, 0.2319010893971509d0, 0.2319010893971509d0/)
    gpoint45(26,:) = (/0.2319010893971509d0, 0.5132800333608811d0, 0.2319010893971509d0/)
    gpoint45(27,:) = (/0.2319010893971509d0, 0.2319010893971509d0, 0.5132800333608811d0/)
    gpoint45(28,:) = (/0.2319010893971509d0, 0.0229177878448171d0, 0.5132800333608811d0/)
    gpoint45(29,:) = (/0.0229177878448171d0, 0.5132800333608811d0, 0.2319010893971509d0/)
    gpoint45(30,:) = (/0.5132800333608811d0, 0.2319010893971509d0, 0.0229177878448171d0/)
    gpoint45(31,:) = (/0.2319010893971509d0, 0.5132800333608811d0, 0.0229177878448171d0/)
    gpoint45(32,:) = (/0.0229177878448171d0, 0.2319010893971509d0, 0.5132800333608811d0/)
    gpoint45(33,:) = (/0.5132800333608811d0, 0.0229177878448171d0, 0.2319010893971509d0/)
    gpoint45(34,:) = (/0.7303134278075384d0, 0.0379700484718286d0, 0.0379700484718286d0/)
    gpoint45(35,:) = (/0.0379700484718286d0, 0.7303134278075384d0, 0.0379700484718286d0/)
    gpoint45(36,:) = (/0.0379700484718286d0, 0.0379700484718286d0, 0.7303134278075384d0/)
    gpoint45(37,:) = (/0.1937464752488044d0, 0.0379700484718286d0, 0.0379700484718286d0/)
    gpoint45(38,:) = (/0.0379700484718286d0, 0.1937464752488044d0, 0.0379700484718286d0/)
    gpoint45(39,:) = (/0.0379700484718286d0, 0.0379700484718286d0, 0.1937464752488044d0/)
    gpoint45(40,:) = (/0.0379700484718286d0, 0.7303134278075384d0, 0.1937464752488044d0/)
    gpoint45(41,:) = (/0.7303134278075384d0, 0.1937464752488044d0, 0.0379700484718286d0/)
    gpoint45(42,:) = (/0.1937464752488044d0, 0.0379700484718286d0, 0.7303134278075384d0/)
    gpoint45(43,:) = (/0.0379700484718286d0, 0.1937464752488044d0, 0.7303134278075384d0/)
    gpoint45(44,:) = (/0.7303134278075384d0, 0.0379700484718286d0, 0.1937464752488044d0/)
    gpoint45(45,:) = (/0.1937464752488044d0, 0.7303134278075384d0, 0.0379700484718286d0/)

    gweight45(1 ) = -0.2359620398477557d0
    gweight45(2 ) =  0.0244878963560562d0
    gweight45(3 ) =  0.0244878963560562d0
    gweight45(4 ) =  0.0244878963560562d0
    gweight45(5 ) =  0.0244878963560562d0
    gweight45(6 ) =  0.0039485206398261d0
    gweight45(7 ) =  0.0039485206398261d0
    gweight45(8 ) =  0.0039485206398261d0
    gweight45(9 ) =  0.0039485206398261d0
    gweight45(10) =  0.0263055529507371d0
    gweight45(11) =  0.0263055529507371d0
    gweight45(12) =  0.0263055529507371d0
    gweight45(13) =  0.0263055529507371d0
    gweight45(14) =  0.0263055529507371d0
    gweight45(15) =  0.0263055529507371d0
    gweight45(16) =  0.0829803830550589d0
    gweight45(17) =  0.0829803830550589d0
    gweight45(18) =  0.0829803830550589d0
    gweight45(19) =  0.0829803830550589d0
    gweight45(20) =  0.0829803830550589d0
    gweight45(21) =  0.0829803830550589d0
    gweight45(22) =  0.0254426245481023d0
    gweight45(23) =  0.0254426245481023d0
    gweight45(24) =  0.0254426245481023d0
    gweight45(25) =  0.0254426245481023d0
    gweight45(26) =  0.0254426245481023d0
    gweight45(27) =  0.0254426245481023d0
    gweight45(28) =  0.0254426245481023d0
    gweight45(29) =  0.0254426245481023d0
    gweight45(30) =  0.0254426245481023d0
    gweight45(31) =  0.0254426245481023d0
    gweight45(32) =  0.0254426245481023d0
    gweight45(33) =  0.0254426245481023d0
    gweight45(34) =  0.0134324384376852d0
    gweight45(35) =  0.0134324384376852d0
    gweight45(36) =  0.0134324384376852d0
    gweight45(37) =  0.0134324384376852d0
    gweight45(38) =  0.0134324384376852d0
    gweight45(39) =  0.0134324384376852d0
    gweight45(40) =  0.0134324384376852d0
    gweight45(41) =  0.0134324384376852d0
    gweight45(42) =  0.0134324384376852d0
    gweight45(43) =  0.0134324384376852d0
    gweight45(44) =  0.0134324384376852d0
    gweight45(45) =  0.0134324384376852d0

  end subroutine gaussian_integral45



subroutine local_nodes()

  localpointp2(1, :) = (/0.0d0,0.0d0,0.0d0/)
  localpointp2(2, :) = (/1.0d0,0.0d0,0.0d0/)
  localpointp2(3, :) = (/0.0d0,1.0d0,0.0d0/)
  localpointp2(4, :) = (/0.0d0,0.0d0,1.0d0/)
  localpointp2(5, :) = (/0.5d0,0.0d0,0.0d0/)
  localpointp2(6, :) = (/0.5d0,0.5d0,0.0d0/)
  localpointp2(7, :) = (/0.0d0,0.5d0,0.0d0/)
  localpointp2(8, :) = (/0.0d0,0.0d0,0.5d0/)
  localpointp2(9, :) = (/0.5d0,0.0d0,0.5d0/)
  localpointp2(10,:) = (/0.0d0,0.5d0,0.5d0/)

end subroutine local_nodes

  subroutine basis_p2(list)

    double precision,dimension(:),intent(in) :: list
    ! double precision,dimension(1:10) :: phi_p2
    
    phi_p2(1) =(1.0d0-list(1)-list(2)-list(3))*(1.0d0-2.0d0*list(1)-2.0d0*list(2)-2.0d0*list(3))  !!! 2(1-x-y-z)(1/2-x-y-z)
    phi_p2(2) =2.0d0*list(1)*(list(1)-0.5d0)                                                      !!! 2x(x-1/2)
    phi_p2(3) =2.0d0*list(2)*(list(2)-0.5d0)                                                      !!! 2y(y-1/2)
    phi_p2(4) =2.0d0*list(3)*(list(3)-0.5d0)                                                      !!! 2z(z-1/2)
    phi_p2(5) =4.0d0*list(1)*(1.0d0-list(1)-list(2)-list(3))                                      !!! 4x(2-x-y-z)
    phi_p2(6) =4.0d0*list(1)*list(2)                                                              !!! 4xy
    phi_p2(7) =4.0d0*list(2)*(1.0d0-list(1)-list(2)-list(3))                                      !!! 4y(2-x-y-z)
    phi_p2(8) =4.0d0*list(3)*(1.0d0-list(1)-list(2)-list(3))                                      !!! 4z(2-x-y-z)
    phi_p2(9) =4.0d0*list(1)*list(3)                                                              !!! 4xz
    phi_p2(10)=4.0d0*list(2)*list(3)                                                              !!! 4yz

  end subroutine basis_p2

subroutine basis_p2_del(list)

    double precision,dimension(:),intent(in) :: list
    ! double precision,dimension(1:10,1:3) :: phi_p2_del

    phi_p2_del(1,1:3) =(/-3.0d0+4.0d0*list(1)+4.0d0*list(2)+4.0d0*list(3),-3.0d0+4.0d0*list(1)+4.0d0*list(2)+4.0d0*list(3),&
                         -3.0d0+4.0d0*list(1)+4.0d0*list(2)+4.0d0*list(3)/)
    phi_p2_del(2,1:3) =(/4.0d0*list(1)-1.0d0,0.0d0,0.0d0/)
    phi_p2_del(3,1:3) =(/0.0d0,4.0d0*list(2)-1.0d0,0.0d0/)
    phi_p2_del(4,1:3) =(/0.0d0,0.0d0,4.0d0*list(3)-1.0d0/)
    phi_p2_del(5,1:3) =(/4.0d0*(1.0d0-2.0d0*list(1)-list(2)-list(3)),-4.0d0*list(1),-4.0d0*list(1)/)
    phi_p2_del(6,1:3) =(/4.0d0*list(2),4.0d0*list(1),0.0d0/)
    phi_p2_del(7,1:3) =(/-4.0d0*list(2),4.0d0*(1.0d0-list(1)-2.0d0*list(2)-list(3)),-4.0d0*list(2)/)
    phi_p2_del(8,1:3) =(/-4.0d0*list(3),-4.0d0*list(3),4.0d0*(1.0d0-list(1)-list(2)-2.0d0*list(3))/)
    phi_p2_del(9,1:3) =(/4.0d0*list(3),0.0d0,4.0d0*list(1)/)
    phi_p2_del(10,1:3)=(/0.0d0,4.0d0*list(3),4.0d0*list(2)/)
    
end subroutine basis_p2_del

subroutine basis_p3(list)

  double precision,dimension(:) :: list
  ! double precision,dimension(1:20) :: phi_p3
  
phi_p3(1) =9.0d0/2.0d0*(1.0d0-list(1)-list(2)-list(3))*(2.0d0/3.0d0-list(1)-list(2)-list(3))*(1.0d0/3.0d0-list(1)-list(2)-list(3))   !!! 9/2*(1-x-y-z)(2/3-x-y-z)(1/3-x-y-z)
phi_p3(2) =9.0d0/2.0d0*list(1)*(list(1)-1.0d0/3.0d0)*(list(1)-2.0d0/3.0d0)                                                           !!! 9/2*x(x-1/3)(x-2/3)
phi_p3(3) =9.0d0/2.0d0*list(2)*(list(2)-1.0d0/3.0d0)*(list(2)-2.0d0/3.0d0)                                                           !!! 9/2*y(y-1/3)(y-2/3)
phi_p3(4) =9.0d0/2.0d0*list(3)*(list(3)-1.0d0/3.0d0)*(list(3)-2.0d0/3.0d0)                                                           !!! 9/2*z(z-1/3)(z-2/3)
phi_p3(5) =27.0d0/2.0d0*list(1)*(1.0d0-list(1)-list(2)-list(3))*(2.0d0/3.0d0-list(1)-list(2)-list(3))                                !!! 27/2*x(1-x-y-z)(2/3-x-y-z)
phi_p3(6) =27.0d0/2.0d0*list(1)*(1.0d0-list(1)-list(2)-list(3))*(list(1)-1.0d0/3.0d0)                                                !!! 27/2*x(1-x-y-z)(x-1/3)
phi_p3(7) =27.0d0/2.0d0*list(1)*(list(1)-1.0d0/3.0d0)*list(2)                                                                        !!! 27/2*x(x-1/3)y
phi_p3(8) =27.0d0/2.0d0*list(1)*list(2)*(list(2)-1.0d0/3.0d0)                                                                        !!! 27/2*xy(y-1/3)
phi_p3(9) =27.0d0/2.0d0*list(2)*(1.0d0-list(1)-list(2)-list(3))*(list(2)-1.0d0/3.0d0)                                                !!! 27/2*y(1-x-y-z)(y-1/3)
phi_p3(10)=27.0d0/2.0d0*list(2)*(1.0d0-list(1)-list(2)-list(3))*(2.0d0/3.0d0-list(1)-list(2)-list(3))                                !!! 27/2*y(1-x-y-z)(2/3-x-y-z)
phi_p3(11)=27.0d0/2.0d0*list(3)*(1.0d0-list(1)-list(2)-list(3))*(2.0d0/3.0d0-list(1)-list(2)-list(3))                                !!! 27/2*z(1-x-y-z)(2/3-x-y-z)
phi_p3(12)=27.0d0/2.0d0*list(3)*(1.0d0-list(1)-list(2)-list(3))*(list(3)-1.0d0/3.0d0)                                                !!! 27/2*z(1-x-y-z)(z-1/3)
phi_p3(13)=27.0d0/2.0d0*list(1)*list(3)*(list(1)-1.0d0/3.0d0)                                                                        !!! 27/2*xz(x-1/3)
phi_p3(14)=27.0d0/2.0d0*list(1)*list(3)*(list(3)-1.0d0/3.0d0)                                                                        !!! 27/2*xz(z-1/3)
phi_p3(15)=27.0d0/2.0d0*list(2)*list(3)*(list(2)-1.0d0/3.0d0)                                                                        !!! 27/2*yz(y-1/3)
phi_p3(16)=27.0d0/2.0d0*list(2)*list(3)*(list(3)-1.0d0/3.0d0)                                                                        !!! 27/2*xz(z-1/3)
phi_p3(17)=27.0d0*list(1)*list(2)*list(3)                                                                                            !!! 27*xyz
phi_p3(18)=27.0d0*list(2)*list(3)*(1.0d0-list(1)-list(2)-list(3))                                                                    !!! 27*yz(1-x-y-z)
phi_p3(19)=27.0d0*list(1)*list(3)*(1.0d0-list(1)-list(2)-list(3))                                                                    !!! 27*xz(1-x-y-z)
phi_p3(20)=27.0d0*list(1)*list(2)*(1.0d0-list(1)-list(2)-list(3))                                                                    !!! 27*xy(1-x-y-z)

end subroutine basis_p3



subroutine basis_p3_del(list)

  double precision,dimension(:) :: list
  ! double precision,dimension(1:20,1:3) :: phi_p3_del

phi_p3_del(1,1:3) =(/-(9.0d0*(list(1)+list(2)+list(3))*(3.0d0*(list(1)+list(2)+list(3))-4.0d0)+11.0d0)/2.0d0,&
                -(9.0d0*(list(1)+list(2)+list(3))*(3.0d0*(list(1)+list(2)+list(3))-4.0d0)+11.0d0)/2.0d0,&
                -(9.0d0*(list(1)+list(2)+list(3))*(3.0d0*(list(1)+list(2)+list(3))-4.0d0)+11.0d0)/2.0d0/)

phi_p3_del(2,1:3) =(/(27.0d0*list(1)**2-18.0d0*list(1)+2.0d0)/2.0d0,0.0d0,0.0d0/)
phi_p3_del(3,1:3) =(/0.0d0,(27.0d0*list(2)**2-18.0d0*list(2)+2.0d0)/2.0d0,0.0d0/)
phi_p3_del(4,1:3) =(/0.0d0,0.0d0,(27.0d0*list(3)**2-18.0d0*list(3)+2.0d0)/2.0d0/)

phi_p3_del(5,1:3) =(/(81.0d0*list(1)**2+(108.0d0*list(2)+108.0d0*list(3)-90.0d0)*list(1)+27.0d0*list(2)**2+27.0d0*list(3)**2+&
54.0d0*list(2)*list(3)-45.0d0*list(2)-45.0d0*list(3)+18.0d0)/2.0d0,1.0d0/2.0d0*9.0d0*list(1)*(6.0d0*(list(1)+list(2)+list(3))-&
5.0d0),1.0d0/2.0d0*9.0d0*list(1)*(6.0d0*(list(1)+list(2)+list(3))-5.0d0)/)

phi_p3_del(6,1:3) =(/-9.0d0/2.0d0*(9.0d0*list(1)**2+(6.0d0*list(2)+6.0d0*list(3)-8.0d0)*list(1)-list(2)-list(3)+1.0d0),-27.0d0/&
2.0d0*list(1)*(list(1)-1.0d0/3.0d0),-27.0d0/2.0d0*list(1)*(list(1)-1.0d0/3.0d0)/)

phi_p3_del(7,1:3) =(/(6.0d0*list(1)-1.0d0)*9.0d0*list(2)/2.0d0,27.0d0/2.0d0*list(1)*(list(1)-1.0d0/3.0d0),0.0d0/)
phi_p3_del(8,1:3) =(/27.0d0/2.0d0*list(2)*(list(2)-1.0d0/3.0d0),(6.0d0*list(2)-1.0d0)*9.0d0*list(1)/2.0d0,0.0d0/)

phi_p3_del(9,1:3) =(/-27.0d0/2.0d0*list(2)*(list(2)-1.0d0/3.0d0),-9.0d0/2.0d0*(9.0d0*list(2)**2+(6.0d0*list(1)+6.0d0*list(3)-&
8.0d0)*list(2)-list(1)-list(3)+1.0d0),-27.0d0/2.0d0*list(2)*(list(2)-1.0d0/3.0d0)/)

phi_p3_del(10,1:3)=(/1.0d0/2.0d0*9.0d0*list(2)*(6.0d0*(list(1)+list(2)+list(3))-5.0d0),(81.0d0*list(2)**2+(108.0d0*list(1)+&
108.0d0*list(3)-90.0d0)*list(2)+27.0d0*list(1)**2+27.0d0*list(3)**2+54.0d0*list(1)*list(3)-45.0d0*list(1)-45.0d0*list(3)+18.0d0)&
/2.0d0,1.0d0/2.0d0*9.0d0*list(2)*(6.0d0*(list(1)+list(2)+list(3))-5.0d0)/)

phi_p3_del(11,1:3)=(/1.0d0/2.0d0*9.0d0*list(3)*(6.0d0*(list(1)+list(2)+list(3))-5.0d0),1.0d0/2.0d0*9.0d0*list(3)*(6.0d0*(list(1)&
+list(2)+list(3))-5.0d0),(81.0d0*list(3)**2+(108.0d0*list(2)+108.0d0*list(1)-90.0d0)*list(3)+27.0d0*list(2)**2+27.0d0*list(1)**2+&
54.0d0*list(2)*list(1)-45.0d0*list(2)-45.0d0*list(1)+18.0d0)/2.0d0/)

phi_p3_del(12,1:3)=(/-27.0d0/2.0d0*list(3)*(list(3)-1.0d0/3.0d0),-27.0d0/2.0d0*list(3)*(list(3)-1.0d0/3.0d0),-9.0d0/2.0d0*(9.0d0&
*list(3)**2+(6.0d0*list(2)+6.0d0*list(1)-8.0d0)*list(3)-list(2)-list(1)+1.0d0)/)

phi_p3_del(13,1:3)=(/(6.0d0*list(1)-1.0d0)*9.0d0*list(3)/2.0d0,0.0d0,27.0d0/2.0d0*list(1)*(list(1)-1.0d0/3.0d0)/)
phi_p3_del(14,1:3)=(/27.0d0/2.0d0*list(3)*(list(3)-1.0d0/3.0d0),0.0d0,(6.0d0*list(3)-1.0d0)*9.0d0*list(1)/2.0d0/)
phi_p3_del(15,1:3)=(/0.0d0,(6.0d0*list(2)-1.0d0)*9.0d0*list(3)/2.0d0,27.0d0/2.0d0*list(2)*(list(2)-1.0d0/3.0d0)/)
phi_p3_del(16,1:3)=(/0.0d0,27.0d0/2.0d0*list(3)*(list(3)-1.0d0/3.0d0),(6.0d0*list(3)-1.0d0)*9.0d0*list(2)/2.0d0/)
phi_p3_del(17,1:3)=(/27.0d0*list(2)*list(3),27.0d0*list(1)*list(3),27.0d0*list(1)*list(2)/)

phi_p3_del(18,1:3)=(/-27.0d0*list(2)*list(3),27.0d0*list(3)*((1.0d0-list(1)-list(2)-list(3))-list(2)),&
27.0d0*list(2)*((1.0d0-list(1)-list(2)-list(3))-list(3))/)

phi_p3_del(19,1:3)=(/27.0d0*list(3)*((1.0d0-list(1)-list(2)-list(3))-list(1)),-27.0d0*list(1)*list(3),&
27.0d0*list(1)*((1.0d0-list(1)-list(2)-list(3))-list(3))/)

phi_p3_del(20,1:3)=(/27.0d0*list(2)*((1.0d0-list(1)-list(2)-list(3))-list(1)),27.0d0*list(1)*((1.0d0-&
list(1)-list(2)-list(3))-list(2)),-27.0d0*list(1)*list(2)/)
  
end subroutine basis_p3_del

subroutine mass_mat()

    p1_matrix(1,:) =  (/1.0d0/60.0d0,1.0d0/120.0d0,1.0d0/120.0d0,1.0d0/120.0d0/)
    p1_matrix(2,:) =  (/1.0d0/120.0d0,1.0d0/60.0d0,1.0d0/120.0d0,1.0d0/120.0d0/)
    p1_matrix(3,:) =  (/1.0d0/120.0d0,1.0d0/120,1.0d0/60.0d0,1.0d0/120.0d0/)
    p1_matrix(4,:) =  (/1.0d0/120.0d0,1.0d0/120,1.0d0/120.0d0,1/60.0d0/)


    p2_matrix(1,:)  =  (/1.0d0/420.0d0,1.0d0/2520.0d0,1.0d0/2520.0d0,1.0d0/2520.0d0,-1.0d0/630.0d0,-1.0d0/420.0d0,&
                             -1.0d0/630.0d0,-1.0d0/630.0d0,-1.0d0/420.0d0,-1.0d0/420.0d0/)
    p2_matrix(2,:)  =  (/1.0d0/2520.0d0,1.0d0/420.0d0,1.0d0/2520.0d0,1.0d0/2520.0d0,-1.0d0/630.0d0,-1.0d0/630.0d0,&
                             -1.0d0/420.0d0,-1.0d0/420.0d0,-1.0d0/630.0d0,-1.0d0/420.0d0/)
    p2_matrix(3,:)  =  (/1.0d0/2520.0d0,1.0d0/2520.0d0,1.0d0/420.0d0,1.0d0/2520.0d0,-1.0d0/420.0d0,-1.0d0/630.0d0,&
                             -1.0d0/630.0d0,-1.0d0/420.0d0,-1.0d0/420.0d0,-1.0d0/630.0d0/)
    p2_matrix(4,:)  =  (/1.0d0/2520.0d0,1.0d0/2520.0d0,1.0d0/2520.0d0,1.0d0/420.0d0,-1.0d0/420.0d0,-1.0d0/420.0d0,&
                             -1.0d0/420.0d0,-1.0d0/630.0d0,-1.0d0/630.0d0,-1.0d0/630.0d0/)
    p2_matrix(5,:)  =  (/-1.0d0/630.0d0,-1.0d0/630.0d0,-1.0d0/420.0d0,-1.0d0/420.0d0,4.0d0/315.0d0,2.0d0/315.0d0,&
                               2.0d0/315.0d0,2.0d0/315.0d0,2.0d0/315.0d0,1.0d0/315.0d0/)
    p2_matrix(6,:)  =  (/-1.0d0/420.0d0,-1.0d0/630.0d0,-1.0d0/630.0d0,-1.0d0/420.0d0,2.0d0/315.0d0,4.0d0/315.0d0,&
                               2.0d0/315.0d0,1.0d0/315.0d0,2.0d0/315.0d0,2.0d0/315.0d0/)
    p2_matrix(7,:)  =  (/-1.0d0/630.0d0,-1.0d0/420.0d0,-1.0d0/630.0d0,-1.0d0/420.0d0,2.0d0/315.0d0,2.0d0/315.0d0,&
                               4.0d0/315.0d0,2.0d0/315.0d0,1.0d0/315.0d0,2.0d0/315.0d0/)
    p2_matrix(8,:)  =  (/-1.0d0/630.0d0,-1.0d0/420.0d0,-1.0d0/420.0d0,-1.0d0/630.0d0,2.0d0/315.0d0,1.0d0/315.0d0,&
                               2.0d0/315.0d0,4.0d0/315.0d0,2.0d0/315.0d0,2.0d0/315.0d0/)
    p2_matrix(9,:)  =  (/-1.0d0/420.0d0,-1.0d0/630.0d0,-1.0d0/420.0d0,-1.0d0/630.0d0,2.0d0/315.0d0,2.0d0/315.0d0,&
                               1.0d0/315.0d0,2.0d0/315.0d0,4.0d0/315.0d0,2.0d0/315.0d0/)
    p2_matrix(10,:) =  (/-1.0d0/420.0d0,-1.0d0/420.0d0,-1.0d0/630.0d0,-1.0d0/630.0d0,1.0d0/315.0d0,2.0d0/315.0d0,&
                               2.0d0/315.0d0,2.0d0/315.0d0,2.0d0/315.0d0,4.0d0/315.0d0/)


    p3_matrix(1,:) =(/ 1.0d0/1680.0d0,  1.0d0/13440.0d0, 1.0d0/13440.0d0, 1.0d0/13440.0d0, -1.0d0/2240.0d0,  1.0d0/4480.0d0,&
    1.0d0/8960.0d0,   1.0d0/8960.0d0,  1.0d0/4480.0d0, -1.0d0/2240.0d0, -1.0d0/2240.0d0,  1.0d0/4480.0d0,  1.0d0/8960.0d0, &
      1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  3.0d0/2240.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0 /)
    p3_matrix(2,:) =(/ 1.0d0/13440.0d0, 1.0d0/1680.0d0,  1.0d0/13440.0d0, 1.0d0/13440.0d0,  1.0d0/4480.0d0, -1.0d0/2240.0d0,&
    -1.0d0/2240.0d0,   1.0d0/4480.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0, -1.0d0/2240.0d0, &
      1.0d0/4480.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  3.0d0/4480.0d0,  3.0d0/2240.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0 /)
    p3_matrix(3,:) =(/ 1.0d0/13440.0d0, 1.0d0/13440.0d0, 1.0d0/1680.0d0,  1.0d0/13440.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,&
    1.0d0/4480.0d0,  -1.0d0/2240.0d0, -1.0d0/2240.0d0,  1.0d0/4480.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0, &
      1.0d0/8960.0d0, -1.0d0/2240.0d0,  1.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/2240.0d0,  3.0d0/4480.0d0 /)
    p3_matrix(4,:) =(/ 1.0d0/13440.0d0, 1.0d0/13440.0d0, 1.0d0/13440.0d0, 1.0d0/1680.0d0,   1.0d0/8960.0d0,  1.0d0/8960.0d0,&
    1.0d0/8960.0d0,   1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/4480.0d0, -1.0d0/2240.0d0,  1.0d0/4480.0d0, &
    -1.0d0/2240.0d0,  1.0d0/4480.0d0, -1.0d0/2240.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/2240.0d0 /)
    p3_matrix(5,:) =(/-1.0d0/2240.0d0,  1.0d0/4480.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,   9.0d0/2240.0d0, -9.0d0/4480.0d0,&
    -9.0d0/8960.0d0,   0.0d0,          -9.0d0/8960.0d0,  9.0d0/4480.0d0,  9.0d0/4480.0d0, -9.0d0/8960.0d0, -9.0d0/8960.0d0, &
      0.0d0,           0.0d0,           0.0d0,          -9.0d0/4480.0d0,  0.0d0,           0.0d0,           0.0d0          /)
    p3_matrix(6,:) =(/ 1.0d0/4480.0d0, -1.0d0/2240.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  -9.0d0/4480.0d0,  9.0d0/2240.0d0,&
    9.0d0/4480.0d0,  -9.0d0/8960.0d0,  0.0d0,          -9.0d0/8960.0d0, -9.0d0/8960.0d0,  0.0d0,           9.0d0/4480.0d0, &
    -9.0d0/8960.0d0,  0.0d0,           0.0d0,           0.0d0,          -9.0d0/4480.0d0,  0.0d0,           0.0d0          /)
    p3_matrix(7,:) =(/ 1.0d0/8960.0d0, -1.0d0/2240.0d0,  1.0d0/4480.0d0,  1.0d0/8960.0d0,  -9.0d0/8960.0d0,  9.0d0/4480.0d0,&
    9.0d0/2240.0d0,  -9.0d0/4480.0d0, -9.0d0/8960.0d0,  0.0d0,           0.0d0,           0.0d0,           9.0d0/4480.0d0, &
    -9.0d0/8960.0d0, -9.0d0/8960.0d0,  0.0d0,           0.0d0,          -9.0d0/4480.0d0,  0.0d0,           0.0d0          /)
    p3_matrix(8,:) =(/ 1.0d0/8960.0d0,  1.0d0/4480.0d0, -1.0d0/2240.0d0,  1.0d0/8960.0d0,   0.0d0,          -9.0d0/8960.0d0,&
    -9.0d0/4480.0d0,   9.0d0/2240.0d0,  9.0d0/4480.0d0, -9.0d0/8960.0d0,  0.0d0,           0.0d0,          -9.0d0/8960.0d0, &
      0.0d0,           9.0d0/4480.0d0, -9.0d0/8960.0d0,  0.0d0,           0.0d0,          -9.0d0/4480.0d0,  0.0d0          /)
    p3_matrix(9,:) =(/ 1.0d0/4480.0d0,  1.0d0/8960.0d0, -1.0d0/2240.0d0,  1.0d0/8960.0d0,  -9.0d0/8960.0d0,   0.0d0,         &
    -9.0d0/8960.0d0,   9.0d0/4480.0d0,  9.0d0/2240.0d0, -9.0d0/4480.0d0, -9.0d0/8960.0d0,  0.0d0,           0.0d0,          &
      0.0d0,           9.0d0/4480.0d0, -9.0d0/8960.0d0,  0.0d0,           0.0d0,          -9.0d0/4480.0d0,  0.0d0          /)
    p3_matrix(10,:)=(/-1.0d0/2240.0d0,  1.0d0/8960.0d0,  1.0d0/4480.0d0,  1.0d0/8960.0d0,   9.0d0/4480.0d0, -9.0d0/8960.0d0,&
    0.0d0,           -9.0d0/8960.0d0, -9.0d0/4480.0d0,  9.0d0/2240.0d0,  9.0d0/4480.0d0, -9.0d0/8960.0d0,  0.0d0,          &
      0.0d0,          -9.0d0/8960.0d0,  0.0d0,          -9.0d0/4480.0d0,  0.0d0,           0.0d0,           0.0d0          /)

    p3_matrix(11,:)=(/-1.0d0/2240.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/4480.0d0,   9.0d0/4480.0d0, -9.0d0/8960.0d0,&
    0.0d0,            0.0d0,          -9.0d0/8960.0d0,  9.0d0/4480.0d0,  9.0d0/2240.0d0, -9.0d0/4480.0d0,  0.0d0,          &
    -9.0d0/8960.0d0,  0.0d0,          -9.0d0/8960.0d0, -9.0d0/4480.0d0,  0.0d0,           0.0d0,           0.0d0          /)
    p3_matrix(12,:)=(/ 1.0d0/4480.0d0,  1.0d0/8960.0d0,  1.0d0/8960.0d0, -1.0d0/2240.0d0,  -9.0d0/8960.0d0,  0.0d0,         &
    0.0d0,            0.0d0,           0.0d0,          -9.0d0/8960.0d0, -9.0d0/4480.0d0,  9.0d0/2240.0d0, -9.0d0/8960.0d0, &
      9.0d0/4480.0d0, -9.0d0/8960.0d0,  9.0d0/4480.0d0,  0.0d0,           0.0d0,           0.0d0,          -9.0d0/4480.0d0 /)
    p3_matrix(13,:)=(/ 1.0d0/8960.0d0, -1.0d0/2240.0d0,  1.0d0/8960.0d0,  1.0d0/4480.0d0,  -9.0d0/8960.0d0,  9.0d0/4480.0d0,&
    9.0d0/4480.0d0,  -9.0d0/8960.0d0,  0.0d0,           0.0d0,           0.0d0,          -9.0d0/8960.0d0,  9.0d0/2240.0d0, &
    -9.0d0/4480.0d0,  0.0d0,          -9.0d0/8960.0d0,  0.0d0,          -9.0d0/4480.0d0,  0.0d0,           0.0d0          /)
    p3_matrix(14,:)=(/ 1.0d0/8960.0d0,  1.0d0/4480.0d0,  1.0d0/8960.0d0, -1.0d0/2240.0d0,   0.0d0,          -9.0d0/8960.0d0,&
    -9.0d0/8960.0d0,   0.0d0,           0.0d0,           0.0d0,          -9.0d0/8960.0d0,  9.0d0/4480.0d0, -9.0d0/4480.0d0, &
      9.0d0/2240.0d0, -9.0d0/8960.0d0,  9.0d0/4480.0d0,  0.0d0,           0.0d0,           0.0d0,          -9.0d0/4480.0d0 /)
    p3_matrix(15,:)=(/ 1.0d0/8960.0d0,  1.0d0/8960.0d0, -1.0d0/2240.0d0,  1.0d0/4480.0d0,   0.0d0,           0.0d0,         &
    -9.0d0/8960.0d0,   9.0d0/4480.0d0,  9.0d0/4480.0d0, -9.0d0/8960.0d0,  0.0d0,          -9.0d0/8960.0d0,  0.0d0,          &
    -9.0d0/8960.0d0,  9.0d0/2240.0d0, -9.0d0/4480.0d0,  0.0d0,           0.0d0,          -9.0d0/4480.0d0,  0.0d0          /)
    p3_matrix(16,:)=(/ 1.0d0/8960.0d0,  1.0d0/8960.0d0,  1.0d0/4480.0d0, -1.0d0/2240.0d0,   0.0d0,           0.0d0,         &
    0.0d0,           -9.0d0/8960.0d0, -9.0d0/8960.0d0,  0.0d0,          -9.0d0/8960.0d0,  9.0d0/4480.0d0, -9.0d0/8960.0d0, &
      9.0d0/4480.0d0, -9.0d0/4480.0d0,  9.0d0/2240.0d0,  0.0d0,           0.0d0,           0.0d0,          -9.0d0/4480.0d0 /)
    p3_matrix(17,:)=(/ 3.0d0/2240.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,  -9.0d0/4480.0d0,  0.0d0,         &
    0.0d0,            0.0d0,           0.0d0,          -9.0d0/4480.0d0, -9.0d0/4480.0d0,  0.0d0,           0.0d0,          &
      0.0d0,           0.0d0,           0.0d0,           9.0d0/560.0d0,   9.0d0/1120.0d0,  9.0d0/1120.0d0,  9.0d0/1120.0d0 /)
    p3_matrix(18,:)=(/ 3.0d0/4480.0d0,  3.0d0/2240.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,   0.0d0,          -9.0d0/4480.0d0,&
    -9.0d0/4480.0d0,   0.0d0,           0.0d0,           0.0d0,           0.0d0,           0.0d0,          -9.0d0/4480.0d0, &
      0.0d0,           0.0d0,           0.0d0,           9.0d0/1120.0d0,  9.0d0/560.0d0,   9.0d0/1120.0d0,  9.0d0/1120.0d0 /)
    p3_matrix(19,:)=(/ 3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/2240.0d0,  3.0d0/4480.0d0,   0.0d0,           0.0d0,         &
    0.0d0,           -9.0d0/4480.0d0, -9.0d0/4480.0d0,  0.0d0,           0.0d0,           0.0d0,           0.0d0,          &
      0.0d0,          -9.0d0/4480.0d0,  0.0d0,           9.0d0/1120.0d0,  9.0d0/1120.0d0,  9.0d0/560.0d0,   9.0d0/1120.0d0 /)
    p3_matrix(20,:)=(/ 3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/4480.0d0,  3.0d0/2240.0d0,   0.0d0,           0.0d0,         &
    0.0d0,            0.0d0,           0.0d0,           0.0d0,           0.0d0,          -9.0d0/4480.0d0,  0.0d0,          &
    -9.0d0/4480.0d0,  0.0d0,          -9.0d0/4480.0d0,  9.0d0/1120.0d0,  9.0d0/1120.0d0,  9.0d0/1120.0d0,  9.0d0/560.0d0  /)
 

end subroutine mass_mat




function interpo_gaussian(e1,e2,e3,e4,g,f1,f2,f3,f4) result(f)

    integer,intent(in) :: e1,e2,e3,e4
    double precision,intent(in) :: f1,f2,f3,f4
    double precision,dimension(1:3),intent(in) :: g
    double precision :: V,V1,V2,V3,V4,f
  
  
    call jacobian(e1,e2,e3,e4)
    V=Jdet
    call jacobian1(g(1),g(2),g(3),point(e2,1),point(e2,2),point(e2,3)&
                  ,point(e3,1),point(e3,2),point(e3,3),point(e4,1),point(e4,2),point(e4,3))
    V1=Jdet1
    call jacobian1(point(e1,1),point(e1,2),point(e1,3),g(1),g(2),g(3)&
                  ,point(e3,1),point(e3,2),point(e3,3),point(e4,1),point(e4,2),point(e4,3))
    V2=Jdet1
    call jacobian1(point(e1,1),point(e1,2),point(e1,3),point(e2,1),point(e2,2),point(e2,3)&
                  ,g(1),g(2),g(3),point(e4,1),point(e4,2),point(e4,3))
    V3=Jdet1
    call jacobian1(point(e1,1),point(e1,2),point(e1,3),point(e2,1),point(e2,2),point(e2,3)&
                  ,point(e3,1),point(e3,2),point(e3,3),g(1),g(2),g(3))
    V4=Jdet1
  
  
    f=f1*V1/V+f2*V2/V+f3*V3/V+f4*V4/V
  
  
  end function interpo_gaussian


  subroutine psi_interpolation0(Ne,Nquadrature,states,degree)
    integer,intent(in) :: Ne,Nquadrature,states,degree
    integer :: i,ii,jj,kk
    double precision,dimension(:),allocatable :: location!,basis_function
    double precision,dimension(:,:),allocatable :: f
    double precision :: f0


    if (degree==2) then

    allocate(location(1:3))
    allocate(f(1:10,1:states))
    ! allocate(basis_function(1:10))

    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        psi_g(i,ii)%location=location

        call basis_p2(gpoint(ii,:))

        f(1 ,:)=psi(ele(i,1 ),1:states)
        f(2 ,:)=psi(ele(i,2 ),1:states)
        f(3 ,:)=psi(ele(i,3 ),1:states)
        f(4 ,:)=psi(ele(i,4 ),1:states)
        f(5 ,:)=psi(ele(i,5 ),1:states)
        f(6 ,:)=psi(ele(i,6 ),1:states)
        f(7 ,:)=psi(ele(i,7 ),1:states)
        f(8 ,:)=psi(ele(i,8 ),1:states)
        f(9 ,:)=psi(ele(i,9 ),1:states)
        f(10,:)=psi(ele(i,10),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,10
            f0=f0+f(kk,jj)*phi_p2(kk)
          enddo
          psi_g(i,ii)%state_value(jj)=f0
        enddo
      
      enddo
    enddo

    deallocate(location)
    deallocate(f)
    ! deallocate(basis_function)


    else if (degree==3) then

      allocate(location(1:3))
      allocate(f(1:20,1:states))
      ! allocate(basis_function(1:10))
  
      do i=1,Ne
        call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
        call gaussian_integral
  
        do ii=1,Nquadrature
          location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
  
          psi_g(i,ii)%location=location
  
          call basis_p3(gpoint(ii,:))
  
          f(1 ,:)=psi(ele(i,1 ),1:states)
          f(2 ,:)=psi(ele(i,2 ),1:states)
          f(3 ,:)=psi(ele(i,3 ),1:states)
          f(4 ,:)=psi(ele(i,4 ),1:states)
          f(5 ,:)=psi(ele(i,5 ),1:states)
          f(6 ,:)=psi(ele(i,6 ),1:states)
          f(7 ,:)=psi(ele(i,7 ),1:states)
          f(8 ,:)=psi(ele(i,8 ),1:states)
          f(9 ,:)=psi(ele(i,9 ),1:states)
          f(10,:)=psi(ele(i,10),1:states)
          f(11,:)=psi(ele(i,11),1:states)
          f(12,:)=psi(ele(i,12),1:states)
          f(13,:)=psi(ele(i,13),1:states)
          f(14,:)=psi(ele(i,14),1:states)
          f(15,:)=psi(ele(i,15),1:states)
          f(16,:)=psi(ele(i,16),1:states)
          f(17,:)=psi(ele(i,17),1:states)
          f(18,:)=psi(ele(i,18),1:states)
          f(19,:)=psi(ele(i,19),1:states)
          f(20,:)=psi(ele(i,20),1:states)
  
          
          do jj=1,states
            f0=0.0d0
            do kk=1,20
              f0=f0+f(kk,jj)*phi_p3(kk)
            enddo
            psi_g(i,ii)%state_value(jj)=f0
          enddo
        
        enddo
      enddo
  
      deallocate(location)
      deallocate(f)
      ! deallocate(basis_function)

    end if

  end subroutine psi_interpolation0

  subroutine psi_interpolation(Ne,Nquadrature,states,type,degree)
    integer,intent(in) :: Ne,Nquadrature,states,type,degree
    integer :: i,ii,jj,kk
    double precision,dimension(:),allocatable :: location!,basis_function
    double precision,dimension(:,:),allocatable :: f
    double precision :: f0


    if (degree==2) then

    allocate(location(1:3))
    allocate(f(1:10,1:states))
    ! allocate(basis_function(1:10))


    if (type==0) then

    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        psi_g(i,ii)%location=location

        call basis_p2(gpoint(ii,:))

        f(1 ,:)=psi(ele(i,1 ),1:states)
        f(2 ,:)=psi(ele(i,2 ),1:states)
        f(3 ,:)=psi(ele(i,3 ),1:states)
        f(4 ,:)=psi(ele(i,4 ),1:states)
        f(5 ,:)=psi(ele(i,5 ),1:states)
        f(6 ,:)=psi(ele(i,6 ),1:states)
        f(7 ,:)=psi(ele(i,7 ),1:states)
        f(8 ,:)=psi(ele(i,8 ),1:states)
        f(9 ,:)=psi(ele(i,9 ),1:states)
        f(10,:)=psi(ele(i,10),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,10
            f0=f0+f(kk,jj)*phi_p2(kk)
          enddo
          psi_g(i,ii)%state_value(jj)=f0
        enddo
      
      enddo
    enddo
  
  else if (type==1) then
    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        point_g((i-1)*Nquadrature+ii,:)=location

        call basis_p2(gpoint(ii,:))

        f(1 ,:)=psi(ele(i,1 ),1:states)
        f(2 ,:)=psi(ele(i,2 ),1:states)
        f(3 ,:)=psi(ele(i,3 ),1:states)
        f(4 ,:)=psi(ele(i,4 ),1:states)
        f(5 ,:)=psi(ele(i,5 ),1:states)
        f(6 ,:)=psi(ele(i,6 ),1:states)
        f(7 ,:)=psi(ele(i,7 ),1:states)
        f(8 ,:)=psi(ele(i,8 ),1:states)
        f(9 ,:)=psi(ele(i,9 ),1:states)
        f(10,:)=psi(ele(i,10),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,10
            f0=f0+f(kk,jj)*phi_p2(kk)
          enddo
          psi_point_g((i-1)*Nquadrature+ii,jj)=f0
        enddo
      
      enddo
    enddo

  end if

  deallocate(location)
  deallocate(f)
  ! deallocate(basis_function)

else if (degree==3) then


  allocate(location(1:3))
    allocate(f(1:20,1:states))
    ! allocate(basis_function(1:10))


    if (type==0) then

    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        psi_g(i,ii)%location=location

        call basis_p3(gpoint(ii,:))

        f(1 ,:)=psi(ele(i,1 ),1:states)
        f(2 ,:)=psi(ele(i,2 ),1:states)
        f(3 ,:)=psi(ele(i,3 ),1:states)
        f(4 ,:)=psi(ele(i,4 ),1:states)
        f(5 ,:)=psi(ele(i,5 ),1:states)
        f(6 ,:)=psi(ele(i,6 ),1:states)
        f(7 ,:)=psi(ele(i,7 ),1:states)
        f(8 ,:)=psi(ele(i,8 ),1:states)
        f(9 ,:)=psi(ele(i,9 ),1:states)
        f(10,:)=psi(ele(i,10),1:states)
        f(11,:)=psi(ele(i,11),1:states)
        f(12,:)=psi(ele(i,12),1:states)
        f(13,:)=psi(ele(i,13),1:states)
        f(14,:)=psi(ele(i,14),1:states)
        f(15,:)=psi(ele(i,15),1:states)
        f(16,:)=psi(ele(i,16),1:states)
        f(17,:)=psi(ele(i,17),1:states)
        f(18,:)=psi(ele(i,18),1:states)
        f(19,:)=psi(ele(i,19),1:states)
        f(20,:)=psi(ele(i,20),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,20
            f0=f0+f(kk,jj)*phi_p3(kk)
          enddo
          psi_g(i,ii)%state_value(jj)=f0
        enddo
      
      enddo
    enddo
  
  else if (type==1) then
    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        point_g((i-1)*Nquadrature+ii,:)=location

        call basis_p3(gpoint(ii,:))

        f(1 ,:)=psi(ele(i,1 ),1:states)
        f(2 ,:)=psi(ele(i,2 ),1:states)
        f(3 ,:)=psi(ele(i,3 ),1:states)
        f(4 ,:)=psi(ele(i,4 ),1:states)
        f(5 ,:)=psi(ele(i,5 ),1:states)
        f(6 ,:)=psi(ele(i,6 ),1:states)
        f(7 ,:)=psi(ele(i,7 ),1:states)
        f(8 ,:)=psi(ele(i,8 ),1:states)
        f(9 ,:)=psi(ele(i,9 ),1:states)
        f(10,:)=psi(ele(i,10),1:states)
        f(11,:)=psi(ele(i,11),1:states)
        f(12,:)=psi(ele(i,12),1:states)
        f(13,:)=psi(ele(i,13),1:states)
        f(14,:)=psi(ele(i,14),1:states)
        f(15,:)=psi(ele(i,15),1:states)
        f(16,:)=psi(ele(i,16),1:states)
        f(17,:)=psi(ele(i,17),1:states)
        f(18,:)=psi(ele(i,18),1:states)
        f(19,:)=psi(ele(i,19),1:states)
        f(20,:)=psi(ele(i,20),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,20
            f0=f0+f(kk,jj)*phi_p3(kk)
          enddo
          psi_point_g((i-1)*Nquadrature+ii,jj)=f0
        enddo
      
      enddo
    enddo

  end if

  deallocate(location)
  deallocate(f)
  ! deallocate(basis_function)

end if


  end subroutine psi_interpolation

  subroutine psi_interpolation_HF(Ne,Nquadrature,states,type)
    integer,intent(in) :: Ne,Nquadrature,states,type
    integer :: i,ii,jj,kk
    double precision,dimension(:),allocatable :: location,basis_function
    double precision,dimension(:,:),allocatable :: f
    double precision :: f0


    allocate(location(1:3))
    allocate(f(1:10,1:states))
    allocate(basis_function(1:10))


    if (type==0) then

    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        psi_g(i,ii)%location=location

        call basis_p2(gpoint(ii,:))

        f(1 ,:)=psi_i(ele(i,1 ),1:states)
        f(2 ,:)=psi_i(ele(i,2 ),1:states)
        f(3 ,:)=psi_i(ele(i,3 ),1:states)
        f(4 ,:)=psi_i(ele(i,4 ),1:states)
        f(5 ,:)=psi_i(ele(i,5 ),1:states)
        f(6 ,:)=psi_i(ele(i,6 ),1:states)
        f(7 ,:)=psi_i(ele(i,7 ),1:states)
        f(8 ,:)=psi_i(ele(i,8 ),1:states)
        f(9 ,:)=psi_i(ele(i,9 ),1:states)
        f(10,:)=psi_i(ele(i,10),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,10
            f0=f0+f(kk,jj)*phi_p2(kk)
          enddo
          psi_g(i,ii)%state_value(jj)=f0
        enddo
      
      enddo
    enddo
  
  else if (type==1) then
    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        point_g((i-1)*Nquadrature+ii,:)=location

        call basis_p2(gpoint(ii,:))

        f(1 ,:)=psii(ele(i,1 ),1:states)
        f(2 ,:)=psii(ele(i,2 ),1:states)
        f(3 ,:)=psii(ele(i,3 ),1:states)
        f(4 ,:)=psii(ele(i,4 ),1:states)
        f(5 ,:)=psii(ele(i,5 ),1:states)
        f(6 ,:)=psii(ele(i,6 ),1:states)
        f(7 ,:)=psii(ele(i,7 ),1:states)
        f(8 ,:)=psii(ele(i,8 ),1:states)
        f(9 ,:)=psii(ele(i,9 ),1:states)
        f(10,:)=psii(ele(i,10),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,10
            f0=f0+f(kk,jj)*phi_p2(kk)
          enddo
          psi_point_g((i-1)*Nquadrature+ii,jj)=f0
        enddo
      
      enddo
    enddo

  end if



  end subroutine psi_interpolation_HF

  ! subroutine psi_gradient_interpolation(Ne,Nquadrature,states,type)
  !   integer,intent(in) :: Ne,Nquadrature,states,type
  !   integer :: i,ii,jj,kk
  !   double precision,dimension(:),allocatable :: location,basis_function
  !   double precision,dimension(:,:),allocatable :: f
  !   double precision :: f0


  !   allocate(location(1:3))
  !   allocate(f(1:10,1:states))
  !   allocate(basis_function(1:10))
  

  !   do i=1,Ne
  !     call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
  !     call gaussian_integral

  !     do ii=1,Nquadrature
  !       location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

  !       point_g((i-1)*Nquadrature+ii,:)=location

  !       basis_function = basis_p2(gpoint(ii,:))

  !       f(1 ,:)=psi(ele(i,1 ),1:states)
  !       f(2 ,:)=psi(ele(i,2 ),1:states)
  !       f(3 ,:)=psi(ele(i,3 ),1:states)
  !       f(4 ,:)=psi(ele(i,4 ),1:states)
  !       f(5 ,:)=psi(ele(i,5 ),1:states)
  !       f(6 ,:)=psi(ele(i,6 ),1:states)
  !       f(7 ,:)=psi(ele(i,7 ),1:states)
  !       f(8 ,:)=psi(ele(i,8 ),1:states)
  !       f(9 ,:)=psi(ele(i,9 ),1:states)
  !       f(10,:)=psi(ele(i,10),1:states)

        
  !       do jj=1,states
  !         f0=0.0d0
  !         do kk=1,10
  !           f0=f0+f(kk,jj)*basis_function(kk)
  !         enddo
  !         psi_point_g((i-1)*Nquadrature+ii,jj)=f0
  !       enddo
      
  !     enddo
  !   enddo





  ! end subroutine psi_gradient_interpolation


  subroutine psi_interpolation_GW(Ne,Nquadrature,states)
    integer,intent(in) :: Ne,Nquadrature,states
    integer :: i,ii,jj,kk
    double precision,dimension(:),allocatable :: location!,basis_function
    double precision,dimension(:,:),allocatable :: f
    double precision :: f0


    allocate(location(1:3))
    allocate(f(1:10,1:states))
    ! allocate(basis_function(1:10))

    do i=1,Ne
      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call gaussian_integral

      do ii=1,Nquadrature
        ! location = matmul(J0,gpoint(ii,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)

        ! psi_g(i,ii)%location=location

        call basis_p2(gpoint(ii,:))

        f(1 ,:)=psi_GW(ele(i,1 ),1:states)
        f(2 ,:)=psi_GW(ele(i,2 ),1:states)
        f(3 ,:)=psi_GW(ele(i,3 ),1:states)
        f(4 ,:)=psi_GW(ele(i,4 ),1:states)
        f(5 ,:)=psi_GW(ele(i,5 ),1:states)
        f(6 ,:)=psi_GW(ele(i,6 ),1:states)
        f(7 ,:)=psi_GW(ele(i,7 ),1:states)
        f(8 ,:)=psi_GW(ele(i,8 ),1:states)
        f(9 ,:)=psi_GW(ele(i,9 ),1:states)
        f(10,:)=psi_GW(ele(i,10),1:states)

        
        do jj=1,states
          f0=0.0d0
          do kk=1,10
            f0=f0+f(kk,jj)*phi_p2(kk)
          enddo
          ! psi_g(i,ii)%state_value(jj)=f0
          psi_GW_g((i-1)*Nquadrature+ii,jj)=f0
        enddo
      
      enddo
    enddo

    deallocate(location)
    deallocate(f)
    ! deallocate(basis_function)

  end subroutine psi_interpolation_GW




  subroutine construct_volumegweight(Ne,Nquadrature)

    integer,intent(in) :: Ne,Nquadrature
    integer :: i,j

    volumegweight=0.0d0

    call gaussian_integral

    do i=1,Ne
        do j=1,Nquadrature
            volumegweight((i-1)*Nquadrature+j)=volume(i)*gweight(j)
        enddo
    enddo

end subroutine construct_volumegweight



subroutine construct_NnToNg(Ne,Nquadrature,Nlocal)

    integer,intent(in) :: Ne,Nquadrature,Nlocal
    integer :: i,j,k,l,ii,neigh_size
    complex(kind=(kind(1.0d0))) :: x,y
    ! double precision, dimension(:), allocatable :: basis_function
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: xy_AB

    call gaussian_integral

    ! allocate(basis_function(10))

    y=cmplx(0.0d0,0.0d0,8)

    if (Nlocal==10) then

    do i=1,Ne
        do j=1,Nquadrature
            call basis_p2(gpoint(j,:))
            do k=1,Nlocal
                x=cmplx(phi_p2(k),0.0d0,8)
                call neigh_NnToNg((i-1)*Nquadrature+j)%insert(ele(i,k))
                call neigh_NnToNg((i-1)*Nquadrature+j)%insertAB(ele(i,k),x,y)
            enddo
        enddo
    enddo


    neigh_size = 0
    do i=1,Ne*Nquadrature
        neigh_size=neigh_size+neigh_NnToNg(i)%size
    enddo


    allocate(xy_AB(1:2))

    allocate(NnToNg(neigh_size))
    allocate(NnToNg_IA(Ne*Nquadrature+1))
    allocate(NnToNg_JA(neigh_size))
    allocate(NnToNg_IA_pntrb(Ne*Nquadrature))
    allocate(NnToNg_IA_pntre(Ne*Nquadrature))

    ii=0
    do i=1,Ne*Nquadrature
        NnToNg_IA(i)=ii+1
        do l=1,neigh_NnToNg(i)%size
            ii=ii+1
            NnToNg_JA(ii)=neigh_NnToNg(i)%get(l)
            xy_AB=neigh_NnToNg(i)%getab(l)
            NnToNg(ii)=real(xy_AB(1))
        enddo
    enddo

    NnToNg_IA(Ne*Nquadrature+1)=neigh_size+1

    NnToNg_IA_pntrb=NnToNg_IA(1:Ne*Nquadrature)
    NnToNg_IA_pntre=NnToNg_IA(2:Ne*Nquadrature+1)


  else if (Nlocal==20) then


    do i=1,Ne
      do j=1,Nquadrature
          call basis_p3(gpoint(j,:))
          do k=1,Nlocal
              x=cmplx(phi_p3(k),0.0d0,8)
              call neigh_NnToNg((i-1)*Nquadrature+j)%insert(ele(i,k))
              call neigh_NnToNg((i-1)*Nquadrature+j)%insertAB(ele(i,k),x,y)
          enddo
      enddo
  enddo


  neigh_size = 0
  do i=1,Ne*Nquadrature
      neigh_size=neigh_size+neigh_NnToNg(i)%size
  enddo


  allocate(xy_AB(1:2))

  allocate(NnToNg(neigh_size))
  allocate(NnToNg_IA(Ne*Nquadrature+1))
  allocate(NnToNg_JA(neigh_size))
  allocate(NnToNg_IA_pntrb(Ne*Nquadrature))
  allocate(NnToNg_IA_pntre(Ne*Nquadrature))

  ii=0
  do i=1,Ne*Nquadrature
      NnToNg_IA(i)=ii+1
      do l=1,neigh_NnToNg(i)%size
          ii=ii+1
          NnToNg_JA(ii)=neigh_NnToNg(i)%get(l)
          xy_AB=neigh_NnToNg(i)%getab(l)
          NnToNg(ii)=real(xy_AB(1))
      enddo
  enddo

  NnToNg_IA(Ne*Nquadrature+1)=neigh_size+1

  NnToNg_IA_pntrb=NnToNg_IA(1:Ne*Nquadrature)
  NnToNg_IA_pntre=NnToNg_IA(2:Ne*Nquadrature+1)

end if !!!! end p2/p3 selection


    deallocate(xy_AB)
    ! deallocate(basis_function)



end subroutine construct_NnToNg




  subroutine set_linkedlistAB(Ne,Nlocal,Nquadrature)

    integer,intent(in) :: Ne,Nlocal,Nquadrature
    integer :: i,m,n,k
    complex(kind=(kind(1.0d0))) :: x,y
    logical :: mtest
    ! double precision, dimension(:,:),allocatable :: del_m,del_n
    ! double precision, dimension(:), allocatable :: basis_m,basis_n
    double precision :: A_temp,V_temp,ik,k2, x_real,y_real,x_imag,y_imag



    ! allocate(del_m(1:Nlocal,1:3))
    ! allocate(del_n(1:Nlocal,1:3))
    ! allocate(basis_m(1:Nlocal))
    ! allocate(basis_n(1:Nlocal))


  allocate(volume(1:Ne))

    if (Nlocal==10) then
    
  do i=1,Ne

    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    call mass_mat
    call gaussian_integral
    
    volume(i)=(Jdet)/6.0d0

    do m=1,Nlocal
        do n=1,Nlocal
            A_temp=0
            V_temp=0
            ik=0.0d0
            k2=0.0d0
            do k=1,Nquadrature
                call basis_p2_del(gpoint(k,:))
                ! call basis_p2_del(gpoint(k,:))
                ! call basis_p2(gpoint(k,:))
                ! call basis_p2(gpoint(k,:))
                A_temp = A_temp + dot_product(matmul(Jit,phi_p2_del(m,:)),matmul(Jit,phi_p2_del(n,:)))*gweight(k)
                ! V_temp0 = V_temp0 + potential9990(i0,k0)*basis_m0(m0)*basis_n0(n0)*gweight(k0)
            enddo

            x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
            y_real=Jdet*p2_matrix(m,n)

            x_imag=0.0d0
            y_imag=0.0d0
            
            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)
            

            
            mtest=.true.
            if (neigh_AB(ele(i,m))%find(ele(i,n)))  mtest=.false.
            if (mtest) call neigh_AB(ele(i,m))%insert(ele(i,n))
            call neigh_AB(ele(i,m))%insertAB(ele(i,n),x,y)


        enddo
    enddo

    

  enddo

else if (Nlocal==20) then

  do i=1,Ne

    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    call mass_mat
    call gaussian_integral
    
    volume(i)=(Jdet)/6.0d0

    do m=1,Nlocal
        do n=1,Nlocal
            A_temp=0
            V_temp=0
            ik=0.0d0
            k2=0.0d0
            do k=1,Nquadrature
                call basis_p3_del(gpoint(k,:))
                ! call basis_p2_del(gpoint(k,:))
                ! call basis_p2(gpoint(k,:))
                ! call basis_p2(gpoint(k,:))
                A_temp = A_temp + dot_product(matmul(Jit,phi_p3_del(m,:)),matmul(Jit,phi_p3_del(n,:)))*gweight(k)
                ! V_temp0 = V_temp0 + potential9990(i0,k0)*basis_m0(m0)*basis_n0(n0)*gweight(k0)
            enddo

            x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
            y_real=Jdet*p3_matrix(m,n)

            x_imag=0.0d0
            y_imag=0.0d0
            
            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)
            
            mtest=.true.
            if (neigh_AB(ele(i,m))%find(ele(i,n)))  mtest=.false.
            if (mtest) call neigh_AB(ele(i,m))%insert(ele(i,n))
            call neigh_AB(ele(i,m))%insertAB(ele(i,n),x,y)

            ! print *,i,m,n

        enddo
    enddo

  enddo


end if !!! p2/p3 selection

! print *,'toto'

  ! deallocate(del_m)
  ! deallocate(del_n)
  ! deallocate(basis_m)
  ! deallocate(basis_n)


end subroutine set_linkedlistAB



subroutine set_linkedlistC(Ne,Nlocal,Nquadrature,orbital)

  integer,intent(in) :: Ne,Nlocal,Nquadrature,orbital
  integer :: i,m,n,k
  complex(kind=(kind(1.0d0))) :: x,y
  logical :: mtest
  ! double precision, dimension(:,:),allocatable :: del_m,del_n
  ! double precision, dimension(:), allocatable :: basis_m,basis_n
  double precision :: A_temp,C_temp,ik,k2, x_real,y_real,x_imag,y_imag



  ! allocate(del_m(1:Nlocal,1:3))
  ! allocate(del_n(1:Nlocal,1:3))
  ! allocate(basis_m(1:Nlocal))
  ! allocate(basis_n(1:Nlocal))


! allocate(volume(1:Ne))


do i=1,Ne

  call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
  call mass_mat
  call gaussian_integral
  
  ! volume(i)=(Jdet)/6.0d0

  do m=1,Nlocal
      do n=1,Nlocal
          A_temp=0.0d0
          C_temp=0.0d0
          ik=0.0d0
          k2=0.0d0
          do k=1,Nquadrature
              ! call basis_p2_del(gpoint(k,:))
              ! call basis_p2_del(gpoint(k,:))
              ! call basis_p2(gpoint(k,:))
              call basis_p2(gpoint(k,:))
              ! A_temp = A_temp + dot_product(matmul(Jit,del_m(m,:)),matmul(Jit,del_n(n,:)))*gweight(k)
              ! V_temp0 = V_temp0 + potential9990(i0,k0)*basis_m0(m0)*basis_n0(n0)*gweight(k0)
    
              C_temp = C_temp + psi_GW_g((i-1)*Nquadrature+k,orbital)*phi_p2(m)*phi_p2(n)*gweight(k)
          enddo

          x_real=C_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
          y_real=0.0d0!Jdet*p2_matrix(m,n)

          x_imag=0.0d0
          y_imag=0.0d0
          
          x=cmplx(x_real,x_imag,8)
          y=cmplx(y_real,y_imag,8)
          

          
          mtest=.true.
          if (neigh_C(ele(i,m))%find(ele(i,n)))  mtest=.false.
          if (mtest) call neigh_C(ele(i,m))%insert(ele(i,n))
          call neigh_C(ele(i,m))%insertAB(ele(i,n),x,y)


      enddo
  enddo

enddo

! deallocate(del_m)
! deallocate(del_n)
! deallocate(basis_m)
! deallocate(basis_n)



end subroutine set_linkedlistC


subroutine set_linkedlistAB_efem(Nn,Ne,Nlocal,Nquadrature,Nquadrature_efem,at_c,r0,t,Z)

  integer,intent(in) :: Nn,Ne,Nlocal,Nquadrature,Nquadrature_efem
  double precision, intent(in) :: r0,t,Z
  double precision,dimension(1:3),intent(in) :: at_c
  integer :: i,m,n,k,ii,kk
  complex(kind=(kind(1.0d0))) :: x,y
  logical :: mtest
  ! double precision, dimension(:,:),allocatable :: del_m,del_n
  ! double precision, dimension(:), allocatable :: basis_m,basis_n
  double precision :: A_temp,V_temp,ik,k2, x_real,y_real,x_imag,y_imag
  double precision,dimension(:),allocatable :: location,gweight_temp
  double precision,dimension(:,:),allocatable :: gpoint_temp
  
  
  allocate(location(3))

  if (Nquadrature_efem==45) then
    call gaussian_integral45
    allocate(gweight_temp(45))
    allocate(gpoint_temp(45,3))
    gweight_temp(:)=gweight45(:)
    gpoint_temp(:,:)=gpoint45(:,:)
  else if (Nquadrature_efem==31) then
    call gaussian_integral31
    allocate(gweight_temp(31))
    allocate(gpoint_temp(31,3))
    gweight_temp(:)=gweight31(:)
    gpoint_temp(:,:)=gpoint31(:,:)
  else if (Nquadrature_efem==24) then
    call gaussian_integral24
    allocate(gweight_temp(24))
    allocate(gpoint_temp(24,3))
    gweight_temp(:)=gweight24(:)
    gpoint_temp(:,:)=gpoint24(:,:)
  else if (Nquadrature_efem==15) then
    call gaussian_integral15
    allocate(gweight_temp(15))
    allocate(gpoint_temp(15,3))
    gweight_temp(:)=gweight15(:)
    gpoint_temp(:,:)=gpoint15(:,:)
  else if (Nquadrature_efem==14) then
    call gaussian_integral14
    allocate(gweight_temp(14))
    allocate(gpoint_temp(14,3))
    gweight_temp(:)=gweight14(:)
    gpoint_temp(:,:)=gpoint14(:,:)
  else if (Nquadrature_efem==11) then
    call gaussian_integral
    allocate(gweight_temp(11))
    allocate(gpoint_temp(11,3))
    gweight_temp(:)=gweight(:)
    gpoint_temp(:,:)=gpoint(:,:)
  end if


  ! allocate(del_m(1:Nlocal,1:3))
  ! allocate(del_n(1:Nlocal,1:3))
  ! allocate(basis_m(1:Nlocal))
  ! allocate(basis_n(1:Nlocal))


allocate(volume(1:Ne))

  if (Nlocal==10) then
  
do i=1,Ne

  kk=0
  do ii=1,Nlocal
    if (ele(i,ii)==9) kk=kk+1
  enddo

  if (kk==0) then

  call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
  call mass_mat
  call gaussian_integral
  
  volume(i)=(Jdet)/6.0d0

  do m=1,Nlocal
      do n=1,Nlocal
          A_temp=0
          V_temp=0
          ik=0.0d0
          k2=0.0d0
          do k=1,Nquadrature
              call basis_p2_del(gpoint(k,:))
              ! call basis_p2_del(gpoint(k,:))
              ! call basis_p2(gpoint(k,:))
              ! call basis_p2(gpoint(k,:))
              A_temp = A_temp + dot_product(matmul(Jit,phi_p2_del(m,:)),matmul(Jit,phi_p2_del(n,:)))*gweight(k)
              ! V_temp0 = V_temp0 + potential9990(i0,k0)*basis_m0(m0)*basis_n0(n0)*gweight(k0)
          enddo

          x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
          y_real=Jdet*p2_matrix(m,n)

          x_imag=0.0d0
          y_imag=0.0d0
          
          x=cmplx(x_real,x_imag,8)
          y=cmplx(y_real,y_imag,8)
          

          
          mtest=.true.
          if (neigh_AB(ele(i,m))%find(ele(i,n)))  mtest=.false.
          if (mtest) call neigh_AB(ele(i,m))%insert(ele(i,n))
          call neigh_AB(ele(i,m))%insertAB(ele(i,n),x,y)


      enddo
  enddo

  else 

    ! print *,i

      call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
      call mass_mat
      call gaussian_integral
      ! call gaussian_integral45
      
      volume(i)=(Jdet)/6.0d0

      do m=1,Nlocal
          do n=1,Nlocal
              A_temp=0.0d0
              do k=1,Nquadrature_efem
                  call basis_p2_del(gpoint_temp(k,:))
                  A_temp = A_temp + dot_product(matmul(Jit,phi_p2_del(m,:)),matmul(Jit,phi_p2_del(n,:)))*gweight_temp(k)
              enddo

              x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
              y_real=Jdet*p2_matrix(m,n)
              x_imag=0.0d0
              y_imag=0.0d0
              x=cmplx(x_real,x_imag,8)
              y=cmplx(y_real,y_imag,8)
              mtest=.true.
              if (neigh_AB(ele(i,m))%find(ele(i,n)))  mtest=.false.
              if (mtest) call neigh_AB(ele(i,m))%insert(ele(i,n))
              call neigh_AB(ele(i,m))%insertAB(ele(i,n),x,y)
          enddo
      enddo

      do m=1,Nlocal
          do n=1,Nlocal

              A_temp=0.0d0
              do k=1,Nquadrature_efem
                  call basis_p2_del(gpoint_temp(k,:)) 
                  location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                  ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                  ! call efem_basis_grad(location,at_c,Z,r0) 
                  call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                  call efem_basis_grad(location,at_c,Z,r0,t) 
                  A_temp = A_temp + dot_product(matmul(Jit,phi_p2_del(m,:)*efem_phi(1)+phi_p2(m)*efem_phi_grad(1,:))&
                                              ,matmul(Jit,phi_p2_del(n,:)*efem_phi(1)+phi_p2(n)*efem_phi_grad(1,:)))*gweight_temp(k)
              enddo

              x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0

              
              V_temp=0.d00
              do k=1,Nquadrature_efem
                call basis_p2(gpoint_temp(k,:))
                location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                V_temp=V_temp+efem_phi(1)*phi_p2(m)*efem_phi(1)*phi_p2(n)*gweight_temp(k)
              enddo

              y_real=Jdet*V_temp/6.0d0 !!!!!!! new basis 


              x_imag=0.0d0
              y_imag=0.0d0
              x=cmplx(x_real,x_imag,8)
              y=cmplx(y_real,y_imag,8)
              mtest=.true.
              if (neigh_AB(enriched_ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
              if (mtest) call neigh_AB(enriched_ele(i,m))%insert(enriched_ele(i,n))
              call neigh_AB(enriched_ele(i,m))%insertAB(enriched_ele(i,n),x,y)
          enddo
      enddo



      do m=1,Nlocal
        do n=1,Nlocal

            A_temp=0.0d0
            do k=1,Nquadrature_efem
                call basis_p2_del(gpoint_temp(k,:)) 
                location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                ! call efem_basis_grad(location,at_c,Z,r0) 
                call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                call efem_basis_grad(location,at_c,Z,r0,t) 
                A_temp = A_temp + dot_product(matmul(Jit,phi_p2_del(m,:)*efem_phi(1)+phi_p2(m)*efem_phi_grad(1,:))&
                                             ,matmul(Jit,phi_p2_del(n,:)))*gweight_temp(k)
            enddo

            x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0

            
            V_temp=0.d00
            do k=1,Nquadrature_efem
              call basis_p2(gpoint_temp(k,:))
              location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
              ! call efem_psi_1s_truncated(location,at_c,Z,r0)
              call efem_psi_1s_truncated(location,at_c,Z,r0,t)
              V_temp=V_temp+efem_phi(1)*phi_p2(m)*phi_p2(n)*gweight_temp(k)
            enddo

            y_real=Jdet*V_temp/6.0d0 !!!!!!! new basis 


            x_imag=0.0d0
            y_imag=0.0d0
            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)
            mtest=.true.
            if (neigh_AB(enriched_ele(i,m))%find(ele(i,n)))  mtest=.false.
            if (mtest) call neigh_AB(enriched_ele(i,m))%insert(ele(i,n))
            call neigh_AB(enriched_ele(i,m))%insertAB(ele(i,n),x,y)
        enddo
    enddo


    do m=1,Nlocal
      do n=1,Nlocal

          A_temp=0.0d0
          do k=1,Nquadrature_efem
              call basis_p2_del(gpoint_temp(k,:)) 
              location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
              ! call efem_psi_1s_truncated(location,at_c,Z,r0)
              ! call efem_basis_grad(location,at_c,Z,r0) 
              call efem_psi_1s_truncated(location,at_c,Z,r0,t)
              call efem_basis_grad(location,at_c,Z,r0,t) 
              A_temp = A_temp + dot_product(matmul(Jit,phi_p2_del(m,:))&
                                           ,matmul(Jit,phi_p2_del(n,:)*efem_phi(1)+phi_p2(n)*efem_phi_grad(1,:)))*gweight_temp(k)
          enddo

          x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0

          
          V_temp=0.d00
          do k=1,Nquadrature_efem
            call basis_p2(gpoint_temp(k,:))
            location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
            ! call efem_psi_1s_truncated(location,at_c,Z,r0)
            call efem_psi_1s_truncated(location,at_c,Z,r0,t)
            V_temp=V_temp+phi_p2(m)*efem_phi(1)*phi_p2(n)*gweight_temp(k)
          enddo

          y_real=Jdet*V_temp/6.0d0 !!!!!!! new basis 


          x_imag=0.0d0
          y_imag=0.0d0
          x=cmplx(x_real,x_imag,8)
          y=cmplx(y_real,y_imag,8)
          mtest=.true.
          if (neigh_AB(ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
          if (mtest) call neigh_AB(ele(i,m))%insert(enriched_ele(i,n))
          call neigh_AB(ele(i,m))%insertAB(enriched_ele(i,n),x,y)
      enddo
  enddo


  end if

enddo

! print *,'done'




else if (Nlocal==20) then

  do i=1,Ne

    kk=0
    do ii=1,Nlocal
      if (ele(i,ii)==9) kk=kk+1
    enddo
  
    if (kk==0) then
  
    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    call mass_mat
    call gaussian_integral
    
    volume(i)=(Jdet)/6.0d0
  
    do m=1,Nlocal
        do n=1,Nlocal
            A_temp=0
            V_temp=0
            ik=0.0d0
            k2=0.0d0
            do k=1,Nquadrature
                call basis_p3_del(gpoint(k,:))
                ! call basis_p2_del(gpoint(k,:))
                ! call basis_p2(gpoint(k,:))
                ! call basis_p2(gpoint(k,:))
                A_temp = A_temp + dot_product(matmul(Jit,phi_p3_del(m,:)),matmul(Jit,phi_p3_del(n,:)))*gweight(k)
                ! V_temp0 = V_temp0 + potential9990(i0,k0)*basis_m0(m0)*basis_n0(n0)*gweight(k0)
            enddo
  
            x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
            y_real=Jdet*p3_matrix(m,n)
  
            x_imag=0.0d0
            y_imag=0.0d0
            
            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)
            
  
            
            mtest=.true.
            if (neigh_AB(ele(i,m))%find(ele(i,n)))  mtest=.false.
            if (mtest) call neigh_AB(ele(i,m))%insert(ele(i,n))
            call neigh_AB(ele(i,m))%insertAB(ele(i,n),x,y)
  
  
        enddo
    enddo
  
    else 
  
      ! print *,i
  
        call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
        call mass_mat
        call gaussian_integral
        call gaussian_integral45
        
        volume(i)=(Jdet)/6.0d0
  
        do m=1,Nlocal
            do n=1,Nlocal
                A_temp=0.0d0
                do k=1,Nquadrature_efem
                    call basis_p3_del(gpoint_temp(k,:))
                    A_temp = A_temp + dot_product(matmul(Jit,phi_p3_del(m,:)),matmul(Jit,phi_p3_del(n,:)))*gweight_temp(k)
                enddo
  
                x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                y_real=Jdet*p3_matrix(m,n)
                x_imag=0.0d0
                y_imag=0.0d0
                x=cmplx(x_real,x_imag,8)
                y=cmplx(y_real,y_imag,8)
                mtest=.true.
                if (neigh_AB(ele(i,m))%find(ele(i,n)))  mtest=.false.
                if (mtest) call neigh_AB(ele(i,m))%insert(ele(i,n))
                call neigh_AB(ele(i,m))%insertAB(ele(i,n),x,y)
            enddo
        enddo
  
        do m=1,Nlocal
            do n=1,Nlocal
  
                A_temp=0.0d0
                do k=1,Nquadrature_efem
                    call basis_p3_del(gpoint_temp(k,:)) 
                    location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                    ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                    ! call efem_basis_grad(location,at_c,Z,r0) 
                    call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                    call efem_basis_grad(location,at_c,Z,r0,t) 
                    A_temp = A_temp + dot_product(matmul(Jit,phi_p3_del(m,:)*efem_phi(1)+phi_p3(m)*efem_phi_grad(1,:))&
                                              ,matmul(Jit,phi_p3_del(n,:)*efem_phi(1)+phi_p3(n)*efem_phi_grad(1,:)))*gweight_temp(k)
                enddo
  
                x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
  
                
                V_temp=0.d00
                do k=1,Nquadrature_efem
                  call basis_p3(gpoint_temp(k,:))
                  location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                  ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                  call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                  V_temp=V_temp+efem_phi(1)*phi_p3(m)*efem_phi(1)*phi_p3(n)*gweight_temp(k)
                enddo
  
                y_real=Jdet*V_temp/6.0d0 !!!!!!! new basis 
  
  
                x_imag=0.0d0
                y_imag=0.0d0
                x=cmplx(x_real,x_imag,8)
                y=cmplx(y_real,y_imag,8)
                mtest=.true.
                if (neigh_AB(enriched_ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
                if (mtest) call neigh_AB(enriched_ele(i,m))%insert(enriched_ele(i,n))
                call neigh_AB(enriched_ele(i,m))%insertAB(enriched_ele(i,n),x,y)
            enddo
        enddo
  
  
  
        do m=1,Nlocal
          do n=1,Nlocal
  
              A_temp=0.0d0
              do k=1,Nquadrature_efem
                  call basis_p3_del(gpoint_temp(k,:)) 
                  location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                  ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                  ! call efem_basis_grad(location,at_c,Z,r0)
                  call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                  call efem_basis_grad(location,at_c,Z,r0,t) 
                  A_temp = A_temp + dot_product(matmul(Jit,phi_p3_del(m,:)*efem_phi(1)+phi_p3(m)*efem_phi_grad(1,:))&
                                               ,matmul(Jit,phi_p3_del(n,:)))*gweight_temp(k)
              enddo
  
              x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
  
              
              V_temp=0.d00
              do k=1,Nquadrature_efem
                call basis_p3(gpoint_temp(k,:))
                location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                V_temp=V_temp+efem_phi(1)*phi_p3(m)*phi_p3(n)*gweight_temp(k)
              enddo
  
              y_real=Jdet*V_temp/6.0d0 !!!!!!! new basis 
  
  
              x_imag=0.0d0
              y_imag=0.0d0
              x=cmplx(x_real,x_imag,8)
              y=cmplx(y_real,y_imag,8)
              mtest=.true.
              if (neigh_AB(enriched_ele(i,m))%find(ele(i,n)))  mtest=.false.
              if (mtest) call neigh_AB(enriched_ele(i,m))%insert(ele(i,n))
              call neigh_AB(enriched_ele(i,m))%insertAB(ele(i,n),x,y)
          enddo
      enddo
  
  
      do m=1,Nlocal
        do n=1,Nlocal
  
            A_temp=0.0d0
            do k=1,Nquadrature_efem
                call basis_p3_del(gpoint_temp(k,:)) 
                location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                ! call efem_basis_grad(location,at_c,Z,r0) 
                call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                call efem_basis_grad(location,at_c,Z,r0,t) 
                A_temp = A_temp + dot_product(matmul(Jit,phi_p3_del(m,:))&
                                             ,matmul(Jit,phi_p3_del(n,:)*efem_phi(1)+phi_p3(n)*efem_phi_grad(1,:)))*gweight_temp(k)
            enddo
  
            x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
  
            
            V_temp=0.d00
            do k=1,Nquadrature_efem
              call basis_p3(gpoint_temp(k,:))
              location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
              ! call efem_psi_1s_truncated(location,at_c,Z,r0)
              call efem_psi_1s_truncated(location,at_c,Z,r0,t)
              V_temp=V_temp+phi_p3(m)*efem_phi(1)*phi_p3(n)*gweight_temp(k)
            enddo
  
            y_real=Jdet*V_temp/6.0d0 !!!!!!! new basis 
  
  
            x_imag=0.0d0
            y_imag=0.0d0
            x=cmplx(x_real,x_imag,8)
            y=cmplx(y_real,y_imag,8)
            mtest=.true.
            if (neigh_AB(ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
            if (mtest) call neigh_AB(ele(i,m))%insert(enriched_ele(i,n))
            call neigh_AB(ele(i,m))%insertAB(enriched_ele(i,n),x,y)
        enddo
    enddo
  
  
    end if
  
  enddo

end if !!! p2/p3 selection


deallocate(location)
deallocate(gweight_temp)
deallocate(gpoint_temp)


! print *,'toto'

! deallocate(del_m)
! deallocate(del_n)
! deallocate(basis_m)
! deallocate(basis_n)


end subroutine set_linkedlistAB_efem








! subroutine efem_basis_grad(xyz,location,Z,r0,t)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: Z,r0,t
!   double precision :: r,pi,f
!   double precision,dimension(:),allocatable :: local_location

!   pi=4.0d0*atan(1.0d0)

!   call efem_psi_1s_grad(xyz,location,Z,r0,t)
!   call truncationfunction_h_grad(xyz,location,r0,t)

!   allocate(local_location(3))

!   local_location = (/0.0d0,0.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   if (r<r0) then
!     efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)
!   else
!     efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)*truncationfunction_h(xyz,location,local_location,r0,t)&
!                           +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(1,1:3)
!   end if


!   local_location = (/1.0d0,0.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   if (r<r0) then
!     efem_phi_grad(2,1:3) = efem_1s_grad(2,1:3)
!   else
!     efem_phi_grad(2,1:3) = efem_1s_grad(2,1:3)*truncationfunction_h(xyz,location,local_location,r0,t)&
!                           +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(2,1:3)
!   end if

!   local_location = (/0.0d0,1.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   if (r<r0) then
!     efem_phi_grad(3,1:3) = efem_1s_grad(3,1:3)
!   else
!     efem_phi_grad(3,1:3) = efem_1s_grad(3,1:3)*truncationfunction_h(xyz,location,local_location,r0,t)&
!                           +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(3,1:3)
!   end if

!   local_location = (/0.0d0,0.0d0,1.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   if (r<r0) then
!     efem_phi_grad(4,1:3) = efem_1s_grad(4,1:3)
!   else
!     efem_phi_grad(4,1:3) = efem_1s_grad(4,1:3)*truncationfunction_h(xyz,location,local_location,r0,t)&
!                           +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(4,1:3)
!   end if



! end subroutine efem_basis_grad






! subroutine efem_psi_1s_grad(xyz,location,Z,r0,t)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: Z,r0,t
!   double precision :: r,pi,f
!   double precision,dimension(:),allocatable :: local_location

!   pi=4.0d0*atan(1.0d0)

!   allocate(local_location(3))

!   local_location = (/0.0d0,0.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   efem_1s_grad(1,1:3)=(/Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(1)+local_location(1)-xyz(1))/r,&
!                         Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(2)+local_location(2)-xyz(2))/r,&
!                         Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(3)+local_location(3)-xyz(3))/r/)

!   local_location = (/1.0d0,0.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)
  
!   efem_1s_grad(2,1:3)=(/Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(1)+local_location(1)-xyz(1))/r,&
!                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(2)+local_location(2)-xyz(2))/r,&
!                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(3)+local_location(3)-xyz(3))/r/)

!   local_location = (/0.0d0,1.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   efem_1s_grad(3,1:3)=(/Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(1)+local_location(1)-xyz(1))/r,&
!                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(2)+local_location(2)-xyz(2))/r,&
!                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(3)+local_location(3)-xyz(3))/r/)

!   local_location = (/0.0d0,0.0d0,1.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   efem_1s_grad(4,1:3)=(/Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(1)+local_location(1)-xyz(1))/r,&
!                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(2)+local_location(2)-xyz(2))/r,&
!                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(3)+local_location(3)-xyz(3))/r/)

! end subroutine efem_psi_1s_grad


! subroutine truncationfunction_h_grad(xyz,location,r0,t)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: r0,t
!   double precision :: r,r_tilde,u1,u2
!   double precision,dimension(:),allocatable :: local_location

!   allocate(local_location(3))


!     local_location(1:3) = (/0.0d0,0.0d0,0.0d0/)

!     r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                                 +(xyz(3)-location(3)-local_location(3))**2)

!     r_tilde = 1.0d0-t*(r-r0)/r0

!     u1 = truncationfunction_u(r_tilde)
!     u2 = truncationfunction_u(1.0d0-r_tilde)

!     h_grad(1,1:3) = (/u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(1)+local_location(1)-xyz(1))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(2)+local_location(2)-xyz(2))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(3)+local_location(3)-xyz(3))/r/)

  

!     local_location(1:3) = (/1.0d0,0.0d0,0.0d0/)
    
!     r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                                 +(xyz(3)-location(3)-local_location(3))**2)

!     r_tilde = 1.0d0-t*(r-r0)/r0

!     u1 = truncationfunction_u(r_tilde)
!     u2 = truncationfunction_u(1.0d0-r_tilde)

!     h_grad(2,1:3) = (/u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(1)+local_location(1)-xyz(1))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(2)+local_location(2)-xyz(2))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(3)+local_location(3)-xyz(3))/r/)

  

!     local_location(1:3) = (/0.0d0,1.0d0,0.0d0/)
    
!     r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                                 +(xyz(3)-location(3)-local_location(3))**2)

!     r_tilde = 1.0d0-t*(r-r0)/r0

!     u1 = truncationfunction_u(r_tilde)
!     u2 = truncationfunction_u(1.0d0-r_tilde)

!     h_grad(3,1:3) = (/u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(1)+local_location(1)-xyz(1))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(2)+local_location(2)-xyz(2))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(3)+local_location(3)-xyz(3))/r/)



!     local_location(1:3) = (/0.0d0,0.0d0,1.0d0/)
    
!     r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                                 +(xyz(3)-location(3)-local_location(3))**2)

!     r_tilde = 1.0d0-t*(r-r0)/r0

!     u1 = truncationfunction_u(r_tilde)
!     u2 = truncationfunction_u(1.0d0-r_tilde)

!     h_grad(4,1:3) = (/u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(1)+local_location(1)-xyz(1))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(2)+local_location(2)-xyz(2))/r,&
!                       u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
!                       *(location(3)+local_location(3)-xyz(3))/r/)



! end subroutine truncationfunction_h_grad




! subroutine efem_psi_1s_truncated(xyz,location,Z,r0,t)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: Z,r0,t
!   double precision :: r,pi,f
!   double precision,dimension(:),allocatable :: local_location

!   pi=4.0d0*atan(1.0d0)

!   allocate(local_location(3))

!   local_location = (/0.0d0,0.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)
!   efem_phi(1) = 1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*truncationfunction_h(xyz,location,local_location,r0,t)

!   local_location = (/1.0d0,0.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)
!   efem_phi(2) = 1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*truncationfunction_h(xyz,location,local_location,r0,t)

!   local_location = (/0.0d0,1.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)
!   efem_phi(3) = 1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*truncationfunction_h(xyz,location,local_location,r0,t)

!   local_location = (/0.0d0,0.0d0,1.0d0/)
!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)
!   efem_phi(4) = 1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*truncationfunction_h(xyz,location,local_location,r0,t)

!   deallocate(local_location)

! end subroutine efem_psi_1s_truncated

! function truncationfunction_h(xyz,location,local_location,r0,t) result(h)
!   double precision,dimension(:),intent(in) :: xyz,location,local_location
!   double precision,intent(in)  :: r0,t
!   double precision :: r,r_tilde,h,u1,u2

!   r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!                               +(xyz(3)-location(3)-local_location(3))**2)

!   r_tilde = 1.0d0-t*(r-r0)/r0

!   u1 = truncationfunction_u(r_tilde)
!   u2 = truncationfunction_u(1.0d0-r_tilde)

!   h = u1/(u1+u2)

! end function truncationfunction_h

! function truncationfunction_u(r_tilde) result(u)
!   double precision,intent(in)  :: r_tilde
!   double precision :: u

!   if ((r_tilde<0.0d0).or.(r_tilde==0.0d0)) then
!     u = 0.0d0
!   else
!     u = exp(-1.0d0/r_tilde)
!   end if

! end function truncationfunction_u




! subroutine efem_basis_grad(xyz,location,Z,r0,t)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: Z,r0,t
!   double precision :: r,pi,f
!   ! double precision,dimension(:),allocatable :: local_location

!   pi=4.0d0*atan(1.0d0)

!   call efem_psi_1s_grad(xyz,location,Z,r0,t)
!   call truncationfunction_h_grad(xyz,location,r0,t)

!   ! allocate(local_location(3))

!   ! local_location = (/0.0d0,0.0d0,0.0d0/)
!   r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
!                               +(xyz(3)-location(3))**2)

!   if (r<r0) then
!     efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)
!   else
!     efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)*truncationfunction_h(xyz,location,r0,t)&
!                           +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(1,1:3)
!   end if


!   ! local_location = (/1.0d0,0.0d0,0.0d0/)
!   ! r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!   !                             +(xyz(3)-location(3)-local_location(3))**2)

!   ! if (r<r0) then
!   !   efem_phi_grad(2,1:3) = efem_1s_grad(2,1:3)
!   ! else
!   !   efem_phi_grad(2,1:3) = efem_1s_grad(2,1:3)*truncationfunction_h(xyz,location,local_location,r0,t)&
!   !                         +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(2,1:3)
!   ! end if

!   ! local_location = (/0.0d0,1.0d0,0.0d0/)
!   ! r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!   !                             +(xyz(3)-location(3)-local_location(3))**2)

!   ! if (r<r0) then
!   !   efem_phi_grad(3,1:3) = efem_1s_grad(3,1:3)
!   ! else
!   !   efem_phi_grad(3,1:3) = efem_1s_grad(3,1:3)*truncationfunction_h(xyz,location,local_location,r0,t)&
!   !                         +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(3,1:3)
!   ! end if

!   ! local_location = (/0.0d0,0.0d0,1.0d0/)
!   ! r = sqrt((xyz(1)-location(1)-local_location(1))**2+(xyz(2)-location(2)-local_location(2))**2&
!   !                             +(xyz(3)-location(3)-local_location(3))**2)

!   ! if (r<r0) then
!   !   efem_phi_grad(4,1:3) = efem_1s_grad(4,1:3)
!   ! else
!   !   efem_phi_grad(4,1:3) = efem_1s_grad(4,1:3)*truncationfunction_h(xyz,location,local_location,r0,t)&
!   !                         +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(4,1:3)
!   ! end if



! end subroutine efem_basis_grad
















!!!!!! Extended Finite Element function truncation function 

! !!!!!!!!!!! J.E.Pask paper 
! !!!!!!!!!!! https://www.sciencedirect.com/science/article/pii/S2352431616302048
! !!!!!!!!!!! https://ars.els-cdn.com/content/image/1-s2.0-S2352431616302048-mmc1.pdf

! subroutine efem_psi_1s_truncated(xyz,location,Z,r0)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: Z,r0
!   double precision :: r,pi,f

!   pi=4.0d0*atan(1.0d0)

!   r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
!                               +(xyz(3)-location(3))**2)
!   efem_phi(1) = 1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*truncationfunction_h(xyz,location,r0)

! end subroutine efem_psi_1s_truncated

! subroutine efem_basis_grad(xyz,location,Z,r0)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in) :: Z,r0
!   double precision :: r,pi,f

!   pi=4.0d0*atan(1.0d0)

!   call efem_psi_1s_grad(xyz,location,Z,r0)
!   call truncationfunction_h_grad(xyz,location,r0)

!   r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
!                               +(xyz(3)-location(3))**2)

!   if (r<r0) then
!     efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)
!   else
!     efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)*truncationfunction_h(xyz,location,r0)&
!                           +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(1,1:3)
!   end if

! end subroutine efem_basis_grad

! subroutine efem_psi_1s_grad(xyz,location,Z,r0)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: Z,r0
!   double precision :: r,pi,f

!   pi=4.0d0*atan(1.0d0)

!   r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
!                               +(xyz(3)-location(3))**2)

!   efem_1s_grad(1,1:3)=(/Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(1)-xyz(1))/r,&
!                         Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(2)-xyz(2))/r,&
!                         Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(3)-xyz(3))/r/)

! end subroutine efem_psi_1s_grad

! function truncationfunction_h(xyz,location,r0) result(h)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: r0
!   double precision :: r,h

!   r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
!                               +(xyz(3)-location(3))**2)

!   if ( r < r0 ) then
  
!     h = 1.0d0+20.0d0*r**7.0d0/r0**7.0d0-70.0d0*r**6.0d0/r0**6.0d0+84.0d0*r**5.0d0/r0**5.0d0-35.0d0*r**4.0d0/r0**4.0d0

!   else

!     h = 0.0d0

!   end if
  

! end function truncationfunction_h

! subroutine truncationfunction_h_grad(xyz,location,r0)
!   double precision,dimension(:),intent(in) :: xyz,location
!   double precision,intent(in)  :: r0
!   double precision :: r,temp

!     r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
!                                 +(xyz(3)-location(3))**2)

!     temp = 7.0d0*20.0d0*r**5.0d0/r0**7.0d0-6.0d0*70.0d0*r**4.0d0/r0**6.0d0&
!           +5.0d0*84.0d0*r**3.0d0/r0**5.0d0-4.0d0*35.0d0*r**2.0d0/r0**4.0d0

!     if ( r < r0 ) then 

!       h_grad(1,1:3) = (/temp*(xyz(1)-location(1)),temp*(xyz(2)-location(2)),temp*(xyz(3)-location(3))/)

!     else 

!       h_grad(1,1:3) = (/0.0d0,0.0d0,0.0d0/)

!     end if

! end subroutine truncationfunction_h_grad




!!!!!!!!!!! Michigan paper 
!!!!!!!!!!! https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.035112
!!!!!!!!!!! https://journals.aps.org/prb/abstract/10.1103/PhysRevB.104.085112

subroutine efem_psi_1s_truncated(xyz,location,Z,r0,t)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z,r0,t
  double precision :: r,pi,f

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
                              +(xyz(3)-location(3))**2)
  efem_phi(1) = 1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*truncationfunction_h(xyz,location,r0,t)

end subroutine efem_psi_1s_truncated

subroutine efem_basis_grad(xyz,location,Z,r0,t)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in) :: Z,r0,t
  double precision :: r,pi,f

  pi=4.0d0*atan(1.0d0)

  call efem_psi_1s_grad(xyz,location,Z)
  call truncationfunction_h_grad(xyz,location,r0,t)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
                              +(xyz(3)-location(3))**2)

  if (r<r0) then
    efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)
  else
    efem_phi_grad(1,1:3) = efem_1s_grad(1,1:3)*truncationfunction_h(xyz,location,r0,t)&
                          +1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*h_grad(1,1:3)
  end if

end subroutine efem_basis_grad

subroutine efem_psi_1s_grad(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
                              +(xyz(3)-location(3))**2)

  efem_1s_grad(1,1:3)=(/Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(1)-xyz(1))/r,&
                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(2)-xyz(2))/r,&
                        Z*1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)*(location(3)-xyz(3))/r/)

end subroutine efem_psi_1s_grad

function truncationfunction_h(xyz,location,r0,t) result(h)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: r0,t
  double precision :: r,r_tilde,h,u1,u2

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
                              +(xyz(3)-location(3))**2)

  r_tilde = 1.0d0-t*(r-r0)/r0

  u1 = truncationfunction_u(r_tilde)
  u2 = truncationfunction_u(1.0d0-r_tilde)

  h = u1/(u1+u2)

end function truncationfunction_h

function truncationfunction_u(r_tilde) result(u)
  double precision,intent(in)  :: r_tilde
  double precision :: u

  if ((r_tilde<0.0d0).or.(r_tilde==0.0d0)) then
    u = 0.0d0
  else
    u = exp(-1.0d0/r_tilde)
  end if

end function truncationfunction_u

subroutine truncationfunction_h_grad(xyz,location,r0,t)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: r0,t
  double precision :: r,r_tilde,u1,u2

    r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2&
                                +(xyz(3)-location(3))**2)

    r_tilde = 1.0d0-t*(r-r0)/r0

    u1 = truncationfunction_u(r_tilde)
    u2 = truncationfunction_u(1.0d0-r_tilde)

    h_grad(1,1:3) = (/u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
                      *(location(1)-xyz(1))/r,&
                      u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
                      *(location(2)-xyz(2))/r,&
                      u1*u2/(u1+u2)**2 * (2.0d0*r_tilde**2-2.0d0*r_tilde+1.0d0)/(r_tilde*(1.0d0-r_tilde))**2*(-t/r0)&
                      *(location(3)-xyz(3))/r/)

end subroutine truncationfunction_h_grad








subroutine psi_1s(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  psi_square_100 = 1.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*exp(-Z*r)

end subroutine psi_1s

subroutine psi_2s(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  psi_square_200 = 1.0d0/sqrt(32.0d0*pi)*Z**(3.0d0/2.0d0)*(2.0d0-Z*r)*exp(-Z*r/2.0d0)

end subroutine psi_2s

subroutine psi_2x(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi
  double precision :: sin_theta,cos_phi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  if (r-0.0d0<1E-8) then
      sin_theta = 1.0d0
      cos_phi = 1.0d0
  else
      sin_theta = sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)/r
      cos_phi = (location(1)-xyz(1))/sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)
  end if

  if (sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)-0.0d0<1E-8) cos_phi = 0.0d0
      

  psi_square_21x = 1.0d0/sqrt(32.0d0*pi)*Z**(3.0d0/2.0d0)*(Z*r)*exp(-Z*r/2.0d0)*sin_theta*cos_phi

end subroutine psi_2x

subroutine psi_2y(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi
  double precision :: sin_theta,sin_phi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  if (r-0.0d0<1E-8) then
      sin_theta = 1.0d0
      sin_phi = 1.0d0
  else
      sin_theta = sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)/r
      sin_phi = (location(2)-xyz(2))/sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)
  end if

  if (sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)-0.0d0<1E-8) sin_phi = 0.0d0

  psi_square_21y = 1.0d0/sqrt(32.0d0*pi)*Z**(3.0d0/2.0d0)*(Z*r)*exp(-Z*r/2.0d0)*sin_theta*sin_phi

end subroutine psi_2y

subroutine psi_2z(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi
  double precision :: cos_theta

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  if (r-0.0d0<1E-8) then
      cos_theta = 1.0d0
  else
      cos_theta = (location(3)-xyz(3))/r
  end if

  psi_square_21z = 1.0d0/sqrt(32.0d0*pi)*Z**(3.0d0/2.0d0)*(Z*r)*exp(-Z*r/2.0d0)*cos_theta

end subroutine psi_2z


subroutine psi_3s(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  psi_square_300 = 1.0d0/(81.0d0*sqrt(3.0d0*pi))*Z**(3.0d0/2.0d0)*(27.0d0-18.0d0*(Z*r)+2.0d0*(Z*r)**2)*exp(-Z*r/3.0d0)

end subroutine psi_3s

subroutine psi_3x(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi
  double precision :: sin_theta,cos_phi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  if (r-0.0d0<1E-8) then
      sin_theta = 1.0d0
      cos_phi = 1.0d0
  else
      sin_theta = sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)/r
      cos_phi = (location(1)-xyz(1))/sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)
  end if

  if (sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)-0.0d0<1E-8) cos_phi = 0.0d0
      

  psi_square_31x = 1.0d0/81.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*(6.0d0-Z*r)*(Z*r)*exp(-Z*r/3.0d0)*sin_theta*cos_phi

end subroutine psi_3x

subroutine psi_3y(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi
  double precision :: sin_theta,sin_phi

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  if (r-0.0d0<1E-8) then
      sin_theta = 1.0d0
      sin_phi = 1.0d0
  else
      sin_theta = sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)/r
      sin_phi = (location(2)-xyz(2))/sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)
  end if

  if (sqrt((location(1)-xyz(1))**2+(location(2)-xyz(2))**2)-0.0d0<1E-8) sin_phi = 0.0d0

  psi_square_31y = 1.0d0/81.0d0/sqrt(pi)*Z**(3.0d0/2.0d0)*(6.0d0-Z*r)*(Z*r)*exp(-Z*r/3.0d0)*sin_theta*sin_phi

end subroutine psi_3y

subroutine psi_3z(xyz,location,Z)
  double precision,dimension(:),intent(in) :: xyz,location
  double precision,intent(in)  :: Z
  double precision :: r,pi
  double precision :: cos_theta

  pi=4.0d0*atan(1.0d0)

  r = sqrt((xyz(1)-location(1))**2+(xyz(2)-location(2))**2+(xyz(3)-location(3))**2)

  if (r-0.0d0<1E-8) then
      cos_theta = 1.0d0
  else
      cos_theta = (location(3)-xyz(3))/r
  end if

  psi_square_31z = 1.0d0/81.0d0*sqrt(2.0d0/pi)*Z**(3.0d0/2.0d0)*(6.0d0-Z*r)*(Z*r)*exp(-Z*r/3.0d0)*cos_theta

end subroutine psi_3z


function distance(location1,location2) result(d0)
    
    double precision,dimension(:),intent(in) :: location1,location2
    double precision :: d0

    d0=sqrt((location1(1)-location2(1))**2+(location1(2)-location2(2))**2+(location1(3)-location2(3))**2)

end function distance


subroutine build_quaternionmatrix(Nn)
  integer, intent(in) :: Nn
  integer :: i,l

  do i=1,Nn
      do l=IA(i),IA(i+1)-1
          call quaternion((i-1)*4+1)%insert((JA(l)-1)*4+1)
          call quaternion((i-1)*4+1)%insert((JA(l)-1)*4+2)
          call quaternion((i-1)*4+1)%insert((JA(l)-1)*4+3)
          call quaternion((i-1)*4+1)%insert((JA(l)-1)*4+4)

          call quaternion((i-1)*4+2)%insert((JA(l)-1)*4+1)
          call quaternion((i-1)*4+2)%insert((JA(l)-1)*4+2)
          call quaternion((i-1)*4+2)%insert((JA(l)-1)*4+3)
          call quaternion((i-1)*4+2)%insert((JA(l)-1)*4+4)

          call quaternion((i-1)*4+3)%insert((JA(l)-1)*4+1)
          call quaternion((i-1)*4+3)%insert((JA(l)-1)*4+2)
          call quaternion((i-1)*4+3)%insert((JA(l)-1)*4+3)
          call quaternion((i-1)*4+3)%insert((JA(l)-1)*4+4)
          
          call quaternion((i-1)*4+4)%insert((JA(l)-1)*4+1)
          call quaternion((i-1)*4+4)%insert((JA(l)-1)*4+2)
          call quaternion((i-1)*4+4)%insert((JA(l)-1)*4+3)
          call quaternion((i-1)*4+4)%insert((JA(l)-1)*4+4)
      enddo
  enddo

end subroutine build_quaternionmatrix


end module basisfunctions
