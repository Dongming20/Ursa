module potentials_gw

    use class_linkedlist
    use tools
    use basisfunctions
    use potentials

    implicit none


    
    complex(kind=(kind(1.0d0))),allocatable :: j_imag,j_real,j_zero
    double precision,dimension(:,:),allocatable :: v_kernel,v00,Vhfx_GW,Vhfx_final_GW
    double precision,dimension(:,:),allocatable :: NNmatrix_temp02,NNmatrix_temp03,NNmatrix_temp04
    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: NNmatrix_temp005,NNmatrix_temp,NNmatrix_temp006

    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: bwork1_vectemp1,bwork1_vectemp2,bwork1_vectemp3,bwork1_vectemp4
complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: bwork1_vectemp1_g,bwork1_vectemp2_g,bwork1_vectemp3_g,bwork1_vectemp4_g
    complex(kind=(kind(1.0d0))),dimension(:,:,:),allocatable :: bwork1_backup
    double precision,dimension(:),allocatable :: vectortemp
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: ztemp
    double precision :: vtemp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!! Contour Deformation related variables !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!! scalar !!!!
    double precision, allocatable :: E_guess_min,E_guess_max
    integer, allocatable :: Nguess_complex,N_energy
    double precision, allocatable :: alpha_gw
    double precision, allocatable :: eta_omega1,eta_omega2,eta_time

    !!!! array !!!!
    double precision,dimension(:),allocatable :: E_guess,E_Sigma
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: E_guess_complex,E_Sigma_complex
    double precision,dimension(:),allocatable :: omega_min,omega_max,Xi_min,Xi_max
    double precision,dimension(:,:),allocatable :: omega_range,Xi_range
    ! complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: omega_range,Xi_range

    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: E_guess_complex_HOMO,E_guess_complex_LUMO
    
    
    complex(kind=(kind(1.0d0))),dimension(:,:,:),allocatable :: gw_G_complex,gw_W_complex
    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: IdentityMatrix,gw_Sigma_complex,chi_matrix
    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: gw_Sigma_complex_poles
    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: gw_G_g,gw_W_g,gw_Sigma_g,gw_Sigma_g_poles

    double precision,dimension(:,:),allocatable :: NgNntemp,NnNgtemp,chi_matrix_real,NgNgtemp,NnNgtemp2
    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: zNgNntemp,zNnNgtemp,zNnNgtemp2,zNgNgtemp,zNgNgtemp2

    ! complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: zNgNntemp

    complex(kind=(kind(1.0d0))),dimension(:,:,:),allocatable :: Sigma_c_temp,Sigma_c_temp_g

    complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: SigmaY,Y_primed

    double precision,dimension(:,:),allocatable :: IdentityMatrix_real

    !!!!!! quaternion !!!!!!
    double precision,dimension(:,:),allocatable :: IdenMat_quaternion,G_quaternion

    !!!!!! contour test !!!!!!!
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: contourpoints



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!! Full analytic solution (casida) related variables !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    double precision,dimension(:,:),allocatable :: psi_ia_g_real,psi_n_g_real,psi_n_g_real_occupied!,NgNntemp
    double precision,dimension(:,:),allocatable :: casida_Kx,casida_Kxnj,casida_Kxnj_occupied,casida_Xs,casida_Ys
    double precision,dimension(:),allocatable :: casida_R,casida_omega,casida_R_half,casida_cvomega
    ! complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: casida_Xs,casida_Ys
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: casida_Ec,casida_Ec_HOMO,casida_Ec_LUMO
    double precision,dimension(:,:),allocatable :: NNmatrix_temp12,NNmatrix_temp13



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!     Betha Salpeter Equation related variables     !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    double precision,dimension(:),allocatable :: psi_aa_g_real
    double precision,dimension(:,:),allocatable :: psi_ii_g_real,worktemp0,worktemp00

    double precision,dimension(:,:),allocatable :: gw_W_real

    double precision,dimension(:,:),allocatable :: bse_W,bse_C,bse_H
    double precision,dimension(:),allocatable :: bse_omega
    double precision,dimension(:),allocatable :: E_GW

    double precision,dimension(:,:),allocatable :: bse_Kx

    ! allocate(psi_ia_g_real(Ne*Nquadrature,Nstates*Nempty))
    ! allocate(psi_n_g_real(Ne*Nquadrature,Nempty))
    ! allocate(psi_n_g_real_occupied(Ne*Nquadrature,Nstates))
    ! allocate(NgNntemp(1:Ne*Nquadrature,1:Nstates*Nempty))
    ! ! allocate(NnNgtemp(1:Nn,1:Nstates*Nempty))

    ! allocate(casida_Kx(1:Nstates*Nempty,1:Nstates*Nempty))
    ! allocate(casida_R(1:Nstates*Nempty))
    ! allocate(casida_Xs(1:Nstates*Nempty,1:Nstates*Nempty))
    ! allocate(casida_omega(1:Nstates*Nempty))
    ! allocate(casida_R_half(1:Nstates*Nempty))
    ! ! allocate(casida_V(1:Nstates*Nempty))
    ! allocate(casida_Kxnj(1:Nempty,1:Nstates*Nempty))
    ! allocate(casida_Kxnj_occupied(1:Nstates,1:Nstates*Nempty))
    ! allocate(casida_cvomega(1:Nstates*Nempty))


    ! allocate(NNmatrix_temp12(1:Nn,1:Nstates*Nempty))
    ! allocate(NNmatrix_temp13(1:Nn,1:Nstates*Nempty))

    ! ! allocate(Nguess_complex)
    ! ! allocate(E_guess_complex(1:Nguess_complex))
    ! allocate(casida_Ec(1:Nguess_complex))



    !!!! Matrix inversion using linear solver
    INTEGER          lapack_LWMAX
    ! PARAMETER        ( lapack_LWMAX = 100 )
    INTEGER          lapack_INFO, lapack_LWORK
    DOUBLE PRECISION,dimension(:),allocatable :: lapack_WORK
    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: lapack_WORKc
    Integer,dimension(:),allocatable :: lapack_IPIV

    double precision,dimension(:),allocatable :: snq1,snq2
    integer, allocatable :: ii10,ii11,ii30



    double precision :: conditionnumber

    !!! for PARDISO
    integer, allocatable :: MAXFCT,MNUM,MTYPE,MSGLVL,PHASE,idum,pardiso_info,MTYPE_real
    integer(8),dimension(:),allocatable :: apt_V
    integer,dimension(:),allocatable :: aiparm_V

    integer(8),dimension(:,:),allocatable :: PT
    integer,dimension(:,:),allocatable :: IPARM



    double precision, dimension(:), allocatable :: gpoint_1D,gpoint_1D_3points
    double precision, dimension(:), allocatable :: gweight_1D,gweight_1D_3points
    double precision, dimension(:), allocatable :: gpoint3_1D   
    double precision, dimension(:), allocatable :: gweight3_1D  
    double precision, dimension(:), allocatable :: gpoint20_1D  
    double precision, dimension(:), allocatable :: gweight20_1D 
    double precision, dimension(:), allocatable :: gpoint9_1D   
    double precision, dimension(:), allocatable :: gweight9_1D  
    double precision, dimension(:), allocatable :: gpoint11_1D  
    double precision, dimension(:), allocatable :: gweight11_1D 
    double precision, dimension(:), allocatable :: gpoint7_1D   
    double precision, dimension(:), allocatable :: gweight7_1D  
    double precision, dimension(:), allocatable :: gpoint15_1D  
    double precision, dimension(:), allocatable :: gweight15_1D 
    double precision, dimension(:), allocatable :: gpoint6_1D   
    double precision, dimension(:), allocatable :: gweight6_1D  
    double precision, dimension(:), allocatable :: gpoint8_1D   
    double precision, dimension(:), allocatable :: gweight8_1D  
    double precision, dimension(:), allocatable :: gpoint5_1D   
    double precision, dimension(:), allocatable :: gweight5_1D  
    double precision, dimension(:), allocatable :: gpoint10_1D  
    double precision, dimension(:), allocatable :: gweight10_1D 
    double precision, dimension(:), allocatable :: gpoint16_1D  
    double precision, dimension(:), allocatable :: gweight16_1D 
    double precision, dimension(:), allocatable :: gpoint17_1D  
    double precision, dimension(:), allocatable :: gweight17_1D 
    double precision, dimension(:), allocatable :: gpoint18_1D  
    double precision, dimension(:), allocatable :: gweight18_1D 
    double precision, dimension(:), allocatable :: gpoint19_1D  
    double precision, dimension(:), allocatable :: gweight19_1D 
    double precision, dimension(:), allocatable :: gpoint24_1D  
    double precision, dimension(:), allocatable :: gweight24_1D 
    double precision, dimension(:), allocatable :: gpoint32_1D  
    double precision, dimension(:), allocatable :: gweight32_1D 
    double precision, dimension(:), allocatable :: gpoint40_1D  
    double precision, dimension(:), allocatable :: gweight40_1D 
    double precision, dimension(:), allocatable :: gpoint48_1D  
    double precision, dimension(:), allocatable :: gweight48_1D 

    complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Vhfx_test

    contains

    subroutine construct_v_kernel(Nn,Ne,Nquadrature)

        integer,intent(in) :: Nn,Ne,Nquadrature
        integer :: i,j,k

        do i=1,Nn
            do j=1,Ne
            do k=1,Nquadrature

            v_kernel((j-1)*Nquadrature+k,i)=1.0d0/sqrt((point(i,1)-point_g((j-1)*Nquadrature+k,1))**2&
                                                      +(point(i,2)-point_g((j-1)*Nquadrature+k,2))**2&
                                                      +(point(i,3)-point_g((j-1)*Nquadrature+k,3))**2)

            ! if (isnan(v_kernel((j-1)*Nquadrature+k,i))) stop 'NaN number happened'
            ! if (v_kernel((j-1)*Nquadrature+k,i)-1.0d0==v_kernel((j-1)*Nquadrature+k,i)) stop 'infinity happened'
            enddo
            enddo
        enddo

    end subroutine construct_v_kernel

    subroutine construct_exchangekernel(Nstates,Nn,Ne,Nquadrature)

        integer,intent(in) :: Nstates,Nn,Ne,Nquadrature
        ! double precision, dimension(:,:), intent(in) :: psi_temp
        integer :: k
        double precision, dimension(:,:),allocatable :: A

        allocate(A(size(psi_point_g(:,k)),size(psi(:,k))))
        
        Vhfx_GW(:,:)=0.0d0
        do k=1,Nstates
            ! Vhfx_GW(:,:)=Vhfx_GW(:,:)-outer_product(psi_GW_g(:,k),psi_GW(:,k))
            call outer_product('ge',psi_point_g(:,k),psi(:,k),A)
            Vhfx_GW(:,:)=Vhfx_GW(:,:)-A
            ! Vhfx_GW(1:Ne*Nquadrature,1:Nn)=Vhfx_GW(1:Ne*Nquadrature,1:Nn)-outer_product(psi_point_g(1:Ne*Nquadrature,k),psi(1:Nn,k))
        enddo

        Vhfx_GW=Vhfx_GW*v_kernel

        deallocate(A)

    end subroutine construct_exchangekernel


    subroutine construct_exchangematrix(Nn,Ne,Nquadrature,Nstates)

        integer,intent(in) :: Nn,Ne,Nquadrature,Nstates

        double precision,dimension(:,:), allocatable :: NNmatrix_temp02,NNmatrix_temp03

        allocate(Vhfx_GW(Ne*Nquadrature,Nn))
        allocate(NNmatrix_temp02(Nn,Nn))
        allocate(NNmatrix_temp03(Nn,Nn))

        call construct_exchangekernel(Nstates,Nn,Ne,Nquadrature)


        call mkl_dcsrmm('T',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,Vhfx_GW,Ne*Nquadrature,0.0d0,NNmatrix_temp02,Nn)

        NNmatrix_temp03=transpose(NNmatrix_temp02)

        call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,NNmatrix_temp03,Nn,0.0d0,Vhfx_final_GW,Nn)

        NNmatrix_temp02=transpose(Vhfx_final_GW)
        Vhfx_final_GW=(Vhfx_final_GW+NNmatrix_temp02)/2.0d0

        deallocate(Vhfx_GW)
        deallocate(NNmatrix_temp02)
        deallocate(NNmatrix_temp03)


    end subroutine construct_exchangematrix





    function dense_matrix_inverse(n,A) result(C)
        integer,intent(in) :: n
        complex(kind=(kind(1.0d0))),dimension(n,n), intent(in) :: A
        ! complex(kind=(kind(1.0d0))),dimension(n,n), intent(out) :: B

        integer :: i
        complex(kind=(kind(1.0d0))),dimension(n,n) :: C
        complex(kind=(kind(1.0d0))),dimension(:),allocatable :: work
        integer,dimension(:),allocatable :: ipiv
        integer :: lwork,info

        C = ZZERO
        do i=1,n
            C(i,i)=j_real
        enddo

        lwork=n
        allocate(work(2*lwork))
        allocate(ipiv(n))

        call zsytrf('L',n,A,n,ipiv,work,lwork,info)
        ! call zsytri('L',n,A,n,ipiv,work,info)
        call zsytrs('L',n,n,A,n,ipiv,C,n,info)

    end function dense_matrix_inverse







    subroutine compute_gw_W_CD(Nn,Nstates,Nempty,E_dft,omega,Nquadrature1D_range,Nquadrature1D,nnza,Ne,Nquadrature,Ze,type,hadamard)
        integer,intent(in) :: Nn,Nstates,Nempty,Nquadrature1D_range,Nquadrature1D,nnza,Ne,Nquadrature
        double precision,dimension(1:Nstates),intent(in) :: E_dft
        double precision,dimension(1:Nquadrature1D,1:Nquadrature1D_range),intent(in) :: omega
        ! complex(kind=(kind(1.0d0))),dimension(1:Nquadrature1D,1:Nquadrature1D_range),intent(in) :: omega
        complex(kind=(kind(1.0d0))),intent(in) :: Ze
        character(len=*), intent(in) :: type,hadamard

        integer :: i,ii,k,kk,l,ll
        complex(kind=(kind(1.0d0))) :: valuetemp
        !!!!! pardiso parameters
        integer :: MAXFCT,MNUM,MTYPE,MTYPE_real,MSGLVL,PHASE,idum,pardiso_info
        integer(8),dimension(:),allocatable :: pt_V
        integer,dimension(:),allocatable :: iparm_V

        complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: matrix_temp,znnmatrix,G
        double precision,dimension(:,:),allocatable :: matrix_temp2,dnnmatrix,dnnmatrix2,psi_ia
        double precision,dimension(:,:),allocatable :: NgNntemp,NgNgtemp,NgNgtemp2,NnNgtemp,NnNgtemp2
        complex(kind=(kind(1.0d0))),dimension(:),allocatable :: lapack_work


        gw_W_complex=cmplx(0.0d0,0.0d0,8)

        allocate(lapack_IPIV(1:Nn))

        if (type=='dft') then

        !!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! pardiso parameters !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        MAXFCT=1
        MNUM=1
        MTYPE=6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
        MTYPE_real=1
        MSGLVL=0
        allocate(pt_v(64))
        allocate(iparm_v(64))
        pt_V=0
        call pardisoinit(PT_V,MTYPE,IPARM_V)
        PHASE=13

        allocate(dnnmatrix(Nn,Nn))
        allocate(dnnmatrix2(Nn,Nn))
        allocate(G(Nn,Nn))

        ! allocate(matrix_temp(Nn,Nn))
        ! allocate(znnmatrix(Nn,Nn))

        if (hadamard=='no') then


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do ii=1,Nquadrature1D_range
            do i=1,Nquadrature1D

                dnnmatrix=0.0d0
                ! znnmatrix=ZZERO
                ! NnNgtemp2=0.0d0
    
                ! call pardisoinit(PT_V,MTYPE,IPARM_V)
    
                do k=1,Nstates

                    ! matrix_temp=ZZERO
    
                    ! G_inverse = (E_dft(k)+omega_range(i,ii)*j_imag)*B-H_dft       !!!! imaginary axis
    
                    ! call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,G_inverse,IA,JA,idum,Nn&
                    !                     ,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

                    psaz_dft(1:nnza)=(E_dft(k)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
                    ! psaz_dft(1:nnza)=(E_dft(k)+omega(i,ii)*j_imag+eta_omega1*j_imag+eta_omega2*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
                    ! psaz_dft(1:nnza)=(E_dft(k)+omega(i,ii))*usb(1:nnza)-usa_dft(1:nnza)
call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

! do l=1,Nstates
!     G=G+outer_product(psi(:,l),psi(:,l))*2.0d0*eta_omega2*j_imag/((E_dft(k)+omega(i,ii)*j_imag+eta_omega1*j_imag-E_dft(l))**2&
!     +eta_omega2**2)
! enddo

! matrix_temp=G

!                     psaz_dft(1:nnza)=(E_dft(k)-omega(i,ii)*j_imag+eta_omega1*j_imag+eta_omega2*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
!                     ! psaz_dft(1:nnza)=(E_dft(k)+omega(i,ii))*usb(1:nnza)-usa_dft(1:nnza)
! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

! do l=1,Nstates
!     G=G+outer_product(psi(:,l),psi(:,l))*2.0d0*eta_omega2*j_imag/((E_dft(k)-omega(i,ii)*j_imag+eta_omega1*j_imag-E_dft(l))**2&
!     +eta_omega2**2)
! enddo

! matrix_temp=matrix_temp+G

! print *,G(1,1)
! print *,G(1,2)
! print *,G(2,1)
! print *,G(2,2)

! G=ZZERO

! do l=1,Nn
!     if (l>Nstates) then
!         G=G+outer_product(psi(:,l),psi(:,l))/(E_dft(k)+omega(i,ii)*j_imag+eta_omega1*j_imag-E_dft(l))
!     else
!         G=G+outer_product(psi(:,l),psi(:,l))/(E_dft(k)+omega(i,ii)*j_imag-eta_omega1*j_imag-E_dft(l))
!     end if
! enddo

! print *,G(1,1)
! print *,G(1,2)
! print *,G(2,1)
! print *,G(2,2)

! stop


!     call mkl_dcsrmm('N',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                                     ,dble(G),Nn,0.0d0,NgNntemp,Ne*Nquadrature)

! call mkl_dcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                                     ,transpose(NgNntemp),Nn,0.0d0,NgNgtemp,Ne*Nquadrature)

! NgNgtemp=NgNgtemp*outer_product(psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k))

! call mkl_dcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,4.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                                     ,NgNgtemp,Ne*Nquadrature,0.0d0,NnNgtemp,Nn)


        ! NnNgtemp2=NnNgtemp2+NnNgtemp
                                    
                    call outer_product('ge',psi(:,k),psi(:,k),dnnmatrix2)
                    dnnmatrix=dnnmatrix+4.0d0*dble(G)*dnnmatrix2!outer_product(psi(:,k),psi(:,k))
                    ! znnmatrix=znnmatrix+2.0d0*matrix_temp*outer_product(psi(:,k),psi(:,k))
                
                enddo


            !     print *,znnmatrix(1,1)
            !     print *,znnmatrix(1,2)
            !     print *,znnmatrix(2,1)
            !     print *,znnmatrix(2,2)


            ! allocate(psi_ia(Nn,Nstates*Nempty))
            ! do k=1,Nstates
            !     do kk=1,Nempty
            !         psi_ia(:,(k-1)*Nempty+kk)=psi(:,k)*psi(:,Nstates+kk)
            !     enddo
            ! enddo


            ! znnmatrix=ZZERO

            ! do k=1,Nstates
            !     do kk=1,Nempty

            !         znnmatrix=znnmatrix&
            !             +outer_product(psi_ia(:,(k-1)*Nempty+kk),psi_ia(:,(k-1)*Nempty+kk))&
            !             *(1.0d0/(omega(i,ii)*j_imag-(E_dft(Nstates+kk)-E_dft(k))+j_imag*eta_omega2+j_imag*eta_omega1)&
            !              -1.0d0/(omega(i,ii)*j_imag+(E_dft(Nstates+kk)-E_dft(k))-j_imag*eta_omega2-j_imag*eta_omega1))
            !             ! *(1.0d0/(omega_range(i,ii)*j_imag-(E_GW(Nstates+kk)-E_GW(k)))&
            !             !  -1.0d0/(omega_range(i,ii)*j_imag+(E_GW(Nstates+kk)-E_GW(k))))
                        
                        
            !     enddo
            ! enddo

            ! znnmatrix=2.0d0*znnmatrix

            ! print *,znnmatrix(1,1)
            ! print *,znnmatrix(1,2)
            ! print *,znnmatrix(2,1)
            ! print *,znnmatrix(2,2)

            ! stop


! CALL DGEMM('T','T',Nn,Nn,Ne*Nquadrature,1.0d0,v_kernel,Ne*Nquadrature,NnNgtemp2,Nn,0.0d0,NNmatrix_temp03,Nn)
    
                CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,dnnmatrix,Nn,0.0d0,dnnmatrix2,Nn)
    
                call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,dnnmatrix2,Nn&
                                                                        ,0.0d0,dnnmatrix,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                dnnmatrix2=IdentityMatrix_real-dnnmatrix
                
                dnnmatrix=v00
    
                CALL DGETRF(Nn,Nn,dnnmatrix2,Nn,lapack_IPIV,lapack_INFO)
                CALL DGETRS('N',Nn,Nn,dnnmatrix2,Nn,lapack_IPIV,dnnmatrix,Nn,lapack_INFO)
    
                dnnmatrix2=dnnmatrix-v00

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


                call pardiso(pt_v,MAXFCT,MNUM,MTYPE_real,PHASE,Nn,B,IA,JA,idum,Nn&
                                        ,IPARM_V,MSGLVL,dnnmatrix2,dnnmatrix,pardiso_info)

                ! NNmatrix_temp02=(NNmatrix_temp02+transpose(NNmatrix_temp02))/2.0d0
    
                gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)=dnnmatrix*j_real

                print *,(ii-1)*Nquadrature1D+i

                ! print *,dnnmatrix( 9, 9)
                ! print *,dnnmatrix( 9,10)
                ! print *,dnnmatrix(10, 9)
                ! print *,dnnmatrix(10,10)
    
    
            enddo
        enddo




        else if (hadamard=='yes') then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(NgNntemp(Ne*Nquadrature,Nn))
        allocate(NgNgtemp(Ne*Nquadrature,Ne*Nquadrature))
        allocate(NgNgtemp2(Ne*Nquadrature,Ne*Nquadrature))
        allocate(NnNgtemp(Nn,Ne*Nquadrature))
        allocate(NnNgtemp2(Nn,Ne*Nquadrature))

        ! allocate(dnnmatrix(Nn,Nn))
        ! allocate(dnnmatrix2(Nn,Nn))
        ! allocate(G(Nn,Nn))
        allocate(matrix_temp2(Nn,Nn))
        

        do ii=1,Nquadrature1D_range
            do i=1,Nquadrature1D

                NnNgtemp2=0.0d0
    
                do k=1,Nstates

                    psaz_dft(1:nnza)=(E_dft(k)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)


            call mkl_dcsrmm('N',Ne*Nquadrature,Nn,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                ,dble(G),Nn,0.0d0,NgNntemp,Ne*Nquadrature)

            call mkl_dcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                ,transpose(NgNntemp),Nn,0.0d0,NgNgtemp,Ne*Nquadrature)

            call outer_product('ge',psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k),NgNgtemp2)
            NgNgtemp=NgNgtemp*NgNgtemp2!outer_product(psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k))

            call mkl_dcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,4.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                ,NgNgtemp,Ne*Nquadrature,0.0d0,NnNgtemp,Nn)

                    NnNgtemp2=NnNgtemp2+NnNgtemp
                
                enddo

                CALL DGEMM('T','T',Nn,Nn,Ne*Nquadrature,1.0d0,v_kernel,Ne*Nquadrature,NnNgtemp2,Nn,0.0d0,dnnmatrix2,Nn)

                dnnmatrix=IdentityMatrix_real-dnnmatrix2

                call pardiso(pt_v,MAXFCT,MNUM,MTYPE_real,PHASE,Nn,B,IA,JA,idum,Nn&
                                        ,IPARM_V,MSGLVL,v00,matrix_temp2,pardiso_info)
                
                dnnmatrix2=matrix_temp2
    
                CALL DGETRF(Nn,Nn,dnnmatrix,Nn,lapack_IPIV,lapack_INFO)
                CALL DGETRS('N',Nn,Nn,dnnmatrix,Nn,lapack_IPIV,dnnmatrix2,Nn,lapack_INFO)
    
                dnnmatrix=dnnmatrix2-matrix_temp2
    
                gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)=dnnmatrix*j_real

                print *,(ii-1)*Nquadrature1D+i
    
    
            enddo
        enddo

        deallocate(NgNntemp)
        deallocate(NgNgtemp)
        deallocate(NgNgtemp2)
        deallocate(NnNgtemp)
        deallocate(NnNgtemp2)
        deallocate(matrix_temp2)

        end if



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!       Integral in Hypercomplex plane        !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!         do ii=1,Nquadrature1D_range
!             do i=1,Nquadrature1D

!                 NNmatrix_temp002=(0.0d0,0.0d0)
    
!                 do k=1,Nstates
    

!                     do l=1,Nn
!                         do ll=IA(l),IA(l+1)-1
                    
!                             call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+1,(E_dft(k)*B(ll)-H_dft(ll))*j_real)
!                             call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+2,(-omega(i,ii)*B(ll))*j_real)
!                             ! call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+3,(-aimag(Ze)*B(ll))*j_real)
!                             call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+3,(aimag(Ze)*B(ll))*j_real)
!                             call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+4,(0.0d0,0.0d0))
                    
!                             call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+1,(omega(i,ii)*B(ll))*j_real)
!                             call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+2,(E_dft(k)*B(ll)-H_dft(ll))*j_real)
!                             call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+3,(0.0d0,0.0d0))
!                             ! call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+4,(aimag(Ze)*B(ll))*j_real)
!                             call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+4,(-aimag(Ze)*B(ll))*j_real)
                    
!                             ! call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+1,(aimag(Ze)*B(ll))*j_real)
!                             call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+1,(-aimag(Ze)*B(ll))*j_real)
!                             call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+2,(0.0d0,0.0d0))
!                             call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+3,(E_dft(k)*B(ll)-H_dft(ll))*j_real)
!                             call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+4,(-omega(i,ii)*B(ll))*j_real)
                    
!                             call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+1,(0.0d0,0.0d0))
!                             ! call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+2,(-aimag(Ze)*B(ll))*j_real)
!                             call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+2,(aimag(Ze)*B(ll))*j_real)
!                             call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+3,(omega(i,ii)*B(ll))*j_real)
!                             call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+4,(E_dft(k)*B(ll)-H_dft(ll))*j_real)
                    
!                         enddo
!                     enddo

!                     kk=0
!                     do ll=1,4*Nn
!                         do l=1,quaternion(ll)%size
!                             kk=kk+1
!                             valuetemp=quaternion(ll)%geta(l)
!                             quaternion_A(kk)=dble(valuetemp)   
!                         enddo
!                     enddo

!                     PHASE=12
        
!                 call PARDISO(pt_v,MAXFCT,MNUM,MTYPE_real,PHASE,4*Nn,quaternion_A,quaternion_IA,quaternion_JA,idum,Nn,IPARM_V,MSGLVL&
!                 ,IdenMat_quaternion,G_quaternion,pardiso_info)
                    
!                     PHASE=33
!                     iparm_v(6)=0
!                     iparm_v(13)=1
                    
!                 call PARDISO(pt_v,MAXFCT,MNUM,MTYPE_real,PHASE,4*Nn,quaternion_A,quaternion_IA,quaternion_JA,idum,Nn,IPARM_V,MSGLVL&
!                 ,IdenMat_quaternion,G_quaternion,pardiso_info)
                    
                    
!                     do l=1,Nn
!                         do ll=1,Nn
!                             G(l,ll)=cmplx(G_quaternion((l-1)*4+1,ll),G_quaternion((l-1)*4+3,ll),8)
!                         enddo
!                     enddo

! !                     psaz_dft(1:nnza)=(E_dft(k)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
! ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
                    
!                     NNmatrix_temp002=NNmatrix_temp002+4.0d0*G*outer_product(psi(:,k),psi(:,k))
                
!                 enddo
    
!                 CALL ZGEMM('T','N',Nn,Nn,Nn,ZONE,v00*j_real,Nn,NNmatrix_temp002,Nn,ZZERO,NNmatrix_temp003,Nn)
    
!                 call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B_complex,JA,IA_pntrb,IA_pntre,NNmatrix_temp003,Nn&
!                                                                         ,ZZERO,NNmatrix_temp002,Nn)
    
!                 NNmatrix_temp003=IdentityMatrix-NNmatrix_temp002
!                 NNmatrix_temp002=v00*j_real
    
!                 CALL ZGETRF(Nn,Nn,NNmatrix_temp003,Nn,lapack_IPIV,lapack_INFO)
!                 CALL ZGETRS('N',Nn,Nn,NNmatrix_temp003,Nn,lapack_IPIV,NNmatrix_temp002,Nn,lapack_INFO)
    
!                 NNmatrix_temp003=NNmatrix_temp002-v00

!                 ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


!                 call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,B_complex,IA,JA,idum,Nn&
!                                         ,IPARM_V,MSGLVL,NNmatrix_temp003,NNmatrix_temp002,pardiso_info)
    
!                 gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)=NNmatrix_temp002

!                 print *,(ii-1)*Nquadrature1D+i

!                 ! print *,NNmatrix_temp02( 9, 9)
!                 ! print *,NNmatrix_temp02( 9,10)
!                 ! print *,NNmatrix_temp02(10, 9)
!                 ! print *,NNmatrix_temp02(10,10)
    
    
!             enddo
!         enddo



        deallocate(lapack_IPIV)
        deallocate(pt_v)
        deallocate(iparm_v)

        deallocate(dnnmatrix)
        deallocate(dnnmatrix2)
        deallocate(G)


        else if (type=='hf') then


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(matrix_temp(Nn,Nn))
        allocate(dnnmatrix(Nn,Nn))
        allocate(dnnmatrix2(Nn,Nn))
        allocate(G(Nn,Nn))

        !!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! pardiso parameters !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        MAXFCT=1
        MNUM=1
        MTYPE=6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
        MTYPE_real=1
        MSGLVL=0
        allocate(pt_v(64))
        allocate(iparm_v(64))
        pt_V=0
        call pardisoinit(PT_V,MTYPE,IPARM_V)
        PHASE=13

        lapack_lwork=Nn
        allocate(lapack_work(lapack_lwork))

        do ii=1,Nquadrature1D_range
            do i=1,Nquadrature1D

                dnnmatrix=0.0d0
    
                do k=1,Nstates

!                     psaz_dft(1:nnza)=(E_dft(k)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

                    matrix_temp=(E_dft(k)+omega(i,ii)*j_imag)*S_dense-Hhf_dense
                    G=IdentityMatrix

                    ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
                    ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

                    call zsytrf('U',Nn,matrix_temp,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
                    call zsytrs('U',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

                                    
                    call outer_product('ge',psi(:,k),psi(:,k),dnnmatrix2)
                    dnnmatrix=dnnmatrix+4.0d0*dble(G)*dnnmatrix2!outer_product(psi(:,k),psi(:,k))
                
                enddo
    
                CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,dnnmatrix,Nn,0.0d0,dnnmatrix2,Nn)
    
                call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,dnnmatrix2,Nn&
                                                                        ,0.0d0,dnnmatrix,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                dnnmatrix2=IdentityMatrix_real-dnnmatrix
                
                dnnmatrix=v00
    
                CALL DGETRF(Nn,Nn,dnnmatrix2,Nn,lapack_IPIV,lapack_INFO)
                CALL DGETRS('N',Nn,Nn,dnnmatrix2,Nn,lapack_IPIV,dnnmatrix,Nn,lapack_INFO)
    
                dnnmatrix2=dnnmatrix-v00

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


                call pardiso(pt_v,MAXFCT,MNUM,MTYPE_real,PHASE,Nn,B,IA,JA,idum,Nn&
                                        ,IPARM_V,MSGLVL,dnnmatrix2,dnnmatrix,pardiso_info)

                ! NNmatrix_temp02=(NNmatrix_temp02+transpose(NNmatrix_temp02))/2.0d0
    
                gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)=dnnmatrix*j_real

                print *,(ii-1)*Nquadrature1D+i

                ! print *,dnnmatrix(1,1)
                ! print *,dnnmatrix(1,2)
                ! print *,dnnmatrix(2,1)
                ! print *,dnnmatrix(2,2)
    
    
            enddo
        enddo


        deallocate(lapack_IPIV)
        deallocate(lapack_work)

        deallocate(matrix_temp)
        deallocate(dnnmatrix)
        deallocate(dnnmatrix2)
        deallocate(G)

        ! stop

    end if


    end subroutine compute_gw_W_CD





    subroutine compute_gw_Sigma_CD_imagintegral(Nn,Nstates,Nempty,Ze,E_dft,omega,Xi,Nquadrature1D_range,Nquadrature1D,M0,Y,nnza&
        ,Ne,Nquadrature,type,hadamard)
            integer,intent(in) :: Nn,Nstates,Nempty,Nquadrature1D_range,Nquadrature1D,M0,nnza,Ne,Nquadrature
            complex(kind=(kind(1.0d0))),intent(in) :: Ze
            double precision,dimension(1:Nstates),intent(in) :: E_dft
            double precision,dimension(1:Nquadrature1D,1:Nquadrature1D_range),intent(in) :: omega,Xi
            ! complex(kind=(kind(1.0d0))),dimension(1:Nquadrature1D,1:Nquadrature1D_range),intent(in) :: omega,Xi
            complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0),intent(in) :: Y
            character(len=*), intent(in) :: type,hadamard
    
            integer :: i,ii,k,l,ll
            ! complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0) :: Ytemp1,Ytemp2,Ytemp3,Ytemp4
            complex(kind=(kind(1.0d0))) :: valuetemp
            !!!!! pardiso parameters
            integer :: MAXFCT,MNUM,MTYPE,MTYPE_real,MTYPE_complx,MSGLVL,PHASE,idum,pardiso_info,MTYPE_real_unsy
            integer(8),dimension(:),allocatable :: pt_V
            integer,dimension(:),allocatable :: iparm_V

            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: matrix_temp,G,zNgNntemp,zNgNgtemp,zNgNgtemp2
            double precision,dimension(:,:),allocatable :: matrix_temp2
            complex(kind=(kind(1.0d0))),dimension(:),allocatable :: lapack_work
            integer,dimension(:),allocatable :: ipiv


            

            if (type=='dft') then
    
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! pardiso parameters !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            MAXFCT=1
            MNUM=1
            MTYPE=6!1!6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
            ! MTYPE_real=1
            ! MTYPE_complx=3
            ! MTYPE=1
            MSGLVL=0
            allocate(pt_v(64))
            allocate(iparm_v(64))
            pt_V=0
            call pardisoinit(PT_V,MTYPE,IPARM_V)
            PHASE=13

            allocate(matrix_temp(Nn,Nn))
            allocate(G(Nn,Nn))

            
    
    
            ! SigmaY=(0.0d0,0.0d0)
            ! Ytemp4=(0.0d0,0.0d0)


            if (hadamard=='no') then
    
            gw_Sigma_complex=cmplx(0.0d0,0.0d0,8)
            ! gw_Sigma_g=cmplx(0.0d0,0.0d0,8)


    
        do ii=1,Nquadrature1D_range
            do i=1,Nquadrature1D

    
        psaz_dft(1:nnza)=(Ze+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
                ! psaz_dft(1:nnza)=(Ze+omega(i,ii))*usb(1:nnza)-usa_dft(1:nnza)
call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

! print *,G(1,1)
! print *,G(1,2)
! print *,G(2,1)
! print *,G(2,2)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
        ! NNmatrix_temp=-1.0d0/2.0d0/pi*gw_G_complex*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
        ! gw_Sigma_complex=gw_Sigma_complex+NNmatrix_temp*gweight_1D(i)*(omega_max(ii)-omega_min(ii))!/2.0d0 !!! from ZERO to +infinity only so times 2

        matrix_temp=-1.0d0/2.0d0/pi*G*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
        
        gw_Sigma_complex=gw_Sigma_complex+matrix_temp*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
        *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2


    
    






    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
    !             call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,G,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    !     call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)
    
    !             call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,gw_W_complex(:,:,(ii-1)*Nquadrature1D+i),Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    !     call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp2,Ne*Nquadrature)
    
    ! gw_Sigma_g=gw_Sigma_g-1.0d0/2.0d0/pi*zNgNgtemp*zNgNgtemp2*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
    !             *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2
    
    
    
        enddo
        enddo
    
        

        

        else if (hadamard=='yes') then

            allocate(zNgNntemp(Ne*Nquadrature,Nn))
            allocate(zNgNgtemp(Ne*Nquadrature,Ne*Nquadrature))
            allocate(zNgNgtemp2(Ne*Nquadrature,Ne*Nquadrature))


            ! gw_Sigma_complex=cmplx(0.0d0,0.0d0,8)
            gw_Sigma_g=cmplx(0.0d0,0.0d0,8)


    
        do ii=1,Nquadrature1D_range
            do i=1,Nquadrature1D

    
        psaz_dft(1:nnza)=(dble(Ze)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
                ! psaz_dft(1:nnza)=(Ze+omega(i,ii))*usb(1:nnza)-usa_dft(1:nnza)
call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)


        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
        ! ! NNmatrix_temp=-1.0d0/2.0d0/pi*gw_G_complex*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
        ! ! gw_Sigma_complex=gw_Sigma_complex+NNmatrix_temp*gweight_1D(i)*(omega_max(ii)-omega_min(ii))!/2.0d0 !!! from ZERO to +infinity only so times 2

        ! matrix_temp=-1.0d0/2.0d0/pi*G*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
        
        ! gw_Sigma_complex=gw_Sigma_complex+matrix_temp*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
        ! *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2


    
    






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
                call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,G,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
        call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)
    
                call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,gw_W_complex(:,:,(ii-1)*Nquadrature1D+i),Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
        call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp2,Ne*Nquadrature)
    
    gw_Sigma_g=gw_Sigma_g-1.0d0/2.0d0/pi*zNgNgtemp*zNgNgtemp2*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
                *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2
    
    
    
        enddo
        enddo

        deallocate(zNgNntemp)
        deallocate(zNgNgtemp)
        deallocate(zNgNgtemp2)


    end if

    deallocate(pt_v)
    deallocate(iparm_v)

    deallocate(matrix_temp)
    deallocate(G)


    else if (type=='hf') then

            gw_Sigma_complex=cmplx(0.0d0,0.0d0,8)

            allocate(ipiv(1:Nn))


            allocate(matrix_temp(Nn,Nn))
            allocate(G(Nn,Nn))

            lapack_lwork=3*Nn
            allocate(lapack_work(lapack_lwork))

            do ii=1,Nquadrature1D_range
                do i=1,Nquadrature1D

                    matrix_temp=(dble(Ze)+omega(i,ii)*j_imag)*S_dense-Hhf_dense
                    G=IdentityMatrix

                    call zsytrf('U',Nn,matrix_temp,Nn,ipiv,lapack_work,lapack_lwork,lapack_INFO)
                    call zsytrs('U',Nn,Nn,matrix_temp,Nn,ipiv,G,Nn,lapack_INFO)

                        ! print *,G(1,1)
                        ! print *,G(1,2)
                        ! print *,G(2,1)
                        ! print *,G(2,2)

                    ! print *,S_dense(749,749),Hhf_dense(749,749)


    ! do l=1,Nn
    !     do ll=1,Nn

    !         matrix_temp2((l-1)*4+1,(ll-1)*4+1) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0
    !         matrix_temp2((l-1)*4+1,(ll-1)*4+2) = -omega(i,ii)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+1,(ll-1)*4+3) = -aimag(Ze)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+1,(ll-1)*4+4) = 0.0d0

    !         matrix_temp2((l-1)*4+2,(ll-1)*4+1) = omega(i,ii)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+2,(ll-1)*4+2) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0
    !         matrix_temp2((l-1)*4+2,(ll-1)*4+3) = 0.0d0
    !         matrix_temp2((l-1)*4+2,(ll-1)*4+4) = aimag(Ze)*S_dense(l,ll)

    !         matrix_temp2((l-1)*4+3,(ll-1)*4+1) = aimag(Ze)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+3,(ll-1)*4+2) = 0.0d0
    !         matrix_temp2((l-1)*4+3,(ll-1)*4+3) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0
    !         matrix_temp2((l-1)*4+3,(ll-1)*4+4) = -omega(i,ii)*S_dense(l,ll)

    !         matrix_temp2((l-1)*4+4,(ll-1)*4+1) = 0.0d0
    !         matrix_temp2((l-1)*4+4,(ll-1)*4+2) = -aimag(Ze)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+4,(ll-1)*4+3) = omega(i,ii)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+4,(ll-1)*4+4) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0

    
    !     enddo
    ! enddo

    ! G_quaternion = IdenMat_quaternion

    ! CALL DGETRF(4*Nn,4*Nn,matrix_temp2,4*Nn,ipiv,lapack_INFO)
    ! CALL DGETRS('N',4*Nn,Nn,matrix_temp2,4*Nn,ipiv,G_quaternion,4*Nn,lapack_INFO)

    ! do l=1,Nn
    !     do ll=1,Nn
    !         G(l,ll)=cmplx(G_quaternion((l-1)*4+1,ll),G_quaternion((l-1)*4+3,ll),8)
    !     enddo
    ! enddo

    
    
    
    
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
            ! NNmatrix_temp=-1.0d0/2.0d0/pi*gw_G_complex*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
            ! gw_Sigma_complex=gw_Sigma_complex+NNmatrix_temp*gweight_1D(i)*(omega_max(ii)-omega_min(ii))!/2.0d0 !!! from ZERO to +infinity only so times 2
    
            matrix_temp=-1.0d0/2.0d0/pi*G*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)

            ! print *,matrix_temp(1,1)
            ! print *,matrix_temp(1,2)
            ! print *,matrix_temp(2,1)
            ! print *,matrix_temp(2,2)

            ! print *,dble(G(1,1)),dble(gw_W_complex(1,1,(ii-1)*Nquadrature1D+i))
            ! print *,dble(G(1,2)),dble(gw_W_complex(1,2,(ii-1)*Nquadrature1D+i))
            ! print *,dble(G(2,1)),dble(gw_W_complex(2,1,(ii-1)*Nquadrature1D+i))
            ! print *,dble(G(2,2)),dble(gw_W_complex(2,2,(ii-1)*Nquadrature1D+i))

        ! NNmatrix_temp=-1.0d0/2.0d0/pi*gw_G_complex*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
        ! gw_Sigma_complex=gw_Sigma_complex+NNmatrix_temp*gweight_1D(i)*(omega_max(ii)-omega_min(ii))!/2.0d0 !!! from ZERO to +infinity only so times 2
            
            gw_Sigma_complex=gw_Sigma_complex+matrix_temp*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
            *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2
        
        
            enddo
            enddo

            ! deallocate(lapack_IPIV)
            

            deallocate(ipiv)
            deallocate(matrix_temp)
            deallocate(lapack_work)
            deallocate(G)

        end if

            ! stop
    
    
        end subroutine compute_gw_Sigma_CD_imagintegral




    subroutine compute_gw_Sigma_CD_imagintegral_2(Nn,Nstates,Nempty,Ze,E_dft,omega,Xi,Nquadrature1D_range,Nquadrature1D,M0,Y,nnza&
            ,Ne,Nquadrature,type,hadamard)
            integer,intent(in) :: Nn,Nstates,Nempty,Nquadrature1D_range,Nquadrature1D,M0,nnza,Ne,Nquadrature
            complex(kind=(kind(1.0d0))),intent(in) :: Ze
            double precision,dimension(1:Nstates),intent(in) :: E_dft
            double precision,dimension(1:Nquadrature1D,1:Nquadrature1D_range),intent(in) :: omega,Xi
            complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0),intent(in) :: Y
            character(len=*), intent(in) :: type,hadamard
    
            integer :: i,ii,k,l,ll
            ! complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0) :: Ytemp1,Ytemp2,Ytemp3,Ytemp4
            complex(kind=(kind(1.0d0))) :: valuetemp
            !!!!! pardiso parameters
            integer :: MAXFCT,MNUM,MTYPE,MTYPE_real,MTYPE_complx,MSGLVL,PHASE,idum,pardiso_info,MTYPE_real_unsy
            integer(8),dimension(:),allocatable :: pt_V
            integer,dimension(:),allocatable :: iparm_V

            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: matrix_temp,G,zNgNntemp,zNgNgtemp,zNgNgtemp2
            double precision,dimension(:,:),allocatable :: matrix_temp2
            complex(kind=(kind(1.0d0))),dimension(:),allocatable :: lapack_work
            integer,dimension(:),allocatable :: ipiv


            


            if (type=='dft') then

                if (hadamard=='no') then

                gw_Sigma_complex=cmplx(0.0d0,0.0d0,8)
        
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! pardiso parameters !!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                MAXFCT=1
                MNUM=1
                MTYPE=1!6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
                ! MTYPE_real=1
                ! MTYPE_complx=3
                ! MTYPE=1
                MSGLVL=0
                allocate(pt_v(64))
                allocate(iparm_v(64))
                pt_V=0
                call pardisoinit(PT_V,MTYPE,IPARM_V)
                PHASE=13

                allocate(matrix_temp(Nn,Nn))
                allocate(G(Nn,Nn))
        
        
                ! SigmaY=(0.0d0,0.0d0)
                ! Ytemp4=(0.0d0,0.0d0)
        
                
                ! gw_Sigma_g=cmplx(0.0d0,0.0d0,8)
        
            do ii=1,Nquadrature1D_range
                do i=1,Nquadrature1D
        
        
    !         psaz_dft(1:nnza)=(dble(Ze)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    ! NNmatrix_temp005 = G
    
    ! psaz_dft(1:nnza)=(conjg(Ze)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    ! NNmatrix_temp005 = NNmatrix_temp005 + G

    ! NNmatrix_temp005 = NNmatrix_temp005/2.0d0
        
        do l=1,Nn
            do ll=IA(l),IA(l+1)-1
        
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+1,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+2,(-omega(i,ii)*B(ll))*j_real)
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+3,(-aimag(Ze)*B(ll))*j_real)
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+4,(0.0d0,0.0d0))
        
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+1,(omega(i,ii)*B(ll))*j_real)
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+2,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+3,(0.0d0,0.0d0))
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+4,(aimag(Ze)*B(ll))*j_real)
        
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+1,(aimag(Ze)*B(ll))*j_real)
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+2,(0.0d0,0.0d0))
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+3,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+4,(-omega(i,ii)*B(ll))*j_real)
        
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+1,(0.0d0,0.0d0))
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+2,(-aimag(Ze)*B(ll))*j_real)
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+3,(omega(i,ii)*B(ll))*j_real)
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+4,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
        
            enddo
        enddo
    
        k=0
        do ll=1,4*Nn
            do l=1,quaternion(ll)%size
                k=k+1
                valuetemp=quaternion(ll)%geta(l)
                quaternion_A(k)=dble(valuetemp)   
            enddo
        enddo
        
        ! do l=1,Nn
        !     do ll=IA(l),IA(l+1)-1
        !         if (JA(ll)==l) IdenMat_quaternion((l-1)*4+3,l)=aimag(Ze)*B(ll)
        !     enddo
        ! enddo
        
        ! PHASE=12
        
        ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,4*Nn,quaternion_A,quaternion_IA,quaternion_JA,idum,Nn,IPARM_V,MSGLVL&
        ! ,IdenMat_quaternion,G_quaternion,pardiso_info)
        
        PHASE=13!33
        iparm_v(6)=0
        iparm_v(13)=1
        
        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,4*Nn,quaternion_A,quaternion_IA,quaternion_JA,idum,Nn,IPARM_V,MSGLVL&
        ,IdenMat_quaternion,G_quaternion,pardiso_info)
        
        
        do l=1,Nn
            do ll=1,Nn
                G(l,ll)=cmplx(G_quaternion((l-1)*4+1,ll),G_quaternion((l-1)*4+3,ll),8)
            enddo
        enddo
    
        ! do l=1,Nn
        !     do ll=1,Nn
        !         gw_G_complex(l,ll,(ii-1)*Nquadrature1D+i)=cmplx(G_quaternion((l-1)*4+1,ll),G_quaternion((l-1)*4+3,ll),8)
        !     enddo
        ! enddo
    
        ! gw_G_complex(l,ll,(ii-1)*Nquadrature1D+i)=gw_G_complex(l,ll,(ii-1)*Nquadrature1D+i)*(-1.0d0/2.0d0/pi)&
        ! *gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))*exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2
    
        ! print *,G(1,1)
        ! print *,G(1,2)
        ! print *,G(2,1)
        ! print *,G(2,2)
    
        ! stop


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
            ! NNmatrix_temp=-1.0d0/2.0d0/pi*gw_G_complex*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
            ! gw_Sigma_complex=gw_Sigma_complex+NNmatrix_temp*gweight_1D(i)*(omega_max(ii)-omega_min(ii))!/2.0d0 !!! from ZERO to +infinity only so times 2
    
        matrix_temp=-1.0d0/2.0d0/pi*G*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
            
        gw_Sigma_complex=gw_Sigma_complex+matrix_temp*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
        *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2
    
    






    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
    !             call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,G,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    !     call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)
    
    !             call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,gw_W_complex(:,:,(ii-1)*Nquadrature1D+i),Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    !     call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                                 ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp2,Ne*Nquadrature)
    
    ! gw_Sigma_g=gw_Sigma_g-1.0d0/2.0d0/pi*zNgNgtemp*zNgNgtemp2*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
    !             *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2

        enddo
    enddo

    deallocate(pt_v)
    deallocate(iparm_v)
    deallocate(matrix_temp)
    deallocate(G)


    else if (hadamard=='yes') then


        gw_Sigma_g=cmplx(0.0d0,0.0d0,8)

        allocate(zNgNntemp(Ne*Nquadrature,Nn))
        allocate(zNgNgtemp(Ne*Nquadrature,Ne*Nquadrature))
        allocate(zNgNgtemp2(Ne*Nquadrature,Ne*Nquadrature))
        
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! pardiso parameters !!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!
                MAXFCT=1
                MNUM=1
                MTYPE=1!6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
                ! MTYPE_real=1
                ! MTYPE_complx=3
                ! MTYPE=1
                MSGLVL=0
                allocate(pt_v(64))
                allocate(iparm_v(64))
                pt_V=0
                call pardisoinit(PT_V,MTYPE,IPARM_V)
                PHASE=13

                allocate(matrix_temp(Nn,Nn))
                allocate(G(Nn,Nn))
        
        
                ! SigmaY=(0.0d0,0.0d0)
                ! Ytemp4=(0.0d0,0.0d0)
        
                
                ! gw_Sigma_g=cmplx(0.0d0,0.0d0,8)
        
            do ii=1,Nquadrature1D_range
                do i=1,Nquadrature1D
        
        
    !         psaz_dft(1:nnza)=(dble(Ze)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    ! NNmatrix_temp005 = G
    
    ! psaz_dft(1:nnza)=(conjg(Ze)+omega(i,ii)*j_imag)*usb(1:nnza)-usa_dft(1:nnza)
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    ! NNmatrix_temp005 = NNmatrix_temp005 + G

    ! NNmatrix_temp005 = NNmatrix_temp005/2.0d0
        
        do l=1,Nn
            do ll=IA(l),IA(l+1)-1
        
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+1,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+2,(-omega(i,ii)*B(ll))*j_real)
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+3,(-aimag(Ze)*B(ll))*j_real)
                call quaternion((l-1)*4+1)%overwriteA((JA(ll)-1)*4+4,(0.0d0,0.0d0))
        
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+1,(omega(i,ii)*B(ll))*j_real)
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+2,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+3,(0.0d0,0.0d0))
                call quaternion((l-1)*4+2)%overwriteA((JA(ll)-1)*4+4,(aimag(Ze)*B(ll))*j_real)
        
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+1,(aimag(Ze)*B(ll))*j_real)
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+2,(0.0d0,0.0d0))
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+3,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
                call quaternion((l-1)*4+3)%overwriteA((JA(ll)-1)*4+4,(-omega(i,ii)*B(ll))*j_real)
        
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+1,(0.0d0,0.0d0))
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+2,(-aimag(Ze)*B(ll))*j_real)
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+3,(omega(i,ii)*B(ll))*j_real)
                call quaternion((l-1)*4+4)%overwriteA((JA(ll)-1)*4+4,(dble(Ze)*B(ll)-H_dft(ll))*j_real)
        
            enddo
        enddo
    
        k=0
        do ll=1,4*Nn
            do l=1,quaternion(ll)%size
                k=k+1
                valuetemp=quaternion(ll)%geta(l)
                quaternion_A(k)=dble(valuetemp)   
            enddo
        enddo
        
        ! do l=1,Nn
        !     do ll=IA(l),IA(l+1)-1
        !         if (JA(ll)==l) IdenMat_quaternion((l-1)*4+3,l)=aimag(Ze)*B(ll)
        !     enddo
        ! enddo
        
        ! PHASE=12
        
        ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,4*Nn,quaternion_A,quaternion_IA,quaternion_JA,idum,Nn,IPARM_V,MSGLVL&
        ! ,IdenMat_quaternion,G_quaternion,pardiso_info)
        
        PHASE=13!33
        iparm_v(6)=0
        iparm_v(13)=1
        
        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,4*Nn,quaternion_A,quaternion_IA,quaternion_JA,idum,Nn,IPARM_V,MSGLVL&
        ,IdenMat_quaternion,G_quaternion,pardiso_info)
        
        
        do l=1,Nn
            do ll=1,Nn
                G(l,ll)=cmplx(G_quaternion((l-1)*4+1,ll),G_quaternion((l-1)*4+3,ll),8)
            enddo
        enddo
    
        ! do l=1,Nn
        !     do ll=1,Nn
        !         gw_G_complex(l,ll,(ii-1)*Nquadrature1D+i)=cmplx(G_quaternion((l-1)*4+1,ll),G_quaternion((l-1)*4+3,ll),8)
        !     enddo
        ! enddo
    
        ! gw_G_complex(l,ll,(ii-1)*Nquadrature1D+i)=gw_G_complex(l,ll,(ii-1)*Nquadrature1D+i)*(-1.0d0/2.0d0/pi)&
        ! *gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))*exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2
    
        ! print *,G(1,1)
        ! print *,G(1,2)
        ! print *,G(2,1)
        ! print *,G(2,2)
    
        ! stop


        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !     !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
        !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !     ! NNmatrix_temp=-1.0d0/2.0d0/pi*gw_G_complex*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
        !     ! gw_Sigma_complex=gw_Sigma_complex+NNmatrix_temp*gweight_1D(i)*(omega_max(ii)-omega_min(ii))!/2.0d0 !!! from ZERO to +infinity only so times 2
    
        ! matrix_temp=-1.0d0/2.0d0/pi*G*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
            
        ! gw_Sigma_complex=gw_Sigma_complex+matrix_temp*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
        ! *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2
    
    






        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
        call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        ,G,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)

call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)

call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        ,gw_W_complex(:,:,(ii-1)*Nquadrature1D+i),Nn,ZZERO,zNgNntemp,Ne*Nquadrature)

call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp2,Ne*Nquadrature)

gw_Sigma_g=gw_Sigma_g-1.0d0/2.0d0/pi*zNgNgtemp*zNgNgtemp2*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
*exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2

        enddo
    enddo

    deallocate(zNgNntemp)
    deallocate(zNgNgtemp)
    deallocate(zNgNgtemp2)

    deallocate(pt_v)
    deallocate(iparm_v)
    deallocate(matrix_temp)
    deallocate(G)


    end if !!! hadamard


    else if (type=='hf') then

        gw_Sigma_complex=cmplx(0.0d0,0.0d0,8)


        allocate(ipiv(1:4*Nn))

        allocate(matrix_temp2(4*Nn,4*Nn))
        allocate(matrix_temp(Nn,Nn))
        allocate(G(Nn,Nn))


    

    do ii=1,Nquadrature1D_range
        do i=1,Nquadrature1D

        

    ! do l=1,Nn
    !     do ll=1,Nn

    !         matrix_temp2((l-1)*4+1,(ll-1)*4+1) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0
    !         matrix_temp2((l-1)*4+1,(ll-1)*4+2) = -omega(i,ii)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+1,(ll-1)*4+3) = -aimag(Ze)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+1,(ll-1)*4+4) = 0.0d0

    !         matrix_temp2((l-1)*4+2,(ll-1)*4+1) = omega(i,ii)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+2,(ll-1)*4+2) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0
    !         matrix_temp2((l-1)*4+2,(ll-1)*4+3) = 0.0d0
    !         matrix_temp2((l-1)*4+2,(ll-1)*4+4) = aimag(Ze)*S_dense(l,ll)

    !         matrix_temp2((l-1)*4+3,(ll-1)*4+1) = aimag(Ze)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+3,(ll-1)*4+2) = 0.0d0
    !         matrix_temp2((l-1)*4+3,(ll-1)*4+3) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0
    !         matrix_temp2((l-1)*4+3,(ll-1)*4+4) = -omega(i,ii)*S_dense(l,ll)

    !         matrix_temp2((l-1)*4+4,(ll-1)*4+1) = 0.0d0
    !         matrix_temp2((l-1)*4+4,(ll-1)*4+2) = -aimag(Ze)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+4,(ll-1)*4+3) = omega(i,ii)*S_dense(l,ll)
    !         matrix_temp2((l-1)*4+4,(ll-1)*4+4) = dble(Ze)*S_dense(l,ll)-(Hhf_dense(l,ll)+Hhf_dense(ll,l))/2.0d0

    
    !     enddo
    ! enddo

            do l=1,Nn
                do ll=1,Nn

                    if (l>ll) then

                        matrix_temp2((l-1)*4+1,(ll-1)*4+1) = dble(Ze)*S_dense(ll,l)-Hhf_dense(ll,l)
                        matrix_temp2((l-1)*4+1,(ll-1)*4+2) = -omega(i,ii)*S_dense(ll,l)
                        matrix_temp2((l-1)*4+1,(ll-1)*4+3) = -aimag(Ze)*S_dense(ll,l)
                        matrix_temp2((l-1)*4+1,(ll-1)*4+4) = 0.0d0
            
                        matrix_temp2((l-1)*4+2,(ll-1)*4+1) = omega(i,ii)*S_dense(ll,l)
                        matrix_temp2((l-1)*4+2,(ll-1)*4+2) = dble(Ze)*S_dense(ll,l)-Hhf_dense(ll,l)
                        matrix_temp2((l-1)*4+2,(ll-1)*4+3) = 0.0d0
                        matrix_temp2((l-1)*4+2,(ll-1)*4+4) = aimag(Ze)*S_dense(ll,l)
            
                        matrix_temp2((l-1)*4+3,(ll-1)*4+1) = aimag(Ze)*S_dense(ll,l)
                        matrix_temp2((l-1)*4+3,(ll-1)*4+2) = 0.0d0
                        matrix_temp2((l-1)*4+3,(ll-1)*4+3) = dble(Ze)*S_dense(ll,l)-Hhf_dense(ll,l)
                        matrix_temp2((l-1)*4+3,(ll-1)*4+4) = -omega(i,ii)*S_dense(ll,l)
            
                        matrix_temp2((l-1)*4+4,(ll-1)*4+1) = 0.0d0
                        matrix_temp2((l-1)*4+4,(ll-1)*4+2) = -aimag(Ze)*S_dense(ll,l)
                        matrix_temp2((l-1)*4+4,(ll-1)*4+3) = omega(i,ii)*S_dense(ll,l)
                        matrix_temp2((l-1)*4+4,(ll-1)*4+4) = dble(Ze)*S_dense(ll,l)-Hhf_dense(ll,l)

                    else

                        matrix_temp2((l-1)*4+1,(ll-1)*4+1) = dble(Ze)*S_dense(l,ll)-Hhf_dense(l,ll)
                        matrix_temp2((l-1)*4+1,(ll-1)*4+2) = -omega(i,ii)*S_dense(l,ll)
                        matrix_temp2((l-1)*4+1,(ll-1)*4+3) = -aimag(Ze)*S_dense(l,ll)
                        matrix_temp2((l-1)*4+1,(ll-1)*4+4) = 0.0d0
            
                        matrix_temp2((l-1)*4+2,(ll-1)*4+1) = omega(i,ii)*S_dense(l,ll)
                        matrix_temp2((l-1)*4+2,(ll-1)*4+2) = dble(Ze)*S_dense(l,ll)-Hhf_dense(l,ll)
                        matrix_temp2((l-1)*4+2,(ll-1)*4+3) = 0.0d0
                        matrix_temp2((l-1)*4+2,(ll-1)*4+4) = aimag(Ze)*S_dense(l,ll)
            
                        matrix_temp2((l-1)*4+3,(ll-1)*4+1) = aimag(Ze)*S_dense(l,ll)
                        matrix_temp2((l-1)*4+3,(ll-1)*4+2) = 0.0d0
                        matrix_temp2((l-1)*4+3,(ll-1)*4+3) = dble(Ze)*S_dense(l,ll)-Hhf_dense(l,ll)
                        matrix_temp2((l-1)*4+3,(ll-1)*4+4) = -omega(i,ii)*S_dense(l,ll)
            
                        matrix_temp2((l-1)*4+4,(ll-1)*4+1) = 0.0d0
                        matrix_temp2((l-1)*4+4,(ll-1)*4+2) = -aimag(Ze)*S_dense(l,ll)
                        matrix_temp2((l-1)*4+4,(ll-1)*4+3) = omega(i,ii)*S_dense(l,ll)
                        matrix_temp2((l-1)*4+4,(ll-1)*4+4) = dble(Ze)*S_dense(l,ll)-Hhf_dense(l,ll)

                    end if
        
                enddo
            enddo

            

    G_quaternion = IdenMat_quaternion

    CALL DGETRF(4*Nn,4*Nn,matrix_temp2,4*Nn,ipiv,lapack_INFO)
    CALL DGETRS('N',4*Nn,Nn,matrix_temp2,4*Nn,ipiv,G_quaternion,4*Nn,lapack_INFO)

    do l=1,Nn
        do ll=1,Nn
            G(l,ll)=cmplx(G_quaternion((l-1)*4+1,ll),G_quaternion((l-1)*4+3,ll),8)
        enddo
    enddo

    ! print *,G(1,1)
    ! print *,G(1,2)
    ! print *,G(2,1)
    ! print *,G(2,2)

    ! stop
    
    
    
    
    
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
            ! NNmatrix_temp=-1.0d0/2.0d0/pi*gw_G_complex*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
            ! gw_Sigma_complex=gw_Sigma_complex+NNmatrix_temp*gweight_1D(i)*(omega_max(ii)-omega_min(ii))!/2.0d0 !!! from ZERO to +infinity only so times 2
    
            matrix_temp=-1.0d0/2.0d0/pi*G*gw_W_complex(:,:,(ii-1)*Nquadrature1D+i)
            
            gw_Sigma_complex=gw_Sigma_complex+matrix_temp*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
            *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2
        
        
    
    
    
    
    
    
        !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !     !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
        !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
        
        !             call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        !                                                 ,G,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
        
        !     call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        !                                                 ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)
        
        !             call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        !                                                 ,gw_W_complex(:,:,(ii-1)*Nquadrature1D+i),Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
        
        !     call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        !                                                 ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp2,Ne*Nquadrature)
        
        ! gw_Sigma_g=gw_Sigma_g-1.0d0/2.0d0/pi*zNgNgtemp*zNgNgtemp2*gweight_1D(i)*(Xi_max(ii)-Xi_min(ii))&!/2.0d0 !!! from ZERO to +infinity only so times 2
        !             *exp(alpha_gw*Xi_range(i,ii)/(1.0d0-Xi_range(i,ii)))*alpha_gw/(1.0d0-Xi_range(i,ii))**2!*1.0d0/(1.0d0-Xi_range(i,ii))**2
        
        
    
        
        
            enddo
            enddo

        deallocate(ipiv)
        deallocate(matrix_temp2)
        deallocate(matrix_temp)
        deallocate(G)

        end if
        
            
        
        
        end subroutine compute_gw_Sigma_CD_imagintegral_2






        subroutine compute_gw_Sigma_CD_poles_mtx(Nn,Nstates,Nempty,Ze,E_dft,M00,M0,nnza,Ne,Nquadrature,type,hadamard)
            integer,intent(in) :: Nn,Nstates,Nempty,M00,M0,nnza,Ne,Nquadrature
            complex(kind=(kind(1.0d0))),intent(in) :: Ze
            double precision,dimension(1:Nstates+10),intent(in) :: E_dft
            ! complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0),intent(in) :: Y
            character(len=*), intent(in) :: type,hadamard
    
            integer :: i,ii,k,ll
            ! complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0) :: Ytemp1,Ytemp2,Ytemp3,Ytemp4,Ytemp5,Ytemp6,Ytemp7
            double precision :: norm
            !!!!! pardiso parameters
            integer :: MAXFCT,MNUM,MTYPE,MTYPE_real,MTYPE_complx,MSGLVL,PHASE,idum,pardiso_info
            integer(8),dimension(:),allocatable :: pt_V
            integer,dimension(:),allocatable :: iparm_V

            integer :: lapack_INFO

            double precision,dimension(:,:),allocatable :: matrix_temp,matrix_temp2,matrix_temp3,G_real,matrix_temp0,matrix_temp1
            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: G
            complex(kind=(kind(1.0d0))),dimension(:),allocatable :: lapack_work

            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: NNmatrix_temp005,NNmatrix_temp006,NNmatrix_temp
            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: zNnNgtemp2,zNgNntemp,zNgNgtemp,zNnNgtemp
            double precision,dimension(:,:),allocatable :: NgNgtemp


            

            
            ! gw_Sigma_g_poles=ZZERO


            allocate(lapack_IPIV(1:Nn))

            if (type=='dft') then
    
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! pardiso parameters !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            MAXFCT=1
            MNUM=1
            MTYPE=6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
            MTYPE_real=1
            MTYPE_complx=3
            MSGLVL=0
            allocate(pt_v(64))
            allocate(iparm_v(64))
            pt_V=0
            call pardisoinit(PT_V,MTYPE,IPARM_V)
            PHASE=13

            allocate(matrix_temp0(Nn,Nn))
            allocate(matrix_temp1(Nn,Nn))
            allocate(matrix_temp2(Nn,Nn))
            allocate(G(Nn,Nn))

    

            if (hadamard=='no') then



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
            
            gw_Sigma_complex_poles=ZZERO
    
            do i=1,Nstates
    
                if (E_dft(i)-dble(Ze)>0.0d0) then
        
    
                    matrix_temp0=0.0d0
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp1=0.0d0

        
        
                        psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        matrix_temp1=matrix_temp1+dble(G)
        
                        psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        matrix_temp1=matrix_temp1+dble(G)


                        call outer_product('ge',psi(:,k),psi(:,k),matrix_temp2)
                        matrix_temp0=matrix_temp0+2.0d0*matrix_temp1*matrix_temp2!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,matrix_temp0,Nn,0.0d0,matrix_temp1,Nn)
    
                call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,matrix_temp1,Nn&
                                                                            ,0.0d0,matrix_temp0,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp1=IdentityMatrix_real-matrix_temp0
                
                matrix_temp0=v00
    
                CALL DGETRF(Nn,Nn,matrix_temp1,Nn,lapack_IPIV,lapack_INFO)
                ! CALL DGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL DGETRS('N',Nn,Nn,matrix_temp1,Nn,lapack_IPIV,matrix_temp0,Nn,lapack_INFO)
    
                matrix_temp0=matrix_temp0-v00

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp0*j_real&
                ,G,pardiso_info)  !!!!! G here is meaningless, only for a temp complex matrix

                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp2)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles-G*matrix_temp2!outer_product(psi(:,i),psi(:,i))



                end if
        
            enddo


            do i=Nstates+1,M00
    
                if (E_dft(i)-dble(Ze)<0.0d0) then


                    matrix_temp0=0.0d0
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp1=0.0d0

        
        
                        psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        matrix_temp1=matrix_temp1+dble(G)
        
                        psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        matrix_temp1=matrix_temp1+dble(G)


                        call outer_product('ge',psi(:,k),psi(:,k),matrix_temp2)
                        matrix_temp0=matrix_temp0+2.0d0*matrix_temp1*matrix_temp2!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,matrix_temp0,Nn,0.0d0,matrix_temp1,Nn)
    
                call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,matrix_temp1,Nn&
                                                                            ,0.0d0,matrix_temp0,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp1=IdentityMatrix_real-matrix_temp0
                
                matrix_temp0=v00
    
                CALL DGETRF(Nn,Nn,matrix_temp1,Nn,lapack_IPIV,lapack_INFO)
                ! CALL DGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL DGETRS('N',Nn,Nn,matrix_temp1,Nn,lapack_IPIV,matrix_temp0,Nn,lapack_INFO)
    
                matrix_temp0=matrix_temp0-v00

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp0*j_real&
                ,G,pardiso_info)  !!!!! G here is meaningless, only for a temp matrix

                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp2)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles+G*matrix_temp2!outer_product(psi(:,i),psi(:,i))


                end if
        
            enddo


            else if (hadamard=='yes') then


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            allocate(NNmatrix_temp005(Nn,Nn))
            allocate(zNnNgtemp2(Nn,Ne*Nquadrature))
            allocate(NNmatrix_temp(Nn,Nn))
            allocate(zNgNntemp(Ne*Nquadrature,Nn))
            allocate(zNgNgtemp(Ne*Nquadrature,Ne*Nquadrature))
            allocate(zNnNgtemp(Nn,Ne*Nquadrature))
            allocate(NNmatrix_temp006(Nn,Nn))
            allocate(NgNgtemp(Ne*Nquadrature,Ne*Nquadrature))


            ! gw_Sigma_complex_poles=ZZERO
            gw_Sigma_g_poles=ZZERO
    
            do i=1,Nstates
    
                if (E_dft(i)-dble(Ze)>0.0d0) then
        
    
                    NNmatrix_temp005=cmplx(0.0d0,0.0d0,8)
                    zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        NNmatrix_temp=cmplx(0.0d0,0.0d0,8)

                        psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        NNmatrix_temp=NNmatrix_temp+dble(G)
        
                        psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        NNmatrix_temp=NNmatrix_temp+dble(G)
    
    
        call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                        ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                        ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)

    call outer_product('ge',psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k),NgNgtemp)
    zNgNgtemp=zNgNgtemp*NgNgtemp!outer_product(psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k))
    
    call mkl_zcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                        ,zNgNgtemp,Ne*Nquadrature,ZZERO,zNnNgtemp,Nn)
    
    !     call mkl_zcsrmm('T',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                     ,transpose(zNnNgtemp),Ne*Nquadrature,ZZERO,chi_matrix,Nn)
    
    !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,chi_matrix&
    !     ,NNmatrix_temp,pardiso_info)
    
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,transpose(NNmatrix_temp)&
    !     ,chi_matrix,pardiso_info)
    
        ! NNmatrix_temp005=NNmatrix_temp005+2.0d0*chi_matrix!+4.0d0*G*outer_product(psi(:,k),psi(:,k))

                                        zNnNgtemp2=zNnNgtemp2+2.0d0*zNnNgtemp
        
        
                    enddo
        
        
                    ! CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,NNmatrix_temp02,Nn,0.0d0,NNmatrix_temp03,Nn)
                    ! CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,NNmatrix_temp005,Nn,(0.0d0,0.0d0),NNmatrix_temp,Nn)

            CALL ZGEMM('T','T',Nn,Nn,Ne*Nquadrature,ZONE,v_kernel*j_real,Ne*Nquadrature,zNnNgtemp2,Nn,ZZERO,NNmatrix_temp005,Nn)
    
                    ! call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,NNmatrix_temp03,Nn&
                    !                                                         ,0.0d0,NNmatrix_temp02,Nn)
                    ! call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,NNmatrix_temp,Nn&
                    !                                                         ,(0.0d0,0.0d0),NNmatrix_temp005,Nn)
        
                    NNmatrix_temp=IdentityMatrix-NNmatrix_temp005

        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,v00*j_real&
            ,NNmatrix_temp006,pardiso_info)
        
                    ! Ytemp2=Ytemp3
                    NNmatrix_temp005=NNmatrix_temp006
        
                    CALL ZGETRF(Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,lapack_INFO)
                    ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                    CALL ZGETRS('N',Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,NNmatrix_temp005,Nn,lapack_INFO)
    
                    ! print *,Ytemp2(1:3,1)
        
                    ! Ytemp2=Ytemp2-Ytemp3
                    NNmatrix_temp=NNmatrix_temp005-NNmatrix_temp006
        
        ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp2,Ytemp4,pardiso_info)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,NNmatrix_temp005&
            ! ,NNmatrix_temp,pardiso_info)


            call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
        call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)

            call outer_product('ge',psi_point_g(:,i),psi_point_g(:,i),NgNgtemp)
            zNgNgtemp=zNgNgtemp*NgNgtemp!outer_product(psi_point_g(:,i),psi_point_g(:,i))


                gw_Sigma_g_poles=gw_Sigma_g_poles-zNgNgtemp
    
    
                ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
                !                                     ,NNmatrix_temp,Nn,ZZERO,NNmatrix_temp005,Ne*Nquadrature)
    
                ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
                !                                     ,transpose(NNmatrix_temp005),Nn,ZZERO,NNmatrix_temp,Ne*Nquadrature)
    
    
                ! gw_Sigma_complex_poles=gw_Sigma_complex_poles-NNmatrix_temp
    
            end if
    
            enddo
    
            ! call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp7,Nn,(0.0d0,0.0d0),Ytemp5,Nn)
    
            ! Y_primed=Ytemp5
    
            ! go to 2222
    
            ! Ytemp7=(0.0d0,0.0d0)
        do i=Nstates+1,M00
    
            if (E_dft(i)-dble(Ze)<0.0d0) then
    
    
                NNmatrix_temp005=cmplx(0.0d0,0.0d0,8)
                    zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        NNmatrix_temp=cmplx(0.0d0,0.0d0,8)


                        psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        NNmatrix_temp=NNmatrix_temp+dble(G)
        
                        psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        NNmatrix_temp=NNmatrix_temp+dble(G)
    
    
        call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                        ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                        ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)
    
    call outer_product('ge',psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k),NgNgtemp)
    zNgNgtemp=zNgNgtemp*NgNgtemp!outer_product(psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k))
    
    call mkl_zcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                        ,zNgNgtemp,Ne*Nquadrature,ZZERO,zNnNgtemp,Nn)
    
    !     call mkl_zcsrmm('T',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !                                     ,transpose(zNnNgtemp),Ne*Nquadrature,ZZERO,chi_matrix,Nn)
    
    !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,chi_matrix&
    !     ,NNmatrix_temp,pardiso_info)
    
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,transpose(NNmatrix_temp)&
    !     ,chi_matrix,pardiso_info)
    
        ! NNmatrix_temp005=NNmatrix_temp005+2.0d0*chi_matrix!+4.0d0*G*outer_product(psi(:,k),psi(:,k))

                                        zNnNgtemp2=zNnNgtemp2+2.0d0*zNnNgtemp
        
        
                    enddo
        
        
                    ! CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,NNmatrix_temp02,Nn,0.0d0,NNmatrix_temp03,Nn)
                    ! CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,NNmatrix_temp005,Nn,(0.0d0,0.0d0),NNmatrix_temp,Nn)

            CALL ZGEMM('T','T',Nn,Nn,Ne*Nquadrature,ZONE,v_kernel*j_real,Ne*Nquadrature,zNnNgtemp2,Nn,ZZERO,NNmatrix_temp005,Nn)
    
                    ! call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,NNmatrix_temp03,Nn&
                    !                                                         ,0.0d0,NNmatrix_temp02,Nn)
                    ! call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,NNmatrix_temp,Nn&
                    !                                                         ,(0.0d0,0.0d0),NNmatrix_temp005,Nn)
        
                    NNmatrix_temp=IdentityMatrix-NNmatrix_temp005

        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,v00*j_real&
            ,NNmatrix_temp006,pardiso_info)
        
                    ! Ytemp2=Ytemp3
                    NNmatrix_temp005=NNmatrix_temp006
        
                    CALL ZGETRF(Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,lapack_INFO)
                    ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                    CALL ZGETRS('N',Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,NNmatrix_temp005,Nn,lapack_INFO)
    
                    ! print *,Ytemp2(1:3,1)
        
                    ! Ytemp2=Ytemp2-Ytemp3
                    NNmatrix_temp=NNmatrix_temp005-NNmatrix_temp006
        
        ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp2,Ytemp4,pardiso_info)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,NNmatrix_temp005&
            ! ,NNmatrix_temp,pardiso_info)


            call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
        call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                    ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)
    
    

            call outer_product('ge',psi_point_g(:,i),psi_point_g(:,i),NgNgtemp)
            zNgNgtemp=zNgNgtemp*NgNgtemp!outer_product(psi_point_g(:,i),psi_point_g(:,i))


                gw_Sigma_g_poles=gw_Sigma_g_poles+zNgNgtemp
    
    
                ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
                !                                     ,NNmatrix_temp,Nn,ZZERO,NNmatrix_temp005,Ne*Nquadrature)
    
                ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
                !                                     ,transpose(NNmatrix_temp005),Nn,ZZERO,NNmatrix_temp,Ne*Nquadrature)
    
    
                ! gw_Sigma_complex_poles=gw_Sigma_complex_poles-NNmatrix_temp
    
    
            end if
    
            enddo

            deallocate(NNmatrix_temp005)
            deallocate(zNnNgtemp2)
            deallocate(NNmatrix_temp)
            deallocate(zNgNntemp)
            deallocate(zNgNgtemp)
            deallocate(zNnNgtemp)
            deallocate(NNmatrix_temp006)
            deallocate(NgNgtemp)


            end if  !!! hadamard 










            
    
    
            deallocate(lapack_IPIV)
            deallocate(pt_v)
            deallocate(iparm_v)

            deallocate(matrix_temp0)
            deallocate(matrix_temp1)
            deallocate(matrix_temp2)
            deallocate(G)

            





            else if (type=='hf') then


            !!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! pardiso parameters !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            MAXFCT=1
            MNUM=1
            MTYPE=6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
            MTYPE_real=1
            MTYPE_complx=3
            MSGLVL=0
            allocate(pt_v(64))
            allocate(iparm_v(64))
            pt_V=0
            call pardisoinit(PT_V,MTYPE,IPARM_V)
            PHASE=13


            ! allocate(matrix_temp0(Nn,Nn))
            allocate(matrix_temp(Nn,Nn))
            allocate(matrix_temp2(Nn,Nn))
            allocate(matrix_temp3(Nn,Nn))
            allocate(G_real(Nn,Nn))
            allocate(G(Nn,Nn))
            
            lapack_lwork=Nn
            allocate(lapack_work(lapack_lwork))
            

            gw_Sigma_complex_poles=ZZERO
            ! gw_Sigma_g_poles=ZZERO
    
            do i=1,Nstates
    
                if (E_dft(i)-dble(Ze)>0.0d0) then
        
    
                    matrix_temp3=0.0d0
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp2=0.0d0

        
        
            !             psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

            matrix_temp=(E_dft(k)+(E_dft(i)-dble(Ze)))*S_dense-Hhf_dense
            G_real=IdentityMatrix_real

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call dsytrf('U',Nn,matrix_temp,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call dsytrs('U',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G_real,Nn,lapack_INFO)

            ! print *,'toto'

        
                        matrix_temp2=matrix_temp2+G_real
        
            !             psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

            matrix_temp=(E_dft(k)-(E_dft(i)-dble(Ze)))*S_dense-Hhf_dense
            G_real=IdentityMatrix_real

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call dsytrf('U',Nn,matrix_temp,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call dsytrs('U',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G_real,Nn,lapack_INFO)

            ! print *,'toto'

        
                        matrix_temp2=matrix_temp2+G_real


            call outer_product('ge',psi(:,k),psi(:,k),matrix_temp)
            matrix_temp3=matrix_temp3+2.0d0*matrix_temp2*matrix_temp!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,matrix_temp3,Nn,0.0d0,matrix_temp2,Nn)
    
                call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,matrix_temp2,Nn&
                                                                            ,0.0d0,matrix_temp3,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp2=IdentityMatrix_real-matrix_temp3
                
                matrix_temp3=v00
    
                CALL DGETRF(Nn,Nn,matrix_temp2,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL DGETRS('N',Nn,Nn,matrix_temp2,Nn,lapack_IPIV,matrix_temp3,Nn,lapack_INFO)
    
                matrix_temp3=matrix_temp3-v00

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp3*j_real&
                ,G,pardiso_info)

                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles-G*matrix_temp!outer_product(psi(:,i),psi(:,i))

                end if

        
            enddo


            do i=Nstates+1,M00
    
                if (E_dft(i)-dble(Ze)<0.0d0) then


                    matrix_temp3=0.0d0
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp2=0.0d0

        
        
            !             psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

            matrix_temp=(E_dft(k)+(E_dft(i)-dble(Ze)))*S_dense-Hhf_dense
            G_real=IdentityMatrix_real

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call dsytrf('U',Nn,matrix_temp,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call dsytrs('U',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G_real,Nn,lapack_INFO)

        
                        matrix_temp2=matrix_temp2+G_real
        
            !             psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

            matrix_temp=(E_dft(k)-(E_dft(i)-dble(Ze)))*S_dense-Hhf_dense
            G_real=IdentityMatrix_real

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call dsytrf('U',Nn,matrix_temp,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call dsytrs('U',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G_real,Nn,lapack_INFO)

        
                        matrix_temp2=matrix_temp2+G_real


            call outer_product('ge',psi(:,k),psi(:,k),matrix_temp)
            matrix_temp3=matrix_temp3+2.0d0*matrix_temp2*matrix_temp!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,matrix_temp3,Nn,0.0d0,matrix_temp2,Nn)
    
                call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,matrix_temp2,Nn&
                                                                            ,0.0d0,matrix_temp3,Nn)


            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp2=IdentityMatrix_real-matrix_temp3
                
                matrix_temp3=v00
    
                CALL DGETRF(Nn,Nn,matrix_temp2,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL DGETRS('N',Nn,Nn,matrix_temp2,Nn,lapack_IPIV,matrix_temp3,Nn,lapack_INFO)
    
                matrix_temp3=matrix_temp3-v00

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp3*j_real&
                ,G,pardiso_info)

                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles+G*matrix_temp!outer_product(psi(:,i),psi(:,i))*j_real


                end if
        
            enddo

            ! deallocate(matrix_temp0)
            deallocate(matrix_temp)
            deallocate(matrix_temp2)
            deallocate(matrix_temp3)
            deallocate(G_real)
            deallocate(G)
            deallocate(lapack_work)

            deallocate(lapack_IPIV)
            deallocate(pt_v)
            deallocate(iparm_v)



        end if

    
    
        end subroutine compute_gw_Sigma_CD_poles_mtx



        subroutine compute_gw_Sigma_CD_poles_mtx_2(Nn,Nstates,Nempty,Ze,E_dft,M00,M0,nnza,Ne,Nquadrature,type,hadamard)
            integer,intent(in) :: Nn,Nstates,Nempty,M00,M0,nnza,Ne,Nquadrature
            complex(kind=(kind(1.0d0))),intent(in) :: Ze
            double precision,dimension(1:Nstates+10),intent(in) :: E_dft
            ! complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0),intent(in) :: Y
            character(len=*), intent(in) :: type,hadamard
    
            integer :: i,ii,k,ll
            ! complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0) :: Ytemp1,Ytemp2,Ytemp3,Ytemp4,Ytemp5,Ytemp6,Ytemp7
            double precision :: norm
            !!!!! pardiso parameters
            integer :: MAXFCT,MNUM,MTYPE,MTYPE_real,MTYPE_complx,MSGLVL,PHASE,idum,pardiso_info
            integer(8),dimension(:),allocatable :: pt_V
            integer,dimension(:),allocatable :: iparm_V

            integer :: lapack_INFO

            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: matrix_temp,matrix_temp1,matrix_temp2,G
            complex(kind=(kind(1.0d0))),dimension(:),allocatable :: lapack_work
            double precision,dimension(:,:),allocatable :: matrix_temp_real,NgNgtemp_real


            ! double precision,dimension(:,:),allocatable :: matrix_temp,matrix_temp2,matrix_temp3,G_real,matrix_temp0,matrix_temp1
            ! complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: G
            ! complex(kind=(kind(1.0d0))),dimension(:),allocatable :: lapack_work

            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: NNmatrix_temp005,NNmatrix_temp006,NNmatrix_temp
            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: zNnNgtemp2,zNgNntemp,zNgNgtemp,zNnNgtemp
            

            allocate(lapack_IPIV(1:Nn))

            if (type=='dft') then
    
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! pardiso parameters !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            MAXFCT=1
            MNUM=1
            MTYPE=6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
            MTYPE_real=1
            MTYPE_complx=3
            MSGLVL=0
            allocate(pt_v(64))
            allocate(iparm_v(64))
            pt_V=0
            call pardisoinit(PT_V,MTYPE,IPARM_V)
            PHASE=13
    
                allocate(matrix_temp(Nn,Nn))
                allocate(matrix_temp1(Nn,Nn))
                allocate(matrix_temp_real(Nn,Nn))
                allocate(G(Nn,Nn))


                if (hadamard=='no') then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Approximated integral without Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! allocate(matrix_temp(Nn,Nn))
            
            ! lapack_lwork=Nn
            ! allocate(lapack_work(lapack_lwork))


            gw_Sigma_complex_poles=ZZERO
            
            ! gw_Sigma_g_poles=ZZERO
    
            do i=1,Nstates
    
                if (E_dft(i)-dble(Ze)>0.0d0) then
                ! if ((E_dft(i)-dble(Ze)<0.0d0).and.(aimag(Ze)>0.0d0)) then

                    ! print *,'toto'
        
    
                    matrix_temp=cmplx(0.0d0,0.0d0,8)
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp1=cmplx(0.0d0,0.0d0,8)

        
        
                        psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                    matrix_temp1=matrix_temp1+G
        
                        psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                    matrix_temp1=matrix_temp1+G

                    call outer_product('ge',psi(:,k),psi(:,k),matrix_temp_real)
                    matrix_temp=matrix_temp+2.0d0*matrix_temp1*matrix_temp_real!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,matrix_temp,Nn,(0.0d0,0.0d0),matrix_temp1,Nn)
    
                call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,matrix_temp1,Nn&
                                                                            ,(0.0d0,0.0d0),matrix_temp,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp1=IdentityMatrix-matrix_temp
                
                matrix_temp=v00*j_real
    
                CALL ZGETRF(Nn,Nn,matrix_temp1,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL ZGETRS('N',Nn,Nn,matrix_temp1,Nn,lapack_IPIV,matrix_temp,Nn,lapack_INFO)
    
                matrix_temp=matrix_temp-v00*j_real

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


                call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp&
                ,matrix_temp1,pardiso_info)


                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp_real)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles-matrix_temp1*matrix_temp_real!outer_product(psi(:,i),psi(:,i))


                end if
        
            enddo


            do i=Nstates+1,M00
    
                if (E_dft(i)-dble(Ze)<0.0d0) then


                    matrix_temp=cmplx(0.0d0,0.0d0,8)
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp1=cmplx(0.0d0,0.0d0,8)

        
        
                        psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        matrix_temp1=matrix_temp1+G
        
                        psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
                        matrix_temp1=matrix_temp1+G


           


                        call outer_product('ge',psi(:,k),psi(:,k),matrix_temp_real)
                        matrix_temp=matrix_temp+2.0d0*matrix_temp1*matrix_temp_real!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,matrix_temp,Nn,(0.0d0,0.0d0),matrix_temp1,Nn)
    
                call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,matrix_temp1,Nn&
                                                                            ,(0.0d0,0.0d0),matrix_temp,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp1=IdentityMatrix-matrix_temp
                
                matrix_temp=v00*j_real
    
                CALL ZGETRF(Nn,Nn,matrix_temp1,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL ZGETRS('N',Nn,Nn,matrix_temp1,Nn,lapack_IPIV,matrix_temp,Nn,lapack_INFO)
    
                matrix_temp=matrix_temp-v00*j_real

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


                call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp&
                ,matrix_temp1,pardiso_info)

                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp_real)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles+matrix_temp1*matrix_temp_real!outer_product(psi(:,i),psi(:,i))


                end if
        
            enddo
    
    
    

            else if (hadamard=='yes') then


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!    Fully integral using Hadamard product    !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        allocate(NNmatrix_temp005(Nn,Nn))
        allocate(zNnNgtemp2(Nn,Ne*Nquadrature))
        allocate(NNmatrix_temp(Nn,Nn))
        allocate(zNgNntemp(Ne*Nquadrature,Nn))
        allocate(zNgNgtemp(Ne*Nquadrature,Ne*Nquadrature))
        allocate(zNnNgtemp(Nn,Ne*Nquadrature))
        allocate(NNmatrix_temp006(Nn,Nn))
        allocate(NgNgtemp_real(Ne*Nquadrature,Ne*Nquadrature))


        ! gw_Sigma_complex_poles=ZZERO
        gw_Sigma_g_poles=ZZERO

        do i=1,Nstates

            if (E_dft(i)-dble(Ze)>0.0d0) then
    

                NNmatrix_temp005=cmplx(0.0d0,0.0d0,8)
                zNnNgtemp2=ZZERO
    
                do k=1,Nstates
    
                    NNmatrix_temp=cmplx(0.0d0,0.0d0,8)

                    psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
                    NNmatrix_temp=NNmatrix_temp+dble(G)
    
                    psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
                    NNmatrix_temp=NNmatrix_temp+dble(G)


    call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                    ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)

call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                    ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)

call outer_product('ge',psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k),NgNgtemp_real)
zNgNgtemp=zNgNgtemp*NgNgtemp_real!outer_product(psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k))

call mkl_zcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                    ,zNgNgtemp,Ne*Nquadrature,ZZERO,zNnNgtemp,Nn)

!     call mkl_zcsrmm('T',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                                     ,transpose(zNnNgtemp),Ne*Nquadrature,ZZERO,chi_matrix,Nn)

!     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,chi_matrix&
!     ,NNmatrix_temp,pardiso_info)

! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,transpose(NNmatrix_temp)&
!     ,chi_matrix,pardiso_info)

    ! NNmatrix_temp005=NNmatrix_temp005+2.0d0*chi_matrix!+4.0d0*G*outer_product(psi(:,k),psi(:,k))

                                    zNnNgtemp2=zNnNgtemp2+2.0d0*zNnNgtemp
    
    
                enddo
    
    
                ! CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,NNmatrix_temp02,Nn,0.0d0,NNmatrix_temp03,Nn)
                ! CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,NNmatrix_temp005,Nn,(0.0d0,0.0d0),NNmatrix_temp,Nn)

        CALL ZGEMM('T','T',Nn,Nn,Ne*Nquadrature,ZONE,v_kernel*j_real,Ne*Nquadrature,zNnNgtemp2,Nn,ZZERO,NNmatrix_temp005,Nn)

                ! call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,NNmatrix_temp03,Nn&
                !                                                         ,0.0d0,NNmatrix_temp02,Nn)
                ! call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,NNmatrix_temp,Nn&
                !                                                         ,(0.0d0,0.0d0),NNmatrix_temp005,Nn)
    
                NNmatrix_temp=IdentityMatrix-NNmatrix_temp005

    call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,v00*j_real&
        ,NNmatrix_temp006,pardiso_info)
    
                ! Ytemp2=Ytemp3
                NNmatrix_temp005=NNmatrix_temp006
    
                CALL ZGETRF(Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL ZGETRS('N',Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,NNmatrix_temp005,Nn,lapack_INFO)

                ! print *,Ytemp2(1:3,1)
    
                ! Ytemp2=Ytemp2-Ytemp3
                NNmatrix_temp=NNmatrix_temp005-NNmatrix_temp006
    
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp2,Ytemp4,pardiso_info)
        ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,NNmatrix_temp005&
        ! ,NNmatrix_temp,pardiso_info)


        call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)

    call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)

        call outer_product('ge',psi_point_g(:,i),psi_point_g(:,i),NgNgtemp_real)
        zNgNgtemp=zNgNgtemp*NgNgtemp_real!outer_product(psi_point_g(:,i),psi_point_g(:,i))


            gw_Sigma_g_poles=gw_Sigma_g_poles-zNgNgtemp


            ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
            !                                     ,NNmatrix_temp,Nn,ZZERO,NNmatrix_temp005,Ne*Nquadrature)

            ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
            !                                     ,transpose(NNmatrix_temp005),Nn,ZZERO,NNmatrix_temp,Ne*Nquadrature)


            ! gw_Sigma_complex_poles=gw_Sigma_complex_poles-NNmatrix_temp

        end if

        enddo

        ! call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp7,Nn,(0.0d0,0.0d0),Ytemp5,Nn)

        ! Y_primed=Ytemp5

        ! go to 2222

        ! Ytemp7=(0.0d0,0.0d0)
    do i=Nstates+1,M00

        if (E_dft(i)-dble(Ze)<0.0d0) then


            NNmatrix_temp005=cmplx(0.0d0,0.0d0,8)
                zNnNgtemp2=ZZERO
    
                do k=1,Nstates
    
                    NNmatrix_temp=cmplx(0.0d0,0.0d0,8)


                    psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
                    NNmatrix_temp=NNmatrix_temp+dble(G)
    
                    psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
        call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
                    NNmatrix_temp=NNmatrix_temp+dble(G)


    call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                    ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)

call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                    ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)

call outer_product('ge',psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k),NgNgtemp_real)
zNgNgtemp=zNgNgtemp*NgNgtemp_real!outer_product(psi_point_g(:,k)*volumegweight(:),psi_point_g(:,k))

call mkl_zcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                    ,zNgNgtemp,Ne*Nquadrature,ZZERO,zNnNgtemp,Nn)

!     call mkl_zcsrmm('T',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
!                                     ,transpose(zNnNgtemp),Ne*Nquadrature,ZZERO,chi_matrix,Nn)

!     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,chi_matrix&
!     ,NNmatrix_temp,pardiso_info)

! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,transpose(NNmatrix_temp)&
!     ,chi_matrix,pardiso_info)

    ! NNmatrix_temp005=NNmatrix_temp005+2.0d0*chi_matrix!+4.0d0*G*outer_product(psi(:,k),psi(:,k))

                                    zNnNgtemp2=zNnNgtemp2+2.0d0*zNnNgtemp
    
    
                enddo
    
    
                ! CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,NNmatrix_temp02,Nn,0.0d0,NNmatrix_temp03,Nn)
                ! CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,NNmatrix_temp005,Nn,(0.0d0,0.0d0),NNmatrix_temp,Nn)

        CALL ZGEMM('T','T',Nn,Nn,Ne*Nquadrature,ZONE,v_kernel*j_real,Ne*Nquadrature,zNnNgtemp2,Nn,ZZERO,NNmatrix_temp005,Nn)

                ! call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,NNmatrix_temp03,Nn&
                !                                                         ,0.0d0,NNmatrix_temp02,Nn)
                ! call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,NNmatrix_temp,Nn&
                !                                                         ,(0.0d0,0.0d0),NNmatrix_temp005,Nn)
    
                NNmatrix_temp=IdentityMatrix-NNmatrix_temp005

    call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,v00*j_real&
        ,NNmatrix_temp006,pardiso_info)
    
                ! Ytemp2=Ytemp3
                NNmatrix_temp005=NNmatrix_temp006
    
                CALL ZGETRF(Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL ZGETRS('N',Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,NNmatrix_temp005,Nn,lapack_INFO)

                ! print *,Ytemp2(1:3,1)
    
                ! Ytemp2=Ytemp2-Ytemp3
                NNmatrix_temp=NNmatrix_temp005-NNmatrix_temp006
    
    ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp2,Ytemp4,pardiso_info)
        ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,NNmatrix_temp005&
        ! ,NNmatrix_temp,pardiso_info)


        call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)

    call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                ,transpose(zNgNntemp),Nn,ZZERO,zNgNgtemp,Ne*Nquadrature)



        call outer_product('ge',psi_point_g(:,i),psi_point_g(:,i),NgNgtemp_real)
        zNgNgtemp=zNgNgtemp*NgNgtemp_real!outer_product(psi_point_g(:,i),psi_point_g(:,i))


            gw_Sigma_g_poles=gw_Sigma_g_poles+zNgNgtemp


            ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
            !                                     ,NNmatrix_temp,Nn,ZZERO,NNmatrix_temp005,Ne*Nquadrature)

            ! call mkl_zcsrmm('N',Nn,Nn,Nn,ZONE,matdescrb,B*j_real,JA,IA_pntrb,IA_pntre&
            !                                     ,transpose(NNmatrix_temp005),Nn,ZZERO,NNmatrix_temp,Ne*Nquadrature)


            ! gw_Sigma_complex_poles=gw_Sigma_complex_poles-NNmatrix_temp


        end if

        enddo

        deallocate(NNmatrix_temp005)
        deallocate(zNnNgtemp2)
        deallocate(NNmatrix_temp)
        deallocate(zNgNntemp)
        deallocate(zNgNgtemp)
        deallocate(zNnNgtemp)
        deallocate(NNmatrix_temp006)
        deallocate(NgNgtemp_real)

            
            end if !!! hadamard
    
    
            deallocate(lapack_IPIV)
            deallocate(pt_v)
            deallocate(iparm_v)

            deallocate(matrix_temp)
            deallocate(matrix_temp1)
            deallocate(matrix_temp_real)
            deallocate(G)

    else if (type=='hf') then

        !!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! pardiso parameters !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        MAXFCT=1
        MNUM=1
        MTYPE=6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
        MTYPE_real=1
        MTYPE_complx=3
        MSGLVL=0
        allocate(pt_v(64))
        allocate(iparm_v(64))
        pt_V=0
        call pardisoinit(PT_V,MTYPE,IPARM_V)
        PHASE=13

            allocate(matrix_temp(Nn,Nn))
            allocate(matrix_temp1(Nn,Nn))
            allocate(matrix_temp2(Nn,Nn))
            allocate(matrix_temp_real(Nn,Nn))
            allocate(G(Nn,Nn))
            
            lapack_lwork=Nn
            allocate(lapack_work(lapack_lwork))


    
            gw_Sigma_complex_poles=ZZERO
            ! gw_Sigma_g_poles=ZZERO
    
            do i=1,Nstates
    
                if (E_dft(i)-dble(Ze)>0.0d0) then
                ! if ((E_dft(i)-dble(Ze)<0.0d0).and.(aimag(Ze)>0.0d0)) then

                    ! print *,'toto'
        
    
                    matrix_temp1=cmplx(0.0d0,0.0d0,8)
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp=cmplx(0.0d0,0.0d0,8)

        
        
            !             psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
            !             NNmatrix_temp=NNmatrix_temp+G
        
            !             psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
            !             NNmatrix_temp=NNmatrix_temp+G


            matrix_temp2=(E_dft(k)+(E_dft(i)-Ze))*S_dense-Hhf_dense
            G=IdentityMatrix

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call zsytrf('U',Nn,matrix_temp2,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call zsytrs('U',Nn,Nn,matrix_temp2,Nn,lapack_IPIV,G,Nn,lapack_INFO)

        
                        matrix_temp=matrix_temp+G
        
            !             psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

            matrix_temp2=(E_dft(k)-(E_dft(i)-Ze))*S_dense-Hhf_dense
            G=IdentityMatrix

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call zsytrf('U',Nn,matrix_temp2,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call zsytrs('U',Nn,Nn,matrix_temp2,Nn,lapack_IPIV,G,Nn,lapack_INFO)

        
                        matrix_temp=matrix_temp+G


                        call outer_product('ge',psi(:,k),psi(:,k),matrix_temp_real)
                        matrix_temp1=matrix_temp1+2.0d0*matrix_temp*matrix_temp_real!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,matrix_temp1,Nn,(0.0d0,0.0d0),matrix_temp,Nn)
    
                call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,matrix_temp,Nn&
                                                                            ,(0.0d0,0.0d0),matrix_temp1,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp=IdentityMatrix-matrix_temp1
                
                matrix_temp1=v00*j_real
    
                CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,matrix_temp1,Nn,lapack_INFO)
    
                matrix_temp1=matrix_temp1-v00*j_real

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


                call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp1&
                ,matrix_temp,pardiso_info)

                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp_real)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles-matrix_temp*matrix_temp_real!outer_product(psi(:,i),psi(:,i))
                


                end if
        
            enddo


            do i=Nstates+1,M00
    
                if (E_dft(i)-dble(Ze)<0.0d0) then


                    matrix_temp1=cmplx(0.0d0,0.0d0,8)
                    ! zNnNgtemp2=ZZERO
        
                    do k=1,Nstates
        
                        matrix_temp=cmplx(0.0d0,0.0d0,8)

        
        
            !             psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
            !             NNmatrix_temp=NNmatrix_temp+G
        
            !             psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
            !             NNmatrix_temp=NNmatrix_temp+G


            matrix_temp2=(E_dft(k)+(E_dft(i)-Ze))*S_dense-Hhf_dense
            G=IdentityMatrix

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call zsytrf('L',Nn,matrix_temp2,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call zsytrs('L',Nn,Nn,matrix_temp2,Nn,lapack_IPIV,G,Nn,lapack_INFO)

        
                        matrix_temp=matrix_temp+G
        
            !             psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
            ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)

            matrix_temp2=(E_dft(k)-(E_dft(i)-Ze))*S_dense-Hhf_dense
            G=IdentityMatrix

            ! CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
            ! CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,G,Nn,lapack_INFO)

            call zsytrf('L',Nn,matrix_temp2,Nn,lapack_IPIV,lapack_work,lapack_lwork,lapack_INFO)
            call zsytrs('L',Nn,Nn,matrix_temp2,Nn,lapack_IPIV,G,Nn,lapack_INFO)

        
                        matrix_temp=matrix_temp+G


                        call outer_product('ge',psi(:,k),psi(:,k),matrix_temp_real)
                        matrix_temp1=matrix_temp1+2.0d0*matrix_temp*matrix_temp_real!outer_product(psi(:,k),psi(:,k))



                    enddo



                CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,matrix_temp1,Nn,(0.0d0,0.0d0),matrix_temp,Nn)
    
                call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,matrix_temp,Nn&
                                                                            ,(0.0d0,0.0d0),matrix_temp1,Nn)

            ! NNmatrix_temp02=NNmatrix_temp03

                matrix_temp=IdentityMatrix-matrix_temp1
                
                matrix_temp1=v00*j_real
    
                CALL ZGETRF(Nn,Nn,matrix_temp,Nn,lapack_IPIV,lapack_INFO)
                ! CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
                CALL ZGETRS('N',Nn,Nn,matrix_temp,Nn,lapack_IPIV,matrix_temp1,Nn,lapack_INFO)
    
                matrix_temp1=matrix_temp1-v00*j_real

                ! call pardisoinit(PT_V,MTYPE_real,IPARM_V)


                call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,matrix_temp1&
                ,matrix_temp,pardiso_info)

                call outer_product('ge',psi(:,i),psi(:,i),matrix_temp_real)
                gw_Sigma_complex_poles=gw_Sigma_complex_poles+matrix_temp*matrix_temp_real!outer_product(psi(:,i),psi(:,i))


                end if
        
            enddo

            deallocate(matrix_temp)
            deallocate(lapack_work)
            deallocate(lapack_IPIV)
            deallocate(pt_v)
            deallocate(iparm_v)

            deallocate(matrix_temp1)
            deallocate(matrix_temp2)
            deallocate(matrix_temp_real)
            deallocate(G)

        end if
    
    
        end subroutine compute_gw_Sigma_CD_poles_mtx_2











    !     subroutine compute_gw_Sigma_CD_poles(Nn,Nstates,Nempty,Ze,E_dft,M00,M0,Y,nnza,Ne,Nquadrature)
    !         integer,intent(in) :: Nn,Nstates,Nempty,M00,M0,nnza,Ne,Nquadrature
    !         complex(kind=(kind(1.0d0))),intent(in) :: Ze
    !         double precision,dimension(1:M00),intent(in) :: E_dft
    !         complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0),intent(in) :: Y
    
    !         integer :: i,ii,k,ll
    !         complex(kind=(kind(1.0d0))),dimension(1:Nn,1:M0) :: Ytemp1,Ytemp2,Ytemp3,Ytemp4,Ytemp5,Ytemp6,Ytemp7
    !         double precision :: norm
    !         !!!!! pardiso parameters
    !         integer :: MAXFCT,MNUM,MTYPE,MTYPE_real,MTYPE_complx,MSGLVL,PHASE,idum,pardiso_info
    !         integer(8),dimension(:),allocatable :: pt_V
    !         integer,dimension(:),allocatable :: iparm_V

    !         allocate(lapack_IPIV(1:Nn))
    
    !         !!!!!!!!!!!!!!!!!!!!!!!!!!
    !         !!! pardiso parameters !!!
    !         !!!!!!!!!!!!!!!!!!!!!!!!!!
    !         MAXFCT=1
    !         MNUM=1
    !         MTYPE=6 !!!! 3 for complex and structurally symmetric matrix, 1 for real and structurally symmetric matrix
    !         MTYPE_real=1
    !         MTYPE_complx=3
    !         MSGLVL=0
    !         allocate(pt_v(64))
    !         allocate(iparm_v(64))
    !         pt_V=0
    !         call pardisoinit(PT_V,MTYPE,IPARM_V)
    !         PHASE=13
    
    
    !         call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,Y,Nn,(0.0d0,0.0d0),Ytemp1,Nn)
    
            
    !         ! Y_primed=(0.0d0,0.0d0)
    
    !         Ytemp7=(0.0d0,0.0d0)
    !         ! gw_Sigma_g_poles=(0.0d0,0.0d0)
    
    !         do i=1,Nstates
    
    !             if (E_dft(i)-dble(Ze)>0.0d0) then
        
    !                 ! print *,E(i),E_guess(m)
    
    !                 do ii=1,M0
    !                     Ytemp2(:,ii)=Ytemp1(:,ii)*psi(:,i)
    !                 enddo
        
    !                 CALL ZGEMM('N','N',Nn,M0,Nn,(1.0d0,0.0d0),v00*j_real,Nn,Ytemp2,Nn,(0.0d0,0.0d0),Ytemp3,Nn)
        
    
    
    
    !                 NNmatrix_temp005=cmplx(0.0d0,0.0d0,8)
        
    !                 do k=1,Nstates
        
    !                     NNmatrix_temp=cmplx(0.0d0,0.0d0,8)
        
    !                     ! G_inverse = (E_dft(k)+(E_dft(i)-Ze))*B-H_dft
        
    !                     ! call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,G_inverse,IA,JA,idum,Nn&
    !                     !                     ,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
    !                     psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !         call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
    !                     NNmatrix_temp=NNmatrix_temp+G
        
    !                     ! G_inverse = (E_dft(k)-(E_dft(i)-Ze))*B-H_dft
        
    !                     ! call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,G_inverse,IA,JA,idum,Nn&
    !                     !                     ,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
    !                     psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !         call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
        
    !                     NNmatrix_temp=NNmatrix_temp+G
    
    
    
    
    
    
    ! !     call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ! !                                     ,NNmatrix_temp,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    ! ! call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ! !                                     ,transpose(zNgNntemp),Nn,ZZERO,G_g_complex,Ne*Nquadrature)
    
    ! !     G_g_complex=G_g_complex*outer_product(psi_point_g(:,k),psi_point_g(:,k))*outer_product(volumegweight,volumegweight)
    
    ! ! call mkl_zcsrmm('T',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ! !                                     ,G_g_complex,Ne*Nquadrature,ZZERO,zNnNgtemp,Nn)
    
    ! !     call mkl_zcsrmm('T',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    ! !                                     ,transpose(zNnNgtemp),Ne*Nquadrature,ZZERO,chi_matrix_complex,Nn)
    
    ! !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,chi_matrix_complex&
    ! !     ,NNmatrix_temp003,pardiso_info)
    
    ! !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,transpose(NNmatrix_temp003)&
    ! !     ,chi_matrix_complex,pardiso_info)
    
    
    ! !     NNmatrix_temp002=NNmatrix_temp002+2.0d0*chi_matrix_complex!+4.0d0*G*outer_product(psi(:,k),psi(:,k))
        
        
    !                     NNmatrix_temp005=NNmatrix_temp005+2.0d0*NNmatrix_temp*outer_product(psi(:,k),psi(:,k))
    
        
    !                 enddo
        
        
    !                 ! CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,NNmatrix_temp02,Nn,0.0d0,NNmatrix_temp03,Nn)
    !                 CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,NNmatrix_temp005,Nn,(0.0d0,0.0d0),NNmatrix_temp,Nn)
    
    !                 ! call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,NNmatrix_temp03,Nn&
    !                 !                                                         ,0.0d0,NNmatrix_temp02,Nn)
    !                 call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,NNmatrix_temp,Nn&
    !                                                                         ,(0.0d0,0.0d0),NNmatrix_temp005,Nn)
        
    !                 NNmatrix_temp=IdentityMatrix-NNmatrix_temp005
        
    !                 Ytemp2=Ytemp3
    !                 ! NNmatrix_temp002=v00*j_real
        
    !                 CALL ZGETRF(Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,lapack_INFO)
    !                 CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
    !                 ! CALL ZGETRS('N',Nn,Nn,NNmatrix_temp003,Nn,lapack_IPIV,NNmatrix_temp002,Nn,lapack_INFO)
    
    !                 ! print *,Ytemp2(1:3,1)
        
    !                 Ytemp2=Ytemp2-Ytemp3
    !                 ! NNmatrix_temp002=NNmatrix_temp002-v00*j_real
        
    !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp2,Ytemp4,pardiso_info)
    !         ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,NNmatrix_temp002&
    !         ! ,NNmatrix_temp003,pardiso_info)
    
    
    !                 do ii=1,M0
    !                     Ytemp4(:,ii)=-Ytemp4(:,ii)*psi(:,i)
    !                 enddo
    
    !         ! call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp4,Nn,(0.0d0,0.0d0),Ytemp5,Nn)
    
    !                 ! print *,Ytemp5(1:3,1)
    
    
    !     !     call mkl_zcsrmm('N',Ne*Nquadrature,Nn,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !     !                                             ,NNmatrix_temp003,Nn,ZZERO,zNgNntemp,Ne*Nquadrature)
    
    !     ! call mkl_zcsrmm('N',Ne*Nquadrature,Ne*Nquadrature,Nn,ZONE,matdescra,NnToNg*j_real,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
    !     !                                             ,transpose(zNgNntemp),Nn,ZZERO,gw_W_g,Ne*Nquadrature)
    
    
    !         ! gw_W_g=gw_W_g*outer_product(psi_point_g(:,i),psi_point_g(:,i))
    
    !         ! gw_Sigma_g_poles=gw_Sigma_g_poles-gw_W_g*outer_product(psi_point_g(:,i),psi_point_g(:,i))
        
        
    !                 Ytemp7=Ytemp7+Ytemp4
        
    
    
    
    
    !     !         Ytemp2=Ytemp3
    !     !         ll=0
    
    !     !         whileloop1: do while (.true.)
    
    !     !             ! if (ll>0) Ytemp2=Ytemp5
    
    !     !             Ytemp5=(0.0d0,0.0d0)
    
    !     !             do k=1,Nstates
    
    !     !                 do ii=1,M0
    !     !                     Ytemp4(:,ii)=Ytemp2(:,ii)*psi(:,k)
    !     !                 enddo
    
    !     !                 psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !     ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp4,Ytemp6,pardiso_info)
    
    !     !                 do ii=1,M0
    !     !                     Ytemp6(:,ii)=Ytemp6(:,ii)*psi(:,k)
    !     !                 enddo
    
    !     !                 Ytemp5=Ytemp5+Ytemp6
    
    
    !     !                 psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !     ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp4,Ytemp6,pardiso_info)
    
    !     !                 do ii=1,M0
    !     !                     Ytemp6(:,ii)=Ytemp6(:,ii)*psi(:,k)
    !     !                 enddo
    
    !     !                 Ytemp5=Ytemp5+Ytemp6
    
    !     !             enddo
    
    !     !             Ytemp5=2.0d0*Ytemp5
    
    !     !             CALL ZGEMM('T','N',Nn,M0,Nn,(1.0d0,0.0d0),v00*j_real,Nn,Ytemp5,Nn,(0.0d0,0.0d0),Ytemp4,Nn)
    
    
    !     !     call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp4,Nn,(0.0d0,0.0d0),Ytemp5,Nn)
    
    !     !             Ytemp5=Ytemp5+Ytemp3
    
    !     !             norm=0.0d0
    !     !             do ii=1,M0
    !     !                 norm=norm+norm2(abs(Ytemp5(:,ii)-Ytemp2(:,ii)))
    !     !             enddo
    
    !     !             ! if (norm<10**(LOG10(minval(abs(Ytemp5(:,:))))-(LOG10(dble(Nn*M0))-4.0d0))) then
    !     !             if (norm<10**(LOG10(minval(abs(Ytemp5(:,:)))))) then
    !     !                 exit whileloop1
    !     !             end if
    
    !     !             Ytemp2=Ytemp5
    
    !     !             ll=ll+1
    
    !     !         enddo whileloop1
    
    !     !         print *,'ll',ll
    
    !     !         ! print *,Ytemp5(1:3,1)
    
    !     !         Ytemp5=Ytemp5-Ytemp3
    
    !     !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp5,Ytemp4,pardiso_info)
    
    !     !         do ii=1,M0
    !     !             Ytemp4(:,ii)=Ytemp4(:,ii)*psi(:,i)
    !     !         enddo
    
    
    !     !         Ytemp7=Ytemp7-Ytemp4
    
    !         end if
    
    !         enddo
    
    !         ! call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp7,Nn,(0.0d0,0.0d0),Ytemp5,Nn)
    
    !         ! Y_primed=Ytemp5
    
    !         ! go to 2222
    
    !         ! Ytemp7=(0.0d0,0.0d0)
    !     do i=Nstates+1,M00
    
    !         if (E_dft(i)-dble(Ze)<0.0d0) then
    
    !             ! print *,E(i),E_guess(m)
    
    !             do ii=1,M0
    !                 Ytemp2(:,ii)=Ytemp1(:,ii)*psi(:,i)
    !             enddo
    
    !             CALL ZGEMM('N','N',Nn,M0,Nn,(1.0d0,0.0d0),v00*j_real,Nn,Ytemp2,Nn,(0.0d0,0.0d0),Ytemp3,Nn)
    
    
    
    
    !             NNmatrix_temp005=cmplx(0.0d0,0.0d0,8)
    
    !             do k=1,Nstates
    
    !                 NNmatrix_temp=cmplx(0.0d0,0.0d0,8)
    
    !                 ! G_inverse = (E_dft(k)+(E_dft(i)-Ze))*B-H_dft
    
    !                 ! call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,G_inverse,IA,JA,idum,Nn&
    !                 !                     ,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    !                 psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !         call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    !                 NNmatrix_temp=NNmatrix_temp+G
    
    !                 ! G_inverse = (E_dft(k)-(E_dft(i)-Ze))*B-H_dft
    
    !                 ! call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,G_inverse,IA,JA,idum,Nn&
    !                 !                     ,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    !                 psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !         call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,Nn,IPARM_V,MSGLVL,IdentityMatrix,G,pardiso_info)
    
    !                 NNmatrix_temp=NNmatrix_temp+G
    
    
    !                 NNmatrix_temp005=NNmatrix_temp005+2.0d0*NNmatrix_temp*outer_product(psi(:,k),psi(:,k))
    
    !             enddo
    
    
    !             ! CALL DGEMM('T','N',Nn,Nn,Nn,1.0d0,v00,Nn,NNmatrix_temp02,Nn,0.0d0,NNmatrix_temp03,Nn)
    !             CALL ZGEMM('T','N',Nn,Nn,Nn,(1.0d0,0.0d0),v00*j_real,Nn,NNmatrix_temp005,Nn,(0.0d0,0.0d0),NNmatrix_temp,Nn)
    
    !             ! call mkl_dcsrmm('N',Nn,Nn,Nn,1.0d0,matdescrb,B,JA,IA_pntrb,IA_pntre,NNmatrix_temp03,Nn&
    !             !                                                         ,0.0d0,NNmatrix_temp02,Nn)
    !             call mkl_zcsrmm('N',Nn,Nn,Nn,(1.0d0,0.0d0),matdescrb,B*j_real,JA,IA_pntrb,IA_pntre,NNmatrix_temp,Nn&
    !                                                                     ,(0.0d0,0.0d0),NNmatrix_temp005,Nn)
    
    !             NNmatrix_temp=IdentityMatrix-NNmatrix_temp005
    
    !             Ytemp2=Ytemp3
    
    !             CALL ZGETRF(Nn,Nn,NNmatrix_temp,Nn,lapack_IPIV,lapack_INFO)
    !             CALL ZGETRS('N',Nn,M0,NNmatrix_temp,Nn,lapack_IPIV,Ytemp2,Nn,lapack_INFO)
    
    !             Ytemp2=Ytemp2-Ytemp3
    
    !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp2,Ytemp4,pardiso_info)
    
    
    !             do ii=1,M0
    !                 Ytemp4(:,ii)=Ytemp4(:,ii)*psi(:,i)
    !             enddo
    
    !             ! call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp4,Nn,(0.0d0,0.0d0),Ytemp5,Nn)
    
    
    !             Ytemp7=Ytemp7+Ytemp4
    
    
    
    
    
    !     !         Ytemp2=Ytemp3
    !     !         ll=0
    
    !     !         whileloop2: do while (.true.)
    
    !     !             Ytemp5=(0.0d0,0.0d0)
    
    !     !             do k=1,Nstates
    
    !     !                 do ii=1,M0
    !     !                     Ytemp4(:,ii)=Ytemp2(:,ii)*psi(:,k)
    !     !                 enddo
    
    !     !                 psaz_dft(1:nnza)=(E_dft(k)+(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !     ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp4,Ytemp6,pardiso_info)
    
    !     !                 do ii=1,M0
    !     !                     Ytemp6(:,ii)=Ytemp6(:,ii)*psi(:,k)
    !     !                 enddo
    
    !     !                 Ytemp5=Ytemp5+Ytemp6
    
    
    !     !                 psaz_dft(1:nnza)=(E_dft(k)-(E_dft(i)-Ze))*usb(1:nnza)-usa_dft(1:nnza)
    !     ! call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,psaz_dft(1),uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp4,Ytemp6,pardiso_info)
    
    !     !                 do ii=1,M0
    !     !                     Ytemp6(:,ii)=Ytemp6(:,ii)*psi(:,k)
    !     !                 enddo
    
    !     !                 Ytemp5=Ytemp5+Ytemp6
    
    !     !             enddo
    
    !     !             Ytemp5=2.0d0*Ytemp5
    
    !     !             CALL ZGEMM('T','N',Nn,M0,Nn,(1.0d0,0.0d0),v00*j_real,Nn,Ytemp5,Nn,(0.0d0,0.0d0),Ytemp4,Nn)
    
    
    !     !     call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp4,Nn,(0.0d0,0.0d0),Ytemp5,Nn)
    
    
    !     !             Ytemp5=Ytemp5+Ytemp3
    
    !     !             norm=0.0d0
    !     !             do ii=1,M0
    !     !                 norm=norm+norm2(abs(Ytemp5(:,ii)-Ytemp2(:,ii)))
    !     !             enddo
    
    !     !             ! if (norm<10**(LOG10(minval(abs(Ytemp5(:,:))))-(LOG10(dble(Nn*M0))-4.0d0))) then
    !     !             if (norm<10**(LOG10(minval(abs(Ytemp5(:,:)))))) then
    !     !                 exit whileloop2
    !     !             end if
    
    !     !             Ytemp2=Ytemp5
    
    !     !             ll=ll+1
    
    !     !         enddo whileloop2
    
    !     !         print *,'ll',ll
    
    !     !         Ytemp5=Ytemp5-Ytemp3
    
    !     !     call PARDISO(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,usb(1:nnza)*j_real,uisa,ujsa,idum,M0,IPARM_V,MSGLVL,Ytemp5,Ytemp6,pardiso_info)
    
    !     !         do ii=1,M0
    !     !             Ytemp6(:,ii)=Ytemp6(:,ii)*psi(:,i)
    !     !         enddo
    
    
    !     !         Ytemp7=Ytemp7+Ytemp6
    
    !         end if
    
    !         enddo
    
    !         ! 2222 continue
    
    !         ! call mkl_zcsrmm('N',Nn,M0,Nn,(1.0d0,0.0d0),matdescrb,B_complex,JA,IA_pntrb,IA_pntre,Ytemp7,Nn,(0.0d0,0.0d0),Ytemp5,Nn)
    
    !         Y_primed=Ytemp7
    
    
    !         deallocate(lapack_IPIV)
    !         deallocate(pt_v)
    !         deallocate(iparm_v)
    
    
    !     end subroutine compute_gw_Sigma_CD_poles





















        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!           casida self-energy            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        subroutine compute_gw_Casida_Kx(Nn,Ne,Nquadrature,Nstates,Nempty)
            integer,intent(in) :: Nn,Ne,Nquadrature,Nstates,Nempty
            integer :: k
    
            CALL DGEMM('T','N',Nn,Nstates*Nempty,Ne*Nquadrature,1.0d0,v_kernel,Ne*Nquadrature,psi_ia_g_real,Ne*Nquadrature&
            ,0.0d0,NNmatrix_temp12,Nn)
    
            call mkl_dcsrmm('N',Ne*Nquadrature,Nstates*Nempty,Nn,1.0d0,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
                                                        ,NNmatrix_temp12,Nn,0.0d0,NgNntemp,Ne*Nquadrature)
    
            do k=1,Nstates*Nempty
                NgNntemp(:,k)=NgNntemp(:,k)*volumegweight(:)
            enddo
    
            ! CALL DGEMM('N','N',Nstates*Nempty,Nstates*Nempty,Ne*Nquadrature,1.0d0,transpose(psi_ia_g_real),Nstates*Nempty&
            !                                     ,NgNntemp,Ne*Nquadrature,0.0d0,casida_Kx,Nstates*Nempty)
            CALL DGEMM('T','N',Nstates*Nempty,Nstates*Nempty,Ne*Nquadrature,1.0d0,psi_ia_g_real,Ne*Nquadrature&
                                                ,NgNntemp,Ne*Nquadrature,0.0d0,casida_Kx,Nstates*Nempty)
    
    
        end subroutine compute_gw_Casida_Kx
    
        subroutine compute_gw_Casida_matrix(Nstates,Nempty,E_casida)
            integer,intent(in) :: Nstates,Nempty
            double precision,dimension(:),intent(in) :: E_casida
            integer :: k,kk
            double precision,dimension(:,:),allocatable :: temp
    
            casida_Kx=4.0d0*casida_Kx
            ! casida_Kx=2.0d0*casida_Kx
            ! casida_Kx=casida_Kx
    
            do k=1,Nstates
                do kk=1,Nempty
                    casida_R((k-1)*Nempty+kk)=E_casida(Nstates+kk)-E_casida(k)
                    casida_R_half((k-1)*Nempty+kk)=sqrt(E_casida(Nstates+kk)-E_casida(k))
                enddo
            enddo
    
            do k=1,Nstates*Nempty
                casida_Kx(k,k)=casida_Kx(k,k)+casida_R(k)
            enddo

            allocate(temp(Nstates*Nempty,Nstates*Nempty))
            call outer_product('ge',casida_R_half,casida_R_half,temp)
    
            casida_Kx=casida_Kx*temp!outer_product(casida_R_half,casida_R_half)
    
    
        end subroutine compute_gw_Casida_matrix
    
        subroutine solve_Casida_matrix(Nstates,Nempty)
    
            integer,intent(in) :: Nstates,Nempty
        
            integer :: DSYEV_lwork,DSYEV_info,i,j
            double precision,dimension(:),allocatable:: DSYEV_work
            double precision,dimension(:,:),allocatable :: AA,BB,vl,vr,A,B
            ! double precision,dimension(:),allocatable :: w,beta,alphar,alphai

            integer,dimension(:),allocatable :: fpm
            integer :: info,M0,M,loop
            double precision :: epsout
            complex(kind=(kind(1.0d0))) :: Emid
            double precision :: r00
            complex(kind=(kind(1.0d0))),dimension(:),allocatable :: omega00
            complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: X00
            double precision,dimension(:),allocatable :: res00

            double precision,dimension(:),allocatable:: temp
            complex(kind=(kind(1.0d0))),dimension(:),allocatable :: ztemp

        
            !!!!! solve eigenvalue problem
            DSYEV_lwork = 3*Nstates*Nempty-1
            allocate(DSYEV_work(DSYEV_lwork))
            CALL DSYEV('V','U',Nstates*Nempty,casida_Kx,Nstates*Nempty,casida_omega,DSYEV_work,DSYEV_lwork,DSYEV_info)
        
            print *,casida_omega(1),casida_omega(Nstates*Nempty)
        
            casida_Xs=casida_Kx
        
            deallocate(DSYEV_work)





!             allocate(AA(2*Nstates*Nempty,2*Nstates*Nempty))
!             allocate(BB(2*Nstates*Nempty,2*Nstates*Nempty))

!             AA(1:Nstates*Nempty,1:Nstates*Nempty)=casida_Kx(1:Nstates*Nempty,1:Nstates*Nempty)
!             do i=1,Nstates*Nempty
!                 AA(i,i)=AA(i,i)+casida_R(i)
!             enddo
!             AA(Nstates*Nempty+1:2*Nstates*Nempty,Nstates*Nempty+1:2*Nstates*Nempty)=AA(1:Nstates*Nempty,1:Nstates*Nempty)
!             AA(1:Nstates*Nempty,Nstates*Nempty+1:2*Nstates*Nempty)=casida_Kx(1:Nstates*Nempty,1:Nstates*Nempty)
!             AA(Nstates*Nempty+1:2*Nstates*Nempty,1:Nstates*Nempty)=casida_Kx(1:Nstates*Nempty,1:Nstates*Nempty)

!             BB=0.0d0
!             do i=1,Nstates*Nempty
!                 BB(i,i)=1.0d0
!             enddo
!             do i=Nstates*Nempty+1,2*Nstates*Nempty
!                 BB(i,i)=-1.0d0
!             enddo


!             M0=2*Nstates*Nempty!2*Nstates*Nempty
!             allocate(omega00(M0),X00(2*Nstates*Nempty,2*M0),res00(2*M0))
!             Emid=cmplx(1.0d7,0.0d0,8)
!             r00=1.0d7

!             allocate(fpm(64))
!             allocate(temp(2*Nstates*Nempty))
!             allocate(ztemp(2*Nstates*Nempty))
        
!             call feastinit(fpm)
!             ! fpm(1)=1
!             fpm(15)=0!2!0!1
!     call dfeast_gegv(2*Nstates*Nempty,AA,2*Nstates*Nempty,BB,2*Nstates*Nempty,fpm,epsout,loop,Emid,r00,M0,omega00,X00,M,res00,info)

! ! call zfeast_sygv('U',2*Nstates*Nempty,AA*ZONE,2*Nstates*Nempty,BB*ZONE,2*Nstates*Nempty,fpm,epsout,loop,Emid,r00,M0,omega00,X00,M&
! ! ,res00,info)

!             print *,info,M

!     !         ! print *,X00(1:4,1)
!     !         ! print *,X00(1:4,M0+1)
!     !         ! do i=1,10
!     !         !     print *,omega00(i),X00(1,i),X00(M0+i)
!     !         ! enddo
!     !         print *,'======'
!     !         do i=1,10
!     !         print *,omega00(i),res00(i)!,res00(M0+i),res00(i)
!     !         enddo
!     !         print *,omega00(Nstates*Nempty)
!     !         print *,'======'
!     !         print *,X00(1:3,1)
!     !         print *,X00(1:3,M0+1)
!     !         print *,'======'

!     !         allocate(ztemp(4))

!     !         ztemp(1)=dot_product(X00(1:Nstates*Nempty,1),X00(1:Nstates*Nempty,1))
!     ! ztemp(1)=ztemp(1)-dot_product(X00(Nstates*Nempty+1:2*Nstates*Nempty,1),X00(Nstates*Nempty+1:2*Nstates*Nempty,1))

!     ! print *,ztemp(1)

!     ! ztemp(1)=dot_product(X00(1:Nstates*Nempty,1),X00(1:Nstates*Nempty,2))
!     ! ztemp(1)=ztemp(1)-dot_product(X00(Nstates*Nempty+1:2*Nstates*Nempty,1),X00(Nstates*Nempty+1:2*Nstates*Nempty,2))

!     ! print *,ztemp(1)


!     !         ! print *,'======'
!     !         ! print *,omega00(1:4)
!     !         ! print *,omega00(M0/2-3:M0/2)
!     !         ! print *,omega00(M0/2+1:M0/2+4)
!     !         ! print *,omega00(M0-3:M0)

!     !         stop

!     !         allocate(A(4,4),B(4,4))

!     !         A(1,:)=(/0.2d0,0.1d0,0.05d0,0.03d0/)
!     !         A(2,:)=(/0.1d0,0.2d0,0.03d0,0.05d0/)
!     !         A(3,:)=(/0.05d0,0.03d0,0.2d0,0.1d0/)
!     !         A(4,:)=(/0.03d0,0.05d0,0.1d0,0.2d0/)

!     !         B(1,:)=(/1.0d0,0.0d0,0.0d0,0.0d0/)
!     !         B(2,:)=(/0.0d0,1.0d0,0.0d0,0.0d0/)
!     !         B(3,:)=(/0.0d0,0.0d0,-1.0d0,0.0d0/)
!     !         B(4,:)=(/0.0d0,0.0d0,0.0d0,-1.0d0/)

!     !         M0=4
!     !         allocate(omega00(M0),X00(4,2*M0),res00(2*M0))
!     !         Emid=cmplx(0.0d0,0.0d0,8)
!     !         r00=1.0d0

!     !         allocate(fpm(64))
!     !         allocate(temp(4))
!     !         allocate(ztemp(4))
        
!     !         call feastinit(fpm)
!     !         ! fpm(1)=1
!     !         fpm(15)=0!1
!     !         call dfeast_gegv(4,A,4,B,4,fpm,epsout,loop,Emid,r00,M0,omega00,X00,M,res00,info)

!     !         print *,info,M

!     !         print *,X00(1:4,1)
!     !         print *,X00(1:4,M+1)

!     !         ztemp(1)=dot_product(X00(1:2,M0+1),X00(1:2,1))
!     ! ztemp(1)=ztemp(1)-dot_product(X00(2+1:2*2,M0+1),X00(2+1:2*2,1))

!     ! print *,ztemp(1)

!     !         call zsymv('U', 4, (1.0d0,0.0d0), B*(1.0d0,0.0d0), 4, X00(:,1), 1, (0.0d0,0.0d0), ztemp, 1)

!     !         print *,dot_product(X00(:,M0+1),ztemp)

!     !         temp(1)=dot_product(dble(X00(1:2,1)),dble(X00(1:2,2)))
!     ! temp(1)=temp(1)-dot_product(dble(X00(2+1:2*2,1)),dble(X00(2+1:2*2,2)))

!     !         print *,temp(1)

!     !         stop





            

!     !         temp(1)=dot_product(dble(X00(1:Nstates*Nempty,M0+1)),dble(X00(1:Nstates*Nempty,1)))
!     ! temp(1)=temp(1)-dot_product(dble(X00(Nstates*Nempty+1:2*Nstates*Nempty,M0+1)),dble(X00(Nstates*Nempty+1:2*Nstates*Nempty,1)))

!     !         print *,temp(1)

!     !         temp(1)=dot_product(dble(X00(1:Nstates*Nempty,M0+1)),dble(X00(1:Nstates*Nempty,2)))
!     ! temp(1)=temp(1)-dot_product(dble(X00(Nstates*Nempty+1:2*Nstates*Nempty,M0+1)),dble(X00(Nstates*Nempty+1:2*Nstates*Nempty,2)))

!     !         print *,temp(1)

!     !         stop





!         !     ztemp=0.0d0

!         !     do i=1,Nstates*Nempty
!         !         ztemp(i)=dot_product(X00(1:Nstates*Nempty,i),X00(1:Nstates*Nempty,i))
!         ! ztemp(i)=ztemp(i)-dot_product(X00(Nstates*Nempty+1:2*Nstates*Nempty,i),X00(Nstates*Nempty+1:2*Nstates*Nempty,i))
!         !     enddo

!         !     ! print *,temp(1)

!         !     do i=1,Nstates*Nempty
!         !         X00(:,i)=X00(:,i)/sqrt(ztemp(i))
!         !     enddo

!             ! casida_Xs=dble(X00(1:Nstates*Nempty,1:Nstates*Nempty))
!             ! casida_Ys=dble(X00(Nstates*Nempty+1:2*Nstates*Nempty,1:Nstates*Nempty))
!             casida_Xs=X00(1:Nstates*Nempty,1:Nstates*Nempty)
!             casida_Ys=X00(Nstates*Nempty+1:2*Nstates*Nempty,1:Nstates*Nempty)
!             casida_omega=dble(omega00(1:Nstates*Nempty))

!             deallocate(AA)
!             deallocate(BB)
!             deallocate(fpm)
!             deallocate(temp)

        
        end subroutine solve_Casida_matrix

    
        subroutine compute_gw_Casida_Kxnj(Nn,Ne,Nquadrature,Nstates,Nempty)
            integer,intent(in) :: Nn,Ne,Nquadrature,Nstates,Nempty
    
            !!!!! for empty states !!!!!
            CALL DGEMM('T','N',Nempty,Nstates*Nempty,Ne*Nquadrature,1.0d0,psi_n_g_real,Ne*Nquadrature&
                                                    ,NgNntemp,Ne*Nquadrature,0.0d0,casida_Kxnj,Nempty)
    
            !!!!! for occupied states !!!!!
            CALL DGEMM('T','N',Nstates,Nstates*Nempty,Ne*Nquadrature,1.0d0,psi_n_g_real_occupied,Ne*Nquadrature&
                                                    ,NgNntemp,Ne*Nquadrature,0.0d0,casida_Kxnj_occupied,Nstates)
    
    
        end subroutine compute_gw_Casida_Kxnj
    
        subroutine compute_gw_casida_Ec_complex(Eguess,E_casida,Nguess,Nstates,Nempty)
            integer,intent(in) :: Nguess,Nstates,Nempty
            complex(kind=(kind(1.0d0))),dimension(:),intent(in) :: Eguess
            double precision,dimension(:),intent(in) :: E_casida
            integer :: i,m,k,ii
            double precision, dimension(:), allocatable :: vector_temp,vector_temp0
            complex(kind=(kind(1.0d0))), dimension(:), allocatable :: vector_temp1,vector_temp2,vector_temp01,vector_temp02
    
            allocate(vector_temp(Nempty))
            allocate(vector_temp0(Nstates))
            allocate(vector_temp1(Nempty))
            allocate(vector_temp2(Nempty))
            allocate(vector_temp01(Nstates))
            allocate(vector_temp02(Nstates))
    
            casida_Ec=(0.0d0,0.0d0)
    
            do m=1,Nguess
    
            ! if (abs(Eguess(m)%im-0.0d0)<1.0d-10) then
    
                ii=0
                do i=1,Nstates*Nempty
    
                    ii=ii+1
    
                    !!!!!!!! for unoccupied states !!!!!!!!
                    call DGEMV('N',Nempty,Nstates*Nempty,1.0d0,casida_Kxnj,Nempty&
                            ,sqrt((1.0d0/sqrt(casida_omega(i))))*casida_R_half(:)*casida_Xs(:,i),1,0.0d0,vector_temp,1)
                    ! call DGEMV('N',Nempty,Nstates*Nempty,1.0d0,casida_Kxnj,Nempty&
                    !             ,(casida_Xs(:,i)+casida_Ys(:,i)),1,0.0d0,vector_temp,1)
                    ! call DGEMV('N',Nempty,Nstates*Nempty,1.0d0,casida_Kxnj,Nempty&
                    !         ,casida_Xs(:,i),1,0.0d0,vector_temp,1)
    
                    do k=1,Nempty
                        vector_temp1(k)=1.0d0/(Eguess(m)-E_casida(Nstates+k)-sqrt(casida_omega(i)))
                        ! vector_temp1(k)=1.0d0/(Eguess(m)-E_casida(Nstates+k)-casida_omega(i))
                    enddo
    
                    vector_temp2(1:Nempty)=(vector_temp(1:Nempty)*j_real)*vector_temp1(1:Nempty)
    
                    casida_Ec(m)=casida_Ec(m)+2.0d0*dot_product((vector_temp(1:Nempty)*j_real),vector_temp2(1:Nempty))
                    ! casida_Ec(m)=casida_Ec(m)+dot_product((vector_temp(1:Nempty)*j_real),vector_temp2(1:Nempty))
    
                    !!!!!!!! for occupied states !!!!!!!!
                    call DGEMV('N',Nstates,Nstates*Nempty,1.0d0,casida_Kxnj_occupied,Nstates&
                            ,sqrt((1.0d0/sqrt(casida_omega(i))))*casida_R_half(:)*casida_Xs(:,i),1,0.0d0,vector_temp0,1)
                    ! call DGEMV('N',Nstates,Nstates*Nempty,1.0d0,casida_Kxnj_occupied,Nstates&
                    !             ,(casida_Xs(:,i)+casida_Ys(:,i)),1,0.0d0,vector_temp0,1)
                    ! call DGEMV('N',Nstates,Nstates*Nempty,1.0d0,casida_Kxnj_occupied,Nstates&
                            ! ,casida_Xs(:,i),1,0.0d0,vector_temp0,1)
    
                    do k=1,Nstates
                        vector_temp01(k)=1.0d0/(Eguess(m)-E_casida(k)+sqrt(casida_omega(i)))
                        ! vector_temp01(k)=1.0d0/(Eguess(m)-E_casida(k)+casida_omega(i))
                    enddo
                    vector_temp02(1:Nstates)=(vector_temp0(1:Nstates)*j_real)*vector_temp01(1:Nstates)
    
                    casida_Ec(m)=casida_Ec(m)+2.0d0*dot_product((vector_temp0(1:Nstates)*j_real),vector_temp02(1:Nstates))
                    ! casida_Ec(m)=casida_Ec(m)+dot_product((vector_temp0(1:Nstates)*j_real),vector_temp02(1:Nstates))
                    
                enddo
    
            enddo
    
            deallocate(vector_temp)
            deallocate(vector_temp0)
            deallocate(vector_temp1)
            deallocate(vector_temp2)
            deallocate(vector_temp01)
            deallocate(vector_temp02)
    
        end subroutine compute_gw_casida_Ec_complex



        subroutine compute_gw_casida_Ec_complex_2(Eguess,E_casida,Nguess,Nstates,Nempty)
            integer,intent(in) :: Nguess,Nstates,Nempty
            complex(kind=(kind(1.0d0))),dimension(:),intent(in) :: Eguess
            double precision,dimension(:),intent(in) :: E_casida
            integer :: i,m,k,ii
            double precision, dimension(:), allocatable :: vector_temp,vector_temp0
            complex(kind=(kind(1.0d0))), dimension(:), allocatable :: vector_temp1,vector_temp2,vector_temp01,vector_temp02
            complex(kind=(kind(1.0d0))), dimension(:), allocatable :: zvector_temp,zvector_temp0
    
            allocate(vector_temp(Nempty))
            allocate(vector_temp0(Nstates))
            allocate(vector_temp1(Nempty))
            allocate(vector_temp2(Nempty))
            allocate(vector_temp01(Nstates))
            allocate(vector_temp02(Nstates))

            allocate(zvector_temp(Nempty))
            allocate(zvector_temp0(Nstates))
    
            casida_Ec=(0.0d0,0.0d0)
    
            do m=1,Nguess
    
            ! if (abs(Eguess(m)%im-0.0d0)<1.0d-10) then
    
                ii=0
                do i=1,Nstates*Nempty
    
                    ii=ii+1
    
                    !!!!!!!! for unoccupied states !!!!!!!!
                    ! call DGEMV('N',Nempty,Nstates*Nempty,1.0d0,casida_Kxnj,Nempty&
                    !         ,sqrt((1.0d0/sqrt(casida_omega(i))))*casida_R_half(:)*casida_Xs(:,i),1,0.0d0,vector_temp,1)
                    ! call DGEMV('N',Nempty,Nstates*Nempty,1.0d0,casida_Kxnj,Nempty&
                    !             ,(casida_Xs(:,i)+casida_Ys(:,i)),1,0.0d0,vector_temp,1)
                    ! call DGEMV('N',Nempty,Nstates*Nempty,1.0d0,casida_Kxnj,Nempty&
                    !         ,casida_Xs(:,i),1,0.0d0,vector_temp,1)
                    call ZGEMV('N',Nempty,Nstates*Nempty,ZONE,casida_Kxnj*j_real,Nempty&
                            ,(casida_Xs(:,i)+casida_Ys(:,i)),1,ZZERO,zvector_temp,1)
    
                    do k=1,Nempty
                        ! vector_temp1(k)=1.0d0/(Eguess(m)-E_casida(Nstates+k)-sqrt(casida_omega(i)))
                        vector_temp1(k)=1.0d0/(Eguess(m)-E_casida(Nstates+k)-casida_omega(i))
                    enddo
    
                    vector_temp2(1:Nempty)=zvector_temp(1:Nempty)*vector_temp1(1:Nempty)
    
                    ! casida_Ec(m)=casida_Ec(m)+2.0d0*dot_product((vector_temp(1:Nempty)*j_real),vector_temp2(1:Nempty))
                    casida_Ec(m)=casida_Ec(m)+dot_product(conjg(zvector_temp(1:Nempty)),vector_temp2(1:Nempty))
    
                    !!!!!!!! for occupied states !!!!!!!!
                    ! call DGEMV('N',Nstates,Nstates*Nempty,1.0d0,casida_Kxnj_occupied,Nstates&
                    !         ,sqrt((1.0d0/sqrt(casida_omega(i))))*casida_R_half(:)*casida_Xs(:,i),1,0.0d0,vector_temp0,1)
                    ! call DGEMV('N',Nstates,Nstates*Nempty,1.0d0,casida_Kxnj_occupied,Nstates&
                    !             ,(casida_Xs(:,i)+casida_Ys(:,i)),1,0.0d0,vector_temp0,1)
                    call ZGEMV('N',Nstates,Nstates*Nempty,ZONE,casida_Kxnj_occupied*j_real,Nstates&
                            ,(casida_Xs(:,i)+casida_Ys(:,i)),1,ZZERO,zvector_temp0,1)
    
                    do k=1,Nstates
                        ! vector_temp01(k)=1.0d0/(Eguess(m)-E_casida(k)+sqrt(casida_omega(i)))
                        vector_temp01(k)=1.0d0/(Eguess(m)-E_casida(k)+casida_omega(i))
                    enddo
                    vector_temp02(1:Nstates)=zvector_temp0(1:Nstates)*vector_temp01(1:Nstates)
    
                    ! casida_Ec(m)=casida_Ec(m)+2.0d0*dot_product((vector_temp0(1:Nstates)*j_real),vector_temp02(1:Nstates))
                    casida_Ec(m)=casida_Ec(m)+dot_product(conjg(zvector_temp0(1:Nstates)),vector_temp02(1:Nstates))
                    
                enddo
    
            enddo
    
            deallocate(vector_temp)
            deallocate(vector_temp0)
            deallocate(vector_temp1)
            deallocate(vector_temp2)
            deallocate(vector_temp01)
            deallocate(vector_temp02)

            deallocate(zvector_temp)
            deallocate(zvector_temp0)
    
        end subroutine compute_gw_casida_Ec_complex_2







































    subroutine gaussian_integral_1D()

        !!!! Gauss-Legendre rule https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_legendre/quadrature_rules_legendre.html

        ! double precision, dimension(1:11,1:3) :: gpo
        ! double precision, dimension(1:11) :: gweight

        
        allocate(gpoint3_1D  (1:3)) 
        allocate(gweight3_1D (1:3)) 
        allocate(gpoint20_1D (1:20)) 
        allocate(gweight20_1D(1:20)) 
        allocate(gpoint9_1D  (1:9)) 
        allocate(gweight9_1D (1:9)) 
        allocate(gpoint11_1D (1:11)) 
        allocate(gweight11_1D(1:11)) 
        allocate(gpoint7_1D  (1:7)) 
        allocate(gweight7_1D (1:7)) 
        allocate(gpoint15_1D (1:15)) 
        allocate(gweight15_1D(1:15)) 
        allocate(gpoint6_1D  (1:6)) 
        allocate(gweight6_1D (1:6)) 
        allocate(gpoint8_1D  (1:8)) 
        allocate(gweight8_1D (1:8)) 
        allocate(gpoint5_1D  (1:5)) 
        allocate(gweight5_1D (1:5)) 
        allocate(gpoint10_1D (1:10)) 
        allocate(gweight10_1D(1:10)) 
        allocate(gpoint16_1D (1:16)) 
        allocate(gweight16_1D(1:16)) 
        allocate(gpoint17_1D (1:17)) 
        allocate(gweight17_1D(1:17)) 
        allocate(gpoint18_1D (1:18)) 
        allocate(gweight18_1D(1:18)) 
        allocate(gpoint19_1D (1:19)) 
        allocate(gweight19_1D(1:19)) 
        allocate(gpoint24_1D (1:24)) 
        allocate(gweight24_1D(1:24)) 
        allocate(gpoint32_1D (1:32)) 
        allocate(gweight32_1D(1:32)) 
        allocate(gpoint40_1D (1:40)) 
        allocate(gweight40_1D(1:40)) 
        allocate(gpoint48_1D (1:48)) 
        allocate(gweight48_1D(1:48)) 

    
        gpoint3_1D(1)   = -sqrt(3.0d0/5.0d0)
        gpoint3_1D(2)   =  0.0d0
        gpoint3_1D(3)   =  sqrt(3.0d0/5.0d0)
    
    
        gweight3_1D(1)   =  5.0d0/9.0d0
        gweight3_1D(2)   =  8.0d0/9.0d0
        gweight3_1D(3)   =  5.0d0/9.0d0


        gweight20_1D(1 ) = 0.1527533871307258d0
        gweight20_1D(2 ) = 0.1527533871307258d0
        gweight20_1D(3 ) = 0.1491729864726037d0
        gweight20_1D(4 ) = 0.1491729864726037d0
        gweight20_1D(5 ) = 0.1420961093183820d0
        gweight20_1D(6 ) = 0.1420961093183820d0
        gweight20_1D(7 ) = 0.1316886384491766d0
        gweight20_1D(8 ) = 0.1316886384491766d0
        gweight20_1D(9 ) = 0.1181945319615184d0
        gweight20_1D(10) = 0.1181945319615184d0
        gweight20_1D(11) = 0.1019301198172404d0
        gweight20_1D(12) = 0.1019301198172404d0
        gweight20_1D(13) = 0.0832767415767048d0
        gweight20_1D(14) = 0.0832767415767048d0
        gweight20_1D(15) = 0.0626720483341091d0
        gweight20_1D(16) = 0.0626720483341091d0
        gweight20_1D(17) = 0.0406014298003869d0
        gweight20_1D(18) = 0.0406014298003869d0
        gweight20_1D(19) = 0.0176140071391521d0
        gweight20_1D(20) = 0.0176140071391521d0

        gpoint20_1D(1 ) = -0.0765265211334973d0
        gpoint20_1D(2 ) =  0.0765265211334973d0
        gpoint20_1D(3 ) = -0.2277858511416451d0
        gpoint20_1D(4 ) =  0.2277858511416451d0
        gpoint20_1D(5 ) = -0.3737060887154195d0
        gpoint20_1D(6 ) =  0.3737060887154195d0
        gpoint20_1D(7 ) = -0.5108670019508271d0
        gpoint20_1D(8 ) =  0.5108670019508271d0
        gpoint20_1D(9 ) = -0.6360536807265150d0
        gpoint20_1D(10) =  0.6360536807265150d0
        gpoint20_1D(11) = -0.7463319064601508d0
        gpoint20_1D(12) =  0.7463319064601508d0
        gpoint20_1D(13) = -0.8391169718222188d0
        gpoint20_1D(14) =  0.8391169718222188d0
        gpoint20_1D(15) = -0.9122344282513259d0
        gpoint20_1D(16) =  0.9122344282513259d0
        gpoint20_1D(17) = -0.9639719272779138d0
        gpoint20_1D(18) =  0.9639719272779138d0
        gpoint20_1D(19) = -0.9931285991850949d0
        gpoint20_1D(20) =  0.9931285991850949d0


        gweight9_1D(1) = 0.3302393550012598d0
        gweight9_1D(2) = 0.1806481606948574d0
        gweight9_1D(3) = 0.1806481606948574d0
        gweight9_1D(4) = 0.0812743883615744d0
        gweight9_1D(5) = 0.0812743883615744d0
        gweight9_1D(6) = 0.3123470770400029d0
        gweight9_1D(7) = 0.3123470770400029d0
        gweight9_1D(8) = 0.2606106964029354d0
        gweight9_1D(9) = 0.2606106964029354d0


        gpoint9_1D(1) = 0.0000000000000000d0
        gpoint9_1D(2) = -0.8360311073266358d0
        gpoint9_1D(3) = 0.8360311073266358d0
        gpoint9_1D(4) = -0.9681602395076261d0
        gpoint9_1D(5) = 0.9681602395076261d0
        gpoint9_1D(6) = -0.3242534234038089d0
        gpoint9_1D(7) = 0.3242534234038089d0
        gpoint9_1D(8) = -0.6133714327005904d0
        gpoint9_1D(9) = 0.6133714327005904d0


        gweight11_1D(1 ) = 0.272925086777901d0
        gweight11_1D(2 ) = 0.262804544510247d0
        gweight11_1D(3 ) = 0.262804544510247d0
        gweight11_1D(4 ) = 0.233193764591991d0
        gweight11_1D(5 ) = 0.233193764591991d0
        gweight11_1D(6 ) = 0.186290210927734d0
        gweight11_1D(7 ) = 0.186290210927734d0
        gweight11_1D(8 ) = 0.125580369464905d0
        gweight11_1D(9 ) = 0.125580369464905d0
        gweight11_1D(10) = 5.566856711617370d-2
        gweight11_1D(11) = 5.566856711617370d-2


        gpoint11_1D(1 ) =  0.000000000000000d0
        gpoint11_1D(2 ) = -0.269543155952345d0
        gpoint11_1D(3 ) =  0.269543155952345d0
        gpoint11_1D(4 ) = -0.519096129206812d0
        gpoint11_1D(5 ) =  0.519096129206812d0
        gpoint11_1D(6 ) = -0.730152005574049d0
        gpoint11_1D(7 ) =  0.730152005574049d0
        gpoint11_1D(8 ) = -0.887062599768095d0
        gpoint11_1D(9 ) =  0.887062599768095d0
        gpoint11_1D(10) = -0.978228658146057d0
        gpoint11_1D(11) =  0.978228658146057d0



        gweight7_1D(1) = 0.4179591836734694d0
        gweight7_1D(2) = 0.3818300505051189d0
        gweight7_1D(3) = 0.3818300505051189d0
        gweight7_1D(4) = 0.2797053914892766d0
        gweight7_1D(5) = 0.2797053914892766d0
        gweight7_1D(6) = 0.1294849661688697d0
        gweight7_1D(7) = 0.1294849661688697d0


        gpoint7_1D(1) = 0.0000000000000000d0
        gpoint7_1D(2) = 0.4058451513773972d0
        gpoint7_1D(3) = -0.4058451513773972d0
        gpoint7_1D(4) = -0.7415311855993945d0
        gpoint7_1D(5) = 0.7415311855993945d0
        gpoint7_1D(6) = -0.9491079123427585d0
        gpoint7_1D(7) = 0.9491079123427585d0



        gweight15_1D(1 ) =  0.202578241925561d0     
        gweight15_1D(2 ) =  0.198431485327112d0
        gweight15_1D(3 ) =  0.198431485327112d0   
        gweight15_1D(4 ) =  0.186161000015562d0   
        gweight15_1D(5 ) =  0.186161000015562d0
        gweight15_1D(6 ) =  0.166269205816994d0   
        gweight15_1D(7 ) =  0.166269205816994d0     
        gweight15_1D(8 ) =  0.139570677926154d0     
        gweight15_1D(9 ) =  0.139570677926154d0     
        gweight15_1D(10) =  0.107159220467172d0     
        gweight15_1D(11) =  0.107159220467172d0     
        gweight15_1D(12) =  7.036604748810810d-2
        gweight15_1D(13) =  7.036604748810810d-2
        gweight15_1D(14) =  3.075324199611730d-2
        gweight15_1D(15) =  3.075324199611730d-2

        gpoint15_1D(1 ) =  0.000000000000000d0   
        gpoint15_1D(2 ) = -0.201194093997434d0     
        gpoint15_1D(3 ) =  0.201194093997434d0      
        gpoint15_1D(4 ) = -0.394151347077563d0     
        gpoint15_1D(5 ) =  0.394151347077563d0       
        gpoint15_1D(6 ) = -0.570972172608539d0     
        gpoint15_1D(7 ) =  0.570972172608539d0     
        gpoint15_1D(8 ) = -0.724417731360170d0     
        gpoint15_1D(9 ) =  0.724417731360170d0     
        gpoint15_1D(10) = -0.848206583410427d0     
        gpoint15_1D(11) =  0.848206583410427d0     
        gpoint15_1D(12) = -0.937273392400706d0     
        gpoint15_1D(13) =  0.937273392400706d0     
        gpoint15_1D(14) = -0.987992518020485d0     
        gpoint15_1D(15) =  0.987992518020485d0

        
        gweight6_1D(1) = 3.60761573048138607569d-01
        gweight6_1D(2) = 4.67913934572691047389d-01
        gweight6_1D(3) = 1.71324492379170345043d-01
        gweight6_1D(4) = 3.60761573048138607569d-01
        gweight6_1D(5) = 4.67913934572691047389d-01
        gweight6_1D(6) = 1.71324492379170345043d-01

        gpoint6_1D(1) = -6.61209386466264513688d-01
        gpoint6_1D(2) = -2.38619186083196908630d-01
        gpoint6_1D(3) = -9.32469514203152027832d-01
        gpoint6_1D(4) =  6.61209386466264513688d-01
        gpoint6_1D(5) =  2.38619186083196908630d-01
        gpoint6_1D(6) =  9.32469514203152027832d-01

        
        gweight8_1D(1) = 3.62683783378361982976d-01
        gweight8_1D(2) = 3.13706645877887287338d-01
        gweight8_1D(3) = 2.22381034453374470546d-01
        gweight8_1D(4) = 1.01228536290376259154d-01
        gweight8_1D(5) = 3.62683783378361982976d-01
        gweight8_1D(6) = 3.13706645877887287338d-01
        gweight8_1D(7) = 2.22381034453374470546d-01
        gweight8_1D(8) = 1.01228536290376259154d-01

        gpoint8_1D(1) = -1.83434642495649804936d-01
        gpoint8_1D(2) = -5.25532409916328985830d-01
        gpoint8_1D(3) = -7.96666477413626739567d-01
        gpoint8_1D(4) = -9.60289856497536231661d-01
        gpoint8_1D(5) =  1.83434642495649804936d-01
        gpoint8_1D(6) =  5.25532409916328985830d-01
        gpoint8_1D(7) =  7.96666477413626739567d-01
        gpoint8_1D(8) =  9.60289856497536231661d-01


        gweight5_1D(1) = 128.0d0/225.0d0
        gweight5_1D(2) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
        gweight5_1D(3) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
        gweight5_1D(4) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
        gweight5_1D(5) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
        
        gpoint5_1D(1) =  0.0d0
        gpoint5_1D(2) = -(1.0d0/3.0d0)*sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))
        gpoint5_1D(3) = -(1.0d0/3.0d0)*sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))
        gpoint5_1D(4) =  (1.0d0/3.0d0)*sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))
        gpoint5_1D(5) =  (1.0d0/3.0d0)*sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))


        gweight10_1D(1 ) = 2.95524224714752870187d-01
        gweight10_1D(2 ) = 2.69266719309996355105d-01
        gweight10_1D(3 ) = 2.19086362515982044000d-01
        gweight10_1D(4 ) = 1.49451349150580593150d-01
        gweight10_1D(5 ) = 6.66713443086881375920d-02
        gweight10_1D(6 ) = 2.95524224714752870187d-01
        gweight10_1D(7 ) = 2.69266719309996355105d-01
        gweight10_1D(8 ) = 2.19086362515982044000d-01
        gweight10_1D(9 ) = 1.49451349150580593150d-01
        gweight10_1D(10) = 6.66713443086881375920d-02

        gpoint10_1D(1 ) = -1.48874338981631210881d-01
        gpoint10_1D(2 ) = -4.33395394129247190794d-01
        gpoint10_1D(3 ) = -6.79409568299024406207d-01
        gpoint10_1D(4 ) = -8.65063366688984510759d-01
        gpoint10_1D(5 ) = -9.73906528517171720066d-01
        gpoint10_1D(6 ) =  1.48874338981631210881d-01
        gpoint10_1D(7 ) =  4.33395394129247190794d-01
        gpoint10_1D(8 ) =  6.79409568299024406207d-01
        gpoint10_1D(9 ) =  8.65063366688984510759d-01
        gpoint10_1D(10) =  9.73906528517171720066d-01


        gweight16_1D(1 ) = 1.89450610455068496287d-01
        gweight16_1D(2 ) = 1.82603415044923588872d-01
        gweight16_1D(3 ) = 1.69156519395002538183d-01
        gweight16_1D(4 ) = 1.49595988816576732080d-01
        gweight16_1D(5 ) = 1.24628971255533872056d-01
        gweight16_1D(6 ) = 9.51585116824927848073d-02
        gweight16_1D(7 ) = 6.22535239386478928628d-02
        gweight16_1D(8 ) = 2.71524594117540948514d-02
        gweight16_1D(9 ) = 1.89450610455068496287d-01
        gweight16_1D(10) = 1.82603415044923588872d-01
        gweight16_1D(11) = 1.69156519395002538183d-01
        gweight16_1D(12) = 1.49595988816576732080d-01
        gweight16_1D(13) = 1.24628971255533872056d-01
        gweight16_1D(14) = 9.51585116824927848073d-02
        gweight16_1D(15) = 6.22535239386478928628d-02
        gweight16_1D(16) = 2.71524594117540948514d-02
        
        gpoint16_1D(1 ) = -9.50125098376374401877d-02
        gpoint16_1D(2 ) = -2.81603550779258913231d-01
        gpoint16_1D(3 ) = -4.58016777657227386350d-01
        gpoint16_1D(4 ) = -6.17876244402643748452d-01
        gpoint16_1D(5 ) = -7.55404408355003033891d-01
        gpoint16_1D(6 ) = -8.65631202387831743866d-01
        gpoint16_1D(7 ) = -9.44575023073232576090d-01
        gpoint16_1D(8 ) = -9.89400934991649932601d-01
        gpoint16_1D(9 ) =  9.50125098376374401877d-02
        gpoint16_1D(10) =  2.81603550779258913231d-01
        gpoint16_1D(11) =  4.58016777657227386350d-01
        gpoint16_1D(12) =  6.17876244402643748452d-01
        gpoint16_1D(13) =  7.55404408355003033891d-01
        gpoint16_1D(14) =  8.65631202387831743866d-01
        gpoint16_1D(15) =  9.44575023073232576090d-01
        gpoint16_1D(16) =  9.89400934991649932601d-01



        gweight17_1D(1 ) = 0.179446470356207d0     
        gweight17_1D(2 ) = 0.176562705366993d0     
        gweight17_1D(3 ) = 0.176562705366993d0     
        gweight17_1D(4 ) = 0.168004102156450d0     
        gweight17_1D(5 ) = 0.168004102156450d0     
        gweight17_1D(6 ) = 0.154045761076810d0     
        gweight17_1D(7 ) = 0.154045761076810d0     
        gweight17_1D(8 ) = 0.135136368468526d0     
        gweight17_1D(9 ) = 0.135136368468526d0     
        gweight17_1D(10) = 0.111883847193404d0     
        gweight17_1D(11) = 0.111883847193404d0     
        gweight17_1D(12) = 8.503614831717921d-2
        gweight17_1D(13) = 8.503614831717921d-2
        gweight17_1D(14) = 5.545952937398720d-2
        gweight17_1D(15) = 5.545952937398720d-2
        gweight17_1D(16) = 2.414830286854790d-2
        gweight17_1D(17) = 2.414830286854790d-2

        gpoint17_1D(1 ) =  0.000000000000000d0
        gpoint17_1D(2 ) = -0.178484181495848d0 
        gpoint17_1D(3 ) =  0.178484181495848d0  
        gpoint17_1D(4 ) = -0.351231763453876d0   
        gpoint17_1D(5 ) =  0.351231763453876d0    
        gpoint17_1D(6 ) = -0.512690537086477d0    
        gpoint17_1D(7 ) =  0.512690537086477d0     
        gpoint17_1D(8 ) = -0.657671159216691d0    
        gpoint17_1D(9 ) =  0.657671159216691d0     
        gpoint17_1D(10) = -0.781514003896801d0    
        gpoint17_1D(11) =  0.781514003896801d0     
        gpoint17_1D(12) = -0.880239153726986d0     
        gpoint17_1D(13) =  0.880239153726986d0     
        gpoint17_1D(14) = -0.950675521768768d0     
        gpoint17_1D(15) =  0.950675521768768d0     
        gpoint17_1D(16) = -0.990575475314417d0     
        gpoint17_1D(17) =  0.990575475314417d0



        gweight18_1D(1 ) = 0.169142382963144d0     
        gweight18_1D(2 ) = 0.169142382963144d0     
        gweight18_1D(3 ) = 0.164276483745833d0    
        gweight18_1D(4 ) = 0.164276483745833d0   
        gweight18_1D(5 ) = 0.154684675126265d0   
        gweight18_1D(6 ) = 0.154684675126265d0   
        gweight18_1D(7 ) = 0.140642914670651d0   
        gweight18_1D(8 ) = 0.140642914670651d0   
        gweight18_1D(9 ) = 0.122555206711479d0   
        gweight18_1D(10) = 0.122555206711479d0   
        gweight18_1D(11) = 0.100942044106287d0   
        gweight18_1D(12) = 0.100942044106287d0    
        gweight18_1D(13) = 7.642573025488909d-2
        gweight18_1D(14) = 7.642573025488909d-2
        gweight18_1D(15) = 4.971454889496980d-2
        gweight18_1D(16) = 4.971454889496980d-2
        gweight18_1D(17) = 2.161601352648330d-2
        gweight18_1D(18) = 2.161601352648330d-2



        gpoint18_1D(1 ) = -8.477501304173531d-2
        gpoint18_1D(2 ) =  8.477501304173531d-2
        gpoint18_1D(3 ) = -0.251886225691505d0    
        gpoint18_1D(4 ) =  0.251886225691505d0   
        gpoint18_1D(5 ) = -0.411751161462843d0   
        gpoint18_1D(6 ) =  0.411751161462843d0   
        gpoint18_1D(7 ) = -0.559770831073948d0   
        gpoint18_1D(8 ) =  0.559770831073948d0   
        gpoint18_1D(9 ) = -0.691687043060353d0   
        gpoint18_1D(10) =  0.691687043060353d0   
        gpoint18_1D(11) = -0.803704958972523d0   
        gpoint18_1D(12) =  0.803704958972523d0     
        gpoint18_1D(13) = -0.892602466497556d0     
        gpoint18_1D(14) =  0.892602466497556d0     
        gpoint18_1D(15) = -0.955823949571398d0     
        gpoint18_1D(16) =  0.955823949571398d0     
        gpoint18_1D(17) = -0.991565168420931d0     
        gpoint18_1D(18) =  0.991565168420931d0



        gweight19_1D(1 ) = 0.161054449848784d0     
        gweight19_1D(2 ) = 0.158968843393954d0     
        gweight19_1D(3 ) = 0.158968843393954d0   
        gweight19_1D(4 ) = 0.152766042065860d0   
        gweight19_1D(5 ) = 0.152766042065860d0   
        gweight19_1D(6 ) = 0.142606702173607d0   
        gweight19_1D(7 ) = 0.142606702173607d0   
        gweight19_1D(8 ) = 0.128753962539336d0   
        gweight19_1D(9 ) = 0.128753962539336d0   
        gweight19_1D(10) = 0.111566645547334d0   
        gweight19_1D(11) = 0.111566645547334d0     
        gweight19_1D(12) = 9.149002162245000d-2
        gweight19_1D(13) = 9.149002162245000d-2
        gweight19_1D(14) = 6.904454273764120d-2
        gweight19_1D(15) = 6.904454273764120d-2
        gweight19_1D(16) = 4.481422676569960d-2
        gweight19_1D(17) = 4.481422676569960d-2
        gweight19_1D(18) = 1.946178822972650d-2
        gweight19_1D(19) = 1.946178822972650d-2


        gpoint19_1D(1 ) =  0.000000000000000d0
        gpoint19_1D(2 ) = -0.160358645640225d0   
        gpoint19_1D(3 ) =  0.160358645640225d0   
        gpoint19_1D(4 ) = -0.316564099963630d0   
        gpoint19_1D(5 ) =  0.316564099963630d0   
        gpoint19_1D(6 ) = -0.464570741375961d0   
        gpoint19_1D(7 ) =  0.464570741375961d0   
        gpoint19_1D(8 ) = -0.600545304661681d0   
        gpoint19_1D(9 ) =  0.600545304661681d0   
        gpoint19_1D(10) = -0.720966177335229d0    
        gpoint19_1D(11) =  0.720966177335229d0     
        gpoint19_1D(12) = -0.822714656537143d0     
        gpoint19_1D(13) =  0.822714656537143d0     
        gpoint19_1D(14) = -0.903155903614818d0     
        gpoint19_1D(15) =  0.903155903614818d0     
        gpoint19_1D(16) = -0.960208152134830d0     
        gpoint19_1D(17) =  0.960208152134830d0     
        gpoint19_1D(18) = -0.992406843843584d0     
        gpoint19_1D(19) =  0.992406843843584d0


        gweight24_1D(1 ) = 1.27938195346752156976d-01
        gweight24_1D(2 ) = 1.25837456346828296117d-01
        gweight24_1D(3 ) = 1.21670472927803391202d-01
        gweight24_1D(4 ) = 1.15505668053725601353d-01
        gweight24_1D(5 ) = 1.07444270115965634785d-01
        gweight24_1D(6 ) = 9.76186521041138882720d-02
        gweight24_1D(7 ) = 8.61901615319532759152d-02
        gweight24_1D(8 ) = 7.33464814110803057346d-02
        gweight24_1D(9 ) = 5.92985849154367807461d-02
        gweight24_1D(10) = 4.42774388174198061695d-02
        gweight24_1D(11) = 2.85313886289336631809d-02
        gweight24_1D(12) = 1.23412297999871995469d-02
        gweight24_1D(13) = 1.27938195346752156976d-01
        gweight24_1D(14) = 1.25837456346828296117d-01
        gweight24_1D(15) = 1.21670472927803391202d-01
        gweight24_1D(16) = 1.15505668053725601353d-01
        gweight24_1D(17) = 1.07444270115965634785d-01
        gweight24_1D(18) = 9.76186521041138882720d-02
        gweight24_1D(19) = 8.61901615319532759152d-02
        gweight24_1D(20) = 7.33464814110803057346d-02
        gweight24_1D(21) = 5.92985849154367807461d-02
        gweight24_1D(22) = 4.42774388174198061695d-02
        gweight24_1D(23) = 2.85313886289336631809d-02
        gweight24_1D(24) = 1.23412297999871995469d-02

        gpoint24_1D(1 ) = -6.40568928626056260827d-02
        gpoint24_1D(2 ) = -1.91118867473616309153d-01
        gpoint24_1D(3 ) = -3.15042679696163374398d-01
        gpoint24_1D(4 ) = -4.33793507626045138478d-01
        gpoint24_1D(5 ) = -5.45421471388839535649d-01
        gpoint24_1D(6 ) = -6.48093651936975569268d-01
        gpoint24_1D(7 ) = -7.40124191578554364260d-01
        gpoint24_1D(8 ) = -8.20001985973902921981d-01
        gpoint24_1D(9 ) = -8.86415527004401034190d-01
        gpoint24_1D(10) = -9.38274552002732758539d-01
        gpoint24_1D(11) = -9.74728555971309498199d-01
        gpoint24_1D(12) = -9.95187219997021360195d-01
        gpoint24_1D(13) =  6.40568928626056260827d-02
        gpoint24_1D(14) =  1.91118867473616309153d-01
        gpoint24_1D(15) =  3.15042679696163374398d-01
        gpoint24_1D(16) =  4.33793507626045138478d-01
        gpoint24_1D(17) =  5.45421471388839535649d-01
        gpoint24_1D(18) =  6.48093651936975569268d-01
        gpoint24_1D(19) =  7.40124191578554364260d-01
        gpoint24_1D(20) =  8.20001985973902921981d-01
        gpoint24_1D(21) =  8.86415527004401034190d-01
        gpoint24_1D(22) =  9.38274552002732758539d-01
        gpoint24_1D(23) =  9.74728555971309498199d-01
        gpoint24_1D(24) =  9.95187219997021360195d-01




        gweight32_1D(1 ) = 9.65400885147278005666d-02
        gweight32_1D(2 ) = 9.56387200792748594185d-02
        gweight32_1D(3 ) = 9.38443990808045656367d-02
        gweight32_1D(4 ) = 9.11738786957638847129d-02
        gweight32_1D(5 ) = 8.76520930044038111450d-02
        gweight32_1D(6 ) = 8.33119242269467552223d-02
        gweight32_1D(7 ) = 7.81938957870703064685d-02
        gweight32_1D(8 ) = 7.23457941088485062287d-02
        gweight32_1D(9 ) = 6.58222227763618468406d-02
        gweight32_1D(10) = 5.86840934785355471448d-02
        gweight32_1D(11) = 5.09980592623761761959d-02
        gweight32_1D(12) = 4.28358980222266806557d-02
        gweight32_1D(13) = 3.42738629130214331033d-02
        gweight32_1D(14) = 2.53920653092620594561d-02
        gweight32_1D(15) = 1.62743947309056706058d-02
        gweight32_1D(16) = 7.01861000947009660028d-03
        gweight32_1D(17) = 9.65400885147278005666d-02
        gweight32_1D(18) = 9.56387200792748594185d-02
        gweight32_1D(19) = 9.38443990808045656367d-02
        gweight32_1D(20) = 9.11738786957638847129d-02
        gweight32_1D(21) = 8.76520930044038111450d-02
        gweight32_1D(22) = 8.33119242269467552223d-02
        gweight32_1D(23) = 7.81938957870703064685d-02
        gweight32_1D(24) = 7.23457941088485062287d-02
        gweight32_1D(25) = 6.58222227763618468406d-02
        gweight32_1D(26) = 5.86840934785355471448d-02
        gweight32_1D(27) = 5.09980592623761761959d-02
        gweight32_1D(28) = 4.28358980222266806557d-02
        gweight32_1D(29) = 3.42738629130214331033d-02
        gweight32_1D(30) = 2.53920653092620594561d-02
        gweight32_1D(31) = 1.62743947309056706058d-02
        gweight32_1D(32) = 7.01861000947009660028d-03



        gpoint32_1D(1 ) = -4.83076656877383162364d-02
        gpoint32_1D(2 ) = -1.44471961582796493484d-01
        gpoint32_1D(3 ) = -2.39287362252137074544d-01
        gpoint32_1D(4 ) = -3.31868602282127649782d-01
        gpoint32_1D(5 ) = -4.21351276130635345353d-01
        gpoint32_1D(6 ) = -5.06899908932229390044d-01
        gpoint32_1D(7 ) = -5.87715757240762329066d-01
        gpoint32_1D(8 ) = -6.63044266930215200960d-01
        gpoint32_1D(9 ) = -7.32182118740289680412d-01
        gpoint32_1D(10) = -7.94483795967942406965d-01
        gpoint32_1D(11) = -8.49367613732569970160d-01
        gpoint32_1D(12) = -8.96321155766052123971d-01
        gpoint32_1D(13) = -9.34906075937739689159d-01
        gpoint32_1D(14) = -9.64762255587506430761d-01
        gpoint32_1D(15) = -9.85611511545268335400d-01
        gpoint32_1D(16) = -9.97263861849481563534d-01
        gpoint32_1D(17) =  4.83076656877383162364d-02
        gpoint32_1D(18) =  1.44471961582796493484d-01
        gpoint32_1D(19) =  2.39287362252137074544d-01
        gpoint32_1D(20) =  3.31868602282127649782d-01
        gpoint32_1D(21) =  4.21351276130635345353d-01
        gpoint32_1D(22) =  5.06899908932229390044d-01
        gpoint32_1D(23) =  5.87715757240762329066d-01
        gpoint32_1D(24) =  6.63044266930215200960d-01
        gpoint32_1D(25) =  7.32182118740289680412d-01
        gpoint32_1D(26) =  7.94483795967942406965d-01
        gpoint32_1D(27) =  8.49367613732569970160d-01
        gpoint32_1D(28) =  8.96321155766052123971d-01
        gpoint32_1D(29) =  9.34906075937739689159d-01
        gpoint32_1D(30) =  9.64762255587506430761d-01
        gpoint32_1D(31) =  9.85611511545268335400d-01
        gpoint32_1D(32) =  9.97263861849481563534d-01



        gweight40_1D(1 ) = 7.75059479784248112668d-02
        gweight40_1D(2 ) = 7.70398181642479655914d-02
        gweight40_1D(3 ) = 7.61103619006262423723d-02
        gweight40_1D(4 ) = 7.47231690579682641980d-02
        gweight40_1D(5 ) = 7.28865823958040590609d-02
        gweight40_1D(6 ) = 7.06116473912867796979d-02
        gweight40_1D(7 ) = 6.79120458152339038265d-02
        gweight40_1D(8 ) = 6.48040134566010380719d-02
        gweight40_1D(9 ) = 6.13062424929289391679d-02
        gweight40_1D(10) = 5.74397690993915513665d-02
        gweight40_1D(11) = 5.32278469839368243566d-02
        gweight40_1D(12) = 4.86958076350722320604d-02
        gweight40_1D(13) = 4.38709081856732719923d-02
        gweight40_1D(14) = 3.87821679744720176413d-02
        gweight40_1D(15) = 3.34601952825478473933d-02
        gweight40_1D(16) = 2.79370069800234010984d-02
        gweight40_1D(17) = 2.22458491941669572615d-02
        gweight40_1D(18) = 1.64210583819078887131d-02
        gweight40_1D(19) = 1.04982845311528136146d-02
        gweight40_1D(20) = 4.52127709853319125846d-03
        gweight40_1D(21) = 7.75059479784248112668d-02
        gweight40_1D(22) = 7.70398181642479655914d-02
        gweight40_1D(23) = 7.61103619006262423723d-02
        gweight40_1D(24) = 7.47231690579682641980d-02
        gweight40_1D(25) = 7.28865823958040590609d-02
        gweight40_1D(26) = 7.06116473912867796979d-02
        gweight40_1D(27) = 6.79120458152339038265d-02
        gweight40_1D(28) = 6.48040134566010380719d-02
        gweight40_1D(29) = 6.13062424929289391679d-02
        gweight40_1D(30) = 5.74397690993915513665d-02
        gweight40_1D(31) = 5.32278469839368243566d-02
        gweight40_1D(32) = 4.86958076350722320604d-02
        gweight40_1D(33) = 4.38709081856732719923d-02
        gweight40_1D(34) = 3.87821679744720176413d-02
        gweight40_1D(35) = 3.34601952825478473933d-02
        gweight40_1D(36) = 2.79370069800234010984d-02
        gweight40_1D(37) = 2.22458491941669572615d-02
        gweight40_1D(38) = 1.64210583819078887131d-02
        gweight40_1D(39) = 1.04982845311528136146d-02
        gweight40_1D(40) = 4.52127709853319125846d-03



        gpoint40_1D(1 ) = -3.87724175060508219329d-02
        gpoint40_1D(2 ) = -1.16084070675255208481d-01
        gpoint40_1D(3 ) = -1.92697580701371099719d-01
        gpoint40_1D(4 ) = -2.68152185007253681152d-01
        gpoint40_1D(5 ) = -3.41994090825758473008d-01
        gpoint40_1D(6 ) = -4.13779204371605001525d-01
        gpoint40_1D(7 ) = -4.83075801686178712903d-01
        gpoint40_1D(8 ) = -5.49467125095128202056d-01
        gpoint40_1D(9 ) = -6.12553889667980237972d-01
        gpoint40_1D(10) = -6.71956684614179548364d-01
        gpoint40_1D(11) = -7.27318255189927103277d-01
        gpoint40_1D(12) = -7.78305651426519387712d-01
        gpoint40_1D(13) = -8.24612230833311663197d-01
        gpoint40_1D(14) = -8.65959503212259503824d-01
        gpoint40_1D(15) = -9.02098806968874296732d-01
        gpoint40_1D(16) = -9.32812808278676533383d-01
        gpoint40_1D(17) = -9.57916819213791655824d-01
        gpoint40_1D(18) = -9.77259949983774262679d-01
        gpoint40_1D(19) = -9.90726238699457006464d-01
        gpoint40_1D(20) = -9.98237709710559200369d-01
        gpoint40_1D(21) =  3.87724175060508219329d-02
        gpoint40_1D(22) =  1.16084070675255208481d-01
        gpoint40_1D(23) =  1.92697580701371099719d-01
        gpoint40_1D(24) =  2.68152185007253681152d-01
        gpoint40_1D(25) =  3.41994090825758473008d-01
        gpoint40_1D(26) =  4.13779204371605001525d-01
        gpoint40_1D(27) =  4.83075801686178712903d-01
        gpoint40_1D(28) =  5.49467125095128202056d-01
        gpoint40_1D(29) =  6.12553889667980237972d-01
        gpoint40_1D(30) =  6.71956684614179548364d-01
        gpoint40_1D(31) =  7.27318255189927103277d-01
        gpoint40_1D(32) =  7.78305651426519387712d-01
        gpoint40_1D(33) =  8.24612230833311663197d-01
        gpoint40_1D(34) =  8.65959503212259503824d-01
        gpoint40_1D(35) =  9.02098806968874296732d-01
        gpoint40_1D(36) =  9.32812808278676533383d-01
        gpoint40_1D(37) =  9.57916819213791655824d-01
        gpoint40_1D(38) =  9.77259949983774262679d-01
        gpoint40_1D(39) =  9.90726238699457006464d-01
        gpoint40_1D(40) =  9.98237709710559200369d-01



        gweight48_1D(1 ) = 6.47376968126839225006d-02
        gweight48_1D(2 ) = 6.44661644359500822082d-02
        gweight48_1D(3 ) = 6.39242385846481866207d-02
        gweight48_1D(4 ) = 6.31141922862540256548d-02
        gweight48_1D(5 ) = 6.20394231598926639029d-02
        gweight48_1D(6 ) = 6.07044391658938800517d-02
        gweight48_1D(7 ) = 5.91148396983956357477d-02
        gweight48_1D(8 ) = 5.72772921004032157044d-02
        gweight48_1D(9 ) = 5.51995036999841628676d-02
        gweight48_1D(10) = 5.28901894851936670964d-02
        gweight48_1D(11) = 5.03590355538544749590d-02
        gweight48_1D(12) = 4.76166584924904748267d-02
        gweight48_1D(13) = 4.46745608566942804201d-02
        gweight48_1D(14) = 4.15450829434647492133d-02
        gweight48_1D(15) = 3.82413510658307063158d-02
        gweight48_1D(16) = 3.47772225647704388909d-02
        gweight48_1D(17) = 3.11672278327980889025d-02
        gweight48_1D(18) = 2.74265097083569482001d-02
        gweight48_1D(19) = 2.35707608393243791410d-02
        gweight48_1D(20) = 1.96161604573555278142d-02
        gweight48_1D(21) = 1.55793157229438487279d-02
        gweight48_1D(22) = 1.14772345792345394895d-02
        gweight48_1D(23) = 7.32755390127626210220d-03
        gweight48_1D(24) = 3.15334605230583863260d-03
        gweight48_1D(25) = 6.47376968126839225006d-02
        gweight48_1D(26) = 6.44661644359500822082d-02
        gweight48_1D(27) = 6.39242385846481866207d-02
        gweight48_1D(28) = 6.31141922862540256548d-02
        gweight48_1D(29) = 6.20394231598926639029d-02
        gweight48_1D(30) = 6.07044391658938800517d-02
        gweight48_1D(31) = 5.91148396983956357477d-02
        gweight48_1D(32) = 5.72772921004032157044d-02
        gweight48_1D(33) = 5.51995036999841628676d-02
        gweight48_1D(34) = 5.28901894851936670964d-02
        gweight48_1D(35) = 5.03590355538544749590d-02
        gweight48_1D(36) = 4.76166584924904748267d-02
        gweight48_1D(37) = 4.46745608566942804201d-02
        gweight48_1D(38) = 4.15450829434647492133d-02
        gweight48_1D(39) = 3.82413510658307063158d-02
        gweight48_1D(40) = 3.47772225647704388909d-02
        gweight48_1D(41) = 3.11672278327980889025d-02
        gweight48_1D(42) = 2.74265097083569482001d-02
        gweight48_1D(43) = 2.35707608393243791410d-02
        gweight48_1D(44) = 1.96161604573555278142d-02
        gweight48_1D(45) = 1.55793157229438487279d-02
        gweight48_1D(46) = 1.14772345792345394895d-02
        gweight48_1D(47) = 7.32755390127626210220d-03
        gweight48_1D(48) = 3.15334605230583863260d-03



        gpoint48_1D(1 ) = -3.23801709628693620343d-02
        gpoint48_1D(2 ) = -9.70046992094626989322d-02
        gpoint48_1D(3 ) = -1.61222356068891718055d-01
        gpoint48_1D(4 ) = -2.24763790394689061224d-01
        gpoint48_1D(5 ) = -2.87362487355455576728d-01
        gpoint48_1D(6 ) = -3.48755886292160738148d-01
        gpoint48_1D(7 ) = -4.08686481990716729925d-01
        gpoint48_1D(8 ) = -4.66902904750958404535d-01
        gpoint48_1D(9 ) = -5.23160974722233033658d-01
        gpoint48_1D(10) = -5.77224726083972703838d-01
        gpoint48_1D(11) = -6.28867396776513624013d-01
        gpoint48_1D(12) = -6.77872379632663905208d-01
        gpoint48_1D(13) = -7.24034130923814654658d-01
        gpoint48_1D(14) = -7.67159032515740339276d-01
        gpoint48_1D(15) = -8.07066204029442627087d-01
        gpoint48_1D(16) = -8.43588261624393530704d-01
        gpoint48_1D(17) = -8.76572020274247885885d-01
        gpoint48_1D(18) = -9.05879136715569672805d-01
        gpoint48_1D(19) = -9.31386690706554333107d-01
        gpoint48_1D(20) = -9.52987703160430860724d-01
        gpoint48_1D(21) = -9.70591592546247250472d-01
        gpoint48_1D(22) = -9.84124583722826857765d-01
        gpoint48_1D(23) = -9.93530172266350757526d-01
        gpoint48_1D(24) = -9.98771007252426118580d-01
        gpoint48_1D(25) =  3.23801709628693620343d-02
        gpoint48_1D(26) =  9.70046992094626989322d-02
        gpoint48_1D(27) =  1.61222356068891718055d-01
        gpoint48_1D(28) =  2.24763790394689061224d-01
        gpoint48_1D(29) =  2.87362487355455576728d-01
        gpoint48_1D(30) =  3.48755886292160738148d-01
        gpoint48_1D(31) =  4.08686481990716729925d-01
        gpoint48_1D(32) =  4.66902904750958404535d-01
        gpoint48_1D(33) =  5.23160974722233033658d-01
        gpoint48_1D(34) =  5.77224726083972703838d-01
        gpoint48_1D(35) =  6.28867396776513624013d-01
        gpoint48_1D(36) =  6.77872379632663905208d-01
        gpoint48_1D(37) =  7.24034130923814654658d-01
        gpoint48_1D(38) =  7.67159032515740339276d-01
        gpoint48_1D(39) =  8.07066204029442627087d-01
        gpoint48_1D(40) =  8.43588261624393530704d-01
        gpoint48_1D(41) =  8.76572020274247885885d-01
        gpoint48_1D(42) =  9.05879136715569672805d-01
        gpoint48_1D(43) =  9.31386690706554333107d-01
        gpoint48_1D(44) =  9.52987703160430860724d-01
        gpoint48_1D(45) =  9.70591592546247250472d-01
        gpoint48_1D(46) =  9.84124583722826857765d-01
        gpoint48_1D(47) =  9.93530172266350757526d-01
        gpoint48_1D(48) =  9.98771007252426118580d-01

    
    
    end subroutine gaussian_integral_1D





    subroutine compute_cond_number(A, n, cond_num)
        implicit none
        integer, intent(in) :: n   ! Dimension of the square matrix A (n x n)
        double precision, intent(in) :: A(n, n)  ! Input matrix
        double precision, intent(out) :: cond_num  ! Condition number result
        double precision :: s(n), superb(n-1)  ! Singular values and auxiliary array for SVD
        double precision, allocatable :: u(:,:), vt(:,:), work(:)
        integer :: info, lwork
        
        ! Allocate arrays for SVD computation
        allocate(u(n,n), vt(n,n))
        
        ! Workspace query to determine optimal work array size
        lwork = -1
        call dgesvd('A', 'A', n, n, A, n, s, u, n, vt, n, work, lwork, info)
        
        ! Allocate the work array with the recommended size
        lwork = int(work(1))
        allocate(work(lwork))
        
        ! Compute the SVD of A (full SVD with all singular values)
        call dgesvd('A', 'A', n, n, A, n, s, u, n, vt, n, work, lwork, info)
        
        ! Check for success
        if (info /= 0) then
            print *, "SVD failed, INFO =", info
            cond_num = -1.0d0
            return
        endif
        
        ! Compute the condition number as the ratio of the largest to the smallest singular values
        cond_num = maxval(s) / minval(s)
        
        ! Deallocate work arrays
        deallocate(u, vt, work)
    end subroutine compute_cond_number


    
















end module potentials_gw
