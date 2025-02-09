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
        
            casida_Xs=casida_Kx
        
            deallocate(DSYEV_work)

        
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
