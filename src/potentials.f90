module potentials

    use tools
    use basisfunctions
    use class_linkedlist
    !use xc_f90_lib_m

    implicit none

    public
    !!!! Potential !!!!
    double precision, dimension(:), allocatable :: Vi,Vh,Vx,Vc,V_total,vh_temp,Vhfx
    CHARACTER(len=:), allocatable :: dummyc
    ! type(psi_gauss),dimension(:,:),allocatable :: Vh_g,Vx_g,Vc_g
    double precision,dimension(:),allocatable :: Vh_g,Vx_g,Vc_g,Vhfx_g,V_total_g
    double precision,dimension(:,:),allocatable :: hf_v_kernel,hf_Vhfx,Vhfx_final
    double precision,dimension(:,:),allocatable :: NNmatrix00

    !!!!! for libxc 
    double precision,dimension(:),allocatable :: vx_pbe_n,vc_pbe_n,vx_g_pbe_n,vc_g_pbe_n
    double precision,dimension(:,:),allocatable :: vx_pbe_g,vc_pbe_g,vx_g_pbe_g,vc_g_pbe_g
    double precision,dimension(:),allocatable :: ex_pbe,ec_pbe,ex_g_pbe,ec_g_pbe
    



    contains


    subroutine exchangepotential(Nn)

        integer,intent(in) :: Nn
    
        Vx(1:Nn) = -(3.0d0/pi)**(1.0d0/3.0d0)*ni(1:Nn)**(1.0d0/3.0d0)
    
    end subroutine exchangepotential
    
    
    
    subroutine correlationpotential(Nn)
    
        integer,intent(in) :: Nn
        double precision :: rs,eps,Ec_density
        integer :: i
    
        do i=1,Nn
    
            eps = 0.0d0
            if (ni(i)==0.0d0) eps=1E-10
    
            rs = (3.0d0/4.0d0/pi/(ni(i)+eps))**(1.0d0/3.0d0)
    
            if (rs<1.0d0) then
    
    Ec_density = ni(i)*(0.0311d0*log(rs)-0.0480d0+0.002d0*rs*log(rs)-0.0116d0*rs)
    Vc(i) = 0.0311d0*log(rs)-0.0480d0-1.0d0/3.0d0*0.0311d0+2.0d0/3.0d0*0.0002d0*rs*log(rs)+1.0d0/3.0d0*(-2.0d0*0.0116d0-0.002d0)*rs
    
            else if (rs>=1.0d0) then
    
    Ec_density = ni(i)*(-0.1423d0)/(1.0d0+1.0529d0*sqrt(rs)+0.3334d0*rs)
    Vc(i) = Ec_density/(ni(i)+eps)*(1.0d0+7.0d0/6.0d0*1.0529d0*sqrt(rs)+4.0d0/3.0d0*0.3334d0*rs)&
                                        /(1.0d0+1.0529d0*sqrt(rs)+0.3334d0*rs)
    
            end if
            
        enddo
    
    
    end subroutine correlationpotential




    subroutine construct_hf_v_kernel(Nn,Ne,Nquadrature)

        integer,intent(in) :: Nn,Ne,Nquadrature
        integer :: i,j,k

        allocate(hf_v_kernel(Nn,Ne*Nquadrature))

        do i=1,Nn
            do j=1,Ne
            do k=1,Nquadrature
            hf_v_kernel(i,(j-1)*Nquadrature+k)=1.0d0/sqrt((point(i,1)-point_g((j-1)*Nquadrature+k,1))**2&
                                                      +(point(i,2)-point_g((j-1)*Nquadrature+k,2))**2&
                                                      +(point(i,3)-point_g((j-1)*Nquadrature+k,3))**2)
            enddo
            enddo
        enddo

        ! print *,point(9,:)
        ! print *,point_g(1,:)

        ! print *,hf_v_kernel(9,1)

    end subroutine construct_hf_v_kernel

    subroutine construct_hf_exchangekernel(Nstates)

        integer,intent(in) :: Nstates
        ! double precision, dimension(:,:), intent(in) :: psi_temp
        integer :: k
        double precision, dimension(:,:),allocatable :: A

        allocate(A(size(psii(:,k)),size(psi_point_g(:,k))))
        
        hf_Vhfx=0.0d0
        do k=1,Nstates
            call outer_product('ge',psii(:,k),psi_point_g(:,k),A)
            hf_Vhfx = hf_Vhfx - A
        enddo

        hf_Vhfx=hf_Vhfx*hf_v_kernel

        deallocate(A)

    end subroutine construct_hf_exchangekernel












    subroutine set_linkedlistV(Ne,Nlocal,Nquadrature,type)

        integer,intent(in) :: Ne,Nlocal,Nquadrature,type
        integer :: i,m,n,k,ii,ll
        complex(kind=(kind(1.0d0))) :: x,y
        logical :: mtest
        ! double precision, dimension(:,:),allocatable :: del_m,del_n
        ! double precision, dimension(:), allocatable :: basis_m,basis_n
        double precision :: A_temp,V_temp,ik,k2,x_real,y_real,x_imag,y_imag,temp
        double precision,dimension(1:3) :: location

        double precision, dimension(1:3) :: temp_vector
    
    
    
        ! allocate(del_m(1:Nlocal,1:3))
        ! allocate(del_n(1:Nlocal,1:3))
        ! allocate(basis_m(1:Nlocal))
        ! allocate(basis_n(1:Nlocal))


        if (Nlocal==10) then
    
    
        if (type==0) then
    
      do i=1,Ne
    
        call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
        call mass_mat
        call gaussian_integral
    
        do m=1,Nlocal
            do n=1,Nlocal
                A_temp=0.0d0
                V_temp=0.0d0
                ik=0.0d0
                k2=0.0d0
                do k=1,Nquadrature
                    ! call basis_p2_del(gpoint(k,:))
                    ! call basis_p2_del(gpoint(k,:))
                    ! call basis_p2(gpoint(k,:))
                    call basis_p2(gpoint(k,:))

                
                temp=0.0d0
                do ii=1,Nlocal
                  temp=temp+V_total(ele(i,ii))*phi_p2(ii)
                enddo
    
                V_temp = V_temp + temp*phi_p2(m)*phi_p2(n)*gweight(k)

                ! V_temp = V_temp + V_total_g((i-1)*Nquadrature+k)*phi_p2(m)*phi_p2(n)*gweight(k)
    
                
    
                enddo
    
                x_real=V_temp*(Jdet)/6.0d0
                y_real=Jdet*p2_matrix(m,n)
    
                x_imag=0.0d0
                y_imag=0.0d0
                
                x=cmplx(x_real,x_imag,8)
                y=cmplx(y_real,y_imag,8)
    
                
                mtest=.true.
                if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
    
     
              enddo
          enddo
      enddo
    
    
    else if (type==1) then
    
        do i=1,Ne
    
            call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
            call mass_mat
            call gaussian_integral
        
            do m=1,Nlocal
                do n=1,Nlocal
                    A_temp=0.0d0
                    V_temp=0.0d0
                    ik=0.0d0
                    k2=0.0d0
                    do k=1,Nquadrature
                        ! call basis_p2_del(gpoint(k,:))
                        ! call basis_p2_del(gpoint(k,:))
                        ! call basis_p2(gpoint(k,:))
                        call basis_p2(gpoint(k,:))
        
                  location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
    
                do ll=1,Nbat
                    
                    V_temp = V_temp + phi_p2(m)*phi_p2(n)*gweight(k)*(-at(ll)%core)/distance(location(:),at(ll)%c(:))

                  
                enddo
                    
        
                    enddo
        
                    x_real=V_temp*(Jdet)/6.0d0
                    y_real=Jdet*p2_matrix(m,n)
        
                    x_imag=0.0d0
                    y_imag=0.0d0
                    
                    x=cmplx(x_real,x_imag,8)
                    y=cmplx(y_real,y_imag,8)
        
                    
                    mtest=.true.
                    if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                    if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                    call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
        
         
                  enddo
              enddo
          enddo


        else if (type==2) then
    
            do i=1,Ne
          
              call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
              call mass_mat
              call gaussian_integral
          
              do m=1,Nlocal
                  do n=1,Nlocal
                      A_temp=0.0d0
                      V_temp=0.0d0
                      ik=0.0d0
                      k2=0.0d0
                      do k=1,Nquadrature
                        call basis_p2_del(gpoint(k,:))
                        ! call basis_p2_del(gpoint(k,:))
                        ! call basis_p2(gpoint(k,:))
                        call basis_p2(gpoint(k,:))
      
                      
                    

        temp=0.0d0
        do ii=1,Nlocal
            temp=temp+V_total(ele(i,ii))*phi_p2(ii)
        enddo

        V_temp = V_temp + phi_p2(m)*phi_p2(n)*gweight(k)*temp
          
        ! V_temp = V_temp + basis_m(m)*basis_n(n)*gweight(k)*(vx_g_pbe_n((i-1)*Nquadrature+k)+vc_g_pbe_n((i-1)*Nquadrature+k))!*2.0d0
        ! V_temp=V_temp+(dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
        !                 +dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
        !                 *gweight(k)

        ! V_temp=V_temp+(dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
        !                 +dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
        !                 *gweight(k)


        temp_vector=0.0d0
        do ii=1,Nlocal
        temp_vector(:)=temp_vector(:)+(vx_pbe_g(ele(i,ii),:)+vc_pbe_g(ele(i,ii),:))*phi_p2(ii)
        enddo

        V_temp=V_temp+(dot_product(temp_vector(:),matmul(Jit,phi_p2_del(m,:)))*phi_p2(n)&
                        +dot_product(temp_vector(:),matmul(Jit,phi_p2_del(n,:)))*phi_p2(m))&
                        *gweight(k)
          
                      
          
                      enddo
          
                      x_real=V_temp*(Jdet)/6.0d0
                      y_real=Jdet*p2_matrix(m,n)
          
                      x_imag=0.0d0
                      y_imag=0.0d0
                      
                      x=cmplx(x_real,x_imag,8)
                      y=cmplx(y_real,y_imag,8)
          
                      
                      mtest=.true.
                      if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                      if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                      call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
          
           
                    enddo
                enddo
            enddo


        else if (type==3) then
    
            do i=1,Ne
          
              call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
              call mass_mat
              call gaussian_integral
          
              do m=1,Nlocal
                  do n=1,Nlocal
                      A_temp=0.0d0
                      V_temp=0.0d0
                      ik=0.0d0
                      k2=0.0d0
                      do k=1,Nquadrature
                        call basis_p2_del(gpoint(k,:))
                        ! call basis_p2_del(gpoint(k,:))
                        ! call basis_p2(gpoint(k,:))
                        call basis_p2(gpoint(k,:))
      
                      
                    

        ! temp=0.0d0
        ! do ii=1,Nlocal
        !     temp=temp+V_total(ele(i,ii))*basis_m(ii)
        ! enddo

        ! V_temp = V_temp + basis_m(m)*basis_n(n)*gweight(k)*temp
          
        ! V_temp = V_temp + basis_m(m)*basis_n(n)*gweight(k)*(vx_g_pbe_n((i-1)*Nquadrature+k)+vc_g_pbe_n((i-1)*Nquadrature+k))!*2.0d0
        ! V_temp=V_temp+(dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
        !                 +dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
        !                 *gweight(k)

        ! V_temp=V_temp+(dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
        !                 +dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
        !                 *gweight(k)


        temp_vector=0.0d0
        do ii=1,Nlocal
        temp_vector(:)=temp_vector(:)+(vx_pbe_g(ele(i,ii),:)+vc_pbe_g(ele(i,ii),:))*phi_p2(ii)
        enddo

        V_temp=V_temp+(dot_product(temp_vector(:),matmul(Jit,phi_p2_del(m,:)))*phi_p2(n)&
                        +dot_product(temp_vector(:),matmul(Jit,phi_p2_del(n,:)))*phi_p2(m))&
                        *gweight(k)
          
                      
          
                      enddo
          
                      x_real=V_temp*(Jdet)/6.0d0
                      y_real=Jdet*p2_matrix(m,n)
          
                      x_imag=0.0d0
                      y_imag=0.0d0
                      
                      x=cmplx(x_real,x_imag,8)
                      y=cmplx(y_real,y_imag,8)
          
                      
                      mtest=.true.
                      if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                      if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                      call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
          
           
                    enddo
                enddo
            enddo
        
        
        end if



    else if (Nlocal==20) then




        if (type==0) then
    
            do i=1,Ne
          
              call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
              call mass_mat
              call gaussian_integral
          
              do m=1,Nlocal
                  do n=1,Nlocal
                      A_temp=0.0d0
                      V_temp=0.0d0
                      ik=0.0d0
                      k2=0.0d0
                      do k=1,Nquadrature
                          ! call basis_p2_del(gpoint(k,:))
                          ! call basis_p2_del(gpoint(k,:))
                          ! call basis_p2(gpoint(k,:))
                          call basis_p3(gpoint(k,:))
      
                      
                      temp=0.0d0
                      do ii=1,Nlocal
                        temp=temp+V_total(ele(i,ii))*phi_p3(ii)
                      enddo
          
                      V_temp = V_temp + temp*phi_p3(m)*phi_p3(n)*gweight(k)

                    !   V_temp = V_temp + V_total_g((i-1)*Nquadrature+k)*phi_p3(m)*phi_p3(n)*gweight(k)
          
                      
          
                      enddo
          
                      x_real=V_temp*(Jdet)/6.0d0
                      y_real=Jdet*p3_matrix(m,n)
          
                      x_imag=0.0d0
                      y_imag=0.0d0
                      
                      x=cmplx(x_real,x_imag,8)
                      y=cmplx(y_real,y_imag,8)
          
                      
                      mtest=.true.
                      if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                      if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                      call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
          
           
                    enddo
                enddo
            enddo
          
          
          else if (type==1) then
          
              do i=1,Ne
          
                  call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
                  call mass_mat
                  call gaussian_integral
              
                  do m=1,Nlocal
                      do n=1,Nlocal
                          A_temp=0.0d0
                          V_temp=0.0d0
                          ik=0.0d0
                          k2=0.0d0
                          do k=1,Nquadrature
                              ! call basis_p2_del(gpoint(k,:))
                              ! call basis_p2_del(gpoint(k,:))
                              ! call basis_p2(gpoint(k,:))
                              call basis_p3(gpoint(k,:))
              
                        location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
          
                      do ll=1,Nbat
                          
                          V_temp = V_temp + phi_p3(m)*phi_p3(n)*gweight(k)*(-at(ll)%core)/distance(location(:),at(ll)%c(:))
      
                        
                      enddo
                          
              
                          enddo
              
                          x_real=V_temp*(Jdet)/6.0d0
                          y_real=Jdet*p3_matrix(m,n)
              
                          x_imag=0.0d0
                          y_imag=0.0d0
                          
                          x=cmplx(x_real,x_imag,8)
                          y=cmplx(y_real,y_imag,8)
              
                          
                          mtest=.true.
                          if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                          if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                          call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
              
               
                        enddo
                    enddo
                enddo
      
      
              else if (type==2) then
          
                  do i=1,Ne
                
                    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
                    call mass_mat
                    call gaussian_integral
                
                    do m=1,Nlocal
                        do n=1,Nlocal
                            A_temp=0.0d0
                            V_temp=0.0d0
                            ik=0.0d0
                            k2=0.0d0
                            do k=1,Nquadrature
                              call basis_p3_del(gpoint(k,:))
                              ! call basis_p2_del(gpoint(k,:))
                              ! call basis_p2(gpoint(k,:))
                              call basis_p3(gpoint(k,:))
            
                            
                          
      
              temp=0.0d0
              do ii=1,Nlocal
                  temp=temp+V_total(ele(i,ii))*phi_p3(ii)
              enddo
      
              V_temp = V_temp + phi_p3(m)*phi_p3(n)*gweight(k)*temp
                
              ! V_temp = V_temp + basis_m(m)*basis_n(n)*gweight(k)*(vx_g_pbe_n((i-1)*Nquadrature+k)+vc_g_pbe_n((i-1)*Nquadrature+k))!*2.0d0
              ! V_temp=V_temp+(dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
              !                 +dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
              !                 *gweight(k)
      
              ! V_temp=V_temp+(dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
              !                 +dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
              !                 *gweight(k)
      
      
              temp_vector=0.0d0
              do ii=1,Nlocal
              temp_vector(:)=temp_vector(:)+(vx_pbe_g(ele(i,ii),:)+vc_pbe_g(ele(i,ii),:))*phi_p3(ii)
              enddo
      
              V_temp=V_temp+(dot_product(temp_vector(:),matmul(Jit,phi_p3_del(m,:)))*phi_p3(n)&
                              +dot_product(temp_vector(:),matmul(Jit,phi_p3_del(n,:)))*phi_p3(m))&
                              *gweight(k)
                
                            
                
                            enddo
                
                            x_real=V_temp*(Jdet)/6.0d0
                            y_real=Jdet*p3_matrix(m,n)
                
                            x_imag=0.0d0
                            y_imag=0.0d0
                            
                            x=cmplx(x_real,x_imag,8)
                            y=cmplx(y_real,y_imag,8)
                
                            
                            mtest=.true.
                            if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                            if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                            call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
                
                 
                          enddo
                      enddo
                  enddo
      
      
              else if (type==3) then
          
                  do i=1,Ne
                
                    call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
                    call mass_mat
                    call gaussian_integral
                
                    do m=1,Nlocal
                        do n=1,Nlocal
                            A_temp=0.0d0
                            V_temp=0.0d0
                            ik=0.0d0
                            k2=0.0d0
                            do k=1,Nquadrature
                              call basis_p3_del(gpoint(k,:))
                              ! call basis_p2_del(gpoint(k,:))
                              ! call basis_p2(gpoint(k,:))
                              call basis_p3(gpoint(k,:))
            
                            
                          
      
              ! temp=0.0d0
              ! do ii=1,Nlocal
              !     temp=temp+V_total(ele(i,ii))*basis_m(ii)
              ! enddo
      
              ! V_temp = V_temp + basis_m(m)*basis_n(n)*gweight(k)*temp
                
              ! V_temp = V_temp + basis_m(m)*basis_n(n)*gweight(k)*(vx_g_pbe_n((i-1)*Nquadrature+k)+vc_g_pbe_n((i-1)*Nquadrature+k))!*2.0d0
              ! V_temp=V_temp+(dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
              !                 +dot_product(vx_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
              !                 *gweight(k)
      
              ! V_temp=V_temp+(dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_m(m,:)))*basis_n(n)&
              !                 +dot_product(vc_g_pbe_g((i-1)*Nquadrature+k,:),matmul(Jit,del_n(n,:)))*basis_m(m))&
              !                 *gweight(k)
      
      
              temp_vector=0.0d0
              do ii=1,Nlocal
              temp_vector(:)=temp_vector(:)+(vx_pbe_g(ele(i,ii),:)+vc_pbe_g(ele(i,ii),:))*phi_p3(ii)
              enddo
      
              V_temp=V_temp+(dot_product(temp_vector(:),matmul(Jit,phi_p3_del(m,:)))*phi_p3(n)&
                              +dot_product(temp_vector(:),matmul(Jit,phi_p3_del(n,:)))*phi_p3(m))&
                              *gweight(k)
                
                            
                
                            enddo
                
                            x_real=V_temp*(Jdet)/6.0d0
                            y_real=Jdet*p3_matrix(m,n)
                
                            x_imag=0.0d0
                            y_imag=0.0d0
                            
                            x=cmplx(x_real,x_imag,8)
                            y=cmplx(y_real,y_imag,8)
                
                            
                            mtest=.true.
                            if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                            if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                            call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
                
                 
                          enddo
                      enddo
                  enddo
              
              
              end if


    end if  !!!!!!!! end selection p2/p3
    
    
    end subroutine set_linkedlistV


    subroutine set_linkedlistV_efem(Nn,Ne,Nlocal,Nquadrature,Nquadrature_efem,at_c,r0,t,Z,type)

        integer,intent(in) :: Nn,Ne,Nlocal,Nquadrature,Nquadrature_efem,type
        double precision, intent(in) :: r0,t,Z
        double precision,dimension(1:3),intent(in) :: at_c
        integer :: i,m,n,k,ii,ll,kk
        complex(kind=(kind(1.0d0))) :: x,y
        logical :: mtest
        ! double precision, dimension(:,:),allocatable :: del_m,del_n
        ! double precision, dimension(:), allocatable :: basis_m,basis_n
        double precision :: A_temp,V_temp,ik,k2,x_real,y_real,x_imag,y_imag,temp
        double precision,dimension(:),allocatable :: location,gweight_temp
        double precision,dimension(:,:),allocatable :: gpoint_temp

        double precision, dimension(1:3) :: temp_vector

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
    
        do m=1,Nlocal
            do n=1,Nlocal
                A_temp=0.0d0
                V_temp=0.0d0
                ik=0.0d0
                k2=0.0d0
                do k=1,Nquadrature

                    call basis_p2(gpoint(k,:))

                
                    location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
    
                    do ll=1,Nbat
                        
                        V_temp = V_temp + phi_p2(m)*phi_p2(n)*gweight(k)*(-at(ll)%core)/distance(location(:),at(ll)%c(:))
    
                      
                    enddo
    
                
    
                enddo
    
                x_real=V_temp*(Jdet)/6.0d0
                y_real=Jdet*p2_matrix(m,n)
    
                x_imag=0.0d0
                y_imag=0.0d0
                
                x=cmplx(x_real,x_imag,8)
                y=cmplx(y_real,y_imag,8)
    
                
                mtest=.true.
                if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
    
     
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
                        call basis_p2_del(gpoint_temp(k,:))
                        location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                        do ll=1,Nbat
                            A_temp = A_temp + phi_p2(m)*phi_p2(n)*gweight_temp(k)*(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                        enddo
                    enddo
    
                    x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                    y_real=0.0d0!Jdet*p2_matrix(m,n)
                    x_imag=0.0d0
                    y_imag=0.0d0
                    x=cmplx(x_real,x_imag,8)
                    y=cmplx(y_real,y_imag,8)
                    mtest=.true.
                    if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                    if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                    call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)

                enddo
            enddo


            do m=1,Nlocal
                do n=1,Nlocal

                    A_temp=0.0d0
                    do k=1,Nquadrature_efem
                        call basis_p2_del(gpoint_temp(k,:)) 
                        location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                        ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                        call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                        A_temp = A_temp + efem_phi(1)*phi_p2(m)*efem_phi(1)*phi_p2(n)*gweight_temp(k)&
                                                                *(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                    enddo
    
                    x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                    y_real=0.0d0!Jdet*V_temp/6.0d0 !!!!!!! new basis 
    
    
                    x_imag=0.0d0
                    y_imag=0.0d0
                    x=cmplx(x_real,x_imag,8)
                    y=cmplx(y_real,y_imag,8)
                    mtest=.true.
                    if (neigh_V(enriched_ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
                    if (mtest) call neigh_V(enriched_ele(i,m))%insert(enriched_ele(i,n))
                    call neigh_V(enriched_ele(i,m))%insertAB(enriched_ele(i,n),x,y)

                enddo
            enddo


            do m=1,Nlocal
                do n=1,Nlocal

                    A_temp=0.0d0
                    do k=1,Nquadrature_efem
                        call basis_p2_del(gpoint_temp(k,:)) 
                        location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                        ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                        call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                        A_temp = A_temp + efem_phi(1)*phi_p2(m)*phi_p2(n)*gweight_temp(k)&
                                                                *(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                    enddo
    
                    x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                    y_real=0.0d0!Jdet*V_temp/6.0d0 !!!!!!! new basis 
    
    
                    x_imag=0.0d0
                    y_imag=0.0d0
                    x=cmplx(x_real,x_imag,8)
                    y=cmplx(y_real,y_imag,8)
                    mtest=.true.
                    if (neigh_V(enriched_ele(i,m))%find(ele(i,n)))  mtest=.false.
                    if (mtest) call neigh_V(enriched_ele(i,m))%insert(ele(i,n))
                    call neigh_V(enriched_ele(i,m))%insertAB(ele(i,n),x,y)

                enddo
            enddo

            do m=1,Nlocal
                do n=1,Nlocal

                    A_temp=0.0d0
                    do k=1,Nquadrature_efem
                        call basis_p2_del(gpoint_temp(k,:)) 
                        location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                        ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                        call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                        A_temp = A_temp + phi_p2(m)*efem_phi(1)*phi_p2(n)*gweight_temp(k)&
                                                                *(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                    enddo
    
                    x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                    y_real=0.0d0!Jdet*V_temp/6.0d0 !!!!!!! new basis 
    
    
                    x_imag=0.0d0
                    y_imag=0.0d0
                    x=cmplx(x_real,x_imag,8)
                    y=cmplx(y_real,y_imag,8)
                    mtest=.true.
                    if (neigh_V(ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
                    if (mtest) call neigh_V(ele(i,m))%insert(enriched_ele(i,n))
                    call neigh_V(ele(i,m))%insertAB(enriched_ele(i,n),x,y)

                enddo
            enddo
        

            end if


      enddo

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
        
            do m=1,Nlocal
                do n=1,Nlocal
                    A_temp=0.0d0
                    V_temp=0.0d0
                    ik=0.0d0
                    k2=0.0d0
                    do k=1,Nquadrature
    
                        call basis_p3(gpoint(k,:))
    
                    
                        location = matmul(J0,gpoint(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
        
                        do ll=1,Nbat
                            
                            V_temp = V_temp + phi_p3(m)*phi_p3(n)*gweight(k)*(-at(ll)%core)/distance(location(:),at(ll)%c(:))
        
                          
                        enddo
        
                    
        
                    enddo
        
                    x_real=V_temp*(Jdet)/6.0d0
                    y_real=Jdet*p3_matrix(m,n)
        
                    x_imag=0.0d0
                    y_imag=0.0d0
                    
                    x=cmplx(x_real,x_imag,8)
                    y=cmplx(y_real,y_imag,8)
        
                    
                    mtest=.true.
                    if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                    if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                    call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
        
         
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
                            call basis_p3(gpoint_temp(k,:))
                            location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                            do ll=1,Nbat
                            A_temp = A_temp + phi_p3(m)*phi_p3(n)*gweight_temp(k)*(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                            enddo
                        enddo
        
                        x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                        y_real=0.0d0!Jdet*p2_matrix(m,n)
                        x_imag=0.0d0
                        y_imag=0.0d0
                        x=cmplx(x_real,x_imag,8)
                        y=cmplx(y_real,y_imag,8)
                        mtest=.true.
                        if (neigh_V(ele(i,m))%find(ele(i,n)))  mtest=.false.
                        if (mtest) call neigh_V(ele(i,m))%insert(ele(i,n))
                        call neigh_V(ele(i,m))%insertAB(ele(i,n),x,y)
    
                    enddo
                enddo
    
    
                do m=1,Nlocal
                    do n=1,Nlocal
    
                        A_temp=0.0d0
                        do k=1,Nquadrature_efem
                            call basis_p3(gpoint_temp(k,:)) 
                            location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                            ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                            call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                            A_temp = A_temp + efem_phi(1)*phi_p3(m)*efem_phi(1)*phi_p3(n)*gweight_temp(k)&
                                                                    *(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                        enddo
        
                        x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                        y_real=0.0d0!Jdet*V_temp/6.0d0 !!!!!!! new basis 
        
        
                        x_imag=0.0d0
                        y_imag=0.0d0
                        x=cmplx(x_real,x_imag,8)
                        y=cmplx(y_real,y_imag,8)
                        mtest=.true.
                        if (neigh_V(enriched_ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
                        if (mtest) call neigh_V(enriched_ele(i,m))%insert(enriched_ele(i,n))
                        call neigh_V(enriched_ele(i,m))%insertAB(enriched_ele(i,n),x,y)
    
                    enddo
                enddo
    
    
                do m=1,Nlocal
                    do n=1,Nlocal
    
                        A_temp=0.0d0
                        do k=1,Nquadrature_efem
                            call basis_p3(gpoint_temp(k,:)) 
                            location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                            ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                            call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                            A_temp = A_temp + efem_phi(1)*phi_p3(m)*phi_p3(n)*gweight_temp(k)&
                                                                    *(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                        enddo
        
                        x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                        y_real=0.0d0!Jdet*V_temp/6.0d0 !!!!!!! new basis 
        
        
                        x_imag=0.0d0
                        y_imag=0.0d0
                        x=cmplx(x_real,x_imag,8)
                        y=cmplx(y_real,y_imag,8)
                        mtest=.true.
                        if (neigh_V(enriched_ele(i,m))%find(ele(i,n)))  mtest=.false.
                        if (mtest) call neigh_V(enriched_ele(i,m))%insert(ele(i,n))
                        call neigh_V(enriched_ele(i,m))%insertAB(ele(i,n),x,y)
    
                    enddo
                enddo
    
                do m=1,Nlocal
                    do n=1,Nlocal
    
                        A_temp=0.0d0
                        do k=1,Nquadrature_efem
                            call basis_p3(gpoint_temp(k,:)) 
                            location = matmul(J0,gpoint_temp(k,:))+(/point(ele(i,1),1),point(ele(i,1),2),point(ele(i,1),3)/)
                            ! call efem_psi_1s_truncated(location,at_c,Z,r0)
                            call efem_psi_1s_truncated(location,at_c,Z,r0,t)
                            A_temp = A_temp + phi_p3(m)*efem_phi(1)*phi_p3(n)*gweight_temp(k)&
                                                                    *(-at(ll)%core)/distance(location(:),at(ll)%c(:))
                        enddo
        
                        x_real=A_temp*(Jdet)/6.0d0!+V_temp0*(Jdet)/6.0d0
                        y_real=0.0d0!Jdet*V_temp/6.0d0 !!!!!!! new basis 
        
        
                        x_imag=0.0d0
                        y_imag=0.0d0
                        x=cmplx(x_real,x_imag,8)
                        y=cmplx(y_real,y_imag,8)
                        mtest=.true.
                        if (neigh_V(ele(i,m))%find(enriched_ele(i,n)))  mtest=.false.
                        if (mtest) call neigh_V(ele(i,m))%insert(enriched_ele(i,n))
                        call neigh_V(ele(i,m))%insertAB(enriched_ele(i,n),x,y)
    
                    enddo
                enddo
            
    
                end if
    
    
          enddo

        end if

        deallocate(location)
        deallocate(gweight_temp)
        deallocate(gpoint_temp)

    end subroutine set_linkedlistV_efem



    function compute_hartreeBC(colored,Ne,Nquadrature,Nstates,spin,type2,it) result(BC_value)

        integer,intent(in) :: colored,Ne,Nquadrature,Nstates,type2,it
        double precision, intent(in) :: spin
        double precision :: BC_value,r0,n_temp
        integer :: i,j,ii

        call gaussian_integral

    if (type2==0) then
    
        BC_value=0.0d0
        do i=1,Ne
            do j=1,Nquadrature
              
                n_temp=0.0d0
                do ii=1,Nstates
                    if (ii<Nstates) then
                        n_temp=n_temp+2.0d0*psi_point_g((i-1)*Nquadrature+j,ii)**2
                    else if (ii==Nstates) then
                        n_temp=n_temp+spin*psi_point_g((i-1)*Nquadrature+j,ii)**2
                    end if
                enddo
                ! do ii=1,Nstates
                !     if (ii<Nstates-it-1) then
                !         n_temp=n_temp+2.0d0*psi_point_g((i-1)*Nquadrature+j,ii)**2
                !     else if (ii>Nstates-it) then
                !         n_temp=n_temp+(2.0d0-spin)*psi_point_g((i-1)*Nquadrature+j,ii)**2
                !     end if
                ! enddo
    
                r0=sqrt((point_g((i-1)*Nquadrature+j,1)-point(colored,1))**2&
                        +(point_g((i-1)*Nquadrature+j,2)-point(colored,2))**2+(point_g((i-1)*Nquadrature+j,3)-point(colored,3))**2)

                BC_value=BC_value+volume(i)*gweight(j)*n_temp*1.0d0/r0
    
            enddo
            
        enddo


    else if (type2==1) then
    
        BC_value=0.0d0
        do i=1,Ne
            do j=1,Nquadrature
                
                n_temp=0.0d0
                ! do ii=1,Nstates
                !     if (ii<Nstates) then
                !         n_temp=n_temp+2.0d0*psi_point_g((i-1)*Nquadrature+j,ii)**2
                !     else if (ii==Nstates) then
                !         n_temp=n_temp+spin*psi_point_g((i-1)*Nquadrature+j,ii)**2
                !     end if
                ! enddo
                do ii=1,Nstates
                    if (ii<Nstates-it-1) then
                        n_temp=n_temp+2.0d0*psi_point_g((i-1)*Nquadrature+j,ii)**2
                    else if (ii>Nstates-it) then
                        n_temp=n_temp+(2.0d0-spin)*psi_point_g((i-1)*Nquadrature+j,ii)**2
                    end if
                enddo
    
                r0=sqrt((point_g((i-1)*Nquadrature+j,1)-point(colored,1))**2&
                        +(point_g((i-1)*Nquadrature+j,2)-point(colored,2))**2+(point_g((i-1)*Nquadrature+j,3)-point(colored,3))**2)

                BC_value=BC_value+volume(i)*gweight(j)*n_temp*1.0d0/r0
    
            enddo
            
        enddo

    else if (type2==2) then

        BC_value=0.0d0
        do i=1,Ne
            do j=1,Nquadrature
              
                n_temp=0.0d0
                do ii=1,Nstates
                    if (ii<Nstates) then
                        n_temp=n_temp+2.0d0*psi_point_g((i-1)*Nquadrature+j,ii)**2
                    else if (ii==Nstates) then
                        if (it<5) then
                            n_temp=n_temp+(2.0d0-0.1d0-it/10.0d0)*psi_point_g((i-1)*Nquadrature+j,ii)**2
                        else
                            n_temp=n_temp+spin*psi_point_g((i-1)*Nquadrature+j,ii)**2
                        end if
                    end if
                enddo
    
                r0=sqrt((point_g((i-1)*Nquadrature+j,1)-point(colored,1))**2&
                        +(point_g((i-1)*Nquadrature+j,2)-point(colored,2))**2+(point_g((i-1)*Nquadrature+j,3)-point(colored,3))**2)

                BC_value=BC_value+volume(i)*gweight(j)*n_temp*1.0d0/r0
    
            enddo
            
        enddo

    end if
    
    end function compute_hartreeBC


    subroutine set_BC(Nn,Nbc,Ne,Nquadrature,Nstates,spin,type,type2,it)

        integer,intent(in) :: Nn,Nbc,Ne,Nquadrature,Nstates,type,type2,it
        double precision, intent(in) :: spin
        ! integer :: i
        integer :: i,j,jj,num1,num2,num3,index01,index02,index03
        logical :: ftest0

        vh_temp=0.0d0

        if (type==0) then

            do i=1,Nbc
                vh_temp(BC(i))=-Vi(BC(i))
            enddo

        else if (type==1) then
            
            do i=1,Nbc
                vh_temp(BC(i))=compute_hartreeBC(BC(i),Ne,Nquadrature,Nstates,spin,type2,it)
            enddo

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! set boundary conditions !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        A00=A

        
        matBC=0.0d0

        do i=1,Nbc
            num3=0
            do j=1,Nn
                num1=0
                num2=0
                ftest0=.true.
                if (neigh_AB(j)%find(BC(i)))  ftest0=.false.
                if (ftest0) then 
                    cycle
                else
                    index01=neigh_AB(j)%findindex(BC(i))
                    index02=neigh_AB(BC(i))%findindex(j)
                end if

                if (j==1) then 
                    num1=num1+index01
                else
                    do jj=1,j-1
                        num1=num1+neigh_AB(jj)%size
                    enddo
                    num1=num1+index01
                end if
                if (j/=BC(i)) matBC(j,i)=A00(num1)
                ! snq(j)=snq(j)-A00(num1)*vh_temp(BC(i))
                A00(num1)=cmplx(0.0d0,0.0d0,8)

                if (BC(i)==1) then
                    num2=num2+index02
                else
                    do jj=1,BC(i)-1
                        num2=num2+neigh_AB(jj)%size
                    enddo
                    num2=num2+index02
                end if
                A00(num2)=cmplx(0.0d0,0.0d0,8)
            enddo
            index03=neigh_AB(BC(i))%findindex(BC(i))
            if (BC(i)==1) then
                num3=num3+index03
            else
                do jj=1,BC(i)-1
                    num3=num3+neigh_AB(jj)%size
                enddo
                num3=num3+index03
            end if
            A00(num3)=cmplx(1.0d0,0.0d0,8)
            ! snq(BC(i))=vh_temp(BC(i))
        enddo


    end subroutine set_BC





    subroutine hartreepotential(Nn,Ne,Nquadrature,Nstates,Nbc,spin,it,rank)

        integer,intent(in) :: Nn,Ne,Nquadrature,Nstates,Nbc,it,rank
        double precision,intent(in) :: spin

        double precision,dimension(:),allocatable :: snq
        double precision,dimension(:,:),allocatable :: snq0,sn_g
        ! complex(kind=(kind(1.0d0))),dimension(:),allocatable :: snq_cmplx
        !complex(kind=(kind(1.0d0))),dimension(:),allocatable :: xy_AB0
        integer :: i!,j,jj,num1,num2,num3,index01,index02,index03
        double precision :: a

        !!! for mkl_dcsrgemv !!!
        character(len=1) :: transa='N'
        !logical :: ftest0
        !!! for pardiso !!!
        integer :: MAXFCT,MNUM,MTYPE,MSGLVL,PHASE,idum,error_hartree
        integer(8),dimension(:),allocatable :: pt_V
        integer,dimension(:),allocatable :: iparm_V
        
        complex(kind=(kind(1.0d0))) :: jc

        jc=cmplx(0.0d0,1.0d0,8)



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!   solve Poisson equation w/ BC    !!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(snq(Nn))
        ! allocate(snq0(Nn,1))
        ! allocate(sn_g(Ne*Nquadrature,1))

        

        ! if (it==0) then
        
            call mkl_dcsrgemv(transa, Nn, B, IA, JA, 4.0d0*pi*ni , snq)

        ! else

        !     sn_g=0.0d0

        !     do i=1,Nstates
        !         sn_g(:,1)=sn_g(:,1)+spin*psi_point_g(:,i)*psi_point_g(:,i)
        !     enddo

        !     call mkl_dcsrmm('T',Ne*Nquadrature,1,Nn,4.0d0*pi,matdescra,NnToNg,NnToNg_JA,NnToNg_IA_pntrb,NnToNg_IA_pntre&
        !     ,sn_g(:,1)*volumegweight,Ne*Nquadrature,0.0d0,snq0(:,1),Nn)

        !     snq=snq0(:,1)

        ! end if

        ! a=0.0d0
        ! do i=1,size(snq)
        !     a=a+snq(i)
        ! enddo

        ! print *,'========'
        ! print *,a/4.0d0/pi
        ! print *,'========'
            

        do i=1,Nbc
            snq(BC(i))=vh_temp(BC(i))
        enddo

        do i=1,Nbc
            snq(:)=snq(:)-matBC(:,i)*vh_temp(BC(i))
        enddo


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!   solve using Pardiso   !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MAXFCT=1
        MNUM=1
        MTYPE=1     ! 1 for real and spd 3 for complex and sdp
        MSGLVL=0     !0- no output, 1- output
        allocate(pt_v(64))
        allocate(iparm_v(64))
        pt_V=0
        call pardisoinit(PT_V,MTYPE,IPARM_V)
        !IPARM(11,:)=0 ! disable scaling
        PHASE=13 ! 13


        call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,A00,IA,JA,idum,1,IPARM_V,MSGLVL,snq,Vh,error_hartree)


        PHASE=-1

        call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,A00,IA,JA,idum,1,IPARM_V,MSGLVL,snq,Vh,error_hartree)

        if ((rank==0).and.(error_hartree/=0)) then
            print *,'PARDISO_error info -----', error_hartree
            stop
        end if

        deallocate(snq)
        ! deallocate(snq0)
        ! deallocate(sn_g)


    end subroutine hartreepotential





    function coulombpotential(location11,location2,d0,Z00) result(coulombpotential010)

        double precision, dimension(1:3) :: location11,location2
        double precision :: coulombpotential010,d0,dd,eps0
        integer :: Z00

        eps0=1E-8

        dd=sqrt((location11(1)-location2(1))**2+(location11(2)-location2(2))**2+(location11(3)-location2(3))**2)
        if (dd>d0) then
            coulombpotential010=0.0d0
        
        else if (dd>eps0) then
            coulombpotential010 = -Z00/dd 
        else if (dd<eps0) then
            coulombpotential010 = -10000.0d0
            ! coulombpotential010 = -30.0d0
        end if

    end function coulombpotential


    function nucleipotential(location,d) result(potential)

        double precision,dimension(1:3),intent(in) :: location
        double precision,intent(in) ::d
        double precision :: potential
        integer :: i
        ! character(len = 8) :: dummyc

        potential=0.0d0

        do i=1,nbat

            potential=potential+coulombpotential(location,at(i)%c(:),d,at(i)%core)
    
        enddo

        !     if (ll==1) then
        !         Z00=Z001
        !         potential=potential&
        !         +coulombpotential(location,cxyz(ll,:),d,Z00)
        !     else
        !         Z00=Z002
        !         potential=potential&
        !         +coulombpotential(location,cxyz(ll,:),d,Z00)
        !     end if

        ! enddo



    end function nucleipotential

    


    subroutine potential_interpolation(Ne,Nquadrature)
        integer,intent(in) :: Ne,Nquadrature
        integer :: i,ii,kk
        double precision,dimension(:),allocatable :: basis_function
        double precision,dimension(:),allocatable :: f
        double precision :: f0
    
    

        allocate(f(1:10))
        ! allocate(basis_function(1:10))
    
        do i=1,Ne

          call gaussian_integral
    
          do ii=1,Nquadrature
    
            call basis_p2(gpoint(ii,:))
    
            !!!!! Vh !!!!!
            f(1 )=Vh(ele(i,1 ))
            f(2 )=Vh(ele(i,2 ))
            f(3 )=Vh(ele(i,3 ))
            f(4 )=Vh(ele(i,4 ))
            f(5 )=Vh(ele(i,5 ))
            f(6 )=Vh(ele(i,6 ))
            f(7 )=Vh(ele(i,7 ))
            f(8 )=Vh(ele(i,8 ))
            f(9 )=Vh(ele(i,9 ))
            f(10)=Vh(ele(i,10))
    
            f0=0.0d0
            do kk=1,10
            f0=f0+f(kk)*phi_p2(kk)
            enddo
            Vh_g((i-1)*Nquadrature+ii)=f0

            !!!!! Vx !!!!!
            f(1 )=Vx(ele(i,1 ))
            f(2 )=Vx(ele(i,2 ))
            f(3 )=Vx(ele(i,3 ))
            f(4 )=Vx(ele(i,4 ))
            f(5 )=Vx(ele(i,5 ))
            f(6 )=Vx(ele(i,6 ))
            f(7 )=Vx(ele(i,7 ))
            f(8 )=Vx(ele(i,8 ))
            f(9 )=Vx(ele(i,9 ))
            f(10)=Vx(ele(i,10))
    
            f0=0.0d0
            do kk=1,10
            f0=f0+f(kk)*phi_p2(kk)
            enddo
            Vx_g((i-1)*Nquadrature+ii)=f0

            !!!!! Vc !!!!!
            f(1 )=Vc(ele(i,1 ))
            f(2 )=Vc(ele(i,2 ))
            f(3 )=Vc(ele(i,3 ))
            f(4 )=Vc(ele(i,4 ))
            f(5 )=Vc(ele(i,5 ))
            f(6 )=Vc(ele(i,6 ))
            f(7 )=Vc(ele(i,7 ))
            f(8 )=Vc(ele(i,8 ))
            f(9 )=Vc(ele(i,9 ))
            f(10)=Vc(ele(i,10))
    
            f0=0.0d0
            do kk=1,10
            f0=f0+f(kk)*phi_p2(kk)
            enddo
            Vc_g((i-1)*Nquadrature+ii)=f0
          
          enddo
        enddo

        deallocate(f)
    
      end subroutine potential_interpolation




      subroutine pbe_potential_interpolation(Ne,Nquadrature,degree)
        integer,intent(in) :: Ne,Nquadrature,degree
        integer :: i,ii,kk
        double precision,dimension(:),allocatable :: basis_function
        double precision,dimension(:),allocatable :: f
        double precision :: f0
    
    
        if (degree==2) then

        allocate(f(1:10))
        ! allocate(basis_function(1:10))
    
        do i=1,Ne

          call gaussian_integral
    
          do ii=1,Nquadrature
    
            call basis_p2(gpoint(ii,:))
    
            !!!!! Vh !!!!!
            f(1 )=Vh(ele(i,1 ))
            f(2 )=Vh(ele(i,2 ))
            f(3 )=Vh(ele(i,3 ))
            f(4 )=Vh(ele(i,4 ))
            f(5 )=Vh(ele(i,5 ))
            f(6 )=Vh(ele(i,6 ))
            f(7 )=Vh(ele(i,7 ))
            f(8 )=Vh(ele(i,8 ))
            f(9 )=Vh(ele(i,9 ))
            f(10)=Vh(ele(i,10))
    
            f0=0.0d0
            do kk=1,10
            f0=f0+f(kk)*phi_p2(kk)
            enddo
            Vh_g((i-1)*Nquadrature+ii)=f0

            !!!!! Vx !!!!!
            f(1 )=vx_pbe_n(ele(i,1 ))
            f(2 )=vx_pbe_n(ele(i,2 ))
            f(3 )=vx_pbe_n(ele(i,3 ))
            f(4 )=vx_pbe_n(ele(i,4 ))
            f(5 )=vx_pbe_n(ele(i,5 ))
            f(6 )=vx_pbe_n(ele(i,6 ))
            f(7 )=vx_pbe_n(ele(i,7 ))
            f(8 )=vx_pbe_n(ele(i,8 ))
            f(9 )=vx_pbe_n(ele(i,9 ))
            f(10)=vx_pbe_n(ele(i,10))
    
            f0=0.0d0
            do kk=1,10
            f0=f0+f(kk)*phi_p2(kk)
            enddo
            Vx_g((i-1)*Nquadrature+ii)=f0

            !!!!! Vc !!!!!
            f(1 )=vc_pbe_n(ele(i,1 ))
            f(2 )=vc_pbe_n(ele(i,2 ))
            f(3 )=vc_pbe_n(ele(i,3 ))
            f(4 )=vc_pbe_n(ele(i,4 ))
            f(5 )=vc_pbe_n(ele(i,5 ))
            f(6 )=vc_pbe_n(ele(i,6 ))
            f(7 )=vc_pbe_n(ele(i,7 ))
            f(8 )=vc_pbe_n(ele(i,8 ))
            f(9 )=vc_pbe_n(ele(i,9 ))
            f(10)=vc_pbe_n(ele(i,10))
    
            f0=0.0d0
            do kk=1,10
            f0=f0+f(kk)*phi_p2(kk)
            enddo
            Vc_g((i-1)*Nquadrature+ii)=f0
          
          enddo
        enddo

        deallocate(f)


    else if (degree==3) then

            allocate(f(1:20))
            ! allocate(basis_function(1:10))
        
            do i=1,Ne
    
              call gaussian_integral
        
              do ii=1,Nquadrature
        
                call basis_p3(gpoint(ii,:))
        
                !!!!! Vh !!!!!
                f(1 )=Vh(ele(i,1 ))
                f(2 )=Vh(ele(i,2 ))
                f(3 )=Vh(ele(i,3 ))
                f(4 )=Vh(ele(i,4 ))
                f(5 )=Vh(ele(i,5 ))
                f(6 )=Vh(ele(i,6 ))
                f(7 )=Vh(ele(i,7 ))
                f(8 )=Vh(ele(i,8 ))
                f(9 )=Vh(ele(i,9 ))
                f(10)=Vh(ele(i,10))
                f(11)=Vh(ele(i,11))
                f(12)=Vh(ele(i,12))
                f(13)=Vh(ele(i,13))
                f(14)=Vh(ele(i,14))
                f(15)=Vh(ele(i,15))
                f(16)=Vh(ele(i,16))
                f(17)=Vh(ele(i,17))
                f(18)=Vh(ele(i,18))
                f(19)=Vh(ele(i,19))
                f(20)=Vh(ele(i,20))
        
                f0=0.0d0
                do kk=1,20
                f0=f0+f(kk)*phi_p3(kk)
                enddo
                Vh_g((i-1)*Nquadrature+ii)=f0
    
                !!!!! Vx !!!!!
                f(1 )=vx_pbe_n(ele(i,1 ))
                f(2 )=vx_pbe_n(ele(i,2 ))
                f(3 )=vx_pbe_n(ele(i,3 ))
                f(4 )=vx_pbe_n(ele(i,4 ))
                f(5 )=vx_pbe_n(ele(i,5 ))
                f(6 )=vx_pbe_n(ele(i,6 ))
                f(7 )=vx_pbe_n(ele(i,7 ))
                f(8 )=vx_pbe_n(ele(i,8 ))
                f(9 )=vx_pbe_n(ele(i,9 ))
                f(10)=vx_pbe_n(ele(i,10))
                f(11)=vx_pbe_n(ele(i,11))
                f(12)=vx_pbe_n(ele(i,12))
                f(13)=vx_pbe_n(ele(i,13))
                f(14)=vx_pbe_n(ele(i,14))
                f(15)=vx_pbe_n(ele(i,15))
                f(16)=vx_pbe_n(ele(i,16))
                f(17)=vx_pbe_n(ele(i,17))
                f(18)=vx_pbe_n(ele(i,18))
                f(19)=vx_pbe_n(ele(i,19))
                f(20)=vx_pbe_n(ele(i,20))
        
                f0=0.0d0
                do kk=1,20
                f0=f0+f(kk)*phi_p3(kk)
                enddo
                Vx_g((i-1)*Nquadrature+ii)=f0
    
                !!!!! Vc !!!!!
                f(1 )=vc_pbe_n(ele(i,1 ))
                f(2 )=vc_pbe_n(ele(i,2 ))
                f(3 )=vc_pbe_n(ele(i,3 ))
                f(4 )=vc_pbe_n(ele(i,4 ))
                f(5 )=vc_pbe_n(ele(i,5 ))
                f(6 )=vc_pbe_n(ele(i,6 ))
                f(7 )=vc_pbe_n(ele(i,7 ))
                f(8 )=vc_pbe_n(ele(i,8 ))
                f(9 )=vc_pbe_n(ele(i,9 ))
                f(10)=vc_pbe_n(ele(i,10))
                f(11)=vc_pbe_n(ele(i,11))
                f(12)=vc_pbe_n(ele(i,12))
                f(13)=vc_pbe_n(ele(i,13))
                f(14)=vc_pbe_n(ele(i,14))
                f(15)=vc_pbe_n(ele(i,15))
                f(16)=vc_pbe_n(ele(i,16))
                f(17)=vc_pbe_n(ele(i,17))
                f(18)=vc_pbe_n(ele(i,18))
                f(19)=vc_pbe_n(ele(i,19))
                f(20)=vc_pbe_n(ele(i,20))
                
        
                f0=0.0d0
                do kk=1,20
                f0=f0+f(kk)*phi_p3(kk)
                enddo
                Vc_g((i-1)*Nquadrature+ii)=f0
              
              enddo
            enddo
    
            deallocate(f)

        end if
    
      end subroutine pbe_potential_interpolation



      subroutine Vhfx_interpolation(Ne,Nquadrature,degree)
        integer,intent(in) :: Ne,Nquadrature,degree
        integer :: i,ii,kk
        double precision,dimension(:),allocatable :: basis_function
        double precision,dimension(:),allocatable :: f
        double precision :: f0
    
    
        if (degree==2) then

        allocate(f(1:10))
        ! allocate(basis_function(1:10))
    
        do i=1,Ne

          call gaussian_integral
    
          do ii=1,Nquadrature
    
            call basis_p2(gpoint(ii,:))
    
            !!!!! Vhfx !!!!!
            f(1 )=Vhfx(ele(i,1 ))
            f(2 )=Vhfx(ele(i,2 ))
            f(3 )=Vhfx(ele(i,3 ))
            f(4 )=Vhfx(ele(i,4 ))
            f(5 )=Vhfx(ele(i,5 ))
            f(6 )=Vhfx(ele(i,6 ))
            f(7 )=Vhfx(ele(i,7 ))
            f(8 )=Vhfx(ele(i,8 ))
            f(9 )=Vhfx(ele(i,9 ))
            f(10)=Vhfx(ele(i,10))
    
            f0=0.0d0
            do kk=1,10
            f0=f0+f(kk)*phi_p2(kk)
            enddo
            Vhfx_g((i-1)*Nquadrature+ii)=f0
          
          enddo
        enddo

        deallocate(f)

        else if (degree==3) then

            allocate(f(1:20))
            ! allocate(basis_function(1:10))
        
            do i=1,Ne
    
              call gaussian_integral
        
              do ii=1,Nquadrature
        
                call basis_p3(gpoint(ii,:))
        
                !!!!! Vhfx !!!!!
                f(1 )=Vhfx(ele(i,1 ))
                f(2 )=Vhfx(ele(i,2 ))
                f(3 )=Vhfx(ele(i,3 ))
                f(4 )=Vhfx(ele(i,4 ))
                f(5 )=Vhfx(ele(i,5 ))
                f(6 )=Vhfx(ele(i,6 ))
                f(7 )=Vhfx(ele(i,7 ))
                f(8 )=Vhfx(ele(i,8 ))
                f(9 )=Vhfx(ele(i,9 ))
                f(10)=Vhfx(ele(i,10))
                f(11)=Vhfx(ele(i,11))
                f(12)=Vhfx(ele(i,12))
                f(13)=Vhfx(ele(i,13))
                f(14)=Vhfx(ele(i,14))
                f(15)=Vhfx(ele(i,15))
                f(16)=Vhfx(ele(i,16))
                f(17)=Vhfx(ele(i,17))
                f(18)=Vhfx(ele(i,18))
                f(19)=Vhfx(ele(i,19))
                f(20)=Vhfx(ele(i,20))
        
                f0=0.0d0
                do kk=1,20
                f0=f0+f(kk)*phi_p3(kk)
                enddo
                Vhfx_g((i-1)*Nquadrature+ii)=f0
              
              enddo
            enddo
    
            deallocate(f)


        end if
    
      end subroutine Vhfx_interpolation



    !   subroutine libxc_exchange_lda(Ne,Nquadrature)
    !     integer,intent(in) :: Ne,Nquadrature
    !     TYPE(xc_f90_func_t) :: xc_func
    !     TYPE(xc_f90_func_info_t) :: xc_info
    !     ! TYPE(xc_f90_pointer_t) :: xc_func
    !     ! TYPE(xc_f90_pointer_t) :: xc_info
    !     integer :: i, vmajor, vminor, vmicro, func_id = 1!!! GGA_LDA_X id is 1
    !     ! double precision,dimension(:),allocatable :: sigma
    !     ! double precision,dimension(:),allocatable :: vsigma
    !     integer(8) :: Ng

    !     Ng=Ne*Nquadrature

    !     vx_g_pbe_n=0.0d0
    !     ex_g_pbe=0.0d0

    !     ! allocate(vsigma(1:Ne*Nquadrature))
    !     ! allocate(sigma(1:Ne*Nquadrature))
    !     ! do i=1,Ng
    !     ! sigma(i) = nq_g_gradient(i,1)**2+nq_g_gradient(i,2)**2+nq_g_gradient(i,3)**2
    !     ! enddo

    !     ! call xc_f90_func_init(xc_func, xc_info, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     ! call xc_f90_func_init(xc_func, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
    !     call xc_f90_lda_vxc(xc_func, Ng, nq_g(1), vx_g_pbe_n(1))
    !     call xc_f90_lda_exc(xc_func, Ng, nq_g(1), ex_g_pbe(1))
    !     call xc_f90_func_end(xc_func)

    !     ! vx_g_pbe_n=vx_g_pbe_n*nq_g+ex_g_pbe

    !     ! do i=1,Ng
    !     ! vx_g_pbe_g(i,1) = vsigma(i)*2.0d0*nq_g_gradient(i,1)
    !     ! vx_g_pbe_g(i,2) = vsigma(i)*2.0d0*nq_g_gradient(i,2)
    !     ! vx_g_pbe_g(i,3) = vsigma(i)*2.0d0*nq_g_gradient(i,3)
    !     ! enddo

    !     ! deallocate(vsigma)
    !     ! deallocate(sigma)


    !   end subroutine libxc_exchange_lda



    !   subroutine libxc_correlation_lda(Ne,Nquadrature)
    !     integer,intent(in) :: Ne,Nquadrature
    !     TYPE(xc_f90_func_t) :: xc_func
    !     TYPE(xc_f90_func_info_t) :: xc_info
    !     ! TYPE(xc_f90_pointer_t) :: xc_func
    !     ! TYPE(xc_f90_pointer_t) :: xc_info
    !     integer :: i, vmajor, vminor, vmicro, func_id = 9 !!! GGA_LDA_C id is 9
    !     ! double precision,dimension(:),allocatable :: sigma
    !     ! double precision,dimension(:),allocatable :: vsigma
    !     integer(8) :: Ng

    !     Ng=Ne*Nquadrature

    !     vc_g_pbe_n=0.0d0
    !     ec_g_pbe=0.0d0

    !     ! allocate(vsigma(1:Ne*Nquadrature))
    !     ! allocate(sigma(1:Ne*Nquadrature))
    !     ! do i=1,Ng
    !     ! sigma(i) = nq_g_gradient(i,1)**2+nq_g_gradient(i,2)**2+nq_g_gradient(i,3)**2
    !     ! enddo

    !     ! call xc_f90_func_init(xc_func, xc_info, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     ! call xc_f90_func_init(xc_func, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
    !     call xc_f90_lda_vxc(xc_func, Ng, nq_g(1), vc_g_pbe_n(1))
    !     call xc_f90_lda_exc(xc_func, Ng, nq_g(1), ec_g_pbe(1))
    !     call xc_f90_func_end(xc_func)

    !     ! vc_g_pbe_n=vc_g_pbe_n*nq_g+ec_g_pbe

    !     ! do i=1,Ng
    !     ! vx_g_pbe_g(i,1) = vsigma(i)*2.0d0*nq_g_gradient(i,1)
    !     ! vx_g_pbe_g(i,2) = vsigma(i)*2.0d0*nq_g_gradient(i,2)
    !     ! vx_g_pbe_g(i,3) = vsigma(i)*2.0d0*nq_g_gradient(i,3)
    !     ! enddo

    !     ! deallocate(vsigma)
    !     ! deallocate(sigma)


    !   end subroutine libxc_correlation_lda






    !   subroutine libxc_exchange_g_pbe(Ne,Nquadrature)
    !     integer,intent(in) :: Ne,Nquadrature
    !     TYPE(xc_f90_func_t) :: xc_func
    !     TYPE(xc_f90_func_info_t) :: xc_info
    !     ! TYPE(xc_f90_pointer_t) :: xc_func
    !     ! TYPE(xc_f90_pointer_t) :: xc_info
    !     integer :: i, vmajor, vminor, vmicro, func_id = 101 !!! GGA_PBE_X id is 101
    !     double precision,dimension(:),allocatable :: sigma
    !     double precision,dimension(:),allocatable :: vsigma
    !     integer(8) :: Ng

    !     Ng=Ne*Nquadrature

    !     allocate(vsigma(1:Ne*Nquadrature))
    !     allocate(sigma(1:Ne*Nquadrature))
    !     vx_g_pbe_n=0.0d0
    !     vsigma=0.0d0
    !     ex_g_pbe=0.0d0

    !     do i=1,Ng
    !     sigma(i) = nq_g_gradient(i,1)**2+nq_g_gradient(i,2)**2+nq_g_gradient(i,3)**2
    !     enddo

    !     ! call xc_f90_func_init(xc_func, xc_info, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     ! call xc_f90_func_init(xc_func, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
    !     call xc_f90_gga_vxc(xc_func, Ng, nq_g(1), sigma(1), vx_g_pbe_n(1), vsigma(1))
    !     call xc_f90_gga_exc(xc_func, Ng, nq_g(1), sigma(1), ex_g_pbe(1))
    !     call xc_f90_func_end(xc_func)

    !     ! vx_g_pbe_n=vx_g_pbe_n*nq_g+ex_g_pbe
    !     ! vx_g_pbe_n=2.0d0*vx_g_pbe_n


    !     do i=1,Ng
    !     vx_g_pbe_g(i,1) = vsigma(i)*2.0d0*nq_g_gradient(i,1)!*nq_g(i)!*2.0d0
    !     vx_g_pbe_g(i,2) = vsigma(i)*2.0d0*nq_g_gradient(i,2)!*nq_g(i)!*2.0d0
    !     vx_g_pbe_g(i,3) = vsigma(i)*2.0d0*nq_g_gradient(i,3)!*nq_g(i)!*2.0d0
    !     enddo

    !     deallocate(vsigma)
    !     deallocate(sigma)


    !   end subroutine libxc_exchange_g_pbe



    !   subroutine libxc_correlation_g_pbe(Ne,Nquadrature)
    !     integer,intent(in) :: Ne,Nquadrature
    !     TYPE(xc_f90_func_t) :: xc_func
    !     TYPE(xc_f90_func_info_t) :: xc_info
    !     ! TYPE(xc_f90_pointer_t) :: xc_func
    !     ! TYPE(xc_f90_pointer_t) :: xc_info
    !     integer :: i, vmajor, vminor, vmicro, func_id = 130 !!! GGA_PBE_C id is 130
    !     double precision,dimension(:),allocatable :: sigma
    !     double precision,dimension(:),allocatable :: vsigma
    !     integer(8) :: Ng

    !     Ng=Ne*Nquadrature

        

    !     allocate(vsigma(1:Ne*Nquadrature))
    !     allocate(sigma(1:Ne*Nquadrature))
    !     vc_g_pbe_n=0.0d0
    !     vsigma=0.0d0
    !     ec_g_pbe=0.0d0
    !     do i=1,Ng
    !     sigma(i) = nq_g_gradient(i,1)**2+nq_g_gradient(i,2)**2+nq_g_gradient(i,3)**2
    !     enddo
    
    !     ! call xc_f90_func_init(xc_func, xc_info, XC_GGA_C_PBE, XC_UNPOLARIZED)
    !     ! call xc_f90_func_init(xc_func, XC_GGA_C_PBE, XC_UNPOLARIZED)
    !     call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
    !     call xc_f90_gga_vxc(xc_func, Ng, nq_g(1), sigma(1), vc_g_pbe_n(1), vsigma(1))
    !     call xc_f90_gga_exc(xc_func, Ng, nq_g(1), sigma(1), ec_g_pbe(1))
    !     call xc_f90_func_end(xc_func)

    !     ! vc_g_pbe_n=vc_g_pbe_n*nq_g+ec_g_pbe
    !     ! vc_g_pbe_n=2.0d0*vc_g_pbe_n


    !     do i=1,Ng
    !     vc_g_pbe_g(i,1) = vsigma(i)*2.0d0*nq_g_gradient(i,1)!*nq_g(i)!*2.0d0
    !     vc_g_pbe_g(i,2) = vsigma(i)*2.0d0*nq_g_gradient(i,2)!*nq_g(i)!*2.0d0
    !     vc_g_pbe_g(i,3) = vsigma(i)*2.0d0*nq_g_gradient(i,3)!*nq_g(i)!*2.0d0
    !     enddo

    !     deallocate(vsigma)
    !     deallocate(sigma)


    !   end subroutine libxc_correlation_g_pbe




    !   subroutine compute_nq_gradient(Nn,Ne,Nlocal,Nstates)

    !     integer,intent(in) :: Nn,Ne,Nlocal,Nstates
    !     integer :: i,m,n,k
    !     ! double precision,dimension(:,:),allocatable :: grad
    !     double precision,dimension(:), allocatable :: gradient_temp

    !     allocate(gradient_temp(1:3))

    !     ! allocate(nq_gradient(1:Nn,1:3))

    !     call local_nodes

    !     if (Nlocal==10) then

    !     do i=1,Ne
    !         call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    !         do m=1,Nlocal
    !             call basis_p2_del(localpointp2(m,:))
    !             nq_gradient(ele(i,m),:)=0.0d0
    !             ! do n=1,Nlocal
    !             !     do k=1,Nstates
    !             !     ! nq_gradient(ele(i,m),:)=nq_gradient(ele(i,m),:)+nq(ele(i,n))*matmul(Jit,grad(m,:))
    !             !     nq_gradient(ele(i,m),:)=nq_gradient(ele(i,m),:)+psi(ele(i,n),k)*matmul(Jit,grad(m,:))
    !             !     enddo
    !             ! enddo

    !             do k=1,Nstates
    !                 gradient_temp=0.0d0
    !             do n=1,Nlocal
    !                 gradient_temp=gradient_temp+matmul(Jit,phi_p2_del(n,:))*psi(ele(i,n),k)
    !             enddo
    
    !             nq_gradient(ele(i,m),:)&
    !                     =nq_gradient(ele(i,m),:)+4.0d0*psi(ele(i,m),k)*gradient_temp(:)
    !             enddo


    !         enddo
    !     enddo

    !     else if (Nlocal==20) then

    !         do i=1,Ne
    !             call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))
    !             do m=1,Nlocal
    !                 call basis_p3_del(localpointp2(m,:))
    !                 nq_gradient(ele(i,m),:)=0.0d0
    !                 ! do n=1,Nlocal
    !                 !     do k=1,Nstates
    !                 !     ! nq_gradient(ele(i,m),:)=nq_gradient(ele(i,m),:)+nq(ele(i,n))*matmul(Jit,grad(m,:))
    !                 !     nq_gradient(ele(i,m),:)=nq_gradient(ele(i,m),:)+psi(ele(i,n),k)*matmul(Jit,grad(m,:))
    !                 !     enddo
    !                 ! enddo
    
    !                 do k=1,Nstates
    !                     gradient_temp=0.0d0
    !                 do n=1,Nlocal
    !                     gradient_temp=gradient_temp+matmul(Jit,phi_p3_del(n,:))*psi(ele(i,n),k)
    !                 enddo
        
    !                 nq_gradient(ele(i,m),:)&
    !                         =nq_gradient(ele(i,m),:)+4.0d0*psi(ele(i,m),k)*gradient_temp(:)
    !                 enddo
    
    
    !             enddo
    !         enddo

    !     end if

    !     deallocate(gradient_temp)


    !   end subroutine compute_nq_gradient

    !   subroutine libxc_exchange_pbe(Nn)
    !     integer,intent(in) :: Nn
    !     TYPE(xc_f90_func_t) :: xc_func
    !     TYPE(xc_f90_func_info_t) :: xc_info
    !     ! TYPE(xc_f90_pointer_t) :: xc_func
    !     ! TYPE(xc_f90_pointer_t) :: xc_info
    !     integer :: i, vmajor, vminor, vmicro, func_id = 101 !!! GGA_PBE_X id is 101
    !     double precision,dimension(:),allocatable :: sigma
    !     double precision,dimension(:),allocatable :: vsigma
    !     integer(8) :: Nn0

    !     Nn0=real(Nn)

    !     allocate(vsigma(1:Nn))
    !     allocate(sigma(1:Nn))
    !     vx_pbe_n=0.0d0
    !     vsigma=0.0d0
    !     ex_pbe=0.0d0

    !     do i=1,Nn
    !     sigma(i) = nq_gradient(i,1)**2+nq_gradient(i,2)**2+nq_gradient(i,3)**2
    !     enddo

    !     ! call xc_f90_func_init(xc_func, xc_info, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     ! call xc_f90_func_init(xc_func, XC_GGA_X_PBE, XC_UNPOLARIZED)
    !     call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
    !     call xc_f90_gga_vxc(xc_func, Nn0, nq(1), sigma(1), vx_pbe_n(1), vsigma(1))
    !     call xc_f90_gga_exc(xc_func, Nn0, nq(1), sigma(1), ex_pbe(1))
    !     call xc_f90_func_end(xc_func)

    !     ! vx_g_pbe_n=vx_g_pbe_n*nq_g+ex_g_pbe
    !     ! vx_g_pbe_n=2.0d0*vx_g_pbe_n


    !     do i=1,Nn
    !     vx_pbe_g(i,1) = vsigma(i)*2.0d0*nq_gradient(i,1)!*nq_g(i)!*2.0d0
    !     vx_pbe_g(i,2) = vsigma(i)*2.0d0*nq_gradient(i,2)!*nq_g(i)!*2.0d0
    !     vx_pbe_g(i,3) = vsigma(i)*2.0d0*nq_gradient(i,3)!*nq_g(i)!*2.0d0
    !     enddo

    !     deallocate(vsigma)
    !     deallocate(sigma)


    !   end subroutine libxc_exchange_pbe


    !   subroutine libxc_correlation_pbe(Nn)
    !     integer,intent(in) :: Nn
    !     TYPE(xc_f90_func_t) :: xc_func
    !     TYPE(xc_f90_func_info_t) :: xc_info
    !     ! TYPE(xc_f90_pointer_t) :: xc_func
    !     ! TYPE(xc_f90_pointer_t) :: xc_info
    !     integer :: i, vmajor, vminor, vmicro, func_id = 130 !!! GGA_PBE_C id is 130
    !     double precision,dimension(:),allocatable :: sigma
    !     double precision,dimension(:),allocatable :: vsigma
    !     integer(8) :: Nn0

    !     Nn0=real(Nn)

        

    !     allocate(vsigma(1:Nn))
    !     allocate(sigma(1:Nn))
    !     vc_pbe_n=0.0d0
    !     vsigma=0.0d0
    !     ec_pbe=0.0d0
    !     do i=1,Nn
    !     sigma(i) = nq_gradient(i,1)**2+nq_gradient(i,2)**2+nq_gradient(i,3)**2
    !     enddo
    
    !     ! call xc_f90_func_init(xc_func, xc_info, XC_GGA_C_PBE, XC_UNPOLARIZED)
    !     ! call xc_f90_func_init(xc_func, XC_GGA_C_PBE, XC_UNPOLARIZED)
    !     call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
    !     call xc_f90_gga_vxc(xc_func, Nn0, nq(1), sigma(1), vc_pbe_n(1), vsigma(1))
    !     call xc_f90_gga_exc(xc_func, Nn0, nq(1), sigma(1), ec_pbe(1))
    !     call xc_f90_func_end(xc_func)

    !     ! vc_g_pbe_n=vc_g_pbe_n*nq_g+ec_g_pbe
    !     ! vc_g_pbe_n=2.0d0*vc_g_pbe_n


    !     do i=1,Nn
    !     vc_pbe_g(i,1) = vsigma(i)*2.0d0*nq_gradient(i,1)!*nq_g(i)!*2.0d0
    !     vc_pbe_g(i,2) = vsigma(i)*2.0d0*nq_gradient(i,2)!*nq_g(i)!*2.0d0
    !     vc_pbe_g(i,3) = vsigma(i)*2.0d0*nq_gradient(i,3)!*nq_g(i)!*2.0d0
    !     enddo

    !     deallocate(vsigma)
    !     deallocate(sigma)


    !   end subroutine libxc_correlation_pbe



    !   subroutine compute_nq_g_gradient(Nn,Ne,Nlocal,Nquadrature,Nstates)

    !     integer,intent(in) :: Nn,Ne,Nlocal,Nquadrature,Nstates
    !     integer :: i,ii,m,k
    !     double precision, dimension(:,:),allocatable :: grad_m
    !     double precision, dimension(:),allocatable :: basis_function,gradient_temp

    !     ! allocate(basis_function(1:Nlocal))
    !     allocate(gradient_temp(1:3))
    !     ! allocate(grad_m(1:Nlocal,1:3))


    !     nq_g_gradient=0.0d0

    !     ! call mkl_dcsrgemv('N',Nn,NnToNg,NnToNg_IA,NnToNg_JA,nq,nq_g)



    !     call gaussian_integral

    !     nq_g=0.0d0

    !     ! print *, nq_g(1)
    !     ! print *, nq_g_gradient(1,:)

    !     do i=1,Ne

    !         call jacobian(ele(i,1),ele(i,2),ele(i,3),ele(i,4))

    !     do ii=1,Nquadrature
            
    !         ! basis_function = basis_p2(gpoint(ii,:))

    !         do k=1,Nstates
    !         nq_g((i-1)*Nquadrature+ii)=nq_g((i-1)*Nquadrature+ii)+2.0d0*psi_point_g((i-1)*Nquadrature+ii,k)**2
    !         enddo

    !         call basis_p2_del(gpoint(ii,:))

    !         do k=1,Nstates
    !             gradient_temp=0.0d0
    !         do m=1,Nlocal
    !             gradient_temp=gradient_temp+matmul(Jit,phi_p2_del(m,:))*psi(ele(i,m),k)
    !         enddo

    !         nq_g_gradient((i-1)*Nquadrature+ii,:)&
    !                 =nq_g_gradient((i-1)*Nquadrature+ii,:)+4.0d0*psi_point_g((i-1)*Nquadrature+ii,k)*gradient_temp(:)
    !         enddo

    !     enddo
    !     enddo

    !     ! print *, nq_g(1)
    !     ! print *, nq_g_gradient(1,:)

    !     deallocate(gradient_temp)


    !   end subroutine compute_nq_g_gradient




      subroutine set_hf_BC(Nn,Nbc,Ne,Nquadrature,Nstates)

        integer,intent(in) :: Nn,Nbc,Ne,Nquadrature,Nstates
        integer :: i,j,jj,num1,num2,num3,index01,index02,index03
        logical :: ftest0

        vh_temp=0.0d0

        ! do i=1,Nbc
        !     vh_temp(BC(i))=-Vi(BC(i))
        ! enddo


        do i=1,Nbc
            vh_temp(BC(i))=poissonBC(point(BC(i),:),Ne,Nquadrature,Nstates)
        enddo


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! set boundary conditions !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        A00=A

        ! allocate(matBC(Nn,Nbc))
        matBC=0.0d0

        do i=1,Nbc
            num3=0
            do j=1,Nn
                num1=0
                num2=0
                ftest0=.true.
                if (neigh_AB(j)%find(BC(i)))  ftest0=.false.
                if (ftest0) then 
                    cycle
                else
                    index01=neigh_AB(j)%findindex(BC(i))
                    index02=neigh_AB(BC(i))%findindex(j)
                end if

                if (j==1) then 
                    num1=num1+index01
                else
                    do jj=1,j-1
                        num1=num1+neigh_AB(jj)%size
                    enddo
                    num1=num1+index01
                end if
                if (j/=BC(i)) matBC(j,i)=A00(num1)
                ! snq(j)=snq(j)-A00(num1)*vh_temp(BC(i))
                A00(num1)=cmplx(0.0d0,0.0d0,8)

                if (BC(i)==1) then
                    num2=num2+index02
                else
                    do jj=1,BC(i)-1
                        num2=num2+neigh_AB(jj)%size
                    enddo
                    num2=num2+index02
                end if
                A00(num2)=cmplx(0.0d0,0.0d0,8)
            enddo
            index03=neigh_AB(BC(i))%findindex(BC(i))
            if (BC(i)==1) then
                num3=num3+index03
            else
                do jj=1,BC(i)-1
                    num3=num3+neigh_AB(jj)%size
                enddo
                num3=num3+index03
            end if
            A00(num3)=cmplx(1.0d0,0.0d0,8)
            ! snq(BC(i))=vh_temp(BC(i))
        enddo


    end subroutine set_hf_BC

    function poissonBC(x,Ne,Nquadrature,state0) result(bc_value)

        double precision,dimension(1:3),intent(in) :: x
        integer,intent(in) :: Ne,Nquadrature,state0
        double precision,dimension(1:3) :: location
        double precision :: bc_value, r0,n_temp, density_check
        double precision :: psi_temp1, psi_temp2, psi_temp3, psi_temp4, psi_temp5!,f1,f2,f3,f4,psi_temp
        double precision,dimension(:),allocatable :: f1,f2,f3,f4
        integer :: i,j,ii
    
        
    
        density_check=0.0d0
        bc_value=0.0d0
        do i=1,Ne
            do j=1,Nquadrature
              
                n_temp=0.0d0
                ! if (rr0==0) then
                !     do ii=1,state0
                !         if (ii<=1) then
                !             n_temp=n_temp+psi_g01(i,j)%state_value(ii)
                !         else
                !             n_temp=n_temp+psi_g01(i,j)%state_value(ii)
                !         end if
                !     enddo
                ! else if (rr0==1) then
                    do ii=1,state0
                        if (ii<=2) then
                            n_temp=n_temp+2.0d0*psi_point_g((i-1)*Nquadrature+j,ii)**2
                        ! else if (ii==3) then
                        !     n_temp=n_temp+psi_g01(i,j)%state_value(ii)**2
                        else
                            n_temp=n_temp+2.0d0*psi_point_g((i-1)*Nquadrature+j,ii)**2
                        end if
                    enddo
                ! else if (rr0==2) then
                !     do ii=1,state0
                !         if (ii<=1) then
                !             n_temp=n_temp+2.0d0*psi_g01(i,j)%state_value(ii)**2
                !         ! else if (ii==7) then
                !         ! ! else
                !         !     n_temp=n_temp+14.0d0/15.0d0*2.0d0*psi_g01(i,j)%state_value(ii)**2
                !         ! else if (ii==8) then
                !             else
                !             n_temp=n_temp+psi_g01(i,j)%state_value(ii)**2
                !         end if
                !     enddo
                ! end if
    
    
                
    
    
                ! r0=sqrt((psi_g01(i,j)%location(1)-x(1))**2+(psi_g01(i,j)%location(2)-x(2))**2+(psi_g01(i,j)%location(3)-x(3))**2)

r0=sqrt((x(1)-point_g((i-1)*Nquadrature+j,1))**2+(x(2)-point_g((i-1)*Nquadrature+j,2))**2+(x(3)-point_g((i-1)*Nquadrature+j,3))**2)

                    bc_value=bc_value+volume(i)*gweight(j)*n_temp*1.0d0/r0
                
            enddo
        enddo
    
    
    end function poissonBC


    subroutine hf_hartreepotential(Nn,Nbc)

        integer,intent(in) :: Nn,Nbc
        ! double precision,dimension(:),intent(in) :: ni

        double precision,dimension(:),allocatable :: snq
        ! complex(8),dimension(:),allocatable :: snq_cmplx
        !complex(8),dimension(:),allocatable :: xy_AB0
        integer :: i!,j,jj,num1,num2,num3,index01,index02,index03

        !!! for mkl_dcsrgemv !!!
        character(len=1) :: transa='N'
        !logical :: ftest0
        !!! for pardiso !!!
        integer :: MAXFCT,MNUM,MTYPE,MSGLVL,PHASE,idum,error_hartree
        integer(8),dimension(:),allocatable :: pt_V
        integer,dimension(:),allocatable :: iparm_V
        
        complex(8) :: jc

        jc=cmplx(0.0d0,1.0d0,8)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!   solve Poisson equation w/ BC    !!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(snq(Nn))
        
        call mkl_dcsrgemv(transa, Nn, B, IA, JA, 4.0d0*pi*ni , snq)

        do i=1,Nbc
            snq(BC(i))=vh_temp(BC(i))
        enddo

        do i=1,Nbc
            snq(:)=snq(:)-matBC(:,i)*vh_temp(BC(i))
        enddo


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!   solve using Pardiso   !!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MAXFCT=1
        MNUM=1
        MTYPE=1     ! 1 for real and spd 3 for complex and sdp
        MSGLVL=0     !0- no output, 1- output
        allocate(pt_v(64))
        allocate(iparm_v(64))
        pt_V=0
        call pardisoinit(PT_V,MTYPE,IPARM_V)
        !IPARM(11,:)=0 ! disable scaling
        PHASE=13 ! 13


        call pardiso(pt_v,MAXFCT,MNUM,MTYPE,PHASE,Nn,A00,IA,JA,idum,1,IPARM_V,MSGLVL,snq,Vh,error_hartree)

        if (error_hartree/=0) then
            print *,'PARDISO_error info -----', error_hartree
            stop
        end if


    end subroutine hf_hartreepotential

      




end module potentials
