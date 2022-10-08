Program potwell 
    IMPLICIT NONE 
    Integer,parameter :: dp=selected_real_kind(14,200)
    real(dp),parameter :: pi=3.14159265358979_dp
    Integer :: n,npw
    real(dp) :: v0,a,b,d,g
    real(dp),allocatable :: kn(:),e(:),h(:,:),work(:)
    real(dp) :: x,dx,norm,prob
    complex(dp) :: f 
    Integer :: i,j,nr,lwork,info 


    write(*,"('Parameters for potential well: V_0, b >',$)")
    read(*,*) v0,b
    if(v0<=0.0_dp .or. b<=0.0_dp) stop 'wrong input parameters'
    write(*,"(' V_0, b = ',2f10.4)") v0, b

10  write(*,"('Parameters for plane waves: a, n >',$)")
    read(*,*) a,n 
    if(a<=b .or. n<=0) stop 'wrong input parameters'
    write(*,"('a,n=',f8.4,i6)") a,n 
    npw=2*n+1
    allocate(kn(npw),e(npw),work(3*npw),h(npw,npw))

    kn(1)=0.0_dp 
    do i=2,npw-1,2
        kn(i)=(i/2)*2.0_dp*pi/a 
        kn(i+1)=-(i/2)*2.0_dp*pi/a
    end do 

    h(:,:)=0.0_dp 
    do i=1,npw
        do j=1,npw 
            if(i==j) then 
                h(i,j) = kn(i)**2-v0/a*b
            else 
                h(i,j)=-v0/a*sin((kn(j)-kn(i))*b/2.0_dp)/(kn(j)-kn(i))*2.0_dp
            end if 
        end do 
    end do 

    lwork=3*npw 
    call dsyev ('V','U',npw,h,npw,e,work,lwork,info)
    if(info /= 0) stop 'H-matrix diagonalization failed'

    write(*,"(' Lowest eigenvalues: ',3f12.6)") (e(i),i=1,3)

    open(7,file='gs_wfc.out',status='unknown',form='formatted')
    dx=0.01_dp
    nr=nint(a/2.0_dp/dx)
    norm=0.d0 
    do i=-nr,nr
        x=dx*i
        f=0.d0 
        do j=1,npw 
            f=f+h(j,i)*exp((0.0,1.0)*kn(j)*x)/sqrt(a)
        end do 
        prob=f*conjg(f)
        norm=norm+prob*dx 
        write(7,*)x,prob,f 
    end do 

    write(*,"(' norm: ',f12.6)") norm 
    close(7)
    deallocate(h,work,e,kn)
    go to 10 

End program potwell