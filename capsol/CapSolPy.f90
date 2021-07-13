! This program calculates the capacitance of a  tip-sample system  with cylindrical symmetry.
! The tip may consist these parts: a spherical cap, a conical shank, a disk cantilever
! The sample is a dielectric with arbitrary thickness (including zero thickness) over a conducting electrode.
! sample input file is made automatically in the first run
! Author: Ali Sadeghi, Uni Basel, 2011

!==========================================

module defcapsol
    ! variables
          implicit none
          integer , parameter :: dp=kind(1.d0)  ! double precision
          real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
          real(dp), parameter :: e0=8.854187817620d-3 !vacuum permittivity in nN/V2
          real(dp), parameter :: pi_e0 = pi*e0 ! to make units for unitless quantities 
    
    ! tip and sample parameters
         !  real(dp)     :: Rtip,RCant,HCone,dCant,theta
         !  real(dp)     :: eps_r, Hsam 
    ! tip-sample separation: loops from d_min to d_max
          real(dp)     :: d_min,d_max,Z 
    
    ! simulation parameters
         !  real(dp)     :: h0,rho_max,Z_max   ! finest grid spacing and simulation box size (truncation lengths)
          real(dp)     :: qn,qm,ql   ! growth factors 
         !  integer      :: n,m,l ! number of grid points in rho and +/-z directions
          real(dp), allocatable  :: hn(:),r(:),hm(:),zm(:) ! grid-spacing and coordinates of each grid point
          real(dp), allocatable  :: u(:,:),g(:,:)  ! potential and its gradient
    
    ! Definition of the tip for Dirichlet BC
          integer , allocatable  :: ProbeBot(:), ProbeTop(:)   ! tip geometry: bottom and top
          real    , allocatable  :: b(:,:)  ! this applies BC: if b(j,i)=1 then it can vary, if =0 Dirichlet BC (fixed) 
          integer      ::  nApex,nCone,nLever,nEdge
          integer      ::  js
          integer      ::  i,j,k
          real(dp)     ::  Energy,gnrm,x,y
          real(dp)     :: Er,Ez,Force  ! electric field and force 
          logical      :: spheretest, leveronlytest ! for testing purpose
          
    contains
    
    subroutine GenerateGrid(n,m,l, h0, rho_max, Z_max, Rtip,theta,HCone,RCant, dCant, eps_r, Hsam, verbose)
    real*8, INTENT(IN) :: h0, rho_max, z_max, Rtip, theta, HCone, RCant, dCant, eps_r, Hsam
    integer, INTENT(IN) :: n, m, l, verbose  !  __ generates  a non-uniform grid  _______________
    integer Nuni ! number of points in the uniform region
    real(dp) q
    ! growth factors are determined according to h0, # of grid-points and box sizes
    ! Alternatively, one can determine # according to known growth factors, h0 and box sizes
    !   or even determine the box sizes according to known h0, growth factors and #
    ! We prefer to get h0 and box sizes as input parametersi. Then the maximum number
    ! of point-grids that works with the users specific machine-memory is used to find the smallest possible growth factors.
    
    
    Nuni=1  ! They used a larger Nuni before...
    do i=0,Nuni-1 
        r(i)=h0*i
        hn(i)=h0
    enddo
    
    r(Nuni)=h0*Nuni
    
    qn=1.012
    if(verbose>1) print*,'-----------Additional information:'
   !  if(verbose>1) print ('(a,$)'), " Suggested number of grid points (n,m+,m-): "
   !  if(verbose>1) print ('(i5,$)'),  nint(log(1_dp+(rho_max/h0-Nuni)*(qn-1_dp))/log(qn) + Nuni)
    ! find the growth factor
    do  qn=1._dp,1.5,1.d-4
      x=h0*(1-qn**(n-Nuni))/(1-qn)
      if(x>=rho_max - r(Nuni)) exit  ! found
    enddo
    
    hn(Nuni)=h0*(rho_max-r(Nuni))/x
    r(Nuni+1)=sum(hn(0:Nuni))
    do i=Nuni+2,n
      hn(i-1)=hn(i-2)*qn
      r(i)=sum(hn(0:i-1))
    enddo
    hn(n)=hn(n-1)*qn
    !		print*,r(n);  do i=0,n; write(333,*)i,r(i),hn(i); enddo
    
    ! Z+ direction
    Nuni=1 
    hm(1:Nuni)=h0
    do j=1,Nuni ; zm(j)=h0*j ; enddo
    
    q=1.010
      !   if(verbose>1) print('(i5,$)'), nint(log(1_dp+(Z_max/h0-Nuni)*(q-1_dp))/log(q)) + Nuni
    
    do qm=1._dp,1.5_dp, 1.d-4
      x=h0*(1-qm**(m-Nuni))/(1-qm)
      if(x>=Z_max-zm(Nuni)) exit
    enddo
    
    hm(Nuni+1)=h0 * (Z_max-zm(Nuni))/x
    zm(Nuni+1)=sum(hm(1:Nuni+1))
    do j=Nuni+2,m   
        hm(j)=hm(j-1)*qm
        zm(j)=sum(hm(1:j))
    enddo
    
    ! Z- direction
    hm(js:0) = h0 ! gap fiiled with h0
    do j=js,0; zm(j)=h0*j; enddo
    
    if (l+js>0) then 
       q=1.02
      !   if(verbose>1)   print'(i5,$)', nint(log(1_dp+(Hsam/h0)*(q-1_dp))/log(q))
       do ql=1._dp,2.0_dp,1.d-4
          x=h0*(1-ql**(l+js))/(1-ql)
          if(x>=Hsam) exit
       enddo
       hm(js)=h0*Hsam/x
    else 
        hm(js)=h0
    endif
    do j=js-1,-l,-1
        hm(j)=hm(j+1)*ql
        zm(j)=-sum(hm(j+1:0))
    enddo
    Z=-zm(js)
    !!	print*,l,m, zm(-l),zm(m);  do j=-l,m; write(333,*)j,zm(j),hm(j); enddo; stop
    !!do nLever=1,n ; if(r(nLever)>=RCant) exit; enddo;
    !!print '(a,F6.3,a)',  "      <<<<<<<<<<<<<<<<<<<<<<<<<<<<  Tip-Sample dist= " ,Z," nm >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    !print '(a,i10,4i5,a,3f8.4)',     " Number of points (TOTAL, n,m+,m-,gap)  "  &
    !       ,(n+1)*(m+l+1) ,n,m,l+js,-js , " ;    Growth factors rho, z+, z- = " ,qn,qm,ql
   !  if(verbose>1) then 
   !     print '(a,i3)',     "    Additional in gap:   "  ,-js  
   !     !print '(a, i8, a, f5.1)',"  >>> N_TOTAL: ",(n+1)*(m+l+1) , '   Approximate memory (GB): ', ((n+1)*(m+l+1))**1.5 * (8+2)*1.e-9
   !     print '(a, i8, a, f5.1, "GB")',"  >>>>>>>>>>>>    N_TOTAL: ",(n+1)*(m+l+1) , &
   !           ' ; Approximate memory needed: ', (n+1)*(m+l+1)*(min(n,m+l)+5)*8.e-9
   !    print '(a,3f7.4)', " Geometic increment factor for mesh spacin (r,z+,z-):" ,qn,qm,ql
   !  endif
    end subroutine GenerateGrid
    
    
    
    !================================================================== 
    
    subroutine SetupProbe(n,m,l, h0, Rtip,theta,HCone,RCant, dCant, eps_r, Hsam, verbose)
    real*8, INTENT(IN) :: h0, Rtip, theta, HCone, RCant, dCant, eps_r, Hsam
    integer, INTENT(IN) :: n, m, l, verbose
      ! Note: to keep the tip-geometry fixed, Cone Height is constant and it is the sample thickness that varies by separation
    real*8 Ra , Rc , Ez2  
    Ra=Rtip*(1._dp-sin(theta))
    Rc=Rtip*cos(theta)
    b(:,:) =1.0 !free-to-vary points
    ! Dirichlet BC for 
    b(-l,:)=0. ! backelectrode
    b(m,:)=0.  ! and top
    b(:,n)=0.  ! and side walls
    
    ! initial guess for potential: 0 everywhere else on tip
    u=0._dp 
    ! tip geometry (u(tip)=1, b(tip)=0)
    if (spheretest) then 
    do i=0, n-1
        x=r(i)
        if (x< Rtip) then  ! sphere
          nApex=i
          y= Rtip - sqrt(Rtip**2 - x**2)   
          do j=0,m ; if (zm(j) >=y ) exit; enddo ; ProbeBot(i)=j
          do j=1,m ; if (Zm(j)>= y+2*  sqrt(rtip**2 - x**2)) exit; enddo;  ProbeTop(i)=j
             b(ProbeBot(i) : ProbeTop(i),i)=0 
             u(ProbeBot(i) : ProbeTop(i),i)=1
       endif
    enddo
    
    elseif (leveronlytest) then 
    do i=0, n-1
        x=r(i)
        if (x< RCant) then  ! Only the cantilever
         nLever=i
         y= HCone  
          do j=0,m ; if (zm(j) >=y )        exit; enddo ; ProbeBot(i)=j
          do j=1,m ; if(Zm(j)>=HCone+dCant) exit; enddo;  ProbeTop(i)=j
             b(ProbeBot(i) : ProbeTop(i),i)=0 
             u(ProbeBot(i) : ProbeTop(i),i)=1
    
    
        elseif (x< RCant+dCant/2 ) then  ! cantelever rounded edge
          nEdge= i    
          y= HCone + dCant/2 - sqrt((dCant/2)**2 - (x-RCant)**2)   
          do j=0,m ; if (zm(j) >=y ) exit; enddo ; ProbeBot(i)=j
          do j=1,m ; if (Zm(j)>= y+2*  sqrt((dCant/2)**2 - (x-RCant)**2)) exit; enddo;  ProbeTop(i)=j
          !do j=1,m ; if(Zm(j)>=HCone+dCant) exit; enddo;  ProbeTop(i)=j
             b(ProbeBot(i) : ProbeTop(i),i)=0 
             u(ProbeBot(i) : ProbeTop(i),i)=1
       endif
    enddo
    
    
    else
     ! if it is not a sphere test, make tip geometry
    do j=1,m ; if(Zm(j)>=HCone+dCant) exit; enddo;  ProbeTop(:)=j
       
    !  make the tip geometry 
    do i=0, n-1
       x=r(i)
       if(x< Rc) then  ! sphere
         nApex=i
         y=Rtip-sqrt(Rtip**2-x**2)
         do j=0,m ; if (zm(j) >=y )        exit; enddo ; ProbeBot(i)=j
         do j=1,m ; if(Zm(j)>=HCone+dCant) exit; enddo;  ProbeTop(i)=j
             b(ProbeBot(i) : ProbeTop(i),i)=0 
             u(ProbeBot(i) : ProbeTop(i),i)=1
       
       elseif (x< (HCone-Ra)*tan(theta)+Rc )  then  ! cone
          nCone=i
          y= (x-Rc)/tan(theta)+Ra
          do j=0,m ; if (zm(j) >=y )        exit; enddo ; ProbeBot(i)=j
          do j=1,m ; if(Zm(j)>=HCone+dCant) exit; enddo;  ProbeTop(i)=j
             b(ProbeBot(i) : ProbeTop(i),i)=0 
             u(ProbeBot(i) : ProbeTop(i),i)=1
    
        elseif (x< RCant) then  ! lever
         nLever=i
         y= HCone  
          do j=0,m ; if (zm(j) >=y )        exit; enddo ; ProbeBot(i)=j
          do j=1,m ; if(Zm(j)>=HCone+dCant) exit; enddo;  ProbeTop(i)=j
             b(ProbeBot(i) : ProbeTop(i),i)=0 
             u(ProbeBot(i) : ProbeTop(i),i)=1
    
    
        elseif (x< RCant+dCant/2 ) then  ! cantelever rounded edge
          nEdge= i    
          y= HCone + dCant/2 - sqrt((dCant/2)**2 - (x-RCant)**2)   
          do j=0,m ; if (zm(j) >=y ) exit; enddo ; ProbeBot(i)=j
          do j=1,m ; if (Zm(j)>= y+2*  sqrt((dCant/2)**2 - (x-RCant)**2)) exit; enddo;  ProbeTop(i)=j
          !do j=1,m ; if(Zm(j)>=HCone+dCant) exit; enddo;  ProbeTop(i)=j
             b(ProbeBot(i) : ProbeTop(i),i)=0 
             u(ProbeBot(i) : ProbeTop(i),i)=1
       endif
                        
    enddo
    endif
    if(z<=d_min+h0/2) then
       open(600,file='ProbeGeometry.dat')
       write(600,*) " # rho, z,   i,  j,  code:  1=sphere, 2=cone, 3=cant, 4= edge"
       Write(600,55)( r(i),zm(ProbeBot(i)),i,ProbeBot(i),1 , i=1        ,nApex )
       Write(600,55)( r(i),zm(ProbeBot(i)),i,ProbeBot(i),2 , i=nApex+1  ,nCone )
       Write(600,55)( r(i),zm(ProbeBot(i)),i,ProbeBot(i),3 , i=nCone+1  ,nLever)
       Write(600,55)( r(i),zm(ProbeBot(i)),i,ProbeBot(i),4 , i=nLever+1 ,nEdge )
       Write(600,55)( r(i),zm(ProbeTop(i)),i,ProbeTop(i),4 , i=nEdge, nLever+1,-1)
       Write(600,55)( r(i),zm(ProbeTop(i)),i,ProbeTop(i),3 , i=nLever,nCone+1 ,-1)
       Write(600,55)( r(i),zm(ProbeTop(i)),i,ProbeTop(i),2 , i=nCone ,nApex+1 ,-1)
       Write(600,55)( r(i),zm(ProbeTop(i)),i,ProbeTop(i),1 , i=nApex      ,1 ,-1)
       55 format(2e15.5 , 3i6)
       close(600)
    endif
    
    if(verbose>1) then 
    If(SphereTest) then
       print'(" *** Test: only sphere, Rtip:",f12.1)', r(nApex)
    elseIf(leveronlyTest) then
       print'(" *** Test: only cantilever, Rdisk:",f12.1)', r(nLever)
    else
       print'(" Radial extentions cap, cone and  disk:", 3f12.1)', r(nApex),r(nCone),r(nLever)
    endif
    print*
    print*, ' Calculating ....'
    print*
    print*
    endif
    end subroutine  SetupProbe ! __ Geometry of Probe   ___________________
    !=====================================
    
   !  subroutine WriteOutputs(k, eps_r)
   !  integer k
   !  real*8, INTENT(IN) :: eps_r
   !  real*8 Er,Ez,Ezbelow,E2,dFz,Fz
   !  50 format (2e13.3,10e20.10)
   !  51 format (2a13,10a20)
   !  print '(a,i3.3)' ," Writing pot & field in fort.?",k
   !  !write (k+1000,*), "#  " , n, m, l , "  #  n , m , l"
   !  !write (k+1000,51), "# rho(nm)  " , "z(nm)  "  ,  "pot(V)", "E2(1/nm)2 " , "E_rho(1/nm)   " , "E_z (1/nm)  " , "grad"
   !  !!do i=1,n-1
   !  !!  do j=1-l,m-1
   !  !do i=0,n
   !  !   do j=-l,m
   !  !      if(i==0) then 
   !  !            Er =-(u(j,i+1)-u(j,i))/hn(i)
   !  !      elseif(i==n) then 
   !  !            Er = -(u(j,i)-u(j,i-1))/hn(i-1)
   !  !      else
   !  !            Er=-(u(j,i+1)-u(j,i-1))/(hn(i)+hn(i-1)) 
   !  !      endif 
   !  !      if(j==-l) then 
   !  !            Ez =-(u(j+1,i)-u(j,i))/hm(j+1)
   !  !      elseif(j==m .or. j==0 ) then 
   !  !            Ez = -(u(j,i)-u(j-1,i))/hm(j)
   !  !      else
   !  !            Ez=-(u(j+1,i)-u(j-1,i))/(hm(j)+hm(j+1))
   !  !      endif
   !  !      E2=Ez**2+Er**2
   !  !!     write(k+1000,50)r(i),zm(j),u(j,i),g(j,i), E2 ! ,Er,Ez
   !  !     write(k+1000,50)r(i),zm(j),u(j,i) , E2 , Er, Ez
   !  !  enddo
   !  !  write(k+1000,*)
   !  !enddo
   !  !close(k+1000)
   !  !do j=-l,m
   !  !enddo
    
   !  ! force
   !  write (k+2000,*), "#  " , n, m, l , "  #  n , m , l"
   !  write(k+2000,51) " #  rho ", "  F_z (nN/V2) " , " dF_z(nN/V2) " , "E2 (1/nm)2 " , "E_rho   ","  E_z    " ,  " F_z /pi eps_0 V2 "
   !  Fz=0._dp
   !  do i=0, max( nApex,nCone,nLever,nEdge)-5
   !    j= ProbeBot(i)     
   !    Er=-(u(j,i+1)-u(j,i))/hn(i)
   !    Ez=-(u(j,i)-u(j-1,i))/hm(j)
   !    E2=Er**2+Ez**2
   !    x=(r(i)+r(i+1))/2._dp
   !    dFz=E2*x*hn(i)
   !    Fz=Fz+dFz
   !   write(k+2000,50)x,Fz*pi_e0,dFz*pi_e0 ,E2 ,Er,Ez , Fz
   !  enddo
   !  close(k+2000) 
   !  Force=Fz
   !  ! Print '(a,f8.4,e)', '# Z, Force' , Z, Force
   !  ! potential and field on the tip axis: (rho=0,z)
   !  write (3000+k,51) "#z(@rho=0) ","    pot.      ", " E_rho  " , " E_z(1/nm)  "
   !  do j=1-l,m-1
   !    Er=-(u(j,1)-u(j,0))/hn(0)
   !    Ez=-(u(j+1,0)-u(j-1,0))/(hm(j)+hm(j+1))
   !  write(3000+k,50) zm(j),u(j,0),Er,Ez
   !  enddo
   !  close(3000+k)
    
   !  ! Potential and E-field on the interface: (rho,z=0)
   !  write (4000+k,51) "  # rho(@z=0) ," , "    pot. " , "   E_rho  " , " E_z (top) " , "E_z(belove)"     
   !  do i=0,n-1 
   !    if(i==0) then 
   !       Er=-(u(js,1)-u(js,0))/hn(0)
   !    else 
   !      Er=-(u(js,i+1)-u(js,i-1))/(hn(i-1)+hn(i))
   !    endif
   !    Ez=-(u(js+1,i)-u(js,i))/hm(1)
   !    Ezbelow=-(u(js,i)-u(-1,i))/hm(0)
   !  write(4000+k,50) r(i),u(js,i),Er,Ez,-Ezbelow , -eps_r*Ezbelow ! two last should be the same! checjk the output
   !  enddo
   !  close(4000+k)
   !  end subroutine WriteOutputs
    
    !===============================================
    subroutine IntegratedForce(n,m,l, Z, Rtip,theta,HCone,RCant, dCant, eps_r, Hsam)
    real*8, INTENT(IN) :: Z, Rtip, theta, HCone, RCant, dCant, eps_r, Hsam
    integer, INTENT(IN) :: n, m, l
    real*8 E2,df
    character*200 str
    Force=0._dp
    write(str,'(a,f6.4)')"Fz.dat", Z/Rtip
    open (500, file=str)
    write(500,*) "# rho, F_z,dF_z, z(rho), E2 (=>sigma2=e0^2*E2), Er, Ez, Ui,j , Ui+1,j , Ui,j-1"
    do i=0,max(nApex,nCone,nLever,nEdge)
      !j=Probebot(i)    
      !Er=-(u(j,i+1)-u(j,i))/hn(i)
      !Ez=-(u(j,i)-u(j-1,i))/hm(j)
      !Er=-(u(j,i+2)-u(j,i))/(hn(i)+hn(i+1))
      !Ez=-(u(j,i)-u(j-2,i))/(hm(j)+hm(j-1))
      j=Probebot(i) -1 ! avoid stepwise boundariy 
      Er=-(u(j,i+1)-u(j,i))/hn(i)
      Ez=-(u(j,i)-u(j-1,i))/hm(j)
      E2=Er**2+Ez**2
      df=.5_dp*(r(i)+r(i+1)) * hn(i) * E2 ! r_bar * dr * E^2 (integrating a volume element)
      Force=Force+df
      write(500,'(18es15.5)'),.5_dp*(r(i)+r(i+1)) ,Force,dF,zm(Probebot(i)),E2,er,ez,u(j,i),u(j,i+1),u(j-1,i)
      !write(500,'(8es15.5)'),1.*i,1.*j,x,u(j,i),u(j,i+1),u(j-1,i)
    enddo
    ! note that in cylinderical coordinates: dS.n=(2pi r dr)/sin t . sin t =2pi r dr
    close(500)
    end subroutine IntegratedForce
    
    
    Subroutine get_time (time)
         integer v(8)
         character*50 time 
         call date_and_time(VALUES=v)
         write ( time, "(' Time: ',i2.2,'.',i2.2,'.',i4.4,i6.2,':',i2.2,':',i2.2)") , v(3),v(2),v(1),v(5:7)
    end subroutine
end module defcapsol
    
    !=====================================================
    subroutine CapSolCyl(n,m,l_js, h0,rho_max,z_max, d_steps, Rtip,theta_deg,HCone,RCant, &
                        dCant, &
                          eps_r,Hsam, Method_in, test, verbose, C_unitless, Hess)   ! Capicatnce Solver with Cylindrical Symmetry 
    use defcapsol
    implicit none 
    integer, INTENT(IN) :: n, m, l_js, test, verbose, d_steps
    real*8, INTENT(IN) :: h0, rho_max, z_max, Rtip, theta_deg
    real*8, INTENT(IN) :: HCone, RCant, dCant, eps_r, Hsam
    character*10, INTENT(IN)         :: Method_in
    real*8, INTENT(OUT) :: C_unitless
    real*8, dimension(0:min(m+ l_js + d_steps + 1, n+1), (n+1)*( l_js + d_steps + 1)), INTENT(OUT) :: Hess ! as general band storage

    character*10 :: Method
   real*8 :: Z_unitless, theta
    integer*8 iter
    real(8) gnrm0 
    integer l, ij, info ,lmn ,kd,kdA,kdB
    real*8 , allocatable :: gp(:,:)
    integer, allocatable :: ipiv(:)
    real*8 , allocatable :: HA1(:),HA2(:),HA3(:)
    real*8 , allocatable :: HB1(:),HB2(:),HB3(:)
    real*8 , allocatable :: Hd_1(:)
    real*8 , allocatable :: aa(:,:) , bb(:,:)
    real*8 , allocatable :: work(:)
    integer, allocatable :: iwork(:)
    real*8               :: HessNorm,rcond
    real t  
    real*8 Alpha, anoise , energy1  
    real*8 arrayC(4),zk
    real*8 , allocatable :: arrayE(:),arrayZ(:),arrayF(:)
    integer iarray, idstep, ifn
    character*25 fn
    logical output
    
     INTEGER*8 pt(64)
     integer  phase
     INTEGER  idum, solver
     INTEGER maxfct, mnum, mtype,  nrhs,  msglvl 
     INTEGER iparm(64)
     REAL*8  dparm(64) 
     INTEGER  , allocatable :: ia(:)
     INTEGER  , allocatable :: ja(:)
     REAL*8   , allocatable :: a(:) 
     REAL*8  waltime1, waltime2, ddum
     character*50  :: str
    
    Alpha=1_dp

    if (test==1) spheretest=.true.
    if (test==2) leveronlytest=.true.
    if(verbose>2) output=.true.
    
    theta=theta_deg*pi/180._dp
    Method=trim(Method_in)
    
    Z=d_steps * h0

    if (method=='NOSOLVE') then
       print*,'No minimization will be done ! ' 
      goto 1000 
    endif

    
    ! j of the surface; j) apex=0, j_surface=-separation
     js=-nint(Z/h0); l=l_js+(-js)
     lmn=(n+1)*(l+m+1) ; kdA=m+l+1;  kdB=n+1 ; kd=min(kdA,kdB)
    
     allocate (hn(0:n),r(0:n),hm(-l:m),zm(-l:m)) 
     call GenerateGrid(n,m,l, h0, rho_max, Z_max, Rtip,theta,HCone,RCant, dCant, eps_r, Hsam, verbose)   
    
     allocate (u(-l:m,0:n),g(-l:m,0:n),b(-l:m,0:n) ) ! uexact , u0
     allocate (ProbeBot(0:n) , ProbeTop(0:n))
     call SetupProbe(n,m,l, h0, Rtip,theta,HCone,RCant, dCant, eps_r, Hsam, verbose)
    
    
     allocate (aa(-l:m,0:n),bb(-l:m,0:n),gp(0:n,-l:m)) 
     allocate (HA1(lmn),HA2(lmn),HA3(lmn), Hd_1(lmn), HB1(lmn),HB2(lmn),HB3(lmn) , ipiv(lmn))
    
   !  if (Method ==  'LAPACK') 
    if (Method == 'PARDISO') allocate (ia(lmn+1),ja(3*lmn-kdA-1),a(3*lmn-kdA-1) )
    
          
    
    !if(abs(Z-d_min)<h0/2)
    ! print '(a,3F12.6,e12.5)'," rmax,Zmax,Hsam(mm),eps_r: " ,r(n)/1e6 , zm(m)/1e6 , (zm(js)-zm(-l))/1e6, eps_r
    
    !  constructing constant cooficients  aa ,bb. => they depend only  on geometry and Z
    do i=0,n-1
    do j=js+1,m
      aa(j,i)=.5_dp*(r(i)+r(i+1)) *hn(i) / hm(j)  ! r_bar * dr / dz
      bb(j,i)=.5_dp*(r(i)+r(i+1)) /hn(i) * hm(j)  ! r_bar * dz / dr
    enddo
    do j=-l+1,js
      aa(j,i)=.5_dp*(r(i)+r(i+1)) *hn(i) / hm(j)*eps_r ! r_bar * dr / dz * eps_r 
      bb(j,i)=.5_dp*(r(i)+r(i+1)) /hn(i) * hm(j)*eps_r
    enddo
    enddo
    
    
    do i=0,n ;do j=-l,m ;  ij=j+l+1+i*kdA
       if(b(j,i)>0._dp) then 
          HA1(ij)= aa(j,i) +bb(j,i)
          if(i>0) HA1(ij)= HA1(ij)+bb(j,i-1)
          if(j<m) HA1(ij)=HA1(ij)+ aa(j+1,i)  
          HA2(ij)=0._dp
          if(j<m) HA2(ij)=-aa(j+1,i)
          HA3(ij)=-bb(j,i)
      else ! on boundry ==> fixed potentional nodes
         HA1(ij)=1._dp  !on diagonal
         if(ij<lmn)    HA2(ij)      = 0._dp ! just belove diagonal
         if(ij>1)      HA2(ij-1)    = 0._dp ! just left of diagonal
         if(ij<lmn-kdA) HA3(ij)     = 0._dp ! far belove diagonal
         if(ij>kdA)     HA3(ij-kdA) = 0._dp ! for left of diagonal
      endif
    enddo; enddo
    
    
    do j=-l,m ;do i=0,n ; ij=i+1+(j+l)*kdB
       if(b(j,i)>0._dp) then ! Variable coefficient...
          HB1(ij)= aa(j,i) + bb(j,i) ! 
          if(i>0) then
             HB1(ij) = HB1(ij) + bb(j,i-1)
          endif
          if(j<m) then
             HB1(ij) = HB1(ij)+ aa(j+1,i)
          endif
          HB3(ij)=0._dp; if(j<m)  HB3(ij)=-aa(j+1,i)
          HB2(ij)=-bb(j,i)
          HB2(ij)=-bb(j,i)
      else ! on boundry ==> fixed potentional nodes
         HB1(ij)=1._dp  !on diagonal
        if(ij<lmn)      HB2(ij)    = 0._dp ! just belove diagonal
        if(ij>1)        HB2(ij-1)  = 0._dp ! just left of diagonal
         if(ij<lmn-kdB) HB3(ij)    = 0._dp ! far belove diagonal
         if(ij>kdB)     HB3(ij-kdB)= 0._dp ! for left of diagonal
      endif
    enddo; enddo
    ! ____ Hessian ___________
    if (Method=='LAPACK') then
       if(kd==kdB) then 
         Hess=0._dp 
         Hess(0,:) = HB1  ! + 1.d-8
         Hess(1,:) = HB2 
         Hess(kd,:)= HB3
       else
         Hess=0._dp 
         Hess(0,:) = HA1  ! + 1.d-8
         Hess(1,:) = HA2 
         Hess(kd,:)= HA3
       endif
       HessNorm=0._dp
       do i=1,lmn   
             x=abs(Hess(0,i))+abs(Hess(1,i))+abs(Hess(kd,i))
             if(i>1) x=x+abs(Hess(1,i-1))
             if(i>kd) x=x+abs(Hess(kd,i-kd))
             if(x>HessNorm) HessNorm=x
       enddo
      !  if(verbose>1) then
      !    print *, "LBOUND is ", LBOUND(Hess), " UBOUND ", UBOUND(Hess)
      !  endif

         print *, 'Matrix factorization  by LAPACK ...' 
         call dpbtrf( 'L', lmn,kd,Hess,kd+1,info)
         if(info /= 0) then
            print*, info , 'Failure in factorization Hessain matrix ';  stop
         endif

   !  print*, ' Estimating condition number ...' 
   !    allocate (work(3*lmn) , iwork(lmn))  ! as general band storage
   !  call dpbcon('L',lmn,kd,Hess,kd+1,HessNorm,rcond,work,iwork,info)
   !    if(info /= 0) then
   !       print*, info , 'Failure in estimating condition number'; stop
   !    else 
   !       print '(a,2e10.1)',' Condition number of Full Hessian= ' , 1._dp/rcond , HessNorm
   !    endif
   !  deallocate (  work , iwork)  ! as general band storage
    
    ! _ Hessian ___________
    
    elseif  (Method == 'PARDISO') then
           ij=0
           do k=1,lmn 
             ij=ij+1; 
             ia(k)=ij ! ia(k) = k, so the first are 
             ja(ij)=k ! Same...
             a(ij)=HA1(k) ! 
           
           if(k<lmn) then 
             ij=ij+1
             ja(ij)=k+1
             a(ij)=HA2(k)
           endif
         
           if(k<= lmn -kdA) then 
           ij=ij+1
           ja(ij)=k+kdA
            a(ij)=HA3(k)
           endif
         enddo
         ia(lmn+1)=ij+1 
         mtype =2  ! 2 for symmetric positive definite mantrix  , -2 for symmetrtic indefinite 
         solver=  0  ! use sparse direct method
         maxfct=1
         mnum=1
         nrhs=1
         call pardisoinit(pt, mtype, iparm)
       !   print*, 'pt is ', pt
          iparm(6) = 1  !  0=> output on x  , 1 => output on b
          ! iparm(3) = 1   !  .. Numbers of Processors ( value of OMP_NUM_THREADS )
          iparm(8)  = 1   ! max numbers of iterative refinement steps
          msglvl    = 1      ! with statistical information
         print *, 'Matrix factorization ... by PARDISO' 
          phase    = 13  ! Analysis, Numerical Factorization, Solve, Iterative Refinent
          call pardiso (pt, maxfct, mnum, mtype, phase, lmn, a, ia, ja,  idum, nrhs, iparm, msglvl, g, gp, info) 
       if(info /= 0) then
          print*, info , 'Failure in factorization a matrix by PARDISO' , info  ;  stop
       endif
    
    ! Since now, a(1:3lmn-kd-1) in factorized form
    
    else 
       print*,'Not recognizing Method, Exit...' ; stop
    endif
    
    
    !______________________________ minimization 
    
    !write(*,'(a10,a30,2a15,a12,a)') "Itration " ,"Energy (nm.nN/V2) "," Gnrm   " , " Gnrm/Gnrm0  ","Del Energy ", " time(min)"
    iter=-1
    energy1=1.d100
    do  
      iter=iter+1 
    
      !__________ calculate energy and its gradient ________________
    energy=0._dp
    do i=0,n-1; do j=-l+1,m
        energy=energy+aa(j,i)*(u(j,i)-u(j-1,i))**2 + bb(j,i)*(u(j,i+1)-u(j,i))**2   
    enddo; enddo
    
    gnrm=0._dp
    i=0; do j=-l+1,m-1 ! on j=m pot is fixed
        g(j,i)=b(j,i)*(                     - u(j-1,i)*aa(j,i) + & 
               u(j,i)*(            aa(j,i) + aa(j+1,i) + bb(j,i)) -u(j+1,i)*aa(j+1,i)-u(j,i+1)*bb(j,i)) 
        gnrm=gnrm +g(j,i)**2
     enddo
    do i=1,n-1; do j=-l+1,m-1  ! j=m ==> boundry
        g(j,i)=b(j,i)* ( -u(j,i-1)*bb(j,i-1) - u(j-1,i)*aa(j,i) + & 
               u(j,i)*(bb(j,i-1) + aa(j,i) + aa(j+1,i) + bb(j,i)) -u(j+1,i)*aa(j+1,i)-u(j,i+1)*bb(j,i)) 
        gnrm=gnrm +g(j,i)**2
    enddo; enddo
    gnrm=sqrt(gnrm)
    
    
    IF (Method=='LAPACK') THEN
    !_____________________________________________________________________________________
      if (iter<2) gnrm0=gnrm
        
       if(kd==kdB) then ! use B
           do i=0,n; do j=-l,m ;  gp(i,j)=g(j,i); enddo;enddo
           call dpbtrs('L',lmn,kd,1,Hess,kd+1,gp, lmn,info)
           do i=0,n; do j=-l,m ;  g(j,i)=gp(i,j); enddo;enddo
       else ! use A
           call dpbtrs('L',lmn,kd,1,Hess,kd+1,g, lmn,info)
       endif
        if(info /=0) then 
             print* , info , "Fail in dgbtrs" ; stop
        endif
        do i=0, n-1   
                do j=-l,m
                      u(j,i)=u(j,i)-Alpha*g(j,i);   ! for Whole Hessian H
               enddo
        enddo
    !_____________________________________________________________________________________
    ELSEIF (Method=='PARDISO') THEN 
      if (iter<2) gnrm0=gnrm
          phase    = 33  ! Solve, Iterative Refinent
        CALL pardiso (pt, maxfct, mnum, mtype, phase, lmn, a, ia, ja,  idum, nrhs, iparm, msglvl, g, gp, info) 
          
       if(info  /= 0) then
          print*, info , 'Failure in Solving a by PARDISO' , info  ;  stop
       endif
             do i=0, n-1
                do j=-l,m
                       u(j,i)=u(j,i)-Alpha*g(j,i); 
                enddo
             enddo
    
    !_____________________________________________________________________________________
    ELSEIF (Method=='CROSS') THEN ! use tridiag. H
      if (iter<100) gnrm0=gnrm
      if(mod(iter,2)/=0) then  ! A for Crossing 
      !if(mod(iter,2)==0) then  ! A for Crossing 
            call dpttrs (lmn, 1, HA1, HA2, g, lmn, info)  ! Now g is not gradient, but is P=H-1 g
      else                     ! B for Crossing
          do i=0,n; do j=-l,m  ;   gp(i,j)=g(j,i);   enddo;enddo
              call dpttrs (lmn, 1, HB1, HB2, gp, lmn, info)  ! Now g is not gradient, but is P=H-1 g
          do i=0,n; do j=-l,m ;  g(j,i)=gp(i,j) ;enddo;enddo
       endif
        
       if(info /=0) then 
             print* , info , 'Failure in dpttrs'  ; stop
       endif
    
        do i=0, n-1  
                do j=-l,m
                       u(j,i)=u(j,i)-Alpha*g(j,i);   ! for tridiagonal H
               enddo
        enddo
    !_____________________________________________________________________________________
         
     
    ENDIF   ! METHOD
    
    if ( mod(iter,1)==0) then
    !      call cpu_time(t) ;  print '(i10,e,e10.3,f10.2)', Iter,energy*pi_e0,gnrm !, t/60
    endif
      if (gnrm <1.d-5*gnrm0) exit ! Converge    
      if (energy1-energy<anoise) exit 
                    
      energy1=energy
    ! print*,method, iter, gnrm,energy
    enddo  ! next iteratin  
    call cpu_time(t)
    
    ! _________End of minimization _____________________________________________
    
    
    call IntegratedForce(n,m,l, Z, Rtip,theta,HCone,RCant, dCant, eps_r, Hsam)
    !write(* ,'(a,f5.2,3i5,2e,f10.2 , i3 )')    ' FINAL Z,n,m,l,Energy & Force, cpu time: , Iter',Z,n,m,l,energy*pi_e0,Force*pi_e0 , t/60 , Iter
   !  write(* ,'(a,f9.4,e15.7 )')    ' ===S/R,C/pi.e0.R:',Z/Rtip,2_dp*energy/Rtip 
   !  write(* ,'(60x,a,f5.2,e15.7,f10.2 )')    ' +++Separation,Capacitance, cpu time: ',Z,2_dp*energy*pi_e0 ,  t
   !  !write(200 ,'(f10.5,2es )')     Z,energy*2*pi_e0  ! C=2U/V2 and V=1
   !  write(201 ,'(f10.5,2es15.7 )')     Z/Rtip,energy/Rtip,  2._dp*energy/Rtip 
    
   !  j=js/2; Ez=-(u(j+1,0)-u(j-1,0))/(hm(j)+hm(j+1))
   !  write(300,'(f10.4,E13.3, f6.2,E16.6,f8.2, E16.6, a )')   &
   !         z/Rtip,Ez*Rtip, z,Ez,z*10._dp ,Ez/27.211396132_dp*0.052917725_dp   , ' # s, E_mid: R,1/R;   nm, 1/nm; angs, Ha/Bohr'  
    
   !  write(301,*) 
   !  write(301,*) "# along tip axis: rho=0, s=",z
   !  !do j=js-3,3
   !  do j=js-2,2
   !        Ez=-(u(j+1,0)-u(j-1,0))/(hm(j)+hm(j+1))
   !        write(301,'(f6.2,E16.6, f10.4,E13.3,a)') (zm(j)-zm(js)) ,Ez ,(zm(j)-zm(js))/Rtip,Ez*Rtip &
   !            , ' # z (nm) , Ez(1/nm)      , z/R, Ez/(V/R) ' ! 
   !  enddo
    Z_unitless = Z/Rtip
   !  E_unitless = energy/Rtip
    C_unitless = 2._dp*energy/Rtip 


    if(Method=='PARDISO') then
       phase = -1 
     CALL pardiso (pt, maxfct, mnum, mtype, phase, lmn, a, ia, ja,  idum, nrhs, iparm, msglvl, g,gp ,  info) 
    endif
    
     deallocate (hn,r,hm,zm) 
     deallocate (u,g,b ) ! uexact , u0
     deallocate (aa,bb,gp) 
     deallocate (ProbeBot,ProbeTop )
     deallocate (HA1,HA2,HA3, Hd_1, HB1,HB2,HB3, ipiv)
    
   !  if (Method ==  'LAPACK') 
    if (Method == 'PARDISO') deallocate (ia,ja,a)
    
    1000 continue ! end of minimization
    end
    