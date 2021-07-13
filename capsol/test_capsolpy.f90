program testcap
    use defcapsol
    implicit none
    integer :: n, m, l_js, verbose, d_steps, test
    real(dp) :: h0, rho_max, Z_max, rtip, theta_deg, hcone, rcant, hsam, dcant, eps_r
    real(dp) :: C_unitless   ! finest grid spacing and simulation box size (truncation lengths)
    
    character*10 method_in

    real(dp), allocatable :: hess(:, :) ! grid-spacing and coordinates of each grid point



    ! Nuni = 5
    n=500
    m=500
    l_js=20
    h0=0.5_dp
    rho_max=1d6
    Z_max=1d6
    d_steps=20
    rtip=20.0
    theta_deg=15.0
    hcone=15000.0
    rcant=15000.0
    dcant=500.0
    eps_r=5.0
    hsam=10.0
    method_in='LAPACK'
    test=0
    verbose=2
    C_unitless = 0.0_dp


    allocate (Hess(0:min(m+ l_js + d_steps + 1, n+1), (n+1)*( l_js + d_steps  + m + 1)))

    CALL CapSolCyl(n,m,l_js, h0,rho_max,z_max, d_steps, rtip,theta_deg, hcone, rcant, &
    dcant, eps_r,hsam, method_in, test, verbose, C_unitless, Hess)

    deallocate (Hess)
end program testcap