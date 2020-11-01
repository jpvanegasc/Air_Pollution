module constants
    implicit none
    integer, parameter :: cube = 64
    integer, parameter :: Lx = cube, Ly = cube, Lz = cube
    integer, parameter :: D = 3, Q = 19
    
    integer, parameter :: size = Lx*Ly*Lz*Q
    integer, parameter :: x_mult = Ly*Lz*Q
    integer, parameter :: y_mult = Lz*Q
    integer, parameter :: z_mult = Q

    real(8), parameter :: tau = 0.55
    real(8), parameter :: Utau = 1.0/tau
    real(8), parameter :: UmUtau = 1.0 - Utau
end module constants