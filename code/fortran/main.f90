program main
    use lb_d3q19
    implicit none
    integer :: t, t_max
    real(8) :: rho0, v0, zero
    t_max = 1
    rho0 = 1.0
    v0 = 10.0
    zero = 0.0

    call initialize(rho0, zero, zero, zero)
    do t=0, t_max
        call collide
        call impose_fields(v0)
        call propagate
    end do
    call print(v0)

end program main