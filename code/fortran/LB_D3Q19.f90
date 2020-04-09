
module LB_D3Q19
    use constants
    implicit none
    real(8), dimension(Q) :: w = (/1.0/3.0, &
        1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, &
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, &
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0/)
    integer, dimension(Q) :: Vx = (/0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0/)
    integer, dimension(Q) :: Vy = (/0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1/)
    integer, dimension(Q) :: Vz = (/0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1/)
    real(8), dimension(size) :: f, f_new
    
    public :: initialize
    public :: collide
    public :: impose_fields
    public :: propagate
    public :: print

contains
     function get_1D(ix, iy, iz) result(oneD)
        implicit none
        integer :: ix, iy, iz, oneD
        oneD = (ix-1)*x_mult + (iy-1)*y_mult + (iz-1)*z_mult
     end function get_1D

     function rho(ix, iy, iz) result(r)
        implicit none
        integer :: ix, iy, iz, i, pos
        real(8) :: r
        r= 0
        pos = get_1D(ix,iy,iz)
        do i = 1, Q
            r = r + f(pos+i)
        end do
     end function rho

     function Jx(ix, iy, iz) result(J)
        implicit none
        integer :: ix, iy, iz, i, pos
        real(8) :: J
        J = 0
        pos = get_1D(ix,iy,iz)
        do i = 1, Q
            J = J + f(pos+i)+Vx(i)
        end do
     end function Jx
     
     function Jy(ix, iy, iz) result(J)
        implicit none
        integer :: ix, iy, iz, i, pos
        real(8) :: J
        J = 0
        pos = get_1D(ix,iy,iz)
        do i = 1, Q
            J = J + f(pos+i)+Vy(i)
        end do
     end function Jy

     function Jz(ix, iy, iz) result(J)
        implicit none
        integer :: ix, iy, iz, i, pos
        real(8) :: J
        J = 0
        pos = get_1D(ix,iy,iz)
        do i = 1, Q
            J = J + f(pos+i)+Vz(i)
        end do
     end function Jz

    function f_eq(rho0, Ux0, Uy0, Uz0, i) result(f)
        implicit none
        real(8) :: rho0, Ux0, Uy0, Uz0, UdotVi, U2, f
        integer :: i
        UdotVi = Ux0*Vx(i) + Uy0*Vy(i) + Uz0*Vz(i)
        U2 = Ux0*Ux0 + Uy0*Uy0 + Uz0*Uz0
        f = rho0*w(i)*(1.0 + 3.0*UdotVi + 4.5*UdotVi*UdotVi - 1.5*U2)
    end function f_eq

    subroutine collide
        implicit none
        integer :: ix, iy, iz, i, pos
        real(8) :: rho0, Ux0, Uy0, Uz0
        do ix=1, Lx
            do iy=1, Ly
                do iz=1, Lz
                    rho0 = rho(ix, iy, iz)
                    Ux0 = Jx(ix, iy, iz)/rho0
                    Uy0 = Jy(ix, iy, iz)/rho0
                    Uz0 = Jz(ix, iy, iz)/rho0
                    pos = get_1D(ix, iy, iz)
                    do i=1, Q
                        f_new(pos+i) = UmUtau*f(pos+i) + Utau*f_eq(rho0, Ux0, Uy0, Uz0, i)
                    end do
                end do
            end do
        end do
        return
    end subroutine collide

    subroutine propagate
        implicit none
        integer :: ix, iy, iz, i, pos, pos_new, x_pos, y_pos, z_pos
        do ix=1, Lx
            do iy=1, Ly
                do iz=1, Lz
                    pos_new = get_1D(ix, iy, iz)
                    do i=1, Q
                        x_pos = ix + Vx(i)
                        y_pos = iy + Vy(i)
                        z_pos = iz + Vz(i)
                        if (x_pos<=1 .or. x_pos>Lx) then
                            x_pos = ix - Vx(i)
                        end if
                        if (y_pos<=1 .or. y_pos>Ly) then
                            y_pos = iy - Vy(i)
                        end if
                        if (z_pos<=1 .or. z_pos>Lz) then
                            z_pos = iz + Vz(i)
                        end if
                        pos = get_1D(x_pos, y_pos, z_pos)
                        f(pos+i) = f_new(pos_new+i)
                    end do
                end do
            end do
        end do
        return
    end subroutine propagate

    subroutine initialize(rho0, Ux0, Uy0, Uz0)
        implicit none
        real(8) :: rho0, Ux0, Uy0, Uz0
        integer :: ix, iy, iz, i, pos
        do ix=1, Lx
            do iy=1, Ly
                do iz=1, Lz
                    pos = get_1D(ix, iy, iz)
                    do i=1, Q
                        f(pos+i) = f_eq(rho0, Ux0, Uy0, Uz0, i)
                    end do
                end do
            end do
        end do
        return
    end subroutine initialize

    subroutine impose_fields(v)
        implicit none
        integer ix, iy, iz, i, pos
        real(8) :: rho0, v, zero
        zero = 0.0
        do ix=1, Lx
            do iy=1, Ly
                do iz=1, Lz
                    if (ix==1) then
                        pos = get_1D(ix, iy, iz)
                        rho0 = rho(ix, iy, iz)
                        do i=1, Q
                            f_new(pos+i) = f_eq(rho0, v/2.0, zero, zero, i)
                        end do
                    end if
                    if ((ix-Lx/2)*(ix-Lx/2) + (iy-Ly/2)*(iy-Ly/2) + (iz-Lz/2)*(iz-Lz/2) <= cube/4) then 
                        pos = get_1D(ix, iy, iz)
                        rho0 = rho(ix, iy, iz)
                        do i=1, Q
                            f_new(pos+i) = f_eq(rho0, zero, zero, zero, i)
                        end do
                    end if
                end do
            end do
        end do
    end subroutine impose_fields

    subroutine print(v)
        implicit none
        real(8) :: v, rho0, Ux0, Uy0, Uz0
        integer :: ix, iy, iz
        do ix=1, Lx, 4
            do iy=1, Ly, 4
                do iz=1, Lz, 4
                    rho0 = rho(ix, iy, iz)
                    Ux0 = 4*(Jx(ix, iy, iz)/(v*rho0))
                    Uy0 = 4*(Jy(ix, iy, iz)/(v*rho0))
                    Uz0 = 4*(Jz(ix, iy, iz)/(v*rho0))
                    write (*,*) ix, ' ', iy, ' ', iz, ' ', Ux0, ' ', Uy0, ' ', Uz0
                end do
            end do
        end do
    end subroutine print
end module LB_D3Q19