module oas_main

  implicit none

contains

  subroutine manipulate_mesh_main(nx, ny, input_mesh, taper, chord, sweep, xshear, &
    span, yshear, dihedral, zshear, twist, symmetry, rotate_x, mesh)

    implicit none

    integer, intent(in) :: nx, ny
    real(kind=8), intent(in) :: input_mesh(nx, ny, 3), taper, chord(ny)
    real(kind=8), intent(in) :: sweep, xshear(ny), span, yshear(ny)
    real(kind=8), intent(in) :: dihedral, zshear(ny), twist(ny)
    logical, intent(in) :: symmetry, rotate_x

    real(kind=8), intent(out) :: mesh(nx, ny, 3)

    real(kind=8) :: le(ny, 3), te(ny, 3), quarter_chord(ny, 3), p180, tan_theta
    real(kind=8) :: dx(ny), y0, rad_twist(ny), rotation_matrix(ny, 3, 3)
    real(kind=8) :: row(ny, 3), out(3), taper_lins(ny), taper_lins_sym((ny+1)/2)
    real(kind=8) :: rad_theta_x(ny), one, dz_qc(ny-1), dy_qc(ny-1), s(ny), new_span
    real(kind=8) :: dz_qc_l((ny-1)/2), dz_qc_r((ny-1)/2), dy_qc_l((ny-1)/2), dy_qc_r((ny-1)/2)
    real(kind=8) :: computed_span
    integer :: ny2, ix, iy, ind

    p180 = 3.14159265358979323846264338 / 180.
    mesh = input_mesh
    one = 1.

    ! Taper
    le = mesh(1, :, :)
    te = mesh(nx, :, :)
    quarter_chord = 0.25 * te + 0.75 * le

    if (symmetry) then
      computed_span = quarter_chord(ny, 2) - quarter_chord(1, 2)

      do iy=1,ny
        taper_lins(iy) = (quarter_chord(iy, 2) - quarter_chord(1, 2)) / computed_span * (1 - taper) + taper
      end do

      do iy=1,ny
        do ix=1,nx
          do ind=1,3
            mesh(ix, iy, ind) = (mesh(ix, iy, ind) - quarter_chord(iy, ind)) * taper_lins(iy) + &
              quarter_chord(iy, ind)
          end do
        end do
      end do

    else

      computed_span = quarter_chord(ny, 2) - quarter_chord(1, 2)
      ny2 = (ny - 1) / 2

      do iy=1,ny2
        dx(iy) = 1 + quarter_chord(iy, 2) / (computed_span / 2) * taper
      end do

      do iy=1,ny2
        dx(iy) = (quarter_chord(iy, 2) - quarter_chord(1, 2)) / (computed_span / 2) * (1 - taper) + taper
      end do

      do iy=ny,ny2+1,-1
        dx(iy) = -(quarter_chord(iy, 2) - quarter_chord(ny, 2)) / (computed_span / 2) * (1 - taper) + taper
      end do

      do iy=1,ny
        do ix=1,nx
          do ind=1,3
            mesh(ix, iy, ind) = (mesh(ix, iy, ind) - quarter_chord(iy, ind)) * dx(iy) + &
              quarter_chord(iy, ind)
          end do
        end do
      end do
    end if

    ! Scale x
    le = mesh(1, :, :)
    te = mesh(nx, :, :)
    quarter_chord = 0.25 * te + 0.75 * le

    do iy=1,ny
      mesh(:, iy, 1) = (mesh(:, iy, 1) - quarter_chord(iy, 1)) * chord(iy) + &
        quarter_chord(iy, 1)
    end do

    ! Sweep
    le = mesh(1, :, :)
    tan_theta = tan(p180 * sweep)

    if (symmetry) then
      y0 = le(ny, 2)
      dx = -(le(:, 2) - y0) * tan_theta
    else
      ny2 = (ny - 1) / 2
      y0 = le(ny2+1, 2)

      dx(:ny2) = -(le(:ny2, 2) - y0) * tan_theta
      dx(ny2+1:) = (le(ny2+1:, 2) - y0) * tan_theta
    end if

    do ix=1,nx
      mesh(ix, :, 1) = mesh(ix, :, 1) + dx
    end do

    ! x shear
    do ix=1,nx
      mesh(ix, :, 1) = mesh(ix, :, 1) + xshear
    end do

    ! Span
    le = mesh(1, :, :)
    te = mesh(nx, :, :)
    quarter_chord = 0.25 * te + 0.75 * le
    new_span = span

    if (symmetry) then
      new_span = span / 2.
    end if

    s = quarter_chord(:, 2) / (quarter_chord(ny, 2) - quarter_chord(1, 2))
    do ix=1,nx
      mesh(ix, :, 2) = s * new_span
    end do

    ! y shear
    do ix=1,nx
      mesh(ix, :, 2) = mesh(ix, :, 2) + yshear
    end do

    ! Dihedral
    le = mesh(1, :, :)
    tan_theta = tan(p180 * dihedral)

    if (symmetry) then
      y0 = le(ny, 2)
      dx = -(le(:, 2) - y0) * tan_theta
    else
      ny2 = (ny - 1) / 2
      y0 = le(ny2+1, 2)

      dx(:ny2) = -(le(:ny2, 2) - y0) * tan_theta
      dx(ny2+1:) = (le(ny2+1:, 2) - y0) * tan_theta
    end if

    do ix=1,nx
      mesh(ix, :, 3) = mesh(ix, :, 3) + dx
    end do

    ! z shear
    do ix=1,nx
      mesh(ix, :, 3) = mesh(ix, :, 3) + zshear
    end do

    ! Rotate
    le = mesh(1, :, :)
    te = mesh(nx, :, :)
    quarter_chord = 0.25 * te + 0.75 * le

    rad_theta_x(:) = 0.
    if (rotate_x) then
      if (symmetry) then
        dz_qc = quarter_chord(:ny-1,3) - quarter_chord(2:,3)
        dy_qc = quarter_chord(:ny-1,2) - quarter_chord(2:,2)
        rad_theta_x(:ny-1) = atan(dz_qc / dy_qc)
      else
        ny2 = (ny - 1) / 2
        dz_qc_l = quarter_chord(:ny2,3) - quarter_chord(2:ny2+1,3)
        dy_qc_l = quarter_chord(:ny2,2) - quarter_chord(2:ny2+1,2)
        rad_theta_x(:ny2) = atan(dz_qc_l / dy_qc_l)
        dz_qc_r = quarter_chord(ny2+2:,3) - quarter_chord(ny2+1:ny-1,3)
        dy_qc_r = quarter_chord(ny2+2:,2) - quarter_chord(ny2+1:ny-1,2)
        rad_theta_x(ny2+2:) = atan(dz_qc_r / dy_qc_r)
      end if
    end if

    rad_twist = twist * p180
    rotation_matrix(:, :, :) = 0.
    rotation_matrix(:, 1, 1) = cos(rad_twist)
    rotation_matrix(:, 1, 3) = sin(rad_twist)
    rotation_matrix(:, 2, 1) = sin(rad_theta_x) * sin(rad_twist)
    rotation_matrix(:, 2, 2) = cos(rad_theta_x)
    rotation_matrix(:, 2, 3) = -sin(rad_theta_x) * cos(rad_twist)
    rotation_matrix(:, 3, 1) = -cos(rad_theta_x) * sin(rad_twist)
    rotation_matrix(:, 3, 2) = sin(rad_theta_x)
    rotation_matrix(:, 3, 3) = cos(rad_theta_x) * cos(rad_twist)

    do ix=1,nx
      row = mesh(ix, :, :)
      do iy=1,ny
        call matmul2(3, 3, 1, rotation_matrix(iy, :, :), row(iy, :) - quarter_chord(iy, :), out)
        mesh(ix, iy, :) = out
      end do
      mesh(ix, :, :) = mesh(ix, :, :) + quarter_chord
    end do

  end subroutine

  subroutine linspace(l, k, n, z)

    implicit none

    !// Argument declarations
    integer, intent(in) :: n
    real(kind=8), dimension(n), intent(out) :: z
    real(kind=8), intent(in) :: l
    real(kind=8), intent(in) :: k

    !// local variables
    integer :: i
    real(kind=8) :: d

    d = (k - l) / (n - 1)
    z(1) = l
    do i = 2, n-1
      z(i) = z(i-1) + d
    end do
    z(1) = l
    z(n) = k
    return

  end subroutine linspace

  subroutine calc_vonmises_main(nodes, r, disp, E, G, x_gl, n, vonmises)

    implicit none

    ! Input
    integer, intent(in) :: n
    real(kind=8), intent(in) :: nodes(n, 3), r(n-1), disp(n, 6)
    real(kind=8), intent(in) :: E, G, x_gl(3)

    ! Output
    real(kind=8), intent(out) :: vonmises(n-1, 2)

    ! Working
    integer :: ielem
    real(kind=8) :: P0(3), P1(3), L, x_loc(3), y_loc(3), z_loc(3), T(3, 3)
    real(kind=8) :: u0(3), r0(3), u1(3), r1(3), sxx0, sxx1, sxt, tmp
    real(kind=8) :: y_raw(3), z_raw(3), r1r0(3), t1(3), t2(3), t3(3), t4(3)
    real(kind=8) :: nodes2(n, 3), r2(n-1), disp2(n, 6), p1p0(3)

    vonmises(:, :) = 0.
    nodes2 = nodes
    r2 = r
    disp2 = disp

    do ielem=1, n-1

      P0 = nodes2(ielem, :)
      P1 = nodes2(ielem+1, :)
      p1p0 = P1 - P0
      call norm(p1p0, L)

      call unit(p1p0, x_loc)
      call cross(x_loc, x_gl, y_raw)
      call unit(y_raw, y_loc)
      call cross(x_loc, y_loc, z_raw)
      call unit(z_raw, z_loc)

      T(1, :) = x_loc
      T(2, :) = y_loc
      T(3, :) = z_loc

      t1 = disp2(ielem, 1:3)
      t2 = disp2(ielem, 4:6)
      t3 = disp2(ielem+1, 1:3)
      t4 = disp2(ielem+1, 4:6)

      call matmul2(3, 3, 1, T, t1, u0)
      call matmul2(3, 3, 1, T, t2, r0)
      call matmul2(3, 3, 1, T, t3, u1)
      call matmul2(3, 3, 1, T, t4, r1)

      r1r0 = r1 - r0
      tmp = (r1r0(2)**2 + r1r0(3)**2)**.5
      sxx0 = E * (u1(1) - u0(1)) / L + E * r2(ielem) / L * tmp
      sxx1 = E * (u0(1) - u1(1)) / L + E * r2(ielem) / L * tmp
      sxt = G * r2(ielem) * (r1r0(1)) / L

      vonmises(ielem, 1) = (sxx0**2 + 3 * sxt**2)**.5
      vonmises(ielem, 2) = (sxx1**2 + 3 * sxt**2)**.5

    end do

  end subroutine


  subroutine transferdisplacements_main(nx, ny, mesh, disp, w, def_mesh)

    implicit none

    ! Input
    integer, intent(in) :: nx, ny
    real(kind=8), intent(in) :: mesh(nx, ny, 3), disp(ny, 6), w

    ! Output
    real(kind=8), intent(out) :: def_mesh(nx, ny, 3)

    ! Working
    integer :: ind, indx
    real(kind=8) :: Smesh(nx, ny, 3), T(3, 3), T_base(3, 3), vec(3)
    real(kind=8) :: sinr(3), cosr(3), r(3), ref_curve(ny, 3)

    ref_curve = (1-w) * mesh(1, :, :) + w * mesh(nx, :, :)

    def_mesh(:, :, :) = 0.

    T_base(:, :) = 0.
    do ind=1,3
      T_base(ind, ind) = -2.
    end do

    do ind=1,nx
      Smesh(ind, :, :) = mesh(ind, :, :) - ref_curve
    end do

    do ind=1,ny
      r = disp(ind, 4:6)
      cosr = cos(r)
      sinr = sin(r)
      T(:, :) = 0.

      T(1, 1) = cosr(3) + cosr(2)
      T(2, 2) = cosr(3) + cosr(1)
      T(3, 3) = cosr(1) + cosr(2)
      T(1, 2) = -sinr(3)
      T(1, 3) = sinr(2)
      T(2, 1) = sinr(3)
      T(2, 3) = -sinr(1)
      T(3, 1) = -sinr(2)
      T(3, 2) = sinr(1)

      T = T + T_base

      do indx=1,nx
        call matmul2(1, 3, 3, Smesh(indx, ind, :), T, vec)
        def_mesh(indx, ind, :) = def_mesh(indx, ind, :) + vec
      end do

      def_mesh(:, ind, 1) = def_mesh(:, ind, 1) + disp(ind, 1)
      def_mesh(:, ind, 2) = def_mesh(:, ind, 2) + disp(ind, 2)
      def_mesh(:, ind, 3) = def_mesh(:, ind, 3) + disp(ind, 3)
    end do

    def_mesh = def_mesh + mesh

  end subroutine


  subroutine assemblestructmtx_main(n, tot_n_fem, nodes, A, J, Iy, Iz, &
    K_a, K_t, K_y, K_z, &
    cons, E, G, x_gl, T, &
    K_elem, Pelem_a, Pelem_t, Pelem_y, Pelem_z, T_elem, &
    const2, const_y, const_z, K)

    implicit none

    ! Input
    integer, intent(in) :: n, cons, tot_n_fem
    real(kind=8), intent(in) :: nodes(tot_n_fem, 3), A(n-1), J(n-1), Iy(n-1), Iz(n-1)
    real(kind=8), intent(in) :: E, G, x_gl(3)
    real(kind=8), intent(inout) :: K_a(2, 2), K_t(2, 2), K_y(4, 4), K_z(4, 4)
    real(kind=8), intent(inout) :: T(3, 3), K_elem(12, 12), T_elem(12, 12)
    real(kind=8), intent(in) :: Pelem_a(2, 12), Pelem_t(2, 12), Pelem_y(4, 12), Pelem_z(4, 12)
    real(kind=8), intent(in) :: const2(2, 2), const_y(4, 4), const_z(4, 4)

    ! Output
    real(kind=8), intent(out) :: K(6*n+6, 6*n+6)

    ! Working
    real(kind=8) :: P0(3), P1(3), x_loc(3), y_loc(3), z_loc(3), x_cross(3), y_cross(3)
    real(kind=8) :: L, EA_L, GJ_L, EIy_L3, EIz_L3, res(12, 12)
    real(kind=8) :: mat12x12(12, 12), mat12x4(12, 4), mat12x2(12, 2)
    integer ::  num_elems, num_nodes, num_cons, ielem, in0, in1, ind, i
    real(kind=8) :: Pelem_a_T(12, 2), Pelem_t_T(12, 2), K_(6*n+6, 6*n+6)
    real(kind=8) :: Pelem_y_T(12, 4), Pelem_z_T(12, 4), T_elem_T(12, 12)


    num_elems = n - 1
    num_nodes = n
    num_cons = 1 ! only 1 con in current spatialbeam code

    K(:, :) = 0.
    do ielem = 1, num_elems ! loop over num elements
      P0 = nodes(ielem, :)
      P1 = nodes(ielem+1, :)

      call unit(P1 - P0, x_loc)
      call cross(x_loc, x_gl, x_cross)
      call unit(x_cross, y_loc)
      call cross(x_loc, y_loc, y_cross)
      call unit(y_cross, z_loc)

      T(1, :) = x_loc
      T(2, :) = y_loc
      T(3, :) = z_loc

      do ind = 1, 4
        T_elem(3*(ind-1)+1:3*(ind-1)+3, 3*(ind-1)+1:3*(ind-1)+3) = T
      end do

      call norm(P1 - P0, L)
      EA_L = E * A(ielem) / L
      GJ_L = G * J(ielem) / L
      EIy_L3 = E * Iy(ielem) / L**3
      EIz_L3 = E * Iz(ielem) / L**3

      K_a(:, :) = EA_L * const2
      K_t(:, :) = GJ_L * const2

      K_y(:, :) = EIy_L3 * const_y
      K_y(2, :) = K_y(2, :) * L
      K_y(4, :) = K_y(4, :) * L
      K_y(:, 2) = K_y(:, 2) * L
      K_y(:, 4) = K_y(:, 4) * L

      K_z(:, :) = EIz_L3 * const_z
      K_z(2, :) = K_z(2, :) * L
      K_z(4, :) = K_z(4, :) * L
      K_z(:, 2) = K_z(:, 2) * L
      K_z(:, 4) = K_z(:, 4) * L

      K_elem(:, :) = 0.
      call transpose2(2, 12, Pelem_a, Pelem_a_T)
      call matmul2(12, 2, 2, Pelem_a_T, K_a, mat12x2)
      call matmul2(12, 2, 12, mat12x2, Pelem_a, res)
      K_elem = K_elem + res

      call transpose2(2, 12, Pelem_t, Pelem_t_T)
      call matmul2(12, 2, 2, Pelem_t_T, K_t, mat12x2)
      call matmul2(12, 2, 12, mat12x2, Pelem_t, res)
      K_elem = K_elem + res

      call transpose2(4, 12, Pelem_y, Pelem_y_T)
      call matmul2(12, 4, 4, Pelem_y_T, K_y, mat12x4)
      call matmul2(12, 4, 12, mat12x4, Pelem_y, res)
      K_elem = K_elem + res

      call transpose2(4, 12, Pelem_z, Pelem_z_T)
      call matmul2(12, 4, 4, Pelem_z_T, K_z, mat12x4)
      call matmul2(12, 4, 12, mat12x4, Pelem_z, res)
      K_elem = K_elem + res

      call transpose2(12, 12, T_elem, T_elem_T)
      call matmul2(12, 12, 12, T_elem_T, K_elem, mat12x12)
      call matmul2(12, 12, 12, mat12x12, T_elem, res)

      in0 = ielem
      in1 = ielem + 1

      K(6*(in0-1)+1:6*(in0-1)+6, 6*(in0-1)+1:6*(in0-1)+6) = &
      K(6*(in0-1)+1:6*(in0-1)+6, 6*(in0-1)+1:6*(in0-1)+6) + res(:6, :6)

      K(6*(in1-1)+1:6*(in1-1)+6, 6*(in0-1)+1:6*(in0-1)+6) = &
      K(6*(in1-1)+1:6*(in1-1)+6, 6*(in0-1)+1:6*(in0-1)+6) + res(7:, :6)

      K(6*(in0-1)+1:6*(in0-1)+6, 6*(in1-1)+1:6*(in1-1)+6) = &
      K(6*(in0-1)+1:6*(in0-1)+6, 6*(in1-1)+1:6*(in1-1)+6) + res(:6, 7:)

      K(6*(in1-1)+1:6*(in1-1)+6, 6*(in1-1)+1:6*(in1-1)+6) = &
      K(6*(in1-1)+1:6*(in1-1)+6, 6*(in1-1)+1:6*(in1-1)+6) + res(7:, 7:)

    end do

    do i = 1, 6
      K(6*num_nodes+i, 6*cons+i) = 10**9.
      K(6*cons+i, 6*num_nodes+i) = 10**9.
    end do

  end subroutine assemblestructmtx_main

  subroutine transpose2(m, n, mtx, new_mtx)

    implicit none

    integer, intent(in) :: m, n
    real(kind=8), intent(in) :: mtx(m, n)
    real(kind=8), intent(out) :: new_mtx(n, m)

    integer :: i, j

    do i=1,m
      do j=1,n
        new_mtx(j, i) = mtx(i, j)
      end do
    end do

  end subroutine

  subroutine matmul2(m, n, p, A, B, C)

    implicit none

    integer, intent(in) :: m, n, p
    real(kind=8), intent(in) :: A(m, n), B(n, p)
    real(kind=8), intent(out) :: C(m, p)

    integer :: i, j, k

    C(:, :) = 0.

    do i=1,m
      do j=1,p
        do k=1,n
          C(i, j) = C(i, j) + A(i, k) * B(k, j)
        end do
      end do
    end do

  end subroutine

  subroutine assembleaeromtx_main(ny, nx, ny_, nx_, alpha, points, bpts, mesh, skip, symmetry, mtx)

    implicit none

    ! Input
    integer, intent(in) :: ny, nx, ny_, nx_
    real(kind=8), intent(in) :: alpha, mesh(nx_, ny_, 3)
    real(kind=8), intent(in) :: points(nx-1, ny-1, 3), bpts(nx_-1, ny_, 3)
    logical, intent(in) :: skip, symmetry

    ! Output
    real(kind=8), intent(out) :: mtx((nx-1)*(ny-1), (nx_-1)*(ny_-1), 3)

    ! Working
    integer :: el_j, el_i, cp_j, cp_i, el_loc_j, el_loc, cp_loc_j, cp_loc
    real(kind=8) :: pi, P(3), A(3), B(3), u(3), C(3), D(3)
    real(kind=8) :: A_sym(3), B_sym(3), C_sym(3), D_sym(3)
    real(kind=8) :: ur2(3), r1(3), r2(3), r1_mag, r2_mag
    real(kind=8) :: ur1(3), bound(3), dot_ur2, dot_ur1
    real(kind=8) :: edges(3), C_te(3), D_te(3), C_te_sym(3), D_te_sym(3)

    pi = 4.d0*atan(1.d0)

    ! Trailing vortices in AVL follow the x-axis; no cos or sin
    ! u(1) = 1.
    ! u(3) = 0.
    u(1) = cos(alpha * pi / 180.)
    u(2) = 0.
    u(3) = sin(alpha * pi / 180.)

    mtx(:, :, :) = 0.

    do el_j = 1, ny_-1 ! spanwise loop through horseshoe elements
      el_loc_j = (el_j - 1) * (nx_ - 1)
      C_te = mesh(nx_, el_j + 1, :)
      D_te = mesh(nx_, el_j + 0, :)

      if (symmetry) then
        C_te_sym = C_te
        D_te_sym = D_te
        C_te_sym(2) = -C_te_sym(2)
        D_te_sym(2) = -D_te_sym(2)
      end if

      do cp_j = 1, ny-1 ! spanwise loop through control points
        cp_loc_j = (cp_j - 1) * (nx - 1)

        do cp_i = 1, nx-1 ! chordwise loop through control points
          cp_loc = cp_i + cp_loc_j
          P = points(cp_i, cp_j, :)

          r1 = P - D_te
          r2 = P - C_te
          call norm(r1, r1_mag)
          call norm(r2, r2_mag)

          call cross(u, r2, ur2)
          call cross(u, r1, ur1)

          edges(:) = 0.
          call dot(u, r2, dot_ur2)
          call dot(u, r1, dot_ur1)
          edges = ur2 / (r2_mag * (r2_mag - dot_ur2))
          edges = edges - ur1 / (r1_mag * (r1_mag - dot_ur1))

          if (symmetry) then
            r1 = P - D_te_sym
            r2 = P - C_te_sym
            call norm(r1, r1_mag)
            call norm(r2, r2_mag)

            call cross(u, r2, ur2)
            call cross(u, r1, ur1)

            call dot(u, r2, dot_ur2)
            call dot(u, r1, dot_ur1)

            edges = edges - ur2 / (r2_mag * (r2_mag - dot_ur2))
            edges = edges + ur1 / (r1_mag * (r1_mag - dot_ur1))
          end if

          do el_i = nx_-1, 1, -1 ! chordwise loop through horseshoe elements
            el_loc = el_i + el_loc_j

            A = bpts(el_i + 0, el_j + 0, :)
            B = bpts(el_i + 0, el_j + 1, :)

            if (el_i .EQ. nx_ - 1) then
              C = C_te
              D = D_te
            else
              C = bpts(el_i + 1, el_j + 1, :)
              D = bpts(el_i + 1, el_j + 0, :)
            end if

            call calc_vorticity(B, C, P, edges)
            call calc_vorticity(D, A, P, edges)

            bound(:) = 0.
            if (symmetry) then
              A_sym = A
              B_sym = B
              C_sym = C
              D_sym = D
              A_sym(2) = -A_sym(2)
              B_sym(2) = -B_sym(2)
              C_sym(2) = -C_sym(2)
              D_sym(2) = -D_sym(2)

              call calc_vorticity(C_sym, B_sym, P, edges)
              call calc_vorticity(A_sym, D_sym, P, edges)
              call calc_vorticity(B_sym, A_sym, P, bound)
            end if



            if (.not. ((skip) .and. (cp_loc .EQ. el_loc))) then
              call calc_vorticity(A, B, P, bound)
            end if

            mtx(cp_loc, el_loc, :) = edges + bound

          end do
        end do

       end do
    end do

  end subroutine assembleaeromtx_main

  subroutine calc_vorticity(A, B, P, out)

    implicit none

    ! Input
    real(kind=8), intent(in) :: A(3), B(3), P(3)

    ! Output
    real(kind=8), intent(inout) :: out(3)

    ! Working
    real(kind=8) :: r1(3), r2(3), r1_mag, r2_mag, r1r2(3), mag_mult, dot_r1r2

    r1 = P - A
    r2 = P - B

    call norm(r1, r1_mag)
    call norm(r2, r2_mag)

    call cross(r1, r2, r1r2)
    mag_mult = r1_mag * r2_mag

    call dot(r1, r2, dot_r1r2)
    out = out + (r1_mag + r2_mag) * r1r2 / (mag_mult * (mag_mult + dot_r1r2))

  end subroutine calc_vorticity



  subroutine biotsavart(A, B, P, inf, rev, out)

    implicit none

    ! Input
    real(kind=8), intent(in) :: A(3), B(3), P(3)
    logical, intent(in) :: inf, rev

    ! Output
    real(kind=8), intent(inout) :: out(3)

    ! Working
    real(kind=8) :: rPA, rPB, rAB, rH
    real(kind=8) :: cosA, cosB, C(3)
    real(kind=8) :: eps, tmp(3), dot_BAPA, dot_BABA, dot_PBAB

    eps = 1e-5

    call norm(A - P, rPA)
    call norm(B - P, rPB)
    call norm(B - A, rAB)

    call dot(B - A, P - A, dot_BAPA)
    call dot(B - A, B - A, dot_BABA)
    call dot(P - B, A - B, dot_PBAB)

    call norm(P - A - dot_BAPA / dot_BABA * (B - A), rH)
    rH = rH + eps
    cosA = dot_BAPA / (rPA * rAB)
    cosB = dot_PBAB / (rPB * rAB)
    call cross(B - P, A - P, C)
    call unit(C, C)

    if (inf) then
       tmp = -C / rH * (cosA + 1)
    else
       tmp = -C / rH * (cosA + cosB)
    end if

    if (rev) then
       tmp = -tmp
    end if

    out = out + tmp

  end subroutine biotsavart

  subroutine forcecalc_main(v, circ, rho, bpts, nx, ny, num_panels, sec_forces)

    implicit none

    real(kind=8), intent(in) :: v(num_panels, 3), circ(num_panels), rho, bpts(nx-1, ny, 3)
    integer, intent(in) :: nx, ny, num_panels

    real(kind=8), intent(out) :: sec_forces(num_panels, 3)

    real(kind=8) :: bound(num_panels, 3), v_cross_bound(num_panels, 3), tmp(3)
    integer :: i, j, k

    do j=1,ny-1
      do i=1,nx-1
        bound((j-1)*(nx-1) + i, :) = bpts(i, j+1, :) - bpts(i, j, :)
      end do
    end do

    do i=1,num_panels
      call cross(v(i, :), bound(i, :), tmp)
      v_cross_bound(i, :) = tmp
    end do

    do i=1,3
      sec_forces(:, i) = rho * circ * v_cross_bound(:, i)
    end do

  end subroutine

  subroutine momentcalc_main(bpts, cg, chords, widths, S_ref, sec_forces, symmetry, nx, ny, M)

    implicit none

    real(kind=8), intent(in) :: bpts(nx-1, ny, 3)
    integer, intent(in) :: nx, ny
    real(kind=8), intent(in) :: cg(3), S_ref
    real(kind=8), intent(in) :: chords(ny), widths(ny-1)
    logical, intent(in) :: symmetry
    real(kind=8), intent(in) :: sec_forces(nx-1, ny-1, 3)

    real(kind=8), intent(out) :: M(3)

    real(kind=8) :: panel_chords(ny-1), MAC, moment(ny-1, 3), tmp(3)
    integer :: i, j, k

    panel_chords = (chords(2:) + chords(:ny-1)) / 2.
    MAC = 1. / S_ref * sum(panel_chords**2 * widths)

    if (symmetry) then
      MAC = MAC * 2
    end if

    moment(:, :) = 0.
    do j=1,ny-1
      do i=1,nx-1
        call cross((bpts(i, j+1, :) + bpts(i, j, :)) / 2. - cg, sec_forces(i, j, :), tmp)
        moment(j, :) = moment(j, :) + tmp
      end do
    end do
    moment = moment / MAC

    if (symmetry) then
      moment(:, 1) = 0.
      moment(:, 2) = moment(:, 2) * 2
      moment(:, 3) = 0.
    end if

    M = 0.
    do j=1,ny-1
      M = M + moment(j, :)
    end do

  end subroutine

  subroutine compute_normals_main(nx, ny, mesh, normals, S_ref)

    implicit none

    real(kind=8), intent(in) :: mesh(nx, ny, 3)
    integer, intent(in) :: nx, ny

    real(kind=8), intent(out) :: normals(nx-1, ny-1, 3), S_ref

    integer :: i, j
    real(kind=8) :: norms(nx, ny), out(3)

    do i=1,nx-1
      do j=1,ny-1
        call cross(mesh(i, j+1, :) - mesh(i+1, j, :), mesh(i, j, :) - mesh(i+1, j+1, :), out)
        normals(i, j, :) = out
        norms(i, j) = sqrt(sum(normals(i, j, :)**2))
        normals(i, j, :) = normals(i, j, :) / norms(i, j)
      end do
    end do
    S_ref = 0.5 * sum(norms)

  end subroutine


! REAL FUNCTIONS

  subroutine unit(v, U)

    implicit none

    real(kind=8), intent(in) :: v(3)
    real(kind=8), intent(out) :: U(3)
    real(kind=8) :: nm

    call norm(v, nm)
    U(1) = v(1) / nm
    U(2) = v(2) / nm
    U(3) = v(3) / nm

  end subroutine unit

  subroutine norm(v, norm_output)

    implicit none

    real(kind=8), intent(in) :: v(3)
    real(kind=8), intent(out) :: norm_output
    real(kind=8) :: dot_prod, v1(3)

    ! Need to create a copy of v here so the reverse mode AD works correctly
    v1 = v
    !norm = sqrt(dot_product(v, v))
    call dot(v, v1, dot_prod)
    norm_output = dot_prod ** 0.5

  end subroutine norm



  subroutine dot(a, b, dot_prod)

    implicit none

    real(kind=8), intent(in) :: a(3), b(3)
    real(kind=8), intent(out) :: dot_prod

    dot_prod = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

  end subroutine dot



  subroutine cross(A, B, C)

    implicit none

    real(kind=8), intent(in) :: A(3), B(3)
    real(kind=8), intent(out) :: C(3)

    C(1) = A(2) * B(3) - A(3) * B(2)
    C(2) = A(3) * B(1) - A(1) * B(3)
    C(3) = A(1) * B(2) - A(2) * B(1)

  end subroutine cross



end module
