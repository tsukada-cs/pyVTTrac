! pyvttrac computational core.
!
! Fortran port of VTTrac.jl v2.0.0 (https://github.com/tsukada-cs/VTTrac.jl).
! A single public entry point, `track`, is exposed: one Fortran call performs
! the entire multi-step, multi-seed tracking run. 
! The `m` (seed) loop is parallelized with OpenMP;
! everything touched inside it is either read-only shared data or a
! thread-private automatic array, so no explicit locking is needed.
!
! Index convention: this module is entirely 1-based (Fortran-native), mirroring
! VTTrac.jl's 1-based (Julia-native) indexing. The 0-based <-> 1-based
! conversion at the package boundary is handled in Python (src/pyvttrac/api.py),
! exactly as v1 pyVTTrac did around its juliacall calls.
!
! Memory layout: `z` is declared `z(nx,ny,nt)` so that a NumPy C-order
! `(nt,ny,nx)` array and this Fortran array alias the same bytes with no
! transpose/copy (x, the innermost sliding-window scan direction, is
! contiguous in memory either way).
module pyvttrac_core
  use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
  private
  public :: track

  ! status codes (must match pyvttrac.status.Status)
  integer, parameter :: ST_OK = 0
  integer, parameter :: ST_TID_START_OUT_OF_RANGE = 1
  integer, parameter :: ST_TEMPLATE_READ_FAILED = 2
  integer, parameter :: ST_LOW_CONTRAST = 3
  integer, parameter :: ST_PEAK_NOT_INSIDE_TEMPLATE = 4
  integer, parameter :: ST_TID_END_OUT_OF_RANGE = 5
  integer, parameter :: ST_SCORE_COMPUTATION_FAILED = 6
  integer, parameter :: ST_PEAK_NOT_FOUND = 7
  integer, parameter :: ST_SCORE_BELOW_THRESHOLD = 8
  integer, parameter :: ST_VELOCITY_CHANGE_TOO_LARGE = 9

  integer, parameter :: METHOD_XCOR = 0
  integer, parameter :: METHOD_NCOV = 1

contains

  !===========================================================================
  ! Public entry point
  !===========================================================================
  subroutine track( &
       z, nx, ny, nt, t, &
       use_zmiss, zmiss, &
       use_mask, visible, min_samples, &
       nsx, nsy, ixhw, iyhw, &
       itstep, nsteps, &
       n, tid0, x0, y0, vx0, vy0, &
       method, subgrid, subgrid_gaus, &
       sth0, sth1, use_vch, vxch, vych, use_peak_th, peak_th, use_cth, cth, &
       fixed_template, &
       fmiss, imiss, nthreads, &
       cnt, status, tid, xo, yo, vxo, vyo, score, &
       want_zss, zss, want_scr, scr_ary)

    ! -- data --------------------------------------------------------------
    integer, intent(in) :: nx, ny, nt
    real(sp), intent(in) :: z(nx, ny, nt)
    real(dp), intent(in) :: t(nt)
    logical, intent(in) :: use_zmiss
    real(sp), intent(in) :: zmiss
    logical, intent(in) :: use_mask
    logical, intent(in) :: visible(nx, ny, nt)
    integer, intent(in) :: min_samples

    ! -- template / search geometry -----------------------------------------
    integer, intent(in) :: nsx, nsy, ixhw, iyhw

    ! -- time stepping -------------------------------------------------------
    integer, intent(in) :: itstep, nsteps

    ! -- seeds ----------------------------------------------------------------
    integer, intent(in) :: n
    integer, intent(in) :: tid0(n)
    real(dp), intent(in) :: x0(n), y0(n), vx0(n), vy0(n)

    ! -- scoring --------------------------------------------------------------
    integer, intent(in) :: method  ! 0=xcor, 1=ncov
    logical, intent(in) :: subgrid, subgrid_gaus

    ! -- screening --------------------------------------------------------------
    real(dp), intent(in) :: sth0, sth1
    logical, intent(in) :: use_vch
    real(dp), intent(in) :: vxch, vych
    logical, intent(in) :: use_peak_th
    real(sp), intent(in) :: peak_th
    logical, intent(in) :: use_cth
    real(sp), intent(in) :: cth
    logical, intent(in) :: fixed_template

    ! -- missing-value sentinels & threading ------------------------------------
    real(dp), intent(in) :: fmiss
    integer, intent(in) :: imiss
    integer, intent(in) :: nthreads

    ! -- outputs ----------------------------------------------------------------
    integer, intent(out) :: cnt(n)
    integer, intent(out) :: status(n)
    integer, intent(out) :: tid(nsteps + 1, n)
    real(dp), intent(out) :: xo(nsteps + 1, n), yo(nsteps + 1, n)
    real(dp), intent(out) :: vxo(nsteps, n), vyo(nsteps, n)
    real(dp), intent(out) :: score(nsteps, n)

    ! -- optional diagnostics -----------------------------------------------------
    logical, intent(in) :: want_zss, want_scr
    real(sp), intent(inout) :: zss(:, :, :, :)      ! (nsx,nsy,nsteps+1,n) if want_zss
    real(sp), intent(inout) :: scr_ary(:, :, :, :)  ! (2ixhw+1,2iyhw+1,nsteps,n) if want_scr

    ! -- locals (thread-private via OpenMP `private` clause below) --------------
    integer :: m, j, tidf, tidl, kc, lc
    real(dp) :: xcur, ycur, vxg, vyg, dtt, xw, yw, vxw, vyw
    real(dp) :: kpi_r, lpi_r, peak_r, sth
    real(sp) :: fmiss_sp
    real(sp) :: tmpl(nsx, nsy), xd(nsx, nsy), sigx
    logical :: vis_tmpl(nsx, nsy)
    real(sp) :: scr(2 * ixhw + 1, 2 * iyhw + 1)
    real(sp) :: tmpl_diag(nsx, nsy)
    logical :: ok, need_read, peak_failed, vchange_failed
    integer :: nvis
    integer :: nthr

    fmiss_sp = real(fmiss, sp)

    ! Resolved locally rather than via omp_set_num_threads(), which would
    ! mutate process-wide OpenMP state and leak into unrelated track() calls
    ! (and other OpenMP-using libraries) for the rest of the process lifetime.
    nthr = 1
#ifdef _OPENMP
    if (nthreads > 0) then
      nthr = nthreads
    else if (nthreads < 0) then
      nthr = omp_get_num_procs()
    else
      nthr = omp_get_max_threads()
    end if
#endif

    ! -- initial record (unconditional, mirrors VTTrac.jl do_tracking) ----------
    do m = 1, n
      if (subgrid) then
        xo(1, m) = x0(m)
        yo(1, m) = y0(m)
      else
        xo(1, m) = real(nint(x0(m)), dp)
        yo(1, m) = real(nint(y0(m)), dp)
      end if
      tid(1, m) = tid0(m)
      cnt(m) = 0
      status(m) = ST_OK
    end do
    xo(2:nsteps + 1, :) = fmiss
    yo(2:nsteps + 1, :) = fmiss
    tid(2:nsteps + 1, :) = imiss
    vxo = fmiss; vyo = fmiss; score = fmiss

    !$omp parallel do default(none) schedule(dynamic) num_threads(nthr) &
    !$omp   shared(z, nx, ny, nt, t, use_zmiss, zmiss, use_mask, visible, min_samples, &
    !$omp          nsx, nsy, ixhw, iyhw, itstep, nsteps, n, tid0, vx0, vy0, &
    !$omp          method, subgrid, subgrid_gaus, sth0, sth1, use_vch, vxch, vych, &
    !$omp          use_peak_th, peak_th, use_cth, cth, fixed_template, imiss, fmiss, &
    !$omp          fmiss_sp, cnt, status, tid, xo, yo, vxo, vyo, score, &
    !$omp          want_zss, zss, want_scr, scr_ary) &
    !$omp   private(m, j, tidf, tidl, kc, lc, xcur, ycur, vxg, vyg, dtt, xw, yw, &
    !$omp           vxw, vyw, kpi_r, lpi_r, peak_r, sth, tmpl, xd, sigx, vis_tmpl, &
    !$omp           scr, tmpl_diag, ok, need_read, peak_failed, vchange_failed, nvis)
    do m = 1, n
      do j = 1, nsteps
        if (status(m) /= ST_OK) exit

        xcur = xo(j, m)
        ycur = yo(j, m)
        tidf = tid0(m) + (j - 1) * itstep
        tidl = tidf + itstep

        if (tidf < 1 .or. tidf > nt) then
          status(m) = ST_TID_START_OUT_OF_RANGE
          exit
        end if

        need_read = (j == 1 .or. .not. fixed_template)
        if (need_read) then
          if (use_mask) then
            call read_template_visible(z, nx, ny, nt, visible, nsx, nsy, tidf, &
                 xcur, ycur, use_zmiss, zmiss, tmpl, vis_tmpl, ok)
            if (ok) then
              nvis = count(vis_tmpl)
              if (nvis < min_samples) ok = .false.
            end if
          else
            call read_template(z, nx, ny, nt, nsx, nsy, tidf, xcur, ycur, &
                 use_zmiss, zmiss, tmpl, ok)
          end if
          if (.not. ok) then
            status(m) = ST_TEMPLATE_READ_FAILED
            exit
          end if
        end if

        if (use_cth) then
          if (maxval(tmpl) - minval(tmpl) < cth) then
            status(m) = ST_LOW_CONTRAST
            exit
          end if
        end if

        if (use_peak_th) then
          call chk_peak_inside(tmpl, nsx, nsy, peak_th, peak_failed)
          if (peak_failed) then
            status(m) = ST_PEAK_NOT_INSIDE_TEMPLATE
            exit
          end if
        end if

        if (want_zss) then
          if (need_read) then
            zss(:, :, j, m) = tmpl
          else
            call read_template(z, nx, ny, nt, nsx, nsy, tidf, xcur, ycur, &
                 use_zmiss, zmiss, tmpl_diag, ok)
            if (.not. ok) then
              status(m) = ST_TEMPLATE_READ_FAILED
              exit
            end if
            zss(:, :, j, m) = tmpl_diag
          end if
        end if

        if (tidl < 1 .or. tidl > nt) then
          status(m) = ST_TID_END_OUT_OF_RANGE
          exit
        end if
        dtt = t(tidl) - t(tidf)

        if (j == 1) then
          vxg = vx0(m)
          vyg = vy0(m)
        else
          vxg = vxo(j - 1, m)
          vyg = vyo(j - 1, m)
        end if
        kc = nint(xcur + vxg * dtt)
        lc = nint(ycur + vyg * dtt)

        if (use_mask) then
          call sliding_score_masked(z, nx, ny, nt, visible, nsx, nsy, tmpl, vis_tmpl, &
               tidl, method, kc - ixhw, kc + ixhw, lc - iyhw, lc + iyhw, &
               use_zmiss, zmiss, min_samples, fmiss_sp, scr, ok)
        else
          call demean(tmpl, nsx, nsy, xd, sigx)
          call sliding_score(z, nx, ny, nt, nsx, nsy, xd, sigx, tidl, method, &
               kc - ixhw, kc + ixhw, lc - iyhw, lc + iyhw, &
               use_zmiss, zmiss, fmiss_sp, scr, ok)
        end if

        if (.not. ok) then
          status(m) = ST_SCORE_COMPUTATION_FAILED
          exit
        end if

        if (want_scr) scr_ary(:, :, j, m) = scr

        call find_score_peak(scr, 2 * ixhw + 1, 2 * iyhw + 1, subgrid, subgrid_gaus, &
             kpi_r, lpi_r, peak_r, ok)
        if (.not. ok) then
          status(m) = ST_PEAK_NOT_FOUND
          exit
        end if

        if (j == 1) then
          sth = sth0
        else
          sth = sth1
        end if
        if (peak_r < sth) then
          status(m) = ST_SCORE_BELOW_THRESHOLD
          exit
        end if

        xw = kpi_r + real(kc, dp) - 1.0_dp - real(ixhw, dp)
        yw = lpi_r + real(lc, dp) - 1.0_dp - real(iyhw, dp)
        vxw = (xw - xo(j, m)) / dtt
        vyw = (yw - yo(j, m)) / dtt

        if (use_vch .and. j > 1) then
          vchange_failed = (abs(vxw - vxo(j - 1, m)) > vxch .or. abs(vyw - vyo(j - 1, m)) > vych)
          if (vchange_failed) then
            if (j == 2) then
              cnt(m) = 0
              tid(2, m) = imiss
              xo(2, m) = fmiss
              yo(2, m) = fmiss
              vxo(1, m) = fmiss
              vyo(1, m) = fmiss
              score(1, m) = fmiss
            end if
            status(m) = ST_VELOCITY_CHANGE_TOO_LARGE
            exit
          end if
        end if

        cnt(m) = j + 1
        tid(j + 1, m) = tidl
        xo(j + 1, m) = xw
        yo(j + 1, m) = yw
        vxo(j, m) = vxw
        vyo(j, m) = vyw
        score(j, m) = peak_r

        if (want_zss .and. j == nsteps) then
          call read_template(z, nx, ny, nt, nsx, nsy, tidl, xw, yw, &
               use_zmiss, zmiss, tmpl_diag, ok)
          if (.not. ok) then
            status(m) = ST_TEMPLATE_READ_FAILED
            exit
          end if
          zss(:, :, j + 1, m) = tmpl_diag
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine track

  !===========================================================================
  ! Template read-out (bilinearly interpolated at possibly non-integer x,y),
  ! with bounds and (optional) missing-value checking. Mirrors
  ! VTTrac.jl's get_zsub_subgrid + get_zsub_view.
  !===========================================================================
  subroutine read_template(z, nx, ny, nt, nsx, nsy, tid, xc, yc, use_zmiss, zmiss, tmpl, ok)
    integer, intent(in) :: nx, ny, nt, nsx, nsy, tid
    real(sp), intent(in) :: z(nx, ny, nt)
    real(dp), intent(in) :: xc, yc
    logical, intent(in) :: use_zmiss
    real(sp), intent(in) :: zmiss
    real(sp), intent(out) :: tmpl(nsx, nsy)
    logical, intent(out) :: ok

    integer :: xi, yi, x0i, y0i, i, k, isx, isy
    real(sp) :: dx0, dx1, dy0, dy1
    real(dp) :: ddx, ddy

    xi = nint(xc); yi = nint(yc)
    ddx = xc - real(xi, dp); ddy = yc - real(yi, dp)
    isx = 0; if (ddx > 0.0_dp) isx = 1; if (ddx < 0.0_dp) isx = -1
    isy = 0; if (ddy > 0.0_dp) isy = 1; if (ddy < 0.0_dp) isy = -1
    dx0 = real(abs(ddx), sp); dx1 = 1.0_sp - dx0
    dy0 = real(abs(ddy), sp); dy1 = 1.0_sp - dy0

    x0i = xi - nsx / 2
    y0i = yi - nsy / 2
    ok = .not. (x0i + min(0, isx) < 1 .or. x0i + nsx - 1 + abs(isx) > nx .or. &
                y0i + min(0, isy) < 1 .or. y0i + nsy - 1 + abs(isy) > ny)
    if (.not. ok) return

    if (use_zmiss) then
      if (any(z(x0i:x0i + nsx - 1 + abs(isx), y0i:y0i + nsy - 1 + abs(isy), tid) == zmiss)) then
        ok = .false.
        return
      end if
    end if

    if (isx == 0 .and. isy == 0) then
      tmpl = z(x0i:x0i + nsx - 1, y0i:y0i + nsy - 1, tid)
    else
      do k = 1, nsy
        do i = 1, nsx
          tmpl(i, k) = dx1 * dy1 * z(x0i + i - 1,      y0i + k - 1,      tid) &
                     + dx0 * dy1 * z(x0i + i - 1 + isx, y0i + k - 1,      tid) &
                     + dx1 * dy0 * z(x0i + i - 1,       y0i + k - 1 + isy, tid) &
                     + dx0 * dy0 * z(x0i + i - 1 + isx, y0i + k - 1 + isy, tid)
        end do
      end do
    end if
  end subroutine read_template

  !===========================================================================
  ! As read_template, but also bilinearly combines + rounds a `visible` mask
  ! (0/1 weighted average, rounded with ties away from zero, then thresholded).
  ! Mirrors VTTrac.jl's get_zsub_visible_subgrid.
  !===========================================================================
  subroutine read_template_visible(z, nx, ny, nt, visible, nsx, nsy, tid, xc, yc, &
       use_zmiss, zmiss, tmpl, vis, ok)
    integer, intent(in) :: nx, ny, nt, nsx, nsy, tid
    real(sp), intent(in) :: z(nx, ny, nt)
    logical, intent(in) :: visible(nx, ny, nt)
    real(dp), intent(in) :: xc, yc
    logical, intent(in) :: use_zmiss
    real(sp), intent(in) :: zmiss
    real(sp), intent(out) :: tmpl(nsx, nsy)
    logical, intent(out) :: vis(nsx, nsy)
    logical, intent(out) :: ok

    integer :: xi, yi, x0i, y0i, i, k, isx, isy
    real(sp) :: dx0, dx1, dy0, dy1, vsum
    real(dp) :: ddx, ddy

    xi = nint(xc); yi = nint(yc)
    ddx = xc - real(xi, dp); ddy = yc - real(yi, dp)
    isx = 0; if (ddx > 0.0_dp) isx = 1; if (ddx < 0.0_dp) isx = -1
    isy = 0; if (ddy > 0.0_dp) isy = 1; if (ddy < 0.0_dp) isy = -1
    dx0 = real(abs(ddx), sp); dx1 = 1.0_sp - dx0
    dy0 = real(abs(ddy), sp); dy1 = 1.0_sp - dy0

    x0i = xi - nsx / 2
    y0i = yi - nsy / 2
    ok = .not. (x0i + min(0, isx) < 1 .or. x0i + nsx - 1 + abs(isx) > nx .or. &
                y0i + min(0, isy) < 1 .or. y0i + nsy - 1 + abs(isy) > ny)
    if (.not. ok) return

    if (use_zmiss) then
      if (any(z(x0i:x0i + nsx - 1 + abs(isx), y0i:y0i + nsy - 1 + abs(isy), tid) == zmiss)) then
        ok = .false.
        return
      end if
    end if

    if (isx == 0 .and. isy == 0) then
      tmpl = z(x0i:x0i + nsx - 1, y0i:y0i + nsy - 1, tid)
      vis = visible(x0i:x0i + nsx - 1, y0i:y0i + nsy - 1, tid)
    else
      do k = 1, nsy
        do i = 1, nsx
          tmpl(i, k) = dx1 * dy1 * z(x0i + i - 1,      y0i + k - 1,      tid) &
                     + dx0 * dy1 * z(x0i + i - 1 + isx, y0i + k - 1,      tid) &
                     + dx1 * dy0 * z(x0i + i - 1,       y0i + k - 1 + isy, tid) &
                     + dx0 * dy0 * z(x0i + i - 1 + isx, y0i + k - 1 + isy, tid)

          vsum = 0.0_sp
          if (visible(x0i + i - 1,      y0i + k - 1,      tid)) vsum = vsum + dx1 * dy1
          if (visible(x0i + i - 1 + isx, y0i + k - 1,      tid)) vsum = vsum + dx0 * dy1
          if (visible(x0i + i - 1,      y0i + k - 1 + isy, tid)) vsum = vsum + dx1 * dy0
          if (visible(x0i + i - 1 + isx, y0i + k - 1 + isy, tid)) vsum = vsum + dx0 * dy0
          vis(i, k) = nint(vsum) /= 0
        end do
      end do
    end if
    ok = any(vis)
  end subroutine read_template_visible

  !===========================================================================
  ! Demean a template and its (population, uncorrected) standard deviation.
  !===========================================================================
  pure subroutine demean(tmpl, nsx, nsy, xd, sigx)
    integer, intent(in) :: nsx, nsy
    real(sp), intent(in) :: tmpl(nsx, nsy)
    real(sp), intent(out) :: xd(nsx, nsy), sigx
    real(sp) :: xm, s
    xm = sum(tmpl) / real(nsx * nsy, sp)
    xd = tmpl - xm
    s = sum(xd * xd) / real(nsx * nsy, sp)
    sigx = sqrt(s)
  end subroutine demean

  !===========================================================================
  ! Sliding cross-correlation (method=0) / normalized-covariance (method=1)
  ! score between a (demeaned) template and the image at `tid`, over the
  ! window offsets [k0,k1] x [l0,l1]. Mirrors VTTrac.jl's sliding_xcor /
  ! sliding_ncov. `ymean` is updated incrementally (add the incoming column,
  ! remove the outgoing one) between adjacent k -- this running update is
  ! numerically fine (unlike the variance), but vyy/xysum are always
  ! recomputed from scratch for the reason given at `demeaned_sums`.
  !
  ! The per-window vyy_sum/xysum accumulation is inlined here rather than
  ! calling `demeaned_sums` with a `z(...)` array-section actual argument:
  ! that section is strided (z's window is not contiguous when nx /= nsx),
  ! so gfortran must copy it into a contiguous temporary on every one of the
  ! (2*ixhw+1)*(2*iyhw+1) calls per template-step. Inlining removed that
  ! copy and closed a ~15% regression against a benchmark gate (Julia v2.0.0: 0.244s; 
  ! this function's un-inlined predecessor: ~0.257s; inlined: matches or beats Julia).
  !===========================================================================
  subroutine sliding_score(z, nx, ny, nt, nsx, nsy, xd, sigx, tid, method, &
       k0, k1, l0, l1, use_zmiss, zmiss, fmiss_sp, scr, ok)
    integer, intent(in) :: nx, ny, nt, nsx, nsy, tid, method, k0, k1, l0, l1
    real(sp), intent(in) :: z(nx, ny, nt), xd(nsx, nsy), sigx, zmiss, fmiss_sp
    logical, intent(in) :: use_zmiss
    real(sp), intent(out) :: scr(k1 - k0 + 1, l1 - l0 + 1)
    logical, intent(out) :: ok

    integer :: nsx2, nsy2, nsxy, kk, ll, lx0, ly0, i, k
    real(sp) :: ymean, vyy_sum, xysum, sigx2, yleft, yright, d

    nsx2 = nsx / 2; nsy2 = nsy / 2; nsxy = nsx * nsy
    sigx2 = sigx * sigx
    ok = .not. (k0 - nsx2 < 1 .or. k1 + nsx2 > nx .or. l0 - nsy2 < 1 .or. l1 + nsy2 > ny)
    if (.not. ok) return

    if (use_zmiss) then
      if (any(z(k0 - nsx2:k1 + nsx2, l0 - nsy2:l1 + nsy2, tid) == zmiss)) then
        ok = .false.
        return
      end if
    end if

    do ll = l0, l1
      ly0 = ll - nsy2
      lx0 = k0 - nsx2

      ymean = 0.0_sp
      do k = 1, nsy
        do i = 1, nsx
          ymean = ymean + z(lx0 + i - 1, ly0 + k - 1, tid)
        end do
      end do
      ymean = ymean / real(nsxy, sp)

      vyy_sum = 0.0_sp; xysum = 0.0_sp
      do k = 1, nsy
        do i = 1, nsx
          d = z(lx0 + i - 1, ly0 + k - 1, tid) - ymean
          vyy_sum = vyy_sum + d * d
          xysum = xysum + xd(i, k) * d
        end do
      end do
      call store_score(method, vyy_sum, xysum, nsxy, sigx, sigx2, fmiss_sp, scr(1, ll - l0 + 1))

      do kk = k0 + 1, k1
        lx0 = kk - nsx2
        yleft = 0.0_sp; yright = 0.0_sp
        do k = 1, nsy
          yleft = yleft + z(lx0 - 1, ly0 + k - 1, tid)
          yright = yright + z(lx0 + nsx - 1, ly0 + k - 1, tid)
        end do
        ymean = ymean + (yright - yleft) / real(nsxy, sp)

        vyy_sum = 0.0_sp; xysum = 0.0_sp
        do k = 1, nsy
          do i = 1, nsx
            d = z(lx0 + i - 1, ly0 + k - 1, tid) - ymean
            vyy_sum = vyy_sum + d * d
            xysum = xysum + xd(i, k) * d
          end do
        end do
        call store_score(method, vyy_sum, xysum, nsxy, sigx, sigx2, fmiss_sp, scr(kk - k0 + 1, ll - l0 + 1))
      end do
    end do
  end subroutine sliding_score

  pure subroutine store_score(method, vyy_sum, xysum, nsxy, sigx, sigx2, fmiss_sp, out_val)
    integer, intent(in) :: method, nsxy
    real(sp), intent(in) :: vyy_sum, xysum, sigx, sigx2, fmiss_sp
    real(sp), intent(out) :: out_val
    real(sp) :: vyy, vxy
    if (method == METHOD_XCOR) then
      vyy = vyy_sum / real(nsxy, sp)
      if (vyy <= 0.0_sp) then
        out_val = fmiss_sp
      else
        vxy = xysum / real(nsxy, sp)
        out_val = vxy / sqrt(vyy) / sigx
      end if
    else
      vxy = xysum / real(nsxy, sp)
      out_val = vxy / sigx2
    end if
  end subroutine store_score

  !===========================================================================
  ! Masked sliding score (xcor: sxy/sqrt(sxx*syy); ncov: sxy/sxx -- the same
  ! cov(x,y)/var(x) normalization as the unmasked path, per VTTrac.jl v2.0.0's
  ! fix: normalization must not depend on whether masking is active). Returns
  ! ok=false if every window has fewer than min_samples jointly-visible
  ! pixels.
  !
  ! Shortcut (mirrors VTTrac.jl's get_score_xcor_with_visible /
  ! get_score_ncov_with_visible): if the *search-time* super-region (spanning
  ! every window position that will be tested) is entirely visible, fall back
  ! to the plain unmasked computation, ignoring the template's own visibility
  ! altogether. This is checked only against the search time `tid`, not the
  ! template's time -- so a template with masked-out pixels can still be
  ! scored without any masking applied, if the search window happens to be
  ! fully visible. Faithfully reproducing this (rather than "improving" it
  ! to always honor template masking) is required for golden-data parity.
  !===========================================================================
  subroutine sliding_score_masked(z, nx, ny, nt, visible, nsx, nsy, tmpl, vis_tmpl, &
       tid, method, k0, k1, l0, l1, use_zmiss, zmiss, min_samples, fmiss_sp, scr, ok)
    integer, intent(in) :: nx, ny, nt, nsx, nsy, tid, method, k0, k1, l0, l1, min_samples
    real(sp), intent(in) :: z(nx, ny, nt), tmpl(nsx, nsy), zmiss, fmiss_sp
    logical, intent(in) :: visible(nx, ny, nt), vis_tmpl(nsx, nsy), use_zmiss
    real(sp), intent(out) :: scr(k1 - k0 + 1, l1 - l0 + 1)
    logical, intent(out) :: ok

    integer :: nsx2, nsy2, kk, ll, lx0, ly0, i, k, cnt_n
    real(sp) :: sx, sy, mx, my, ddx, ddy, sxx, syy, sxy
    real(sp) :: xd_full(nsx, nsy), sigx_full
    logical :: any_ok, jointly_visible

    nsx2 = nsx / 2; nsy2 = nsy / 2
    ok = .not. (k0 - nsx2 < 1 .or. k1 + nsx2 > nx .or. l0 - nsy2 < 1 .or. l1 + nsy2 > ny)
    if (.not. ok) return

    if (use_zmiss) then
      if (any(z(k0 - nsx2:k1 + nsx2, l0 - nsy2:l1 + nsy2, tid) == zmiss)) then
        ok = .false.
        return
      end if
    end if

    if (all(visible(k0 - nsx2:k1 + nsx2, l0 - nsy2:l1 + nsy2, tid))) then
      call demean(tmpl, nsx, nsy, xd_full, sigx_full)
      call sliding_score(z, nx, ny, nt, nsx, nsy, xd_full, sigx_full, tid, method, &
           k0, k1, l0, l1, .false., zmiss, fmiss_sp, scr, ok)
      return
    end if

    ! Inlined `masked_sums` (see the note on `sliding_score` above for why:
    ! avoids a strided-array-section copy-in on every one of the
    ! (2*ixhw+1)*(2*iyhw+1) window positions).
    any_ok = .false.
    do ll = l0, l1
      ly0 = ll - nsy2
      do kk = k0, k1
        lx0 = kk - nsx2

        cnt_n = 0; sx = 0.0_sp; sy = 0.0_sp
        do k = 1, nsy
          do i = 1, nsx
            if (vis_tmpl(i, k) .and. visible(lx0 + i - 1, ly0 + k - 1, tid)) then
              cnt_n = cnt_n + 1
              sx = sx + tmpl(i, k)
              sy = sy + z(lx0 + i - 1, ly0 + k - 1, tid)
            end if
          end do
        end do

        if (cnt_n < min_samples) then
          scr(kk - k0 + 1, ll - l0 + 1) = fmiss_sp
        else
          mx = sx / real(cnt_n, sp)
          my = sy / real(cnt_n, sp)
          sxx = 0.0_sp; syy = 0.0_sp; sxy = 0.0_sp
          do k = 1, nsy
            do i = 1, nsx
              jointly_visible = vis_tmpl(i, k) .and. visible(lx0 + i - 1, ly0 + k - 1, tid)
              if (jointly_visible) then
                ddx = tmpl(i, k) - mx
                ddy = z(lx0 + i - 1, ly0 + k - 1, tid) - my
                sxx = sxx + ddx * ddx
                syy = syy + ddy * ddy
                sxy = sxy + ddx * ddy
              end if
            end do
          end do
          if (method == METHOD_XCOR) then
            scr(kk - k0 + 1, ll - l0 + 1) = sxy / sqrt(sxx * syy)
          else
            scr(kk - k0 + 1, ll - l0 + 1) = sxy / sxx
          end if
          any_ok = .true.
        end if
      end do
    end do
    ok = any_ok
  end subroutine sliding_score_masked

  !===========================================================================
  ! Peak/prominence-inside-template screening. Mirrors
  ! VTTrac.jl's chk_zsub_peak_inside. Requires nsx>=3 and nsy>=3 (validated in
  ! Python before this is ever called with use_peak_th=.true.).
  !===========================================================================
  pure subroutine chk_peak_inside(tmpl, nsx, nsy, th, failed)
    integer, intent(in) :: nsx, nsy
    real(sp), intent(in) :: tmpl(nsx, nsy), th
    logical, intent(out) :: failed
    real(sp) :: side_max, side_min, inner_max, inner_min

    side_max = max(maxval(tmpl(:, 1)), maxval(tmpl(:, nsy)), &
                    maxval(tmpl(1, 2:nsy - 1)), maxval(tmpl(nsx, 2:nsy - 1)))
    side_min = min(minval(tmpl(:, 1)), minval(tmpl(:, nsy)), &
                    minval(tmpl(1, 2:nsy - 1)), minval(tmpl(nsx, 2:nsy - 1)))
    inner_max = maxval(tmpl(2:nsx - 1, 2:nsy - 1))
    inner_min = minval(tmpl(2:nsx - 1, 2:nsy - 1))

    if (inner_max > side_max + th * (inner_max - inner_min) .or. &
        inner_min < side_min - th * (inner_max - inner_min)) then
      failed = .false.
    else
      failed = .true.
    end if
  end subroutine chk_peak_inside

  !===========================================================================
  ! Find the score peak: integer location of the maximum (ties broken by
  ! *last* occurrence in (outer k, inner l) scan order, to match VTTrac.jl's
  ! `>=`-based findlast-like tie-break -- using `>` here would shift
  ! status/count on ties and break golden-data parity), then optional 5-point
  ! subgrid refinement. Mirrors VTTrac.jl's find_score_peak.
  !===========================================================================
  subroutine find_score_peak(scr, nk, nl, do_subgrid, use_gaussian, kpi, lpi, peak, ok)
    integer, intent(in) :: nk, nl
    real(sp), intent(in) :: scr(nk, nl)
    logical, intent(in) :: do_subgrid, use_gaussian
    real(dp), intent(out) :: kpi, lpi, peak
    logical, intent(out) :: ok

    integer :: ii, jj, ipi, jpi
    real(sp) :: best
    real(dp) :: c, ld, rd, bd, td, kp, lp
    logical :: sub_stat

    best = -huge(best)
    ipi = 0; jpi = 0
    do ii = 1, nk
      do jj = 1, nl
        if (.not. ieee_is_nan(scr(ii, jj)) .and. scr(ii, jj) >= best) then
          best = scr(ii, jj)
          ipi = ii; jpi = jj
        end if
      end do
    end do

    if (ipi == 0) then
      ok = .false.
      kpi = 0.0_dp; lpi = 0.0_dp; peak = 0.0_dp
      return
    end if

    ok = .not. (ipi == 1 .or. ipi == nk .or. jpi == 1 .or. jpi == nl)
    if (.not. ok) then
      kpi = 0.0_dp; lpi = 0.0_dp; peak = 0.0_dp
      return
    end if

    kpi = real(ipi, dp)
    lpi = real(jpi, dp)
    peak = real(best, dp)

    if (do_subgrid) then
      c = real(scr(ipi, jpi), dp)
      ld = real(scr(ipi - 1, jpi), dp)
      rd = real(scr(ipi + 1, jpi), dp)
      bd = real(scr(ipi, jpi - 1), dp)
      td = real(scr(ipi, jpi + 1), dp)
      if (use_gaussian) then
        call subgrid_peak_gaussian(c, ld, rd, bd, td, kp, lp, peak, sub_stat)
      else
        call subgrid_peak_epara(c, ld, rd, bd, td, kp, lp, peak, sub_stat)
      end if
      if (sub_stat) then
        ok = .false.
        kpi = 0.0_dp; lpi = 0.0_dp; peak = 0.0_dp
        return
      end if
      kpi = kpi + kp
      lpi = lpi + lp
    end if
  end subroutine find_score_peak

  !===========================================================================
  ! 5-point elliptic-paraboloid subgrid peak. Mirrors
  ! VTTrac.jl's find_subgrid_peak_5pt_epara.
  !===========================================================================
  pure subroutine subgrid_peak_epara(c, l, r, b, t, kp, lp, peak, stat)
    real(dp), intent(in) :: c, l, r, b, t
    real(dp), intent(out) :: kp, lp, peak
    logical, intent(out) :: stat
    real(dp) :: ld, rd, bd, td, p, q

    ld = l - c; rd = r - c; bd = b - c; td = t - c
    stat = .not. (ld <= 0.0_dp .and. rd <= 0.0_dp .and. bd <= 0.0_dp .and. td <= 0.0_dp &
            .and. (ld < 0.0_dp .or. rd < 0.0_dp) .and. (bd < 0.0_dp .or. td < 0.0_dp))
    if (stat) then
      kp = 0.0_dp; lp = 0.0_dp; peak = 0.0_dp
      return
    end if
    p = -(ld + rd) / 2.0_dp
    q = -(bd + td) / 2.0_dp
    kp = (rd - ld) / (4.0_dp * p)
    lp = (td - bd) / (4.0_dp * q)
    peak = c + p * kp * kp + q * lp * lp
  end subroutine subgrid_peak_epara

  !===========================================================================
  ! 5-point Gaussian subgrid peak (log-space elliptic paraboloid). Mirrors
  ! VTTrac.jl's find_subgrid_peak_5pt_gaus. Requires all five inputs positive.
  !===========================================================================
  pure subroutine subgrid_peak_gaussian(c, l, r, b, t, kp, lp, peak, stat)
    real(dp), intent(in) :: c, l, r, b, t
    real(dp), intent(out) :: kp, lp, peak
    logical, intent(out) :: stat
    real(dp) :: cc, lc, rc, bc, tc, p, q

    stat = .not. (c > 0.0_dp .and. l > 0.0_dp .and. r > 0.0_dp .and. b > 0.0_dp .and. t > 0.0_dp)
    if (stat) then
      kp = 0.0_dp; lp = 0.0_dp; peak = 0.0_dp
      return
    end if
    cc = log(c)
    lc = log(l) - cc; rc = log(r) - cc; bc = log(b) - cc; tc = log(t) - cc
    stat = .not. (lc <= 0.0_dp .and. rc <= 0.0_dp .and. bc <= 0.0_dp .and. tc <= 0.0_dp &
            .and. (lc < 0.0_dp .or. rc < 0.0_dp) .and. (bc < 0.0_dp .or. tc < 0.0_dp))
    if (stat) then
      kp = 0.0_dp; lp = 0.0_dp; peak = 0.0_dp
      return
    end if
    p = -(lc + rc) / 2.0_dp
    q = -(bc + tc) / 2.0_dp
    kp = (rc - lc) / (4.0_dp * p)
    lp = (tc - bc) / (4.0_dp * q)
    peak = exp(cc + p * kp * kp + q * lp * lp)
  end subroutine subgrid_peak_gaussian

end module pyvttrac_core
