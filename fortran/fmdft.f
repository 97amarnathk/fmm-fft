C
C     Distributed 1-dimensional discrete Fourier transform using
C     fast multipole method.
C     Peter McCorquodale, January 21, 1997
C
C     First, each MPI process calls mfti to initialize parameters.
C     Then any number of calls may be made to mft to perform the transform
C     using the same length and number of processors.
C

      double complex function dotp(n, x, r)
C
C     Returns the dot product of complex n-vector x with real n-vector r.
C
      implicit none
      integer n, j
      double complex x(n)
      double precision r(n)

      dotp = x(1) * r(1)
      do j = 2, n
         dotp = dotp + x(j) * r(j)
      end do

      return
      end



      double complex function mn3(n, x, y, z, r)
C
C     Returns dot product of [x, y, z] with r for complex n-vectors x, y, z
C     and real 3n-vector r.
C
      implicit none
      integer n, j
      double complex x(n), y(n), z(n)
      double precision r(3*n)

      mn3 = x(1) * r(1)
      do j = 2, n
         mn3 = mn3 + x(j) * r(j)
      end do
      do j = 1, n
         mn3 = mn3 + y(j) * r(n + j)
      end do
      do j = 1, n
         mn3 = mn3 + z(j) * r(n + n + j)
      end do

      return
      end




      integer function log2(n)
C
C     Returns base-2 logarithm of n.
C
      implicit none
      integer       n, m

      m = n
      log2 = 0
      do while ( m .gt. 1)
         m = m / 2
         log2 = log2 + 1
      end do

      return
      end



      subroutine vx(a, x, mat, p, chnods, mm, acu)
C
C     Checked.  Works OK.
C
      implicit none
      integer p, j, k, l
      double precision a, x(p), th, acu(p), dcos, ac, mat(p, p)
      double precision chnods(p), dtan, mm(p)
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      ac = 3d0 * dtan(a)
      do j = 1, p
         th = -dcos(dfloat(2*j-2)/dfloat(2*(p-1))*pi) * (a+a)
         acu(j) = ac / dtan(0.5d0 * th)
      end do

      do k = 1, p
         mm(k) = chnods(k) - acu(k)
         do j = 1, k-1
            mm(k) = mm(k)*(chnods(k)-acu(j))/(chnods(k)-chnods(j))
         end do
         do j = k+1, p
            mm(k) = mm(k)*(chnods(k)-acu(j))/(chnods(k)-chnods(j))
         end do
      end do

      do j = 1, p
         do k = 1, p
            mat(j, k) = mm(j) / (x(k) - acu(j))
            do l = 1, j-1
               mat(j, k) = mat(j, k) * (x(k)-chnods(l)) / (x(k)-acu(l))
            end do
            do l = j+1, p
               mat(j, k) = mat(j, k) * (x(k)-chnods(l)) / (x(k)-acu(l))
            end do
         end do
      end do

      return
      end



      subroutine wx(a, x, xlen, mat, p, chnods, acu)
C     Checked.  Works OK.
      implicit none
      integer p, j, k, l, xlen
      double precision a, x(xlen), th, acu(p), ac, dcos, mat(p, xlen)
      double precision chnods(p), dtan, mm
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      ac = 1d0 / dtan(a)
      do j = 1, p
         th = -dcos(dfloat(2*j-2)/dfloat(2*(p-1))*pi) * (pi-6*a) + pi
         acu(j) = ac * dtan(0.5d0 * th)
         chnods(j) = -dcos(dfloat(2*j-1)/dfloat(2*p)*pi)
      end do

      do j = 1, p
         mm = chnods(j) - acu(j)
         do k = 1, j-1
            mm = mm * (chnods(j) - acu(k)) / (chnods(j) - chnods(k))
         end do
         do k = j+1, p
            mm = mm * (chnods(j) - acu(k)) / (chnods(j) - chnods(k))
         end do

         do k = 1, xlen
            mat(j, k) = mm / (x(k) - acu(j))
            do l = 1, j-1
               mat(j, k) = mat(j, k) * (x(k)-chnods(l)) / (x(k)-acu(l))
            end do
            do l = j+1, p
               mat(j, k) = mat(j, k) * (x(k)-chnods(l)) / (x(k)-acu(l))
            end do
         end do
      end do

      return
      end



      subroutine initg(nieach, s, initch, p, tang, chnods)
C     Checked.  Works OK.
      integer p, nieach, s, j, k
      double precision initch(p, nieach), th, dtan, a, tang(nieach)
      double precision chnods(p), ta3, pi
      parameter (pi = 3.141592653589793238d0)

      a = pi / dfloat(s+s)
      do j = 1, nieach
         th = dfloat(2*j-nieach-1) / dfloat(nieach) * a
         tang(j) = dtan(th)
      end do
      ta3 = 3d0 * dtan(a)
      do j = 1, nieach
         do k = 1, p
            initch(k, j) = (chnods(k) + ta3 * tang(j)) /
     #           (ta3 - chnods(k) * tang(j))
         end do
      end do

      return
      end



      subroutine shftf(lev, dir, mat, p, chnods, x, wkp, mm)
C     Checked.  Works OK.
      implicit none
      integer p, lev, dir, j
      double precision mat(p, p), a, dtan, ta3, td2, chnods(p)
      double precision x(p), dd, wkp(p), mm(p)
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      dd = dfloat(dir)
      a = pi / dfloat(2**lev)
      ta3 = 3d0 * dtan(a)
      td2 = dtan(0.5d0 * a)
      do j = 1, p
         x(j) = 3d0*td2* (chnods(j) + dd*ta3*td2) / 
     #        (ta3 - dd*td2*chnods(j))
      end do
      call vx(0.5d0 * a, x, mat, p, chnods, wkp, mm)

      return
      end



      subroutine flip(lev, shft, mat, p, chnods, x, wkp, mm)
C     Checked.  Works OK.
      implicit none
      integer p, lev, shft, j
      double precision mat(p, p), a, dtan, c, td2, chnods(p)
      double precision x(p), b, wkp(p), mm(p)
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      a = pi / dfloat(2**lev)
      b = dfloat(shft)
      c = dtan(0.5d0 * a)
      td2 = dtan(b * a)
      do j = 1, p
         x(j) = 3d0*c* (1d0 - td2*c*chnods(j)) / (td2 + c*chnods(j))
      end do
      call vx(0.5d0 * a, x, mat, p, chnods, mm, wkp)

      return
      end




      subroutine shftl(lev, dir, mat, p, chnods, x, acu)
C     Checked.  Works OK.
      implicit none
      integer p, lev, dir, j
      double precision mat(p, p), a, dtan, td2, chnods(p)
      double precision x(p), dd, c, acu(p)
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      dd = dfloat(dir)
      a = pi / dfloat(2**lev)
      c = dtan(0.25d0 * a)
      td2 = dtan(0.5d0 * a)
      do j = 1, p
         x(j) = (chnods(j) - dd) / (1d0/c + dd*c*chnods(j)) / td2
      end do
      call wx(0.5d0 * a, x, p, mat, p, chnods, acu)

      return
      end



      subroutine initmm(qr, np, nieach, t, initch, phi, p)
C
C     Multiplies qr by initch to give phi
C
      implicit none
      integer p, np, nieach, t
      integer box, sc, term, j
      double complex qr(nieach, np, t), phi(p, np-1, t)
      double precision initch(p, nieach)

      do box = 1, t
         do sc = 1, np-1
            do term = 1, p
               phi(term, sc, box) = qr(1, sc+1, box) * initch(term, 1)
               do j = 2, nieach
                  phi(term, sc, box) = phi(term, sc, box) + 
     #                 qr(j, sc+1, box) * initch(term, j)
               end do
            end do
         end do
      end do

      return
      end 




      subroutine evlmf(ffr, s, lq, vec, p, w, chnods, acu)
      implicit none
      integer p, s, lq, term
      double precision ffr, vec(p), a, th, x, w(p), acu(p), chnods(p)
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      a = (pi / dfloat(s)) * 0.5d0
      th = ffr * a
      x = dtan(th) / dtan(a)
      call wx(a, x, 1, w, p, chnods, acu)
      do term = 1, p
         vec(term) = w(term) / dfloat(lq)
      end do

      return
      end



      subroutine mpp(vecs, p, len, matpp, ans)
C
C     Multiply complex vecs(p, len) by real matpp(p, p)
C     and return result in ans(p, len).
C
      implicit none
      integer p, len, term, sc, k
      double precision matpp(p, p)
      double complex vecs(p, len), ans(p, len)

      do sc = 1, len
         do term = 1, p
            ans(term, sc) = vecs(1, sc) * matpp(1, term)
            do k = 2, p
               ans(term, sc) = ans(term, sc) + 
     #              vecs(k, sc) * matpp(k, term)
            end do
         end do
      end do

      return
      end



      subroutine cjall(n, x)
C
C     Takes complex conjugate of double complex vector x of length n.
C
      integer n, j
      double complex x(n)

      do j = 1, n
         x(j) = dconjg(x(j))
      end do

      return
      end



      subroutine mft(lq, qr, dir, np, myid, p, nieach, w, v, sz1, szwk)
C
C     integer lq:  length of qr (total length / np)
C     double complex qr:  input/output (share of it)
C
C     integer np:  number of procs
C     integer myid:  proc id
C     integer p:  number of terms in expansions
C     integer nieach:  number of particles in each box
C
C     double precision w:  non-overwritable work array of size at least:
C     sz1 + 8*p^2*log2(lq/nieach) + (3*nieach+2*p)*nieach*np + 2*np*(np+p)
C     where sz1 is size of non-overwritable work array for FFT routine
C
C     double precision v:  overwritable work array of size at least:
C     szwk + 8*p*(np-1)*(lq/nieach/np + (1+2*p)*log2(lq/nieach) + log2(np) + 2)
C          + 12*np*(nieach + 3) + 3*p
C     where szwk is size of overwritable work array for FFT routine
C
      implicit none
      include 'mpif.h'
      integer   p, lq, np, myid, n, nieach, t, s, sz1, szwk, dir
      integer   ia(57), ind, j, log2np
      integer   log2
      double complex qr(lq)
      double precision w(1), v(1)
      integer   mpierr

      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      s = lq / nieach
      t = s / np
      n = log2(s)
      log2np = log2(np)

      ind = 1
      ia(1) = ind
      ind = ind + sz1
      do j = 2, 9
         ia(j) = ind
         ind = ind + p*p*(n-2)
      end do
      ia(10) = ind
      ind = ind + p*p
      ia(11) = ind
      ind = ind + p*nieach
      ia(12) = ind
      ind = ind + p*(np-1)*nieach
      ia(13) = ind
      ind = ind + p*np/2
      ia(14) = ind
      ind = ind + 3*nieach*np/2
      ia(15) = ind
      ind = ind + 3*nieach*(nieach-1)*(np-1)
      ia(16) = ind
      ind = ind + 3*nieach*np/2
      ia(17) = ind
      ind = ind + p
      ia(18) = ind
      ind = ind + 2 * (np-1)*np

      ind = 1
      do j = 19, 20
         ia(j) = ind
         ind = ind + 2 * p*(np-1)*(2*t + log2np - 3)
      end do
      do j = 23, 24
         ia(j) = ind
         ind = ind + 2 * np
      end do
      do j = 25, 27
         ia(j) = ind
         ind = ind + 2 * p*(np-1)
      end do
      do j = 28, 29
         ia(j) = ind
         ind = ind + 2 * ((np-1)*(nieach+n*2*p))
      end do
      do j = 30, 31
         ia(j) = ind
         ind = ind + 2 * ((np-1)*(nieach+n*2*p)+np/2)
      end do
      ia(32) = ind
      ind = ind + 2 * nieach * np
      ia(33) = ind
      ind = ind + 2 * nieach * (np-1)
      do j = 34, 35
         ia(j) = ind
         ind = ind + 2 * np/2
      end do
      do j = 36, 45
         ia(j) = ind
         ind = ind + 2 * p*(np-1)
      end do
      do j = 46, 47
         ia(j) = ind
         ind = ind + 2 * p*(np-1)*2*n
      end do
      ia(48) = ind
      ind = ind + 2*lq
      ia(49) = ind
      ind = ind + 2 * (np-1)
      do j = 50, 51
         ia(j) = ind
         ind = ind + 2 * np/2
      end do
      ia(52) = ind
      ind = ind + nieach
      do j = 53, 55
         ia(j) = ind
         ind = ind + p
      end do
      do j = 56, 57
         ia(j) = ind
         ind = ind + 3 * np
      end do

      call mftint(qr, lq, np, myid, s, p, n, t, nieach, log2np,
     #     sz1, szwk, dir,
     #     w(ia(1)), w(ia(2)), w(ia(3)), w(ia(4)), w(ia(5)),
     #        w(ia(6)), w(ia(7)),
     #     w(ia(8)), w(ia(9)), w(ia(10)), w(ia(11)),
     #        w(ia(12)), w(ia(13)),
     #     w(ia(14)), w(ia(15)), w(ia(16)), w(ia(17)), w(ia(18)),
     #     v(ia(19)), v(ia(20)),
     #        v(ia(23)), v(ia(24)), v(ia(25)), v(ia(26)), v(ia(27)),
     #     v(ia(28)), v(ia(29)), v(ia(30)), v(ia(31)),
     #        v(ia(32)), v(ia(33)), v(ia(34)), v(ia(35)),
     #     v(ia(36)), v(ia(37)), v(ia(38)), v(ia(39)),
     #        v(ia(40)), v(ia(41)), v(ia(42)),
     #     v(ia(43)), v(ia(44)), v(ia(45)), v(ia(46)), v(ia(47)),
     #        v(ia(48)), v(ia(49)), v(ia(50)), v(ia(51)),
     #     v(ia(52)), v(ia(53)), v(ia(54)), v(ia(55)),
     #     v(ia(56)), v(ia(57)))

      return
      end



      subroutine mfti(lq, np, p, nieach, w, wt, sz1, szwk)
C
C     integer lq:  length of qr (total length / np)
C     integer np:  number of procs
C     integer p:  number of terms in expansions
C     integer nieach:  number of particles in each box
C
C     double precision w:  work array of size at least:
C     ...
C     double precision wt:  work array of size at least:
C     ...
C     where n = log2(s), lnp = log2(np), s = lq/nieach, t = s/np
C
      implicit none
      integer   lq, np, s, p, sz1, szwk
      integer   ind, ia(30), j, n, log2, nieach
      ! w(1) indicates that w is a variable size array (Fortran 77) :/
      double precision w(1), wt(1)

      s = lq / nieach
      n = log2(s)
      ind = 1
      ia(1) = ind
C      ind = ind + 4*lq+15
      ind = ind + sz1
      do j = 2, 9
         ia(j) = ind
         ind = ind + p*p*(n-2)
      end do
      ia(10) = ind
      ind = ind + p*p
      ia(11) = ind
      ind = ind + p*nieach
      ia(12) = ind
      ind = ind + p*(np-1)*nieach
      ia(13) = ind
      ind = ind + p*np/2
      ia(14) = ind
      ind = ind + 3*nieach*np/2
      ia(15) = ind
      ind = ind + 3*nieach*(nieach-1)*(np-1)
      ia(16) = ind
      ind = ind + 3*nieach*np/2
      ia(17) = ind
      ind = ind + p
      ia(18) = ind
      ind = ind + 2*(np-1)*np
      ia(19) = ind
      ind = ind + nieach
      do j = 20, 22
         ia(j) = ind
         ind = ind + p
      end do

      call mftii(lq, np, s, p, nieach, n, sz1, szwk,
     #     w(ia(1)), w(ia(2)), w(ia(3)), w(ia(4)), w(ia(5)),
     #     w(ia(6)), w(ia(7)), w(ia(8)), w(ia(9)), w(ia(10)),
     #     w(ia(11)), w(ia(12)), w(ia(13)), w(ia(14)), w(ia(15)),
     #     w(ia(16)), w(ia(17)), w(ia(18)),
     #     w(ia(19)), w(ia(20)), w(ia(21)), w(ia(22)), wt)

      return
      end




      subroutine mftii(lq, np, s, p, nieach, n, sz1, szwk,
     #     wsave, shftfn, shftfp, flip2n, flip2p,
     #     flip3n, flip3p, shftln, shftlp, topflp,
     #     initch, evalm, evalmh, cotsh, cots,
     #     cotprv, chnods, fo,
     #     wknie, wkp, wkp2, wkp3, wt)
      implicit none
      integer lq, np, s, p, nieach, n, sz1, szwk
      integer t, lev, j, sc, sce, k, log2
      double precision wsave(sz1)
      double precision shftfn(p, p, n - 2), shftfp(p, p, n - 2)
      double precision flip2n(p, p, n - 2), flip2p(p, p, n - 2)
      double precision flip3n(p, p, n - 2), flip3p(p, p, n - 2)
      double precision shftln(p, p, n - 2), shftlp(p, p, n - 2)
      double precision topflp(p, p), initch(p, nieach)
      double precision evalm(p, np-1, nieach), evalmh(p, np/2)
      double precision cotsh(3*nieach, np/2)
      double precision cots(3*nieach, nieach-1, np-1)
      double precision cotprv(3*nieach, np/2)
      double precision chnods(p)
      double complex   fo(np-1, np)
      double precision wknie(nieach), wkp(p), wkp2(p), wkp3(p), wt(1)

      double precision dlq, gfrac, ffr, ang
      double precision dcos, dsin, dtan, dfloat
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      if (np .eq. 1) then
         call myffti(lq, wsave, wt, sz1, szwk)

         return
      endif
      
      dlq = dfloat(lq)
C     Total length is np*lq
      t = s / np
      nieach = lq / s
      n = log2(s)

      do j = 1, p
         chnods(j) = -dcos(dfloat(2*j-1)/dfloat(2*p)*pi)
      end do
      call initg(nieach, s, initch, p, wknie, chnods)
      call flip(2, +2, topflp, p, chnods, wkp, wkp2, wkp3)
      do lev = 1, n-2
         call shftf(lev+2, -1, shftfn(1, 1, lev), p,
     #        chnods, wkp, wkp2, wkp3)
         call shftf(lev+2, +1, shftfp(1, 1, lev), p,
     #        chnods, wkp, wkp2, wkp3)
         call flip(lev+2, -2, flip2n(1, 1, lev), p,
     #        chnods, wkp, wkp2, wkp3)
         call flip(lev+2, +2, flip2p(1, 1, lev), p,
     #        chnods, wkp, wkp2, wkp3)
         call flip(lev+2, -3, flip3n(1, 1, lev), p,
     #        chnods, wkp, wkp2, wkp3)
         call flip(lev+2, +3, flip3p(1, 1, lev), p,
     #        chnods, wkp, wkp2, wkp3)
         call shftl(lev+1, -1, shftln(1, 1, lev), p,
     #        chnods, wkp, wkp2)
         call shftl(lev+1, +1, shftlp(1, 1, lev), p,
     #        chnods, wkp, wkp2)
      end do

      do sc = 1, np-1

         gfrac = dfloat(sc) / dfloat(np)

         do j = 1, nieach
            ffr = (dfloat(2*j-1) - dfloat(2*sc)/dfloat(np)) /
     #           dfloat(nieach) - 1d0
            call evlmf(ffr, s, lq, evalm(1, sc, j), p,
     #           wkp, chnods, wkp2)
         end do

         if (sc .ge. np/2) then
            sce = sc - np/2 + 1
            ffr = (dfloat(2*nieach+1) - dfloat(2*sc)/dfloat(np)) /
     #           dfloat(nieach) - 1d0
            call evlmf(ffr, s, lq, evalmh(1, sce), p,
     #           wkp, chnods, wkp2)
            do j = 1-nieach, 2*nieach
               cotsh(j+nieach, sce) = -1d0 / dtan(pi/dlq *
     #              (dfloat(j - 1 - nieach) + gfrac)) / dlq
            end do
         endif

         do k = 2, nieach
            do j = 1-nieach, 2*nieach
               cots(j+nieach, k-1, sc) = -1d0 / dtan(pi/dlq *
     #              (dfloat(j) - dfloat(k) + gfrac)) / dlq
            end do
         end do

         if (sc .le. np/2) then
            do j = 1-nieach, 2*nieach
               cotprv(j+nieach, sc) = -1d0 / dtan(pi/dlq *
     #              (dfloat(j) - 1d0 + gfrac)) / dlq
            end do
         endif

      end do

      call myffti(lq, wsave, wt, sz1, szwk)

      do sc = 1, np-1
         ang = pi * dfloat(sc) / dfloat(np)
         do j = 1, np
            fo(sc, j) = -zexp(dcmplx(0d0,ang*dfloat(1-2*j)))*dsin(ang)
         end do
      end do

      return
      end




      subroutine mftint(qr, lq, np, myid, s, p, n, t, nieach, log2np,
     #     sz1, szwk, dir,
     #     wsave, shftfn, shftfp, flip2n, flip2p, flip3n, flip3p,
     #     shftln, shftlp, topflp, initch, evalm, evalmh,
     #     cotsh, cots, cotprv, chnods, fo,
     #     phi, psi, sump, sumsec, phils, phirs, phr,
     #     packps, packnr, packns, packpr, qrn, qrp, qcpp, packnn,
     #     phiopp, phim2, phip2, phim3, phip3, shpsi, psiprv,
     #     f2n, f2p, f3, phin, phip, psiev, psit, exts, extr,
     #     wknie, wkp, wkp2, wkp3,
     #     prevd, nextd)
C
C     double complex qr:  input/output (share of it)
C     integer lq:  length of qr (total length / np)
C
C     integer np:  number of procs
C     integer myid:  proc id
C     integer s:  total number of boxes
C
      implicit none
      include 'mpif.h'
      integer   p, lq, np, myid, t, nieach, n, log2np, inc, nl
      integer   sc, sce, j, box, lev, npm1, base, s, term, lup
      integer   lpkp, lpkn, pk, sz1, szwk, dir
      integer   mpista(MPI_STATUS_SIZE), mpierr, tag
      integer   obase, intrvl, intvlo, lr
      integer   previd, nextid, cl, cr, oldnl, ind, pnpm1
      integer   log2

      double complex qr(nieach, np, t)

C     Arrays to be set by initialization routine
      double precision wsave(sz1)
      double precision shftfn(p, p, n - 2), shftfp(p, p, n - 2)
      double precision flip2n(p, p, n - 2), flip2p(p, p, n - 2)
      double precision flip3n(p, p, n - 2), flip3p(p, p, n - 2)
      double precision shftln(p, p, n - 2), shftlp(p, p, n - 2)
      double precision topflp(p, p), initch(p, nieach)
      double precision evalm(p, np-1, nieach), evalmh(p, np/2)
      double precision cotsh(3*nieach, np/2)
      double precision cots(3*nieach, nieach-1, np-1)
      double precision cotprv(3*nieach, np/2)
      double precision chnods(p)
      double complex   fo(np-1, np)

C     Other arrays
      double complex phi(p, np-1, 2*t + log2np - 3)
      double complex psi(p, np-1, 2*t + log2np - 3)
      double complex sump(np), sumsec(np)
      double complex phils(p, np-1), phirs(p, np-1)
      double complex phr(p, np-1)
      double complex packps((np-1)*(nieach+n*2*p))
      double complex packnr((np-1)*(nieach+n*2*p))
      double complex packns((np-1)*(nieach+n*2*p)+np/2)
      double complex packpr((np-1)*(nieach+n*2*p)+np/2)
      double complex qrn(nieach, np), qrp(nieach, np-1)
      double complex qcpp(np/2), packnn(np/2)
      double complex phiopp(p, np-1)
      double complex phim2(p, np-1), phip2(p, np-1)
      double complex phim3(p, np-1), phip3(p, np-1)
      double complex shpsi(p, np-1), psiprv(p, np-1)
      double complex f2n(p, np-1), f2p(p, np-1)
      double complex f3(p, np-1)
      double complex phin(p,np-1,2,n), phip(p,np-1,2,n)
      double complex psiev(nieach, t, np)
      double complex psit(np-1)
      double complex exts(np/2), extr(np/2)
      double precision wknie(nieach), wkp(p), wkp2(p), wkp3(p)
      integer   prevd(3, np), nextd(3, np)

      double precision dlq
      double precision pi
      parameter (pi = 3.141592653589793238d0)

      double complex psit1, psit2, psid1, psid2, psid3
      double complex dotp, mn3

C     Synchronize processes
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)

      if (dir .lt. 0) call cjall(lq, qr)

      if (np .eq. 1) then
         call myfft(lq, qr, wsave, psiev, sz1, szwk)
         if (dir .lt. 0) call cjall(lq, qr)
         return
      endif

      dlq = dfloat(lq)
C     Total length is np*lq
      log2np = log2(np)
      npm1 = np-1
      pnpm1 = p*npm1

      nextid = mod(myid + 1, np)
      previd = mod(myid + npm1, np)

      do lev = 1, log2np
         inc = np / (2**lev)
         do j = 1, 3
            nextd(j, lev) = mod(myid + j*inc, np)
            prevd(j, lev) = mod(myid + 3*np - j*inc, np)
         end do
      end do

      do box = 1, t
         ind = 1
         lr = 1
         do sc = 1, np
            do j = 1, nieach
               qrn(ind, lr) = qr(j, sc, box)
               lr = lr + 1
               if (lr .gt. np) then
                  lr = 1
                  ind = ind + 1
               endif
            end do
         end do
         do sc = 1, np
            do j = 1, nieach
               qr(j, sc, box) = qrn(j, sc)
            end do
         end do
      end do

C--------
C Get sum of each section.  Collect results in sumsec(1:npm1).
C--------

      do sc = 1, npm1
         sump(sc) = 0d0

         do box = 1, t
            do j = 1, nieach
               sump(sc) = sump(sc) + qr(j, sc+1, box)
            end do
         end do

      end do

      call MPI_REDUCE(sump, sumsec, npm1, MPI_DOUBLE_COMPLEX,
     #     MPI_SUM, 0, MPI_COMM_WORLD, mpierr)

      call MPI_BCAST(sumsec, npm1, MPI_DOUBLE_COMPLEX, 0,
     #     MPI_COMM_WORLD, mpierr)

      do sc = 1, npm1
         sumsec(sc) = sumsec(sc) * dcmplx(0d0, 1d0 / dlq)
      end do

C--------
C Step 1:  Lowest level
C--------

      base = t + log2np - 3
      call initmm(qr, np, nieach, t, initch, phi(1, 1, base+1), p)

C--------
C Step 2:  Up the tree
C--------

      lup = log2np
      if (np .eq. 2) lup = 2

      nl = t
C     At lower levels, no communication is required.
      do lev = n-1, lup, -1
         obase = base
         nl = nl/2
         base = base - nl

         do box = 1, nl
            call mpp(phi(1, 1, obase + 2*box-1), p, npm1,
     #           shftfn(1, 1, lev-1), phils)
            call mpp(phi(1, 1, obase + 2*box), p, npm1,
     #           shftfp(1, 1, lev-1), phirs)
            do sc = 1, npm1
               do term = 1, p
                  phi(term, sc, base + box) =
     #                 phils(term, sc) + phirs(term, sc)
               end do
            end do
         end do
      end do

C     At higher levels, communication is required.
      do lev = log2np-1, 2, -1
         obase = base
         base = base - 1
         intrvl = np / (2**lev)
         intvlo = intrvl/2
         if (mod(myid, intrvl) .eq. 0) then
            call MPI_RECV(phr, pnpm1, MPI_DOUBLE_COMPLEX,
     #           nextd(1,lev+1), lev, MPI_COMM_WORLD, mpista, mpierr)
            call mpp(phi(1, 1, obase + 1), p, npm1,
     #           shftfn(1, 1, lev-1), phils)
            call mpp(phr, p, npm1, shftfp(1, 1, lev-1), phirs)
            do sc = 1, npm1
               do term = 1, p
                  phi(term, sc, base + 1) =
     #                 phils(term, sc) + phirs(term, sc)
               end do
            end do
         elseif (mod(myid, intrvl) .eq. intvlo) then
            call MPI_SEND(phi(1,1,obase+1), pnpm1, MPI_DOUBLE_COMPLEX,
     #           prevd(1,lev+1), lev, MPI_COMM_WORLD, mpierr)
         endif
      end do   

      tag = log2np

C     base .eq. 0
     
C--------
C Pack:  phi, qr -> packns, packps
C--------

      pk = 0
      do lev = lup + 1, n
         nl = (2**lev) / np
         base = nl + log2np - 3

         do sc = 1, npm1
            do term = 1, p

               pk = pk + 1
               packns(pk) = phi(term, sc, base + nl-1)
               packps(pk) = phi(term, sc, base + 1)

               pk = pk + 1
               packns(pk) = phi(term, sc, base + nl)
               packps(pk) = phi(term, sc, base + 2)

            end do
         end do

      end do

C     Now pk .eq. (n - lup) * npm1 * 2 * p

      do sc = 1, npm1
         do j = 1, nieach
            pk = pk + 1
            packns(pk) = qr(j, sc+1, t)
            packps(pk) = qr(j, sc+1, 1)
         end do
      end do

C     Now pk .eq. (n - lup)*npm1*2*p + npm1*nieach

      lpkp = pk

      if (t .eq. 1) then
         do sc = np/2, npm1
            sce = sc - np/2 + 1
            packnn(sce) = dotp(nieach, qr(1,sc+1,1), cotsh(1,sce))
         end do
      else
         do sc = np/2, npm1
            sce = sc - np/2 + 1
            pk = pk + 1
            packns(pk) = dotp(nieach, qr(1,sc+1,t-1), cotsh(1,sce))
         end do
      endif

      lpkn = pk

C--------
C Data transfer.
C--------

      tag = tag + 1

      call MPI_SENDRECV(packps, lpkp, MPI_DOUBLE_COMPLEX, previd, tag,
     #     packnr, lpkp, MPI_DOUBLE_COMPLEX, nextid, tag,
     #     MPI_COMM_WORLD, mpista, mpierr)

      tag = tag + 1
      call MPI_SENDRECV(packns, lpkn, MPI_DOUBLE_COMPLEX, nextid, tag,
     #     packpr, lpkn, MPI_DOUBLE_COMPLEX, previd, tag,
     #     MPI_COMM_WORLD, mpista, mpierr)

      if (t .eq. 1) then
         tag = tag + 1
         call MPI_SENDRECV(packnn, np/2, 
     #        MPI_DOUBLE_COMPLEX, nextd(2, log2np), tag,
     #        qcpp, np/2,
     #        MPI_DOUBLE_COMPLEX, prevd(2, log2np), tag,
     #        MPI_COMM_WORLD, mpista, mpierr)

      endif

C--------
C Unpack:  packnr -> phin, qrn
C          packpr -> phip, qrp, qcpp
C--------

      pk = 0
      lr = 0
      do lev = lup + 1, n

         lr = lr + 1
         nl = (2**lev) / np
         base = nl + log2np - 3

         do sc = 1, npm1
            do term = 1, p

               pk = pk + 1
               phin(term,sc,1,lr) = packnr(pk)
               phip(term,sc,1,lr) = packpr(pk)

               pk = pk + 1
               phin(term,sc,2,lr) = packnr(pk)
               phip(term,sc,2,lr) = packpr(pk)

            end do
         end do

      end do

C     Now pk .eq. (n - lup) * npm1 * 2 * p

      do sc = 1, npm1
         do j = 1, nieach
            pk = pk + 1
            qrn(j, sc) = packnr(pk)
            qrp(j, sc) = packpr(pk)
         end do
      end do

C     Now pk .eq. (n - lup)*npm1*2*p + nieach*npm1

      if (t .gt. 1) then
         do sc = 1, np/2
            pk = pk + 1
            qcpp(sc) = packpr(pk)
         end do
      endif

C--------
C Step 3-4:  Down the tree
C--------

C-------- Top flip, which always requires communication

      if (np .eq. 2) then
         tag = tag + 1
         nl = 2
         do box = 1, 2
            tag = tag + 1
            call MPI_SENDRECV(phi(1, 1, box), pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextid, tag,
     #           phiopp, pnpm1, MPI_DOUBLE_COMPLEX, previd, tag,
     #           MPI_COMM_WORLD, mpista, mpierr)
            call mpp(phiopp, p, npm1, topflp, psi(1, 1, box))
         end do

      else
         tag = tag + 1
         nl = 1
         intrvl = np/4
         if (mod(myid, intrvl) .eq. 0) then
            call MPI_SENDRECV(phi(1, 1, 1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(1, 1), tag,
     #           phiopp, pnpm1, MPI_DOUBLE_COMPLEX,
     #           prevd(1, 1), tag,
     #           MPI_COMM_WORLD, mpista, mpierr)
            call mpp(phiopp, p, npm1, topflp, psi(1, 1, 1))
         endif
      endif

      base = 0

C-------- Higher levels requiring communication

      do lev = 3, log2np
         obase = base
         base = base + 1
         intvlo = intrvl
         intrvl = intvlo / 2
         tag = tag + 10
         
         if (mod(myid, intvlo) .eq. 0) then

            call MPI_SEND(psi(1, 1, obase+1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(1, lev), tag+1,
     #           MPI_COMM_WORLD, mpierr)

            call MPI_SENDRECV(phi(1, 1, base+1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(2, lev), tag+2,
     #           phim2, pnpm1,
     #           MPI_DOUBLE_COMPLEX, prevd(2, lev), tag+2,
     #           MPI_COMM_WORLD, mpista, mpierr)

            call MPI_SENDRECV(phi(1, 1, base+1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, prevd(2, lev), tag+3,
     #           phip2, pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(2, lev), tag+3,
     #           MPI_COMM_WORLD, mpista, mpierr)

            call MPI_SENDRECV(phi(1, 1, base+1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(3, lev), tag+5,
     #           phip3, pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(3, lev), tag+4,
     #           MPI_COMM_WORLD, mpista, mpierr)

            call mpp(psi(1,1,obase+1),p,npm1,shftlp(1,1,lev-2),shpsi)
            call mpp(phim2, p, npm1, flip2p(1, 1, lev-2), f2p)
            call mpp(phip2, p, npm1, flip2n(1, 1, lev-2), f2n)
            call mpp(phip3, p, npm1, flip3n(1, 1, lev-2), f3)
            
            do sc = 1, npm1
               do term = 1, p
                  psi(term, sc, base+1) = shpsi(term, sc) +
     #                 f2p(term,sc) + f2n(term,sc) + f3(term,sc)
               end do
            end do
            
         elseif (mod(myid, intvlo) .eq. intrvl) then

            call MPI_SENDRECV(phi(1, 1, base+1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, prevd(3, lev), tag+4,
     #           psiprv, pnpm1,
     #           MPI_DOUBLE_COMPLEX, prevd(1, lev), tag+1,
     #           MPI_COMM_WORLD, mpista, mpierr)

            call MPI_RECV(phim3, pnpm1,
     #           MPI_DOUBLE_COMPLEX, prevd(3, lev), tag+5,
     #           MPI_COMM_WORLD, mpista, mpierr)

            call MPI_SENDRECV(phi(1, 1, base+1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(2, lev), tag+6,
     #           phim2, pnpm1,
     #           MPI_DOUBLE_COMPLEX, prevd(2, lev), tag+6,
     #           MPI_COMM_WORLD, mpista, mpierr)

            call MPI_SENDRECV(phi(1, 1, base+1), pnpm1,
     #           MPI_DOUBLE_COMPLEX, prevd(2, lev), tag+7,
     #           phip2, pnpm1,
     #           MPI_DOUBLE_COMPLEX, nextd(2, lev), tag+7,
     #           MPI_COMM_WORLD, mpista, mpierr)

            call mpp(psiprv, p, npm1, shftln(1,1,lev-2), shpsi)
            call mpp(phim2, p, npm1, flip2p(1, 1, lev-2), f2p)
            call mpp(phip2, p, npm1, flip2n(1, 1, lev-2), f2n)
            call mpp(phim3, p, npm1, flip3p(1, 1, lev-2), f3)

            do sc = 1, npm1
               do term = 1, p
                  psi(term, sc, base+1) = shpsi(term, sc) +
     #                 f2p(term,sc) + f2n(term,sc) + f3(term,sc)
               end do
            end do

         endif
      end do


C-------- Lower levels not requiring communication

      lr = 0
      do lev = lup+1, n
         oldnl = nl
         nl = nl*2
         obase = base
         base = base + oldnl
         lr = lr + 1

         do box = 1, oldnl

C     Left child

            cl = base + 2*box-1

            call mpp(psi(1,1,obase+box),p,npm1,shftlp(1,1,lev-2),shpsi)

            if (box .eq. 1) then
               call mpp(phip(1,1,1,lr),p,npm1,flip2p(1,1,lev-2), f2p)
            else
               call mpp(phi(1,1,cl-2),p,npm1,flip2p(1,1,lev-2), f2p)
            endif

            if (box .eq. oldnl) then
               call mpp(phin(1,1,1,lr),p,npm1,flip2n(1,1,lev-2), f2n)
               call mpp(phin(1,1,2,lr),p,npm1,flip3n(1,1,lev-2), f3)
            else
               call mpp(phi(1,1,cl+2),p,npm1,flip2n(1,1,lev-2), f2n)
               call mpp(phi(1,1,cl+3),p,npm1,flip3n(1,1,lev-2), f3)
            endif

            do sc = 1, npm1
               do term = 1, p
                  psi(term, sc, cl) = shpsi(term, sc) +
     #                 f2p(term,sc) + f2n(term,sc) + f3(term,sc)
               end do

            end do

C     Right child

            cr = base + 2*box

            call mpp(psi(1,1,obase+box),p,npm1,shftln(1,1,lev-2),shpsi)

            if (box .eq. 1) then
               call mpp(phip(1,1,1,lr),p,npm1,flip3p(1,1,lev-2), f3)
               call mpp(phip(1,1,2,lr),p,npm1,flip2p(1,1,lev-2), f2p)
            else
               call mpp(phi(1,1,cr-3),p,npm1,flip3p(1,1,lev-2), f3)
               call mpp(phi(1,1,cr-2),p,npm1,flip2p(1,1,lev-2), f2p)
            endif

            if (box .eq. oldnl) then
               call mpp(phin(1,1,2,lr),p,npm1,flip2n(1,1,lev-2), f2n)
            else
               call mpp(phi(1,1,cr+2),p,npm1,flip2n(1,1,lev-2), f2n)
            endif

            do sc = 1, npm1
               do term = 1, p
                  psi(term, sc, cr) = shpsi(term, sc) +
     #                 f3(term,sc) + f2p(term,sc) + f2n(term,sc)
               end do
            end do

         end do
      end do

C--------
C Step 5:  Evaluate local expansions.
C--------

C     Send psi(:, sc, base+t) * evalmh(:, sc) to nextid
      do sc = np/2, npm1
         sce = sc - np/2 + 1
         exts(sce) = dotp(p, psi(1, sc, base+t), evalmh(1, sce))
      end do
      tag = tag + 10

      call MPI_SEND(exts, np/2,
     #     MPI_DOUBLE_COMPLEX, nextid, tag,
     #     MPI_COMM_WORLD, mpierr)

      call MPI_RECV(extr, np/2,
     #     MPI_DOUBLE_COMPLEX, previd, tag,
     #     MPI_COMM_WORLD, mpista, mpierr)

      do sc = 1, np/2 - 1
         do j = 1, nieach
            do box = 1, t
               psiev(j, box, sc) = 
     #              dotp(p, psi(1, sc, base+box), evalm(1, sc, j))
            end do
         end do
      end do

      sc = np/2
         sce = 1
         j = 1
            box = 1
               psit1 = dotp(p, psi(1,sc,base+box), evalm(1,sc,j))
               psiev(j, box, sc) = 0.5d0 * (psit1 + extr(sce))
            do box = 2, t
               psit1 = dotp(p, psi(1,sc,base+box), evalm(1,sc,j))
               psit2 = dotp(p, psi(1,sc,base+box-1),
     #              evalmh(1, sce))
               psiev(j, box, sc) = 0.5d0 * (psit1 + psit2)
            end do
         do j = 2, nieach
            do box = 1, t
               psiev(j, box, sc) =
     #              dotp(p, psi(1,sc,base+box), evalm(1,sc,j))
            end do
         end do

      do sc = np/2 + 1, npm1
         sce = sc - np/2 + 1
         j = 1
            box = 1
               psiev(j, box, sc) = extr(sce)
            do box = 2, t
               psiev(j, box, sc) =
     #              dotp(p, psi(1,sc,base+box-1), evalmh(1,sce))
            end do
         do j = 2, nieach
            do box = 1, t
               psiev(j, box, sc) =
     #              dotp(p, psi(1,sc,base+box), evalm(1,sc,j))
            end do
         end do
      end do

C--------
C Step 6:  Direct evaluation.
C--------

      do sc = 1, np/2 - 1
         j = 1
         if (t .eq. 1) then
            box = 1
            psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box), qrn(1, sc),
     #           cotprv(1, sc))
         else
            box = 1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box), qr(1, sc+1, box+1),
     #           cotprv(1, sc))
            do box = 2, t-1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qr(1,sc+1,box-1), qr(1,sc+1,box), qr(1,sc+1,box+1),
     #              cotprv(1, sc))
            end do
            box = t
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qr(1, sc+1, box-1), qr(1, sc+1, box), qrn(1, sc),
     #           cotprv(1, sc))
         endif

         do j = 2, nieach
            if (t .eq. 1) then
               box = 1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qrp(1, sc), qr(1, sc+1, box), qrn(1, sc),
     #              cots(1, j-1, sc))
            else
            box = 1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box), qr(1, sc+1, box+1),
     #           cots(1, j-1, sc))
            do box = 2, t-1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qr(1,sc+1,box-1), qr(1,sc+1,box), qr(1,sc+1,box+1),
     #              cots(1, j-1, sc))
            end do
            box = t
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qr(1,sc+1,box-1), qr(1,sc+1,box), qrn(1, sc),
     #           cots(1, j-1, sc))
            endif
         end do
      end do

      sc = np/2
         sce = 1
         j = 1
            if (t .eq. 1) then
               box = 1
               psid1 = mn3(nieach,
     #              qrp(1, sc), qr(1, sc+1, box), qrn(1, sc),
     #              cotprv(1, sc))
               psid2 = dotp(nieach,
     #              qrp(1, sc), cotsh(1+nieach, sce))
               psid3 = dotp(nieach,
     #              qr(1, sc+1, box), cotsh(1+2*nieach, sce))
               psiev(j, box, sc) = psiev(j, box, sc) + 
     #              (psid1 + qcpp(sce) + psid2 + psid3) * 0.5d0

            elseif (t .eq. 2) then

            box = 1
               psid1 = mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box), qr(1, sc+1, box+1),
     #           cotprv(1, sc))
               psid2 = dotp(nieach,
     #              qrp(1, sc), cotsh(1+nieach, sce))
               psid3 = dotp(nieach,
     #              qr(1, sc+1, box), cotsh(1+2*nieach, sce))
               psiev(j, box, sc) = psiev(j, box, sc) +
     #              (psid1 + qcpp(sce) + psid2 + psid3) * 0.5d0
            box = 2
               psid1 = mn3(nieach,
     #           qr(1, sc+1, box-1), qr(1, sc+1, box), qrn(1, sc),
     #           cotprv(1, sc))
               psid2 = mn3(nieach,
     #              qrp(1, sc), qr(1, sc+1, box-1), qr(1, sc+1, box),
     #              cotsh(1, sce))
               psiev(j, box, sc) = psiev(j, box, sc) +
     #              (psid1 + psid2) * 0.5d0

            else

            box = 1
               psid1 = mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box), qr(1, sc+1, box+1),
     #           cotprv(1, sc))
               psid2 = dotp(nieach,
     #              qrp(1, sc), cotsh(1+nieach, sce))
               psid3 = dotp(nieach,
     #              qr(1, sc+1, box), cotsh(1+2*nieach, sce))
               psiev(j, box, sc) = psiev(j, box, sc) + 
     #              (psid1 + qcpp(sce) + psid2 + psid3) * 0.5d0
            box = 2
               psid1 = mn3(nieach,
     #           qr(1,sc+1,box-1), qr(1,sc+1,box), qr(1,sc+1,box+1),
     #           cotprv(1, sc))
               psid2 = mn3(nieach,
     #              qrp(1, sc), qr(1, sc+1, box-1), qr(1, sc+1, box),
     #              cotsh(1, sce))
               psiev(j, box, sc) = psiev(j, box, sc) + 
     #              (psid1 + psid2) * 0.5d0
            do box = 3, t-1
               psid1 = mn3(nieach,
     #           qr(1,sc+1,box-1), qr(1,sc+1,box), qr(1,sc+1,box+1),
     #           cotprv(1, sc))
               psid2 = mn3(nieach,
     #              qr(1,sc+1,box-2), qr(1,sc+1,box-1), qr(1,sc+1,box),
     #              cotsh(1, sce))
               psiev(j, box, sc) = psiev(j, box, sc) +
     #              (psid1 + psid2) * 0.5d0
            end do
            box = t
               psid1 = mn3(nieach,
     #           qr(1, sc+1, box-1), qr(1, sc+1, box), qrn(1, sc),
     #           cotprv(1, sc))
               psid2 = mn3(nieach,
     #              qr(1,sc+1,box-2), qr(1,sc+1,box-1), qr(1,sc+1,box),
     #              cotsh(1, sce))
               psiev(j, box, sc) = psiev(j, box, sc) +
     #              (psid1 + psid2) * 0.5d0
            endif

         do j = 2, nieach
            if (t .eq. 1) then
               box = 1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qrp(1, sc), qr(1, sc+1, box), qrn(1, sc),
     #              cots(1, j-1, sc))
            else
            box = 1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box), qr(1, sc+1, box+1),
     #           cots(1, j-1, sc))
            do box = 2, t-1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qr(1,sc+1,box-1), qr(1,sc+1,box), qr(1,sc+1, box+1),
     #              cots(1, j-1, sc))
            end do
            box = t
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qr(1, sc+1, box-1), qr(1, sc+1, box), qrn(1, sc),
     #           cots(1, j-1, sc))
            endif
         end do

      do sc = np/2 + 1, npm1
         sce = sc - np/2 + 1
         j = 1
            if (t .eq. 1) then
               box = 1
               psid1 = dotp(nieach, qrp(1, sc), cotsh(nieach + 1, sce))
               psid2 = dotp(nieach,qr(1,sc+1,box),cotsh(2*nieach+1,sce))
               psiev(j, box, sc) = psiev(j, box, sc) +
     #              qcpp(sce) + psid1 + psid2
            elseif (t .eq. 2) then
               box = 1
                  psid1 = dotp(nieach, qrp(1,sc), cotsh(nieach+1,sce))
                  psid2 = dotp(nieach, qr(1,sc+1,box), 
     #                 cotsh(2*nieach+1,sce))
                  psiev(j, box, sc) = psiev(j, box, sc) +
     #                 qcpp(sce) + psid1 + psid2
               box = 2
                  psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qrp(1, sc), qr(1, sc+1, box-1), qr(1, sc+1, box),
     #              cotsh(1, sce))
            else
            box = 1
               psid1 = dotp(nieach, qrp(1,sc), cotsh(nieach+1,sce))
               psid2 = dotp(nieach, qr(1,sc+1,box),
     #              cotsh(2*nieach+1,sce))
               psiev(j, box, sc) = psiev(j, box, sc) +
     #              qcpp(sce) + psid1 + psid2
            box = 2
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box-1), qr(1, sc+1, box),
     #           cotsh(1, sce))
            do box = 3, t-1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qr(1,sc+1,box-2), qr(1,sc+1,box-1), qr(1,sc+1,box),
     #              cotsh(1, sce))
            end do
            box = t
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qr(1,sc+1,box-2), qr(1,sc+1,box-1), qr(1,sc+1,box),
     #           cotsh(1, sce))
            endif

         do j = 2, nieach
            if (t .eq. 1) then
               box = 1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qrp(1, sc), qr(1, sc+1, box), qrn(1, sc),
     #              cots(1, j-1, sc))
            else
            box = 1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qrp(1, sc), qr(1, sc+1, box), qr(1, sc+1, box+1),
     #           cots(1, j-1, sc))
            do box = 2, t-1
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #              qr(1,sc+1,box-1), qr(1,sc+1,box), qr(1,sc+1,box+1),
     #              cots(1, j-1, sc))
            end do
            box = t
               psiev(j, box, sc) = psiev(j, box, sc) + mn3(nieach,
     #           qr(1, sc+1, box-1), qr(1, sc+1, box), qrn(1, sc),
     #           cots(1, j-1, sc))
            endif
         end do
      end do

      do box = 1, t
         do pk = 1, nieach
            do j = 1, npm1
               psit(j) = psiev(pk, box, j)
            end do
            do j = 1, np
               psiev(pk, box, j) = qr(pk, 1, box)
               do sc = 1, npm1
                  psiev(pk, box, j) = psiev(pk, box, j) +
     #                 fo(sc, j) * (psit(sc) - sumsec(sc))
               end do
            end do
         end do
      end do

      call MPI_ALLTOALL(psiev, t*nieach,
     #     MPI_DOUBLE_COMPLEX, qr, t*nieach,
     #     MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, mpierr)

      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)

      call myfft(lq, qr, wsave, psiev, sz1, szwk)

      if (dir .lt. 0) call cjall(lq, qr)

      return
      end
