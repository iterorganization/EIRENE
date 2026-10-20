      subroutine eirene_jacobi_eigenvalue ( n, lda, a, it_max, ldv, v,
     .                                     d, it_num, rot_num )
!subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser modifications of the classical
!    Jacobi rotation method with threshold pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!    back to fixed format FORTRAN: Detlev Reiter, Jan. 2019
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
!
!    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
!
!    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
!
      USE EIRMOD_PRECISION
      implicit none

      integer ( kind = 4 ) n
      integer ( kind = 4 ) lda
      integer ( kind = 4 ) ldv

!     real ( DP ) a(n,n)
      real ( DP ) a(lda,n)
      real ( DP ) bw(n)
      real ( DP ) c
      real ( DP ) d(n)
      real ( DP ) g
      real ( DP ) gapq
      real ( DP ) h
      integer ( kind = 4 ) i
      integer ( kind = 4 ) it_max
      integer ( kind = 4 ) it_num
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      integer ( kind = 4 ) l
      integer ( kind = 4 ) m
      integer ( kind = 4 ) p
      integer ( kind = 4 ) q
      integer ( kind = 4 ) rot_num
      real ( DP ) s
      real ( DP ) t
      real ( DP ) tau
      real ( DP ) term
      real ( DP ) termp
      real ( DP ) termq
      real ( DP ) theta
      real ( DP ) thresh
!     real ( DP ) v(n,n)
      real ( DP ) v(ldv,n)
      real ( DP ) w(n)
      real ( DP ) zw(n)

      do j = 1, n
       do i = 1, n
         v(i,j) = 0.0D+00
       end do
       v(j,j) = 1.0D+00
      end do

      do i = 1, n
       d(i) = a(i,i)
      end do

      bw(1:n) = d(1:n)
      zw(1:n) = 0.0D+00
      it_num = 0
      rot_num = 0

      do while ( it_num < it_max )

      it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
      thresh = 0.0D+00
      do j = 1, n
        do i = 1, j - 1
          thresh = thresh + a(i,j) ** 2
        end do
      end do

      thresh = sqrt ( thresh ) / real ( 4 * n, DP )

      if ( thresh == 0.0D+00 ) then
        exit
      end if

      do p = 1, n
        do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and.
     .       termp == abs ( d(p) ) .and.
     .       termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 /
     .         ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
      end do

      bw(1:n) = bw(1:n) + zw(1:n)
      d(1:n) = bw(1:n)
      zw(1:n) = 0.0D+00

      end do
!
!  Restore upper triangle of input matrix.
!
      do j = 1, n
        do i = 1, j - 1
          a(i,j) = a(j,i)
        end do
      end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
      do k = 1, n - 1

        m = k

        do l = k + 1, n
          if ( d(l) < d(m) ) then
            m = l
          end if
        end do

      if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

      end if

      end do

      return
      end subroutine eirene_jacobi_eigenvalue
