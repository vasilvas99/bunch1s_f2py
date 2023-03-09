
!c------------------------------------------------------------------------------
!c             M O D E L S:
      subroutine gpmm2(Y, DY, M, bdef, p1, p2)
!f2py intent(in) Y
!f2py intent(in) M
!f2py intent(in) bdef
!f2py intent(in) p1
!f2py intent(in) p2
!f2py intent(out) DY
!c     13.09.02
!c corresponds to MMII, non-dimensionalized!
!c   Last correction of the model:
!c   It is changed to be in analogy with LW2 with two constants - K and U
!c   and two exponents - p and n, previous        versions of the model in the new
!c   frame should be classified in p = 0 and n_old = n_new - 1
!c   (c) V. Tonchev, 04.01.2012
!c)   In this program - after 28.03.2011
         integer M, I
         real*8 Y(*), DY(*), dk(M + 1), du(M + 1), dtem
         real*8 p1, p2
         real*8 bdef
!c calculate inverse of the distances and its third power:
         do i = 2, M
            dtem = Y(i) - Y(i - 1)
            dk(i) = 1.0d0/dtem**p1
            du(i) = 1.0d0/dtem**p2
         end do
         dtem = Y(1) - Y(M) + M*bdef
         dk(1) = 1.0d0/dtem**p1
         dk(M + 1) = dk(1)
         du(1) = 1.0d0/dtem**p2
         du(M + 1) = du(1)
!c
         do I = 1, M
            DY(I) = dk(i + 1) - dk(i) - du(i + 1) + du(i)
         end do
      end subroutine

!c------------------------------------------------------------------------
      subroutine g1smm(Y, DY, M, bdef, p1, p2)
!f2py intent(in) Y
!f2py intent(in) M
!f2py intent(in) bdef
!f2py intent(in) p1
!f2py intent(in) p2
!f2py intent(out) DY
!c     One-sided MM0 model non-dimensionalized!
!c   (c) V. Tonchev, 15.01.2012
!c)   In this program - after 28.03.2011
         INTEGER M, I
         REAL*8 Y(*), DY(*), dk(M + 1), du(M + 1), dtem
         real*8 p1, p2, bdef
!c calculate inverse of the distances and its third power:
         do i = 2, M
            dtem = Y(i) - Y(i - 1)
            dk(i) = dtem
            du(i) = 1.0d0/dtem**p2
         end do
         dtem = Y(1) - Y(M) + M*bdef
         dk(1) = dtem
         du(1) = 1.0d0/dtem**p2
         du(M + 1) = du(1)
!c
         do I = 1, M
            DY(I) = dk(i) - du(i + 1) + du(i)
         end do
      end subroutine

      subroutine g1slw(Y, DY, M, par)
!f2py intent(in) Y
!f2py intent(in) M
!f2py intent(in) par
!f2py intent(out) DY
         !c???????????????????????????????????
         !C Model:Popkov, Krug, PRB 73, 235430(2006)
         integer M
         !c up to 7.02.07 M was declared as integer AFTER the line below!
         real*8 Y(*), DY(*)
         real*8 one, two, three, be, U, bem, bep
         real*8 par(5)

         one = 1.0d0
         two = 2.0d0
         three = par(5) + 1
         be = par(1)
         U = par(2)
         bem = (one - be)/two
         bep = (one + be)/two

         DY(1) = bem*(Y(2) - Y(1)) + bep*(Y(1) - Y(M) + M) + &
                 &U*(3.0d0/((Y(1) - Y(M) + M)**three) - 3.0d0/((Y(2) - Y(1))**three) - &
                 &one/((Y(M) - Y(M - 1))**three) + one/((Y(3) - Y(2))**three))

         DY(2) = bem*(Y(3) - Y(2)) + bep*(Y(2) - Y(1)) + &
                 &U*(3.0d0/((Y(2) - Y(1))**three) - 3.0d0/((Y(3) - Y(2))**three) - &
                 &one/((Y(1) - Y(M) + M)**three) + one/((Y(4) - Y(3))**three))

         do i = 3, M - 2
            DY(i) = bem*(Y(i + 1) - Y(i)) + bep*(Y(i) - Y(i - 1)) + &
                    &U*(3.0d0/((Y(i) - Y(i - 1))**three) - 3.0d0/((Y(i + 1) - Y(i))**three) - &
                    &one/((Y(i - 1) - Y(i - 2))**three) + one/((Y(i + 2) - Y(i + 1))**three))
         end do

         DY(M - 1) = bem*(Y(M) - Y(M - 1)) + bep*(Y(M - 1) - Y(M - 2)) + &
                     &U*(3.0d0/((Y(M - 1) - Y(M - 2))**three) - 3.0d0/((Y(M) - Y(M - 1))**three) - &
                     &one/((Y(M - 2) - Y(M - 3))**three) + one/((Y(1) - Y(M) + M)**three)  &
                     &)

         DY(M) = bem*(Y(1) - Y(M) + M) + bep*(Y(M) - Y(M - 1)) +&
                     &U*(3.0d0/((Y(M) - Y(M - 1))**three) - 3.0d0/((Y(1) - Y(M) + M)**three) -&
                     &one/((Y(M - 1) - Y(M - 2))**three) + one/((Y(2) - Y(1))**three))
      end subroutine
