
! c------------------------------------------------------------------------------
! c             M O D E L S:
subroutine gpmm2(Y, DY, M, bdef, p1, p2)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) bdef
   intent(in) p1
   intent(in) p2
   real*8, intent(out)::DY(M)
   ! c     13.09.02
   ! c corresponds to MMII, non-dimensionalized!
   ! c   Last correction of the model:
   ! c   It is changed to be in analogy with LW2 with two constants - K and U
   ! c   and two exponents - p and n, previous        versions of the model in the new
   ! c   frame should be classified in p = 0 and n_old = n_new - 1
   ! c   (c) V. Tonchev, 04.01.2012
   ! c)   In this program - after 28.03.2011
   integer M, i
   real(8) dk(M + 1), du(M + 1), dtem
   real(8) p1, p2
   real(8) bdef
   ! c calculate inverse of the distances and its third power:
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
   ! c
   do i = 1, M
      DY(i) = dk(i + 1) - dk(i) - du(i + 1) + du(i)
   end do
end subroutine

! c------------------------------------------------------------------------
subroutine g1smm(Y, DY, M, bdef, p1, p2)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) bdef
   intent(in) p1
   intent(in) p2
   real*8, intent(out)::DY(M)
   ! c     One-sided MM0 model non-dimensionalized!
   ! c   (c) V. Tonchev, 15.01.2012
   ! c)   In this program - after 28.03.2011
   integer M, i
   real(8) dk(M + 1), du(M + 1), dtem
   real(8) p1, p2, bdef
   ! c calculate inverse of the distances and its third power:
   do i = 2, M
      dtem = Y(i) - Y(i - 1)
      dk(i) = dtem
      du(i) = 1.0d0/dtem**p2
   end do
   dtem = Y(1) - Y(M) + M*bdef
   dk(1) = dtem
   du(1) = 1.0d0/dtem**p2
   du(M + 1) = du(1)
   ! c
   do i = 1, M
      DY(i) = dk(i) - du(i + 1) + du(i)
   end do
end subroutine

subroutine g1slw(Y, DY, M, par)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) par
   real*8, intent(out)::DY(M)
   ! c???????????????????????????????????
   ! C Model:Popkov, Krug, PRB 73, 235430(2006)
   integer M
   ! c up to 7.02.07 M was declared as integer AFTER the line below!
   real(8) one, two, three, be, U, bem, bep
   real(8) par(5)

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

subroutine gkrug(Y, DY, M, par)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) par
   real*8, intent(out)::DY(M)
!c     Model: Popkov, Krug, PRB 73, 235430 (2006)
   integer M
!c         up to 7.02.07 M was declared as integer AFTER the line below!
   real*8 one, two, three, be, U, bem, bep
   real*8 par(5)

   one = 1.0d0
   two = 2.0d0
   three = par(5) + 1
   be = par(1)
   U = par(2)
   bem = (one - be)/two
   bep = (one + be)/two
!c
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
            &one/((Y(M - 2) - Y(M - 3))**three) + one/((Y(1) - Y(M) + M)**three) &
            &)

   DY(M) = bem*(Y(1) - Y(M) + M) + bep*(Y(M) - Y(M - 1)) + &
         &U*(3.0d0/((Y(M) - Y(M - 1))**three) - 3.0d0/((Y(1) - Y(M) + M)**three) - &
         &one/((Y(M - 1) - Y(M - 2))**three) + one/((Y(2) - Y(1))**three))
end subroutine

subroutine g_pk2(Y, DY, M, par)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) par
   real*8, intent(out)::DY(M)
!C     evolves from the Model of Popkov, Krug, PRB 73, 235430 (2006), LW2
!c     the stabilization part is the same but the destabilization part is the stabilization one with an opposite signe
!c     and evenually different power.
!c
   integer M, i
!c         up to 7.02.07 M was declared as integer AFTER the line below!
   real*8 one, two, U, K, en, ro
   real*8 par(5)

   one = 1.0d0
   two = 2.0d0
   en = par(5) + 1.0d0
   K = par(1)
   U = par(2)
   ro = par(4) + 1.0d0
!c      in the beginning we take the destabilizing power to be ro = 1
!c
!c
!c
   DY(1) = -K*(3.0d0/((Y(1) - Y(M) + M)**ro) - 3.0d0/((Y(2) - Y(1))**ro) - &
         &one/((Y(M) - Y(M - 1))**ro) + one/((Y(3) - Y(2))**ro) &
         &+ U*(3.0d0/((Y(1) - Y(M) + M)**en) - 3.0d0/((Y(2) - Y(1))**en) &
         &- one/((Y(M) - Y(M - 1))**en) + one/((Y(3) - Y(2))**en)))

   DY(2) = -K*(3.0d0/((Y(2) - Y(1))**ro) - 3.0d0/((Y(3) - Y(2))**ro) - &
         &one/((Y(1) - Y(M) + M)**ro) + one/((Y(4) - Y(3))**ro) &
         &+ U*(3.0d0/((Y(2) - Y(1))**en) - 3.0d0/((Y(3) - Y(2))**en) - &
         &one/((Y(1) - Y(M) + M)**en) + one/((Y(4) - Y(3))**en)))

   do i = 3, M - 2
      DY(i) = -K*(3.0d0/((Y(i) - Y(i - 1))**ro) - 3.0d0/((Y(i + 1) - Y(i))**ro) - &
            &one/((Y(i - 1) - Y(i - 2))**ro) + one/((Y(i + 2) - Y(i + 1))**ro)) &
            &+ U*(3.0d0/((Y(i) - Y(i - 1))**en) - 3.0d0/((Y(i + 1) - Y(i))**en) - &
            &one/((Y(i - 1) - Y(i - 2))**en) + one/((Y(i + 2) - Y(i + 1))**en))
   end do

   DY(M - 1) = -K*(3.0d0/((Y(M - 1) - Y(M - 2))**ro) - 3.0d0/((Y(M) - Y(M - 1))**ro) - &
            &one/((Y(M - 2) - Y(M - 3))**ro) + one/((Y(1) - Y(M) + M))**ro) &
            &+ U*(3.0d0/((Y(M - 1) - Y(M - 2))**en) - 3.0d0/((Y(M) - Y(M - 1))**en) &
            &- one/((Y(M - 2) - Y(M - 3))**en) + one/((Y(1) - Y(M) + M)**en))

   DY(M) = -K*(3.0d0/((Y(M) - Y(M - 1))**ro) - 3.0d0/((Y(1) - Y(M) + M)**ro) - &
            &one/((Y(M - 1) - Y(M - 2))**ro) + one/((Y(2) - Y(1))**ro)) &
            &+ U*(3.0d0/((Y(M) - Y(M - 1))**en) - 3.0d0/((Y(1) - Y(M) + M)**en) - &
            &one/((Y(M - 1) - Y(M - 2))**en) + one/((Y(2) - Y(1))**en))

end subroutine

subroutine g_mm0(Y, DY, M, par)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) par
   real*8, intent(out)::DY(M)
!c     13.09.02
!c corresponds to MMI as introduced by VT in 2002 but slightly changed
!c      to meet the analogy with the model of Popkov and Krug thus
!c      it still has 2 parameters - b and U, like in the PK model         plus
!c      the two powers r and n, note the difference - here
!c      in the stabilization part the terrace widths are raised to (n+1)
!c     (c) VT, May 2008, Lexington KY
   integer M, i
   real*8 d(M + 1), d3(M + 1), dtem
   real*8 on
   data on/1.0d0/
   real*8 U, en, b, bp1, bm1
   real*8 par(5)
!c     first we study the case ro = 1.0 as in the PK model
!c        ro=par(4)
   b = par(1)
   bp1 = (b + 1.0d0)/2.0d0
   bm1 = (1.0d0 - b)/2.0d0
   U = par(2)
   en = par(5) + 1.0d0
!c calculate inverse of the distances and its third power:

   do i = 2, M
      dtem = Y(i) - Y(i - 1)
!c         d(i)=dtem**ro        ! so instead of this line we have:
      d(i) = dtem
      d3(i) = on/dtem**en
   end do

   dtem = Y(1) - Y(M) + M
!c        d(1)=dtem**ro  ! here also:
   d(1) = dtem
   d(M + 1) = d(1)
   d3(1) = on/dtem**en
   d3(M + 1) = d3(1)
!c
   do i = 1, M
      DY(i) = bp1*d(i) + bm1*d(i + 1) - U*(d3(i + 1) - d3(i))
   end do
end subroutine

subroutine g_mm1(Y, DY, M, par)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) par
   real*8, intent(out)::DY(M)
!c     13.09.02
!c corresponds to MMI
   real*8 d(M + 1), d3(M + 1), dtem
   integer M, i
   real*8 on
   data on/1.0d0/
   real*8 ro, a3, en, beta
   real*8 par(5)
   ro = par(1)
   beta = par(2)
   a3 = par(3)
   en = par(5)
!c calculate inverse of the distances and its third power:
   do i = 2, M
      dtem = Y(i) - Y(i - 1)
      d(i) = dtem**ro
      d3(i) = on/dtem**en
   end do
   dtem = Y(1) - Y(M) + M
   d(1) = dtem**ro
   d(M + 1) = d(1)
   d3(1) = on/dtem**en
   d3(M + 1) = d3(1)
!c
   do i = 1, M
      DY(i) = d(i) + beta*d(i + 1) - a3*(d3(i + 1) - d3(i))
   end do
end subroutine


subroutine g_lw(Y, DY, M, par)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) par
   real*8, intent(out)::DY(M)

   real*8 F(M + 1), A(M+1), R(M+1)
   integer M, i
   real*8 on
   data on/1.0d0/
   real*8 p, n
   real*8 par(5)
   
   p = par(1)
   n = par(2)
   

   ! calculcate the forces F
   do i = 2, M
      F(i) = 1.0/(Y(i) - Y(i - 1))
   end do

   ! calculate the b.c.
   F(1) = 1.0/(Y(1) - Y(M) + M)
   F(M + 1) = F(1)

   do i = 2, M
      A(i) = F(i+1)**(p+1) - 2*F(i)**(p+1) + F(i-1)**(p+1)
      R(i) = -(F(i+1)**(n+1) - 2*F(i)**(n+1) + F(i-1)**(n+1))
   end do
   
   A(1) = F(2)**(p+1) - 2*F(1)**(p+1) + F(M+1)**(p+1)
   R(1) = -(F(2)**(n+1) - 2*F(1)**(n+1) + F(M+1)**(n+1))

   A(M+1) = F(1)**(p+1) - 2*F(M+1)**(p+1) + F(M)**(p+1)
   R(M+1) = -(F(1)**(n+1) - 2*F(M+1)**(n+1) + F(M)**(n+1))

   do i = 1, M
      DY(i) = A(i) + R(i)
   end do
end subroutine

subroutine gise2(Y, DY, M, par)
   real*8, intent(in)::Y(M)
   intent(in) M
   intent(in) par
   real*8, intent(out)::DY(M)

   integer M, i
   real*8 d(0:M + 1), d2(0:M + 1), d3(0:M + 2), Ybc(0:M + 1)
   real*8 on, tw, three
   data on, tw/1.d0, 2.d0/
   real*8 gama, l0l, dpl, dml, a3, gama3
   real*8 tp1, tp2, tp3, tp4, tp5, tp6
   real*8 di1, di, dn1, dn
   real*8 ft, st, dyi
   real*8 par(5)
   gama = par(1)
   l0l = par(2)
   dpl = par(3)
   dml = par(4)
   three = par(5) + 1.0d0
   a3 = l0l**three
   gama3 = a3*gama
!c calculate distances:
   do i = 1, M
      Ybc(i) = Y(i)
   end do
   Ybc(0) = y(M) - M
   Ybc(M + 1) = Y(1) + M
!c
   do i = 1, M + 1
      dyi = Ybc(i) - Ybc(i - 1)
      d(i) = dyi
      d2(i) = dyi*dyi/tw
      d3(i) = on/dyi**three
   end do
!c
   d(0) = d(M)
   d2(0) = d2(M)
   d3(0) = d3(M)
   d3(M + 2) = d3(2)
!c
   do i = 1, M
      di = d(i)
      di1 = d(i + 1)
      dn = on/(dpl + dml + di)
      dn1 = on/(dpl + dml + di1)
      tp1 = d3(i + 2)
      tp4 = d3(i + 1)
      tp2 = tw*tp4
      tp3 = d3(i)
      tp5 = tw*tp3
      tp6 = d3(i - 1)
      ft = (d2(i + 1) + dml*di1)*dn1
      ft = ft + (d2(i) + dpl*di)*dn
      st = -(tp2 - tp1 - tp3)*dn1 + (tp5 - tp4 - tp6)*dn
      DY(i) = ft + gama3*st
   end do
end subroutine
