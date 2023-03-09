    subroutine Bunch_1s(model)
c**************************************************************************         
c	Integrates 1D equation of step velocity using 4-th order Runge-Kuta  *
c      to obtain step positions for different models: 					 *
c      gcpcm - growth in the so called "C+ - C-" model					 *
c	 g_ise - growth with Inverse Schwoebel effect						 *
c	 giise - growth with Infinite Schwoebel effect (+ evaporation)		 *
c	 g_mm1, g_mm2 - minimal models									  	 *
c
c      gkrug - Model Popkov_Krug_PRB 73, 235430 (2006), 					 *
c              based on Model Liu_Weeks_PRB 57,23 (1998) => LW2					 *
c
c    !Passes many times along a decade with 9 time points per a decade!																		 *
c	(c) Vesselin Tonchev, Sofia-Clermont Ferrand-Sofia, 2001-2012		 *
c
c***************************************************************************
c                     The models are in these subroutines:
      EXTERNAL        !!! ! non-dimensionalized models .. with L0!
c
      INTEGER         M,IERR,I, Iorders, orders,  npasses
	character       model*5, par_in*12
	integer ipasses
	common /ipxx/ ipasses
c	for the purposes of DE13R:
	include 'integration.h'
	include	'arrays.h'
	real*8 wo(99), ho(99)
	real*8 mindo (99), mindlo
	common /orxx/ wo, ho, mindo
c
c     functions that serve the statistics:
	integer ir
	integer iwrite
	common /iwxx/ iwrite
	real*8 bdef, p1, p2
	real*8 T2PRINT, T3PRINT
c 	Bunch/ Terrace Definitions:
	common /bdxx/ bdef
	common /ppxx/ p1,p2
	character*12 w_h,tn_tm,st_tm,rut,minsz,minst,av2,av3
	common /fnlx/ w_h,tn_tm,st_tm,rut 
	common /fn0x/ minsz,minst,av2,av3
	character*9 alngdr,prfdr,trjdr
	common /drxx/ alngdr
	common /pfxx/ prfdr
	common /trxx/ trjdr
	real*8 on
 	common /xx/ on
	data on /1.0d0/
	character ans
	call startatTime
	
d
c	-----------------------------------------------------------
c       Define which model to study!!!!
c	model='gcpcm' ! 'C+ - C-' Model
c      model='gise2' ! Growth with inverse Schwoebel barrier
c 	model='g2pis' ! Two particle model with inverse Schwoebel barrier on precursors
c	model='g_mm0' ! Minimal Model 0 (MM0)
c	model='g_mm1' ! Minimal Model 1 (MMI)
c	model='gpmm2' ! Minimal Model 2 (MMII) ! Non-dimensionalized
c	model='g1smm'
c	 model = 'MC-An'
c      model='g1slw'
c	???
c      model='gkrug' ! Model Popkov_Krug_PRB 73,(2006),Liu_Weeks_PRB 57,23 (1998)
c      model = 'g_pk2' ! is similar to the previous but the second term with an inverted sign
c
c         the check for criticality ought to be revised and rebuild: ierr=crtchk(model)
c     ------------------------------------------------------------
c       Input the values for the integration routine
c	par_in=model//'_in.txt'
c      open(10,file=par_in)
c  	  
c	READ(10,*) ans	! ans in a character variable and ONLY if ans = 'y' the intial postion is
c                        read from a file
c	READ(10,*) bdef
c      READ(10,*) orders
c	READ(10,*) npasses
c      READ(10,*) TPRINT
c      IF (TPRINT.GT.XEND) TPRINT=XEND
c      READ(10,*) M
c	read(10,*) p1   ! the destabilizing power p in gpMM2
c	read(10,*) p2	! the stabilizing power    n in   gpMM2
c      read(10,*)p3	! null
c	read(10,*) p4	! null
c	read(10,*) p5	! null
c      close(10)
 	bdef = 20.0
      orders = 3
	npasses	= 1
      TPRINT = 1
c      IF (TPRINT.GT.XEND) TPRINT=XEND
      M = 1000
c	 read(10,*) p1   ! the destabilizing power p in gpMM2
c	 read(10,*) p2	! the stabilizing power    n in   gpMM2
c      read(10,*)p3	! null
c	  read(10,*) p4	! null
c	  read(10,*) p5	! null
c      close(10)
c	if (model.eq.'gpmm2') then
c	p1 = p1+1.0d0
c	p2 = p2+1.0d0
c	endif
c	if (model.eq.'g1smm') then
c	p2 = p2+1.0d0
c	endif

c	print *, p1, p2
c
c                                          End of Input operations
c   -----------------------------------------------------------------------
c                                       ... prepare the output operations:
	call clrfls
c        and arrays:
	call initar
c	prepare some files:
     	call fn_MS_I
      call fn_MS_II
c	-----------------------------------------------------------------------
c
c


C    -----------------------------------------------------------------------
	do ipasses = 1, npasses
C    PREPARE THE INITIAL RANDOMLY DEVIATED VICINAL STEP CONFIGURATION:
	call in_pos(Yn,M,model,0)
c    ------------------------------------------------------------------------
	XN = 0.0
	XK = TPRINT
      H=TPRINT
	T2PRINT=TPRINT
c
c                     Loop over the time
		iwrite=0       
	do Iorders = 1, orders
	do ir=1,9
c             One more step in time calling integration subroutine
      if (model.eq.'g_mm1') then
c      CALL DE13R(g_mm1,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)  
      elseif (model.eq.'g_mm1') then
      CALL DE13R(gcpcm,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
      endif	   
c      CALL DE13R(gise2,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
c      CALL DE13R(g2pis,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
c      CALL DE13R(g_mm0,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)

c      CALL DE13R(gpmm2,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
c      CALL DE13R(g1smm,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
c      CALL DE13R(gkrug,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
c       CALL DE13R(g_pk2,M,XN,YN,XK,HMIN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
c       CALL DE13R(g_pk2,M,H,YN,Y)
c	call MC_AN(Y, M)
      IF (IERR.NE.0) THEN
      WRITE(*,*) ' REDUCE START STEP AND ULTIMATE STEP!!!'
      STOP
      ENDIF
      WRITE(*,1030) XK,(Y(I),I=1,5)

	iwrite=iwrite+1
c
c     -------------------------------------------------------------
c     Perform the statistics over step coordinates:
	call dstat(Y,M)
c     --------------------------------------------------------------
c	The right end of the integration interval becomes left!:
      XN=XK
      XK=XK+T2PRINT
      DO I=1,M
         YN(I)=Y(I)
      ENDDO
c
c		Write the profile and the trajectories files
	call wrpro(M,Y,model,ir)
 	call wr_trj(Y,model)

cc
c	End loop over time:
	enddo	                   !!! ir
	T2PRINT = T2PRINT * 10
	XK = T2PRINT
	enddo					   !!! orders
	open (77, file = 'MS_I.dat')
	ilc = 0
	xlk = TPRINT
	T3PRINT = TPRINT
	ilr = 0
	do Ilorders = 1, orders
	do ilr=1,9
	ilc = ilc + 1
	hlo = ho(ilc)+1
	wlo = wo(ilc)
	mindlo = mindo (ilc)
      write(77,3077) xlk,hlo,wlo,wlo/hlo,mindlo ! Time, N, W, <lb>, lmin_g
	XLK = XLK + T3PRINT 
	enddo   !!! ilr
	T3PRINT = T3PRINT * 10
	XLK = T3PRINT
	enddo	!!! ilorders
	close (77)
	enddo                      !!! ipasses
c
3077	format(e15.8,1x,4(f12.4,1x))
	call itime(tt)
	itme = 3600*tt(1)+60*tt(2)+tt(3) - itme
	open (37, file = 'time.dat', access = 'append')
	write(37, 38) 'Stopping at',':   ',tt(1),':',tt(2),':',tt(3)
	write(37,*) ' execution time:', itme
	close(37)
      STOP
 1030 FORMAT(1PE12.3,5E12.4/(12X,5E12.4))
      end
*************************
c        Prepare some files:
	subroutine clrfls
	character*20    doscom
	character*9 alngdr
	common /drxx/	alngdr
	character*9 prfdr
	common /pfxx/ prfdr
	character*9 trjdr
	common /trxx/ trjdr

        	alngdr='.\along_\'
	
       prfdr='.\profil\'
c
        trjdr='.\trajct\'
c
	doscom='mkdir '//'along_'
	call system(doscom)		 
	doscom='mkdir '//'profil'
	call system(doscom)
	doscom='mkdir '//'trajct'
	call system(doscom)
	doscom='del '//alngdr//'*.dat'
	call system(doscom)
	doscom= 'del '//trjdr//'*.dat'
	call system(doscom)
	doscom='del '//prfdr//'*.dat'
	call system(doscom)
	return
	end
**********************
	subroutine wr_trj(Y,model)
	integer i
	real*8 Y(*),xk
	character*5  model
	character*12 fname, fnam2
	character*21 tname, tnam2
 	character*9 trjdr
 	common /trxx/ trjdr
	common /txx/ xk		
c
	fname='tr_'//model//'.dat'
	fnam2='tr-'//model//'.dat'
	tname=trjdr//fname
	tnam2=trjdr//fnam2
	open(10,file=tname,access='append')
	open(11,file=tnam2,access='append')
	write(10,1053) xk,(Y(i),i=1,60)
	write(11,1053) xk,(Y(i),i=61,120)
c     da zapiswame i traektoriite na sledwashtite 60 stapala... 
	close(10)
	close(11)
1053  format(e15.8,60(1x,f10.4))
c
	return
	end
*************************************************************
      subroutine startatTime
	call itime(tt)
	open (37, file = 'time.dat')
	write(37, 38) 'Starting at',':   ',tt(1),':',tt(2),':',tt(3)
	write(*, 38) 'Starting at',':   ',tt(1),':',tt(2),':',tt(3)
	close(37)
38    format(A,3(A,I2.2))
c------------------------------------------------------------
	subroutine wrpro(M,Y,model,ia)
	integer M,i,ia
	real*8 Y(*), tdel,tdel1,one
	character*5  model*5, sno*3	, fname*12
	character*9 prfdr
	character*21 tname
	common /pfxx/ prfdr
	one=1.0d0
	if(ia.lt.1000) then
	write(sno,'(i3.3)')ia
	else
	sno='999'
	endif

	fname=model//sno//'.dat'
	tname=prfdr//fname

	if (model(1:1).eq.'g') then
	open(11,file=tname)
	do i=1,M-1
	if (i.gt.1) then
	tdel=y(i)-y(i-1)
	tdel1=y(i+1)-y(i)
	write(11,1033) y(i),M-(i-1),one/tdel
	write(11,1033) y(i),M-i,one/tdel1
	endif
	enddo
	close(11)
	else
	open(11,file=tname)
	do i=1,M
	tdel=y(i)-y(i-1)
 	tdel1=y(i+1)-y(i)
	write(11,1033) y(i),i-1,one/tdel
	write(11,1033) y(i),i,one/tdel1
	enddo
	endif

 1033 format(1x, f14.6, 1x, I6, 1x, f14.6)
	return
	end
*************************************************************
	subroutine fn_MS_II
	character*12 minsz,minst,av2,av3
	common /fn0x/ minsz,minst,av2,av3
	character*9 alngdr
	common /drxx/	alngdr
	character*21 tname
	 	 minsz=	'!mn-sz.dat'
		 av2= 'av_2sz.dat'
	     av3= 'av_3sz.dat'
	tname=alngdr//minsz
	open(10,file=tname)
	write(10,*)
	close(10)
	tname=alngdr//av2
	open(10,file=tname)
	write(10,*)
	close(10)
      tname=alngdr//av3
	open(10,file=tname)
 	write(10,*)
	close(10)
c
	return
	end

      subroutine dstat(Y,M, bdef, ipasses, iwrite)
c      Distance Statistics 13.02
c     Designed for passing many times along decades
c       (C) Vesselin Tonchev, Clermont Ferrand - Sofia, 2001 -  JAN 2013
c	Here was separated all statistics done using the distances between the steps
c     which come as an input array Y(*)
c	Now it is additionally separated into 3 subroutines doing:
c       l_stat   - only "light" statistics  
c       rugo - calculates the roughness (rugosity) of the surface
c       h_stat   - the "heavy" statistics   
c c
c      On input 
c                          Real*8 array Y(*)
c                          Real*8  bdef - bunch /terrace definition  
c						 integer M - number of steps

c      On output - nothing!
	integer M
	 Real*8 Y(*),DY(M)
c	Bunch / Terrace Definitions:
**********************************************
	call MS_I(Y,M, bdef, ipasses, iwrite) ! prev lstat
c	ier=rugo(M, DY)
c	call MS_II(DY,M) MS - II is commented here because we seek MS I only from this code
**********************************************   
	return
	end
***
	subroutine MS_I(Y,M, bdef, ipasses, iwrite)
    ! All the arguments exchanged with MS_I called previously lstat are on INPUT:
c   Y is an 1-D array and contains the coordinates of the 
c   M steps
c   bdef is a key concept - it defines a distance between two neighboring steps 
c   ipasses gives the number of passes through the same time interval already done while
c   iwrite gives the number of the time step during the current pass
c   
c     the way the program writes in files is changed to have the time as a X axis.
c
 	integer M
 	Real*8 Y(*),DY(M) ! Here I am leaving DM - the distance between the steps to be only an internal array
c	Bunch / Terrace Definitions:
      real*8 bdef
c	common /bdxx/ bdef
	integer ipasses
c	common /ipxx/ ipasses
	integer iwrite
c	common /iwxx/ iwrite
c
	real*8 wo(99), ho(99)
	real*8 mindo (99)
	real*8 w, h
	real*8 mind
	common /orxx/ wo, ho, mindo
	real*8  stime
	real*8 dyi,maxbd,mintd,maxtd
c 
	real*8 avtd,tw,ntd
	integer Tno,i
	real*8 one,zero
c     h , w, mind
	data one,zero /	1.0d0, 0.0d0 /
	character*12 w_h,tn_tm,st_tm,rut
	common /fnlx/  w_h,tn_tm,st_tm,rut
	character*21 tname
	real*8 xk
	common /txx/ xk
	character*9 alngdr
	common /drxx/	alngdr
c
c                        ...initializations...
	w=zero
	h=zero
	tw=zero
	ntd=zero
	mind=bdef
	maxbd=zero
	mintd=M*bdef
      maxtd=zero
c

c          Compute the interstep distances, the minimal one and
c          some other quantities that characterize the bunching:
	do i=2,M
      dyi=Y(i)-Y(i-1)
	dy(i)=dyi
	if(dyi.le.bdef) then
c     counting the number of distances in bunches ...
	h=h+one
c	            ...and the total bunch width 
	w=w+dyi

	if (dyi.lt.mind) mind=dyi
	if (dyi.gt.maxbd) maxbd=dyi
	else
	tw=tw+dyi
	ntd=ntd+one
      if (dyi.lt.mintd) mintd=dyi
	if (dyi.gt.maxtd) maxtd=dyi
	endif
	enddo
	dyi=Y(1)-Y(M)+M*bdef
	dy(1)=dyi
	if(dyi.le.bdef) then
	w=w+dyi
	h=h+one
	if (dyi.lt.mind) mind=dyi
	if (dyi.gt.maxbd) maxbd=dyi
	else
	tw=tw+dyi
	ntd=ntd+one
      if (dyi.lt.mintd) mintd=dyi
	if (dyi.gt.maxtd) maxtd=dyi
	endif
	avtd=tw/ntd
	if(maxbd.gt.bdef) then
	 write(*,*) 'Alert:maxbd!'
	stop
	endif
	if(mintd.lt.bdef) then
	 write(*,*) 'Alert:mintd!'
	stop
	endif
c

**********************************************************************************************
c      Counting the terraces and their width:
	
	Tno=0
	if (dy(1).gt.bdef) Tno=1
	do i=2,M
	dyi=dy(i)
	if (dyi.gt.bdef) then
	  if(dy(i-1).lt.bdef) Tno=Tno+1
	endif	
	enddo
c     check if the terrace in the beginning spannes the boundaries
	if (dy(1).gt.bdef.and.dy(M).gt.bdef) Tno=Tno-1
c
	if (Tno.gt.0) then
	stime=dfloat(M)/dfloat(Tno)
	tw=tw/dfloat(Tno)
	h=h/dfloat(Tno)
	w=w/dfloat(Tno)
	h = ((ipasses -  1)*ho(iwrite)+h)/float(ipasses)
	ho(iwrite) = h
	w = ((ipasses -  1)*wo(iwrite)+w)/float(ipasses)
	wo(iwrite) = w
	mind = ((ipasses -  1)*mindo(iwrite)+mind)/float(ipasses)
	mindo(iwrite) = mind
*****
	tname=alngdr//w_h
 	open(10,file=tname,access='append')
      write(10,3077) xk, h+1,w,w/h,mind	  ! Time, N, W, <lb>, lmin_g
	close(10)
	tname=alngdr//st_tm
 	open(11,file=tname,access='append')
c	write(11,*)  
      write(11,3077)xk,stime,tw,avtd,tw/avtd
	close(11)
	tname=alngdr//tn_tm
 	open(11,file=tname,access='append')
c	write(11,*)  
      write(11,*)xk,tno
	close(11)
	else
	 print *, 'terrace alert!'
	endif
3075	format(4(f16.8,1x),e12.4)
3077	format(e15.8,1x,4(f12.4,1x))
c
	return	
      end
c -----------------------------------------------------------------------------------
***** Calculates the surface width (roughness):
      subroutine rugo(M, DY)
	integer M, i
	real*8 DY(M), ru, fre
	character*12 w_h,tn_tm,st_tm,rut
	common /fnlx/  w_h,tn_tm,st_tm,rut
	character*21 tname
	character*9 alngdr
	common /drxx/	alngdr
	real*8 xk
	common /txx/ xk

	ru=0.0d0
	do i=2,M
	ru=ru+(i-1)**2*DY(i)
	enddo
c still to think of the PBC and the part of the 1/M terrace	ru = sqrt(ru/(1.0d0*M))
c
	fre=0.0d0
	do i=1,M
	fre=fre+1.0d0/DY(i)**2
c	fre=fre+1.0d0/DY(i)**en
	enddo
	tname=alngdr//rut
	open(10,file=tname,access='append')
      write(10,3075) xk,ru,fre
	close(10)
3075	format(e12.4,1x,4(f16.8,1x))
	return	   
	end
*****
	subroutine fn_MS_I
	character*12 w_h,tn_tm,st_tm,rut
	common /fnlx/ w_h,tn_tm,st_tm,rut 
	character*21 tname
	character*9 alngdr
	real*8 wo(99), ho(99)
	real*8 mindo (99)
	common /orxx/ wo, ho, mindo
	common /drxx/	alngdr
	do iord=1,99
	 wo(iord) = 0.0d0
	 mindo(iord) = 0.0d0
	 ho(iord) = 0.0d0
	enddo
	 	 w_h='!!w-h__.dat'
		 tn_tm='!tn-tm__.dat'
	rut='!ru-tm__.dat'
		 st_tm='!st-tm__.dat'
	tname=alngdr//w_h  									 
	open(10,file=tname)
	write(10,*)
	close(10)
	tname=alngdr//rut  
	open(10,file=tname)
	write(10,*)
	close(10)
	tname=alngdr//tn_tm
	open(10,file=tname)
	write(10,*)
	close(10)
	tname=alngdr//st_tm
	open(10,file=tname)
	write(10,*)
	close(10)
	return
	end
*
	subroutine MS_II(DY,M)
	integer M
	Real*8 DY(M)
c	Bunch and Terrace Definitions:
      real*8 bdef
	common /bdxx/ bdef
	real*8 xk
	integer one,zero
	data one,zero /1,0/ 
	integer bno
	integer i,inc,bsz(M/2),ibw,bszi
	real*8 dyi,dyip
c
c         Average quantities:
	real*8 bw(M/2),bd(M/2),mindi(M/2),tbw,tmin,
     &fbd(M/2),lbd(M/2),tla,stime
c          c  
	include 'arrays.h'
c         integration time is:
	common /txx/ xk
c	real*8 par(5)
c	common /pxx/ par
	character*12 minsz,minst,av2,av3
	common /fn0x/ minsz,minst,av2,av3
	character*9 alngdr
	common /drxx/	alngdr
	character*21 tname,tname1
c     Counting the bunches and different characteristics
	 Bno=zero
c
	 dyi=dy(1)
      if (dyi.le.bdef) then
	  tla=dyi
	  tmin=dyi
	  tbw=dyi
	  Bno=one
        fbd(1)=dyi
	  inc=2
 	  do Ip=2,M
	   dyip=dy(Ip)
	   if (dyip.le.Bdef) then
		inc=inc+1
          if(dyip.lt.tmin) tmin=dyip
	    tbw=tbw+dyip
	    tla=dyip
		else       
	      goto 3029
	   endif
	   enddo
	 endif
3029  continue
	if (bno.eq.one) then
      bd(1)=tbw/dfloat(inc-1)
	bw(1)=tbw
	mindi(1)=tmin
	bsz(1)=inc
	lbd(1)=tla
	endif

	do i=2,M
	dyi=dy(i)
	if (dyi.le.Bdef) then
	   if(dy(i-1).gt.Bdef) then
	  tla=dyi
	  tmin=dyi
	  tbw=dyi
	  Bno=Bno+one
        fbd(Bno)=dyi
	  inc=2
	  do Ip=i+1,M
	  dyip=dy(Ip)
	   if (dyip.le.Bdef) then
		inc=inc+1
          if(dyip.lt.tmin) tmin=dyip
	    tbw=tbw+dyip
	    tla=dyip
	    else
          goto 3030
	   endif
	  enddo
c
 3030  continue
      bd(bno)=tbw/dfloat(inc-1)
	bw(bno)=tbw
	mindi(bno)=tmin
	bsz(bno)=inc
	lbd(bno)=tla

	endif
	endif
	enddo
c                  To account for the situation when the bunch crosses the
c                   boundary conditions:
	if (dy(1).le.Bdef. and .dy(M).le.Bdef) then
	ibw=bsz(1)+bsz(bno)-1
	bw(1)=bw(1)+bw(bno)
      bd(1)=bw(1)/dfloat(ibw-1)
	mindi(1)=Dmin1(mindi(1),mindi(bno))	
	bsz(1)=ibw
	fbd(1)=fbd(bno)
	bsz(bno)=0
	bd(bno)=0.0d0
 	bw(bno)=0.0d0
	mindi(bno)=0.0d0
  	lbd(bno)=0.0d0
 	bno=bno-1	
	endif
*****

*****	
	if (bno.gt.0) then
	stime=dfloat(M)/dfloat(bno)
	do i=1,bno
	bszi=bsz(i)
	if(bszi.gt.999) then
	write(*,*) 'Alert:Bunch Size exceeds 999'
	goto 1053
	endif
	tmin=mindi(i)
	avmn(bszi)=avmn(bszi)+tmin
	if(tmin.lt.absmin(bszi)) absmin(bszi)=tmin
 	tmin=fbd(i)
	avf(bszi)=avf(bszi)+tmin
	if(tmin.lt.fdmin(bszi)) fdmin(bszi)=tmin
  	tmin=bd(i)
	avbd(bszi)=avbd(bszi)+tmin
	if(tmin.lt.bdmin(bszi)) bdmin(bszi)=tmin
   	tmin=bw(i)
	avbw(bszi,1)=avbw(bszi,1)+tmin
	avbw(bszi,2)=avbw(bszi,2)+one
 	if(tmin.lt.bwmin(bszi)) bwmin(bszi)=tmin
	avl(bszi)=avl(bszi)+lbd(i)
1053  continue
	enddo
****      
	tname=alngdr//minsz
 	open(10,file=tname)
 	do i=2,300
	 dyi=absmin(i)
	 if(dyi.le.bdef)write(10,2036) i,dyi,fdmin(i),
     & bwmin(i), bdmin(i)
	enddo
	close(10)
******
      tname=alngdr//av2
	tname1=alngdr//av3
	open(10,file=tname)
	open(11,file=tname1)
	do i=2,300
	 dyi=avbw(i,2)
	 if(dyi.gt.zero) then
	 write(10,2036) i,avbw(i,1)/dyi,avbd(i)/dyi
     	 write(11,2036)i,avmn(i)/dyi,avf(i)/dyi,
     &	 avl(i)/dyi
	 endif
	enddo
	close(10)

	
	

	endif
********************************************************
c

3036  format(e14.7,1x,f12.5,1x,i3,1x,4(f7.4,1x),e12.6)
1036  format(e12.6,1x,4(f8.4,1x),i3)
2036  format(i3,1x,4(f8.4,1x),e15.7)
1037  format(e12.6,1x,2(f8.4,1x),i3)
1038  format(e12.6,1x,3(f8.4,1x),i3)
	return
	end

*******************************************************************************


	subroutine initar
	real*8 bdef
	common /bdxx/ bdef
	real*8 bwmin(999),bdmin(999),avbw(999,2),avmn(999),
     & avbd(999),avf(999),avl(999),absmin(999),
     &fdmin(999)
	common /mxx/ absmin,fdmin,bwmin,bdmin,avbw,avmn,
     & avbd,avf,avl
c
	 integer i
	real*8 zero
	zero=0.0d0
	do i=1,300
 	 absmin(i)=bdef
	 fdmin(i)=bdef
	 bdmin(i)=bdef
	 bwmin(i)=bdef*i
	enddo
	do i=1,300
	 avbw(i,1)=zero
	 avbw(i,2)=zero
 	 avmn(i)=zero
 	 avbd(i)=zero
	 avf(i)=zero
	 avl(i)=zero
	enddo
	return
	end

C ---------------------------------------------------------------------

c                   Prepare the initial profile deviated from the vicinal one
c                    randomly using the built-in random number generator
	subroutine in_pos(Yn,M,model,ans)
	integer M, Ifake
	character ans
      real*8 Yn(M), bdef
	common /bdxx/ bdef
	real fake, fake1
	integer itm(3), is1
	REAL*4          RN
      REAL*8          DPOSI
	character*5     model
	character*12 fname
	character*9 prfdr
	character*21 tname
	common /pfxx/ prfdr
C      PREPARE THE INITIAL CONDITIONS WITH THE RANDOM NUMBER GENERATOR:'
C        DPOSI IS THE MAXIMAL DEVIATION OF THE INITIAL POSITION
C        in both directions with respect to THE EQUILIBRIUM ONE:
	if (ans.ne.'y') then 
        DPOSI=0.025d0*bdef
c        reseed the rn generator
      CALL itime(itm)
	is1=(itm(1)+5)*(itm(2)+6)+(itm(2)+2)*(itm(3)+4)
	&+(itm(1)+3)*(itm(3)+8)
	rn=rand(is1)
c                                       give initial values for the step coodinates
      DO I=1,M
      YN(I)=(I-1)*bdef-dposi*(1.0d0-2.0d0*rand(0))
c 	
      END DO
	else
	
	open(11,file='initial_profil.dat')
	do i=1,M
	read(11,*)yn(i),Ifake, fake
	read(11,*)fake, Ifake, fake1
	enddo

	endif
 1033 format(1x, f14.6, 1x, I6, 1x, f14.6)
c
	return 
      end
c	**************************
	integer function crtchk(model)
	real*8 par(5),parcr,oupar
	common /pxx/ par
	character*5 model

	if(model.eq.'gcpcm') then
	 parcr=0.222222222d0
	 oupar=par(3)**3*(par(2)+1.0d0)/(1.0d0-par(2))
	if(oupar.ge.parcr) then
	write(*,*) ' Model ',model,' is not'
	write(*,*)' critical!'
	write(*,*)' The critical parameter is:'
	write(*,*)'0.222222222'
	write(*,*)'and you have: ',oupar
	stop 'sorry ...'
	endif
	endif	
	crtchk=1
	return
	end
C------------------------------------------------------------------------------
c             M O D E L S:
      SUBROUTINE gpmm2(Y,DY,M)
c     13.09.02
c corresponds to MMII, non-dimensionalized!
c   Last correction of the model:
c   It is changed to be in analogy with LW2 with two constants - K and U
c   and two exponents - p and n, previous	versions of the model in the new
c   frame should be classified in p = 0 and n_old = n_new - 1
c   (c) V. Tonchev, 04.01.2012
c)   In this program - after 28.03.2011
      INTEGER          M,I 
      REAL*8           Y(*),DY(*),dk(M+1),du(M+1),dtem
	real*8 on
	common /xx/ on
	real*8 p1,p2
	common /ppxx/ p1,p2
	real*8 bdef
	common /bdxx/ bdef 	
c calculate inverse of the distances and its third power:
      do i=2,M
	 dtem=Y(i)-Y(i-1)
	 dk(i)=on/dtem**p1
	 du(i)=on/dtem**p2
	enddo
	dtem=Y(1) - Y(M) + M*bdef
	dk(1)=on/dtem**p1
	dk(M+1)=dk(1)
	du(1) = on/dtem**p2
	du(M+1)= du(1)
c
      DO I=1,M
	DY(I)=dk(i+1)-dk(i)-du(i+1)+du(i)
       END DO
      RETURN
      END
c------------------------------------------------------------------------
      SUBROUTINE g1smm(Y,DY,M)
c     One-sided MM0 model non-dimensionalized!
c   (c) V. Tonchev, 15.01.2012
c)   In this program - after 28.03.2011
      INTEGER          M,I 
      REAL*8           Y(*),DY(*),dk(M+1),du(M+1),dtem
	real*8 on
	common /xx/ on
	real*8 p1,p2
	common /ppxx/ p1,p2
	real*8 bdef
	common /bdxx/ bdef 	
c calculate inverse of the distances and its third power:
      do i=2,M
	 dtem=Y(i)-Y(i-1)
	 dk(i)=dtem
	 du(i)=on/dtem**p2
	enddo
	dtem=Y(1) - Y(M) + M*bdef
	dk(1)=dtem
	du(1) = on/dtem**p2
	du(M+1)= du(1)
c
      DO I=1,M
	DY(I)=dk(i)-du(i+1)+du(i)
       END DO
      RETURN
      END
c-------------------------------------------------------
	SUBROUTINE g1slw(Y,DY,M)
c???????????????????????????????????
C     Model: Popkov, Krug, PRB 73, 235430 (2006)
	integer M
c	 up to 7.02.07 M was declared as integer AFTER the line below!
	real*8 Y(M), DY(M)
	real*8 one, two, three, be, U, bem, bep	 ! three = stepen na wzaimod. 
	real*8 par(5)
	common /pxx/ par				 ! shte chetem ot wunshen fail

	one = 1.0d0
	two = 2.0d0
	three = par(5) + 1 
	be = par(1) 
	U =  par(2)
	bem=(one - be)/two
	bep=(one + be)/two
c
	DY(1) = bem*(Y(2)-Y(1)) + bep*(Y(1)-Y(M)+M) +
     &U*(3.0d0/((Y(1)-Y(M)+M)**three) - 3.0d0/((Y(2)-Y(1))**three) -
     &one/((Y(M)-Y(M-1))**three) + one/((Y(3)-Y(2))**three))

	DY(2) = bem*(Y(3)-Y(2)) + bep*(Y(2)-Y(1)) +
     &U*(3.0d0/((Y(2)-Y(1))**three) - 3.0d0/((Y(3)-Y(2))**three) -
     &one/((Y(1)-Y(M)+M)**three) + one/((Y(4)-Y(3))**three))

      
	DO i = 3, M-2
	DY(i) = bem*(Y(i+1)-Y(i)) + bep*(Y(i)-Y(i-1)) +
     &U*(3.0d0/((Y(i)-Y(i-1))**three) - 3.0d0/((Y(i+1)-Y(i))**three) -
     &one/((Y(i-1)-Y(i-2))**three) + one/((Y(i+2)-Y(i+1))**three) )

	ENDDO

	
      DY(M-1) = bem*(Y(M)-Y(M-1)) + bep*(Y(M-1)-Y(M-2)) +
     &U*(3.0d0/((Y(M-1)-Y(M-2))**three) - 3.0d0/((Y(M)-Y(M-1))**three) -
     &one/((Y(M-2)-Y(M-3))**three) + one/((Y(1)-Y(M)+M)**three) 
     &)

	
      DY(M) = bem*(Y(1)-Y(M)+M) + bep*(Y(M)-Y(M-1)) +
     &U*(3.0d0/((Y(M)-Y(M-1))**three) - 3.0d0/((Y(1)-Y(M)+M)**three) -
     &one/((Y(M-1)-Y(M-2))**three) + one/((Y(2)-Y(1))**three)) 

	return
	end


*************************************************************

	SUBROUTINE gkrug(Y,DY,M)
C     Model: Popkov, Krug, PRB 73, 235430 (2006)
	integer M
c	 up to 7.02.07 M was declared as integer AFTER the line below!
	real*8 Y(M), DY(M)
	real*8 one, two, three, be, U, bem, bep	 ! three = stepen na wzaimod. 
	real*8 par(5)
	common /pxx/ par				 ! shte chetem ot wunshen fail

	one = 1.0d0
	two = 2.0d0
	three = par(5) + 1 
	be = par(1) 
	U =  par(2)
	bem=(one - be)/two
	bep=(one + be)/two
c
	DY(1) = bem*(Y(2)-Y(1)) + bep*(Y(1)-Y(M)+M) +
     &U*(3.0d0/((Y(1)-Y(M)+M)**three) - 3.0d0/((Y(2)-Y(1))**three) -
     &one/((Y(M)-Y(M-1))**three) + one/((Y(3)-Y(2))**three))

	DY(2) = bem*(Y(3)-Y(2)) + bep*(Y(2)-Y(1)) +
     &U*(3.0d0/((Y(2)-Y(1))**three) - 3.0d0/((Y(3)-Y(2))**three) -
     &one/((Y(1)-Y(M)+M)**three) + one/((Y(4)-Y(3))**three))

      
	DO i = 3, M-2
	DY(i) = bem*(Y(i+1)-Y(i)) + bep*(Y(i)-Y(i-1)) +
     &U*(3.0d0/((Y(i)-Y(i-1))**three) - 3.0d0/((Y(i+1)-Y(i))**three) -
     &one/((Y(i-1)-Y(i-2))**three) + one/((Y(i+2)-Y(i+1))**three) )

	ENDDO

	
      DY(M-1) = bem*(Y(M)-Y(M-1)) + bep*(Y(M-1)-Y(M-2)) +
     &U*(3.0d0/((Y(M-1)-Y(M-2))**three) - 3.0d0/((Y(M)-Y(M-1))**three) -
     &one/((Y(M-2)-Y(M-3))**three) + one/((Y(1)-Y(M)+M)**three) 
     &)

	
      DY(M) = bem*(Y(1)-Y(M)+M) + bep*(Y(M)-Y(M-1)) +
     &U*(3.0d0/((Y(M)-Y(M-1))**three) - 3.0d0/((Y(1)-Y(M)+M)**three) -
     &one/((Y(M-1)-Y(M-2))**three) + one/((Y(2)-Y(1))**three)) 

	return
	end


*************************************************************
*********************************************************
	SUBROUTINE g_pk2(Y,DY,M)
C     evolves from the Model of Popkov, Krug, PRB 73, 235430 (2006), LW2
c     the stabilization part is the same but the destabilization part is the stabilization one with an opposite signe
c     and evenually different power.
c
	integer M
c	 up to 7.02.07 M was declared as integer AFTER the line below!
	real*8 Y(M), DY(M)
	real*8 one, two, U, K, en, ro	 ! en = stepen na wzaimod. 
	real*8 par(5)
	common /pxx/ par				 ! shte chetem ot wunshen fail

	one = 1.0d0
	two = 2.0d0
	en = par(5) + 1.0d0 
	K = par(1) 
	U =  par(2)
 	ro = par(4)	+ 1.0d0
c      in the beginning we take the destabilizing power to be ro = 1
c	
c
c
	DY(1) =
     &-K*(3.0d0/((Y(1)-Y(M)+M)**ro) - 3.0d0/((Y(2)-Y(1))**ro) -
     &one/((Y(M)-Y(M-1))**ro) + one/((Y(3)-Y(2))**ro) 
     &+U*(3.0d0/((Y(1)-Y(M)+M)**en)-3.0d0/((Y(2)-Y(1))**en)
     &-one/((Y(M)-Y(M-1))**en) + one/((Y(3)-Y(2))**en)))

	DY(2) = 
	&- K*(3.0d0/((Y(2)-Y(1))**ro) - 3.0d0/((Y(3)-Y(2))**ro) -
     &one/((Y(1)-Y(M)+M)**ro) + one/((Y(4)-Y(3))**ro)
     &+U*(3.0d0/((Y(2)-Y(1))**en) - 3.0d0/((Y(3)-Y(2))**en) -
     &one/((Y(1)-Y(M)+M)**en) + one/((Y(4)-Y(3))**en)))

      
	DO i = 3, M-2
	DY(i) = 
	&-K*(3.0d0/((Y(i)-Y(i-1))**ro) - 3.0d0/((Y(i+1)-Y(i))**ro) -
     &one/((Y(i-1)-Y(i-2))**ro) + one/((Y(i+2)-Y(i+1))**ro))
	&
     &+U*(3.0d0/((Y(i)-Y(i-1))**en) - 3.0d0/((Y(i+1)-Y(i))**en) -
     &one/((Y(i-1)-Y(i-2))**en) + one/((Y(i+2)-Y(i+1))**en))

	ENDDO

	
      DY(M-1) = 
	&-K*(3.0d0/((Y(M-1)-Y(M-2))**ro) - 3.0d0/((Y(M)-Y(M-1))**ro) -
     &one/((Y(M-2)-Y(M-3))**ro) + one/((Y(1)-Y(M)+M))**ro)
     &+U*(3.0d0/((Y(M-1)-Y(M-2))**en) - 3.0d0/((Y(M)-Y(M-1))**en) 
     &-one/((Y(M-2)-Y(M-3))**en) + one/((Y(1)-Y(M)+M)**en))

	
      DY(M) = 
	& - K*(3.0d0/((Y(M)-Y(M-1))**ro) - 3.0d0/((Y(1)-Y(M)+M)**ro) -
     &one/((Y(M-1)-Y(M-2))**ro) + one/((Y(2)-Y(1))**ro)) 
     &+U*(3.0d0/((Y(M)-Y(M-1))**en) - 3.0d0/((Y(1)-Y(M)+M)**en) -
     &one/((Y(M-1)-Y(M-2))**en) + one/((Y(2)-Y(1))**en)) 

	return
	end


      SUBROUTINE g_mm0(Y,DY,M)
c     13.09.02
c corresponds to MMI as introduced by VT in 2002 but slightly changed
c      to meet the analogy with the model of Popkov and Krug thus
c      it still has 2 parameters - b and U, like in the PK model	 plus
c      the two powers r and n, note the difference - here
c      in the stabilization part the terrace widths are raised to (n+1)
c     (c) VT, May 2008, Lexington KY
      INTEGER          M,I
      REAL*8           Y(*),DY(*),d(M+1),d3(M+1),dtem
	real*8 on
	data on / 1.0d0 /
	  	real*8 U,en,b, bp1, bm1
	real*8 par(5)
	common /pxx/ par
c     first we study the case ro = 1.0 as in the PK model
c	ro=par(4)
	b=par(1)
	bp1 = (b+ 1.0d0)/2.0d0
	bm1 = (1.0d0 - b)/2.0d0													  
	U=par(2)
	en=par(5)+1.0d0
c calculate inverse of the distances and its third power:
      do i=2,M
	 dtem=Y(i)-Y(i-1)
c	 d(i)=dtem**ro	! so instead of this line we have:
	 d(i)=dtem
	 d3(i)=on/dtem**en
	enddo
	dtem=Y(1) - Y(M) + M
c	d(1)=dtem**ro  ! here also:
	d(1)=dtem
	d(M+1)=d(1)
	d3(1) = on/dtem**en
	d3(M+1)= d3(1)
c
      DO I=1,M
	DY(I)=bp1*d(i)+bm1*d(i+1)-U*(d3(i+1)-d3(i))
       END DO
      RETURN
      END
*************************************************************
      SUBROUTINE g_mm1(Y,DY,M)
c     13.09.02
c corresponds to MMI
      REAL*8           Y(*),DY(*),d(M+1),d3(M+1),dtem
      INTEGER          M,I
	real*8 on
	data on / 1.0d0 /
	  	real*8 ro,a3,en,beta
	real*8 par(5)
	common /pxx/ par
	ro=par(1)
	beta=par(2)
	a3=par(3)
	en=par(5)
c calculate inverse of the distances and its third power:
      do i=2,M
	 dtem=Y(i)-Y(i-1)
	 d(i)=dtem**ro
	 d3(i)=on/dtem**en
	enddo
	dtem=Y(1) - Y(M) + M
	d(1)=dtem**ro
	d(M+1)=d(1)
	d3(1) = on/dtem**en
	d3(M+1)= d3(1)
c
      DO I=1,M
	DY(I)=d(i)+beta*d(i+1)-a3*(d3(i+1)-d3(i))
       END DO
      RETURN
      END
*************************************************************
      SUBROUTINE gise2(Y,DY,M)
	INTEGER M,I
      REAL*8          Y(*),DY(*),d(0:M+1),d2(0:M+1),d3(0:M+2),Ybc(0:M+1)      
      real*8 on,tw , three
	data on,tw /1.d0, 2.d0 /
	real*8 gama,l0l,dpl,dml,a3,gama3
	real*8 tp1,tp2,tp3,tp4,tp5,tp6
	real*8 di1,di,dn1,dn
      real*8 ft,st,dyi
	real*8 par(5)
	common /pxx/ par
	gama=par(1)
	l0l=par(2)
	dpl=par(3)
	dml=par(4)
	three = par(5) + 1.0d0
	a3=l0l**three
	gama3=a3*gama
c calculate distances:
      do i=1,M
	 Ybc(i)=Y(i)
	enddo
	Ybc(0)=y(M)-M
	Ybc(M+1)=Y(1)+M
c 
      do i=1,M+1
	dyi=Ybc(i)-Ybc(i-1)
	d(i)=dyi
	d2(i)=dyi*dyi/tw
	d3(i)=on/dyi**three
	enddo
c
	d(0) = d(M)
	d2(0)=d2(M)
	d3(0)=d3(M)
	d3(M+2)=d3(2)
c
       DO I=1,M
	  di=d(I)
	  di1=d(I+1)
	  dn=on/(dpl+dml+di)
	  dn1=on/(dpl+dml+di1)
        tp1=d3(I+2)
	  tp4=d3(I+1)
	  tp2=tw*tp4
	  tp3=d3(I)
	  tp5=tw*tp3
	  tp6=d3(I-1)
	  ft=(d2(I+1)+dml*di1)*dn1
	  ft=ft+(d2(I)+dpl*di)*dn
	  st=-(tp2-tp1-tp3)*dn1+(tp5-tp4-tp6)*dn
       DY(I)=ft+gama3*st
       END DO
      RETURN
      END

      SUBROUTINE DE13R (F,M,H,YN,Y)
      !	(F,M,YN,EPS,P,H,Y,YP,DELTY,YR,DY,IERR)
c                DE13R(F,M,H,YN,Y)
      IMPLICIT        NONE
      EXTERNAL        F	!SUBROUTINE FOR CALCULATION OF THE RIGHT-HAND
	INTEGER         M   !NUMBER OF THE EQUATIONS IN THE SYSTEM
c
c                               !SIDE OF THE SYSTEM OF DIFFERENTIAL EQUATIONS
c                               !SUBROUTINE F(X,Y,DY,M)
c                               !X  - INDEPENDENT VARIABLE
c                               !Y  - VECTOR WITH INDEPENDENT VARIABLES
c                               !DY - VECTOR WITH THE RIGHT-HAND SIDES
c                               !M  _ NUMBER OF THE EQUATIONS
      INTEGER         M        !NUMBER OF THE EQUATIONS IN THE SYSTEM
      REAL*8                 !  XN  we dont need an INITIAL VALUE OF THE INDEPENDENT VARIABLE
     1                YN(*),   !INITIAL CONDITIONS, i.e. the coordinates before the integration step
     2                !  XK,      !END OF THE INTEGRATION INTERVAL, IT IS VALID
c                               !XK.GT.XN.OR.XK.EQ.XN.OR.XK.LT.XN
c                               !OF CONTROL FROM ABSOLUTE TO RELATIVE ERROR
c                               !IF ABS(Y).GE.P EPS IS THE ABSOLUTE ERROR
     6                H,       !VALUE OF THE STEP. POSITIVE, WHEN
c                               !XK.GT.XN AND NEGATIVE, WHEN XK.LT.XN
     7                Y(*)     !THE SOLUTION AT VALUE OF THE ARGUMENT XK
      REAL*8          YP(*),   !WORKING VECTORS
     1                DELTY(*),
     2                YR(*),
     3                DY(*)
      INTEGER         IERR     !FLAG FOR ERRORS
                               !IERR=0  - CORRECT SOLUTION.
                               !IERR=65 - SOME OF THE COMPONENTS OF THE
                               !SOLUTION CANNOT BE CALCULATED WITH ACCURACY
                               !EQUAL TO EPS. THE SOLUTION MUST BE REPEATED
                               !WITH NEW VALUES OF THE PARAMETERS HMIN AND H.

      REAL*8          dabs,
     1                B,
     2                BC,
     3                BD(4),
     4                BE(3),
     5                BF(3),
     6                C,
     7                dsign,
     8                U,
     9                V,
     A                Z
      LOGICAL         A,BUL,W
      INTEGER         BB
	integer   IE
      DATA            BD /0.166666666666666667D0,
     1                    0.333333333333333333D0,
     2                    0.333333333333333333D0,
     3                    0.166666666666666667D0/
      DATA            BE /0.5D0,
     1                    0.5D0,
     2                    1.0D0/
      DATA            BF /-1.0D0,
     1                    -1.0D0,
     2                     1.0D0/
	HMIN =    !MINIMAL PERMISSIBLE ABSOLUTE VALUE OF THE
c                               !INTEGRATION STEP
      EPS =     !PERMISSIBLE ERROR FOR CALCULATION OF ALL
c                               !COMPONENTS OF THE SOLUTION
      P =        !BOUNDARY VALUE IN THE SOLUTION FOR TRANSITION

      IERR=0
      B=32.0D0
      C=XN
      DO IE=1,M
         Y(IE)=YN(IE)
      END DO
      IF (XN.EQ.XK) GO TO 150
      A=.FALSE.
      BUL=.FALSE.
c      H=dsign(H,XK-XN)

   20 CONTINUE
      IF (.NOT.A) GO TO 30
      H=U
      GO TO 150

   30 CONTINUE
      V=XK-C
      IF (dabs(V).GT.dabs(H)) GO TO 40
      U=H
      A=.TRUE.
      H=V

   40 CONTINUE
      BC=C
      DO E=1,M
         YP(E)=Y(E)
         YR(E)=Y(E)
      END DO
      CALL F(YR(1),DY(1),M)
      DO E=1,M
         DY(E)=H*DY(E)
         Y(E)=Y(E)+BD(1)*DY(E)
         DELTY(E)=DY(E)
      END DO
      DO BB=1,3
         C=BC+BE(BB)*H
         DO E=1,M
            YR(E)=YP(E)+BE(BB)*DY(E)
         END DO
         CALL F(YR(1),DY(1),M)
         DO E=1,M
            DY(E)=H*DY(E)
            Y(E)=Y(E)+BD(BB+1)*DY(E)
            DELTY(E)=DELTY(E)+BF(BB)*DY(E)
         END DO
      END DO
      C=BC+H
      W=.TRUE.
      DO E=1,M
         Z=dabs(Y(E))
         IF (Z.LT.P) GO TO 100
         DELTY(E)=DELTY(E)/Y(E)
  100    CONTINUE
      END DO
      DO E=1,M
         Z=dabs(DELTY(E))
         IF (Z.GT.EPS) GO TO 120
            IF (.NOT.W) GO TO 110
               IF (B*Z.LT.EPS) GO TO 110
                  W=.FALSE.
  110    CONTINUE
      END DO
      IF (W) H=H+H
      GO TO 20

  120 CONTINUE
      DO E=1,M
         Y(E)=YP(E)
      END DO
      IF (BUL) GO TO 140
         C=C-H
         H=0.5D0*H
         A=.FALSE.
         IF (dabs(H).GE.HMIN) GO TO 40
            H=dsign(HMIN,H)
            BUL=.TRUE.
            GO TO 40

  140 CONTINUE
      IERR=65


  150 RETURN
      END
c -----------------------------------------------------------
	 subroutine MC_An (M, Y)
	INTEGER          M,I 
      REAL*8           Y(M), YT(M), Dt ! ,DY(*),dk(M+1),du(M+1),dtem
	real*8 bdef
	real r
c 	Bunch/ Terrace Definitions:
	common /bdxx/ bdef
	 do i = 1, M

	 enddo
	do i = 2, M
	Dt = Y(i) - Y(i-1)
	r = rand(0)
	if (Dt.lt.bdef) then
	    Dt=Dt*r
	else
	    Dt = Dt/r
	endif
	Y(i) = Y(i-1) + Dt
	enddo

	return
	end
c