C
      SUBROUTINE USDFLD(FIELD, STATEV, PNEWDT, DIRECT, T, CELENT, TIME, DTIME,
     1 CMNAME, ORNAME, NFIELD, NSTATV, NOEL, NPT, LAYER, KSPT, KSTEP, KINC,
     2 NDI, NSHR, COORD, JMAC, JMTYP, MATLAYO, LACCFLA)
c
       INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD), STATEV(NSTATV), DIRECT(3,3), T(3,3), TIME(2)
      DIMENSION ARRAY(15), JARRAY(15), COORD(*), JMAC(*), JMTYP(*)
cc
c Get SP3 from previous increment
      CALL GETVRM('SP', ARRAY, JARRAY, FLGRAY, JCRD, JMAC, JMTYP,
     3 MATLAYO, LACCFLA)
      FIELD(1) = ARRAY(1)
      FIELD(2) = ARRAY(2)
      FIELD(3) = ARRAY(3)
c   
      STATEV(1) = FIELD(1)
      STATEV(2) = FIELD(2)
      STATEV(3) = FIELD(3)
C 
      RETURN
      END
C ***************************************************************C      
C   This uses traditional units for wrc541.
C
      subroutine creep(decra, deswa, statev, serd, ec, esw, p, qtild,
     1  temp, dtemp, predef, dpred, time, dtime, CMNAME, leximp, lend,
     2  coords, nstatv, noel, npt, layer, kspt, kstep, kinc)
c
        include 'ABA_PARAM.INC'
c
        character*80 CMNAME  
c
        dimension decra(5), deswa(5), statev(nstatv), predef(*),
     1  dpred(*), time(2), coords(*),ec(2),a(4), b(4), c(4), d(4),
     2  aa(10), ab(10), ac(10), ecdot(10)  
c     
c things to check:
c     unitf
c     zhours
c     tfact
c     definition of temp
c     decra(1) and decra(5)
ccc      uses wsrf in damage not strain
c  need to define creep, depvar, electrical conductivity, USDFLD 
c        
c    
      sigr = 0.0  
      if (statev(1).gt.sigr) then 
        sigr = statev(1)
        end if
      if (statev(2).gt.sigr) then 
        sigr = statev(2)
        end if
      if (statev(3).gt.sigr) then 
        sigr = statev(3)
        end if
c      if (sigr.lt.0.10) then sigr = 0.10
c      if(sigr.lt.0.0) then 
c      sigr = 0.0010
c      endif
c'      if(qtild.gt.sigr) then
c      sigr = qtild 
c      endif
C Set units factor, 1 for ksi, 1000 for psi, call matprop
      unitf = 1000.0 
      unitf = 1.0 
c
      call matprop(CMNAME,a,b,c,d,emg,wsrf,ntype)
      wsrf = 1.0
c
c	define use of wsrf
c
	sigr = sigr/unitf/wsrf
	qtild = qtild/unitf
 	statev(4) = sigr
c        sigr = statev(4)
c  
c  statev: 1 - 3 MPS
c          4 sigr            
c          5 min  creep damage
c          6 mean creep damage  
c          7 probability - sigr
c          8 
c          9 omega creep damage defined by deltasr, deltacd
c
c **************************** 
c     Default units, hours, F, ksi
C *********************************************
	abszero = 459.76d0
	tempr = temp + abszero
c  hours per year for probability per year	
	zhours = +24*365d0
	zhours = 8000.0d0
c
c  Factor tfact =3600.0 is for use with coupled thermal analysis in seconds
c  Otherwise tfact = 1.0
c
        tfact = 3600.d0
        tfact = 1.0d0
cc  *******************************************
      zpi = 3.1415926535897932384d0
      ze = 2.718281828d0
c******************************************************************
c  Creep analysis using sinh creep rupture model. comment out "IF" options or change names to use min creep strain rate data
c***************************************************************
c 
	acon = a(1)
	qr = a(2)
	expn = a(3)
	sd = a(4)
C***********************************
c        ntype = 1
C*******************************
c	if (ntype.eq.1) then
c	tempfact = 5.0d0/9.0d0
c      tempfact = 1.0
c     temp in F for xlmlife xlmlife now in imperial units      
c	emg = shcr(tempr,qtild,acon,qr,expn,sd)
c	emg = emg*xlmlife(CMNAME,temp,qtild,-2.1)
c	end if
	
c	emg = shcr(tempr,qtild,acon,qr,expn,sd)/tfact
c	BB1 = xmeanlm(CMNAME, qtild)
c	BB1 = BB1*1000.0/(tempr) - almconst
c	emg = emg*10.0d0**BB1
c	if (ntype.eq.3) then
c	emg = shcr(tempr,qtild,acon,qr,expn,sd)/tfact
c	emg = emg/shcr(tempr,qtild,c(1),b(2),b(3),b(4))
c	end if
c      
	decra(1) = dtime*shcr(tempr,qtild,acon,qr,expn,sd)
c	decra(1) = decra(1)*dexp(ec(1)/emg) 
	if(LEXIMP.eq.1) then
	decra(5) = dtime*shdedq(tempr,qtild,acon,qr,expn,sd)
c	decra(5) = decra(5)*dexp(ec(1)/emg)
	end if	
C ****************
c******************************************************************
c  Damage using wrc 541 LM data
c*************************************************************** 
*******          wsrf for HBB-I-14.10C-1       
c        if (sigr.lt.1.0) then sigr = 1.0
        zlife1 = xlmlife(CMNAME,temp,sigr,1.2d0)
        zlife2 = xlmlife(CMNAME,temp,sigr,0.0d0)
        zlife3 = xlmlife(CMNAME,temp,sigr,1.0d0)
        zlife4 = xlmlife(CMNAME,temp,sigr,-1.650)
c
        statev(9) = zlife1 
c        statev(8) = decra(5)
c        statev(9) = decra(1)
        statev(10) = temp
c        end if
c        statev(12) = zlife3
c******************************************************************       
c******************************************************************
c  Damage and probability
c************************************************************
c 
c        if (zlife1.gt.0.001) then
        statev(5) = statev(5) + dtime/tfact/zlife1 
c        end if
c        if (zlife2.gt.0.001) then
        statev(6) = statev(6) + dtime/tfact/zlife2
c        end if
c        end if
c    HBB creep damage
c        statev(8) = statev(8) + dtime/tfact/zlife4
cc
       fact = zlife3/zlife2
       stdev = dlog(statev(6))/dlog(fact)
       tmod = statev(6)*zlife2
       fact1 = bphi(stdev)
c     
c prob per year
       fact2 = dexp(-0.5d0*stdev*stdev)/((2.0d0*zpi)**0.5d0)
       statev(7) = fact2*stdev/(dlog(statev(6))/dlog(ze))
       statev(8) = statev(7)/tmod*zhours/(1.0d0 - fact1)
c       statev(9) = statev(7) 
c       statev(9) = xmean(CMNAME,sigr)
c
c        statev(6) = xlmpa(CMNAME,0)
c        statev(7) = xmean(CMNAME,sigr)   
90    RETURN  
      END
c  ******************************************************************

c ******** abd and ac removed from xmean and xmin
      real*8 function wsrf1(CMNAME,temp)
      implicit real*8(a-h,o-z)
      character*80 CMNAME
      wsrf1 = 1.0d0
90    Return        
      end
*********************************************************************************************   
c *****************************************************************************
      real*8 function wsrf(CMNAME,temp)
      implicit real*8(a-h,o-z)
c     sourcw= Table HBB-I-14.10C-1200 - 1500 F
      character*80 CMNAME
      wsrf = 1.0d0
      IF (CMNAME.EQ.'A800H_U') Then
      wsrf = -0.00055d0*temp+1.360d0
      End if
90    Return               
      End
********************************************************************** 
      real*8 function xmean(CMNAME,stress)
      implicit real*8(a-h,o-z)
c       dimension ab(10), ac(10)
      character*80 CMNAME
c     WRC 541 data
      IF (CMNAME.EQ.'LCS_U') Then
      a0 = 35093.24
      a1 = -3603.7901
      a2 = -1913.659
      a3 = -250
      End if
      IF (CMNAME.EQ.'GR22_U') Then
      a0 = 43946.4
      a1 = -8390.0
      a2 = 0.0
      a3 = 0.0  
      End if  
      IF (CMNAME.EQ.'GR91_U1') Then
      a0 = 63450.000
      a1 = -1380.000
      a2 = -5139.5230
      a3 = 0.0
      End If
      IF (CMNAME.EQ.'GR91_U2') Then
      a0 = 63450.0
      a1 = -1380.0
      a2 = -5139.9
      a3 = 0.0
      End if
      IF (CMNAME.EQ.'GR91_UW') Then
      a0 = 63450.0
      a1 = -1380.0
      a2 = -5139.523
      a3 = 0.0
      End if 
      IF (CMNAME.EQ.'GR9_UB') Then
      a0 = 63450.0
      a1 = -1380.0
      a2 = -5139.523
      a3 = 0.0
      End if
      IF (CMNAME.EQ.'GR9A_U') Then
      a0 = 45062.15d0
      a1 = -5600.73d0
      a2 = -1649.03d0
      a3 = -224.331d0
      End if
      IF (CMNAME.EQ.'GR9B_U') Then
      a0 = 45062.15d0
      a1 = -5600.73d0
      a2 = -1649.03d0
      a3 = -224.331d0
      End if
      IF (CMNAME.EQ.'GR9C_U') Then
      a0 = 45062.15d0
      a1 = -5600.73d0
      a2 = -1649.03d0
      a3 = -224.331d0
      End if
      IF (CMNAME.EQ.'GR9D_U') Then
      a0 = 45062.15d0
      a1 = -5600.73d0
      a2 = -1649.03d0
      a3 = -224.331d0
      End if
      IF (CMNAME.EQ.'GR9E_U') Then
      a0 = 45062.15d0
      a1 = -5600.73d0
      a2 = -1649.03d0
      a3 = -224.331d0
      End if
      IF (CMNAME.EQ.'GR9F_U') Then
      a0 = 45062.15d0
      a1 = -5600.73d0
      a2 = -1649.03d0
      a3 = -224.331d0
      End if
      IF (CMNAME.EQ.'TP304H_U') Then
      a0 = 4.353946d4
      a1 = -9731.8d0
      a2 = 0.0
      a3 = 0.0
      End if
      xlogs = log(stress)/log(10.0d0)
      xmean = a0+a1*xlogs+a2*xlogs*xlogs+a3*xlogs*xlogs*xlogs
90    Return
      End
c ******************************************************************************
      real*8 function xlmlife(CMNAME,temp,stress,sdev)
      implicit real*8(a-h,o-z)
      character*80 CMNAME
      tr = (temp + 459.76)
      xlmlife = xmean(CMNAME,stress)/tr - xlmpa(CMNAME,sdev)
      xlmlife = 10.0d0**xlmlife
90    RETURN 
      End
c *******************************************************************************
      real*8 function xlmpa(CMNAME,sdev)
      implicit real*8(a-h,o-z)
      character*80 CMNAME
c   ccccccccc   End If
      IF (CMNAME.EQ.'LCS_U') Then
      xlmp2 = 18.15
      xlmp1 = 17.7 
      End if
      IF (CMNAME.EQ.'GR22_U') Then
      xlmp2 = 19.565607
      xlmp1 = 18.9181
      End if  
      IF (CMNAME.EQ.'GR91_U1') Then
      xlmp2 = 30.886006
      xlmp1 = 30.36423
      End If
      IF (CMNAME.EQ.'GR91_U2') Then
      xlmp2 = 30.886006
      xlmp1 = 30.36423
      End If 
      IF (CMNAME.EQ.'GR91_UW') Then
      xlmp2 = 30.886006
      xlmp1 = 30.36423
      End if
      IF (CMNAME.EQ.'GR9_UB') Then
      xlmp2 = 30.886006
      xlmp1 = 30.36423
      End if
      IF (CMNAME.EQ.'GR9A_U') Then
      xlmp2 = 20.946d0
      xlmp1 = 20.5d0
      End if
      IF (CMNAME.EQ.'GR9B_U') Then
      xlmp2 = 20.946d0
      xlmp1 = 20.5d0
      End if
      IF (CMNAME.EQ.'GR9C_U') Then
      xlmp2 = 20.946d0
      xlmp1 = 20.5d0
      End if
      IF (CMNAME.EQ.'GR9D_U') Then
      xlmp2 = 20.946d0
      xlmp1 = 20.5d0
      End if
      IF (CMNAME.EQ.'GR9E_U') Then
      xlmp2 = 20.946d0
      xlmp1 = 20.5d0
      End if
      IF (CMNAME.EQ.'GR9F_U') Then
      xlmp2 = 20.946d0
      xlmp1 = 20.5d0
      End if
      IF (CMNAME.EQ.'TP304H_U') Then
      xlmp2 = 16.145903d0
      xlmp1 = 15.52195d0
      End if
      xlmpa = xlmp1 + (xlmp1 - xlmp2)/1.645*sdev
90    RETURN 
      End
c *******************************************************************************
c
      real*8 function acosh(x)
      implicit real*8(a-h,o-z)
c     
      t1 = dexp(x)
      t2 = dexp(-x)
      acosh = (0.5d0*(t1+t2))
 90   return
      end
      
c ******************************************
      real*8 function bphi(x)
      implicit real*8(a-h,o-z)
c ref wikipedia (Bell)      
      zpi = 3.1412d0
      val = (1.0d0-dexp(-2.0d0/zpi*x*x))**0.5d0
      if (x.eq.0.0d0) then 
      sgn = 0.0d0
      END IF  
      if (x.lt.0.0d0) then 
      sgn = -1.0d0
      END IF
      if (x.gt.0.0d0) then
      sgn = 1.0d0
      END IF
      bphi = (1.0d0 + sgn*val)/2.0d0
90    return  
      end
******************************************************************************
c
      real*8 function shcr(tempr,q,ab,qr,expn,sd)
      implicit real*8(a-h,o-z)
c      dimension  ab(3), qr(3), expn(3), saa(3)
	dot = ab*dexp(-qr/tempr)
	t1 = dexp(q/sd)
      t2 = dexp(-q/sd)
      shcr = dot*(0.5d0*(t1-t2))**expn
 90   return
      end 
c *******************************************************************************
c
      real*8 function shdedq(tempr,q,ab,qr,expn,sd)
      implicit real*8(a-h,o-z)
c      dimension  ab(3), qr(3), expn(3), saa(3)
	dot = ab*dexp(-qr/tempr)
	t1 = dexp(q/sd)
      t2 = dexp(-q/sd)
      shdedq = expn*dot*(0.5d0*(t1-t2))**(expn-1.0d0)
      shdedq = shdedq/sd*0.5d0*(t1+t2)
 90   return   
      end       
c ********************************************************************************* 
	subroutine matprop(CMNAME,a,b,c,d,emg,wsrf,ntype)
	include 'ABA_PARAM.INC'
	character*80 CMNAME 
	dimension a(4),b(4),c(4),d(4)  
	emg = 1  
	wsrf = 1.0 
	ntype = 1
c	sinh strain rate, wrc 541 rupture: check GR5W
c      IF (CMNAME.EQ.'LCSCAP_U') GO TO 202
c      IF (CMNAME.EQ.'LCSHEAD_U') GO TO 202
c      IF (CMNAME.EQ.'LCSFL_U') GO TO 202
       IF (CMNAME.EQ.'LCS_U') GO TO 200
c      IF (CMNAME.EQ.'MCSW_U') GO TO 201
       IF (CMNAME.EQ.'GR91_U1') GO TO 100
c       IF (CMNAME.EQ.'GR91_U2') GO TO 100  
c       IF (CMNAME.EQ.'GR91_UW') GO TO 100
c       IF (CMNAME.EQ.'GR9_UB') Go TO 100
       IF (CMNAME.EQ.'GR22_U') GO TO 120
c      IF (CMNAME.EQ.'GR22_U16') GO TO 120
c      IF (CMNAME.EQ.'GR22_U18') GO TO 120
c      IF (CMNAME.EQ.'GR22_U4') GO TO 120
c      IF (CMNAME.EQ.'GR22_U20') GO TO 120
c      IF (CMNAME.EQ.'GR22_U26') GO TO 120
c      IF (CMNAME.EQ.'parent') GO TO 120
c      IF (CMNAME.EQ.'GR22W_U') GO TO 121
c      IF (CMNAME.EQ.'GR5_U') GO TO 100
c      IF (CMNAME.EQ.'GR5W_U') GO TO 101
c       IF (CMNAME.EQ.'GR9A_U') GO TO 109
c       IF (CMNAME.EQ.'GR9B_U') GO TO 109
c       IF (CMNAME.EQ.'GR9C_U') GO TO 109
c       IF (CMNAME.EQ.'GR9D_U') GO TO 109
c       IF (CMNAME.EQ.'GR9E_U') GO TO 109
c       IF (CMNAME.EQ.'GR9F_U') GO TO 109
c       IF (CMNAME.EQ.'GR9_UW') GO TO 109
c      IF (CMNAME.EQ.'GR94_U') GO TO 109
c      IF (CMNAME.EQ.'GR98_U') GO TO 109
c      IF (CMNAME.EQ.'GR912_U') GO TO 109
c      IF (CMNAME.EQ.'GR9W_U') GO TO 300
c      IF (CMNAME.EQ.'GR91_U') GO TO 102
c      IF (CMNAME.EQ.'GR91W_U') GO TO 103
c      IF (CMNAME.EQ.'SA316L_U') GO TO 104
c      IF (CMNAME.EQ.'SA316H_U') GO TO 105
c      IF (CMNAME.EQ.'SA316HW_U') GO TO 106
      IF (CMNAME.EQ.'TP304H_U') GO TO 107
      IF (CMNAME.EQ.'TP304_U') GO TO 107
c      IF (CMNAME.EQ.'TP304SH_U') GO TO 107
c      IF (CMNAME.EQ.'TP304W_U') GO TO 207
c      IF (CMNAME.EQ.'TP304SHW_U') GO TO 207
c      IF (CMNAME.EQ.'TP347H_U') GO TO 108
c      IF (CMNAME.EQ.'A1_U')GO TO 107
c      IF (CMNAME.EQ.'A2_U')GO TO 107
c      IF (CMNAME.EQ.'A3_U')GO TO 107
c      IF (CMNAME.EQ.'A4_U') GO TO 107
c      IF (CMNAME.EQ.'A3W_U')GO TO 207
c      IF (CMNAME.EQ.'CMV_U') GO TO 130
c      IF (CMNAME.EQ.'A800H_U') GO TO 130
c      
c      GO TO 50
c    tp304h 1100 - 1700 F
107    ax = 21.8638d0
       bx = 4.96396d0
       a(1) = 10.0d0**ax
       a(2) = 10.0d0**bx
       a(3) = 7.624d0
       a(4) = 32.311d0
       emg = 0.0
       wsrf = 1.00d0
       GO TO 50 
c      grade 9
109    ax = 26.2939d0
       bx = 31.0244d0
       a(1) = 10.0d0**ax
       a(2) = 10.0d0**bx
       a(3) = 6.341d0
       a(4) = 31.024d0
       wsrf = 1.0d0  
       IF (CMNAME.EQ.'GR9_UW') wsrf = 0.50
c      grade 9 weld
110     ax = 26.2939d0
       bx = 31.0244d0
       a(1) = 10.0d0**ax
       a(2) = 10.0d0**bx
       a(3) = 6.341d0
       a(4) = 31.024d0
       wsrf = 0.9d0
         
130   aconst = 25.1982d0   
      a(1) = 10.0d0**aconst
      bconst = 5.01601d0
      a(2) = 10.0d0**bconst
      a(3) = 7.13d0
      a(4) = 45.685d0
      emg = 0.006d0
      wsrf = 1.0
      GO TO 50 
c      
100   ax = 33.85175d0
      bx = 5.176690d0
       a(1) = 10.0d0**ax
       a(2) = 10.0d0**bx
       a(3) = 6.270d0
       a(4) = 10.251d0
c       a(4) = 6.290
       wsrf = 1.0d0 
c       IF (CMNAME.EQ.'GR91_UW') wsrf = 0.91d0
c       IF (CMNAME.EQ.'GR9_UB') wsrf = 0.91d0
       GO TO 50  
c 
120   ax = 29.320d0
c     min dustile 
c     ax = 29.200d0
c     mean brittle
      ax = 28.497d0   
      bx = 5.0189
      a(3) = 6.851
      a(4) = 43.363
      a(1) = 10.0d0**ax
      a(2) = 10.0d0**bx
      wsrf = 1.0 
      GO TO 50 
c
c    AWG 2701515 CS sdev = 0.369 
200   ax = 19.584d0
      bx = 4.834d0
      a(1) = 10.0d0**ax
      a(2) = 10.0d0**bx
      a(3) = 5.926d0
      a(4) = 36.28d0
      wsrf = 1.0d0
      GO TO 50
c       
50    return
      end
*****************************************************************
