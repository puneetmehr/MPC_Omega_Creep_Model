c*******************************************************************
c Subroutine written using maincode.py on Monday, 08 November 2021 *
c Name - omega_creep_rev0										   *
c Revision - 0 (written by PM on 04/18/2020) 					   *
c*******************************************************************
c Abaqus subroutine to calculate creep strain using Sinh creep     *
c model and damage due to creep based on API 579 Part 10 using MPC *
c Omega and LMP method based on API 530 and WRC 541.			   *
c*******************************************************************
c Originally written by Peter Carter							   *
c Re-written by Puneet Mehra									   *
c*******************************************************************
c 					Information                                    *
c  statev: 1 - 3 Principal Stresses (PS1<PS2<PS3)                  *
c          4 sige (effective stress)                               *
c          5 min  creep damage (WRC541)                            *
c          6 mean creep damage (WRC541)                            *
c          7 min creep damage (API 530)                            *
c          8 avg creep damage (API 530)                            *
c          9 MPC Omega creep damage                                *
c things to check:												   *
c     unitf                                                        *
c     zhours                                                       *
c     tfact                                                        *
c     definition of temp                                           *
c     decra(1) and decra(5)                                        *
c  need to define creep, depvar, electrical conductivity, USDFLD   *
c 		 														   *
c*******************************************************************
      SUBROUTINE USDFLD(FIELD, STATEV, PNEWDT, DIRECT, T, CELENT, TIME,      
     1 DTIME, CMNAME, ORNAME, NFIELD, NSTATV, NOEL, NPT, LAYER, KSPT, KSTEP, 
     2 KINC, NDI, NSHR, COORD, JMAC, JMTYP, MATLAYO, LACCFLA)                
c                                                                            
      INCLUDE 'ABA_PARAM.INC'                                                
c                                                                            
      CHARACTER*80 CMNAME,ORNAME                                             
      CHARACTER*3  FLGRAY(15)                                                
      DIMENSION FIELD(NFIELD), STATEV(NSTATV), DIRECT(3,3), T(3,3), TIME(2)  
      DIMENSION ARRAY(15), JARRAY(15), COORD(*), JMAC(*), JMTYP(*)           
c********************************************************************        
c Get Principal stresses from previous increment					*         
c********************************************************************        
      CALL GETVRM('SP', ARRAY, JARRAY, FLGRAY, JRCD, JMAC, JMTYP,            
     1 MATLAYO, LACCFLA)                                                     
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
c*********************************************************************
c Creep using the Sinh creep model                                   *
c*********************************************************************
      subroutine creep(decra, deswa, statev, serd, ec, esw, p, qtild, 
     1  temp, dtemp, predef, dpred, time, dtime, CMNAME, leximp, lend,
     2  coords, nstatv, noel, npt, layer, kspt, kstep, kinc)          
c                                                                     
        include 'ABA_PARAM.INC'                                       
c                                                                     
        character*80 CMNAME                                           
c                                                                     
        dimension decra(5), deswa(5), statev(nstatv), predef(*),      
     1  dpred(*), time(2), coords(*),ec(2),a(5), b(5), c(2), d(2)     
c                                                                     
		real*8 acon, sig0, expn, wsrf, t1, t2, shcr, shdedq   
		real*8 zlife1,zlife2,zlife3,zlife4,zlife5   
c*********************************************************************
c   Call matprop for Sinh model parameters based on MPC Omega data   *
c*********************************************************************
C Set units factor, 1 for ksi, 1000 for psi                                     
c     unitf = 1000.0d0                                                         
      unitf = 1.0d0                                                             
      call matprop(CMNAME,acon,sig0,expn,wsrf)                                  
	  wsrf = wsrf                                                                 
	  acon = acon                                                                  
	  expn = expn                                                                  
	  sig0 = sig0                                                                  
c                                                                               
c      wsrf = 1.0d0                                                              
c                                                                               
c	define use of wsrf                                                           
c                                                                               
	  sige=qtild                                                                 
	  sige = sige/unitf/wsrf                                                     
	  qtild = qtild/unitf                                                        
	  statev(4)=sige                                                             
c                                                                               
C *********************************************                                 
c     Default units, hours, F, ksi                                              
C *********************************************                                 
c                                                                               
c  Factor tfact =3600.0 is for use with coupled thermal analysis in seconds     
c  Otherwise tfact = 1.0                                                        
c                                                                               
      tfact = 3600.d0                                                         
      tfact = 1.0d0                                                           
c****************************************************************** 
	  ze = 2.718281828d0                                             
	  tempr=temp+460.0d0                                             
c****************************************************************** 
c  Creep analysis using Sinh creep model                           
c****************************************************************** 
c                                                                   
      t1 = dexp(qtild/sig0)
      t2 = dexp(-qtild/sig0)
      shcr = acon*(0.5d0*(t1-t2))**expn
      shdedq = expn*acon*(0.5d0*(t1-t2))**(expn-1.0d0)
      shdedq = shdedq/(sig0*0.5d0*(t1+t2))
	  decra(1) = dtime*shcr
	  if(LEXIMP.eq.1) then 
	  decra(5) = dtime*shdedq 
	  end if                                                                    
c******************************************************************            
c  Creep damage calculations 									  *             
c******************************************************************            
c min properties WRC 541                                                       
      zlife1 = xlmlife(CMNAME,tempr,sige,0.0d0)                              
c avg properties WRC 541                                                       
      zlife2 = xlmlife(CMNAME,tempr,sige,1.0d0)                              
c min properties API 530                                                       
      zlife3 = xlmlife1(CMNAME,tempr,sige,0.0d0)                             
c avg properties API 530                                                       
      zlife4 = xlmlife1(CMNAME,tempr,sige,1.0d0)                             
c MPC Omega                                                                    
      zlife5 = omegalife(CMNAME,tempr,sige,statev(1),statev(2),statev(3))
c min creep damage using WRC 541                                               
      statev(5) = statev(5) + dtime/tfact/zlife1                             
c avg creep damage using WRC 541                                               
      statev(6) = statev(6) + dtime/tfact/zlife2                             
c min creep damage using API 530                                               
      statev(7) = statev(7) + dtime/tfact/zlife3                             
c avg creep damage using API 530                                               
      statev(8) = statev(8) + dtime/tfact/zlife4                             
c creep damage using MPC Omega method                                          
      statev(9) = statev(9) + dtime/tfact/zlife5	                            
90    RETURN                                                                   
      END                                                                      
c *******************************************************************  
c Calculate life based on WRC 541 LMP Method                           
c *******************************************************************  
      real*8 function xlmlife(CMNAME,tempr,stress,flag)                
      implicit real*8(a-h,o-z)                                         
      real*8 lifemin, lifeavg                                         
      character*80 CMNAME                                              
c********************************************************************  
      IF (CMNAME.EQ.'CS_1-25_CR') Then
      cmin = 22.0548d0
      cavg = 21.55d0
      a0 = 46354.38d0
      a1 = -6946.603d0
      a2 = -343.6751d0
      a3 = 0.0d0
      End if                                                           
c********************************************************************  
      IF (CMNAME.EQ.'CS_1-25_CR_WSRF') Then
      cmin = 22.0548d0
      cavg = 21.55d0
      a0 = 46354.38d0
      a1 = -6946.603d0
      a2 = -343.6751d0
      a3 = 0.0d0
      End if                                                           
c********************************************************************  
      xlogs = log(stress)/log(10.0d0)                                  
      xwrc = a0+a1*xlogs+a2*xlogs**2+a3*xlogs**3          
      lifemin = xwrc/tempr - cmin 
      lifemin = 10.0d0**lifemin 
      lifeavg = xwrc/tempr - cavg 
      lifeavg = 10.0d0**lifeavg 
c min properties, flag=0                                               
	  IF (flag.EQ.0.0d0) Then                                           
	    xlmlife=lifemin 
	  ELSE IF (flag.EQ.1.0d0) Then                                           
	    xlmlife=lifeavg 
	  End If
90    RETURN                                                           
      End                                                              
c********************************************************************  
c Calculate life based on the LMP method (API 530)					    
c********************************************************************  
c                                                                      
      real*8 function xlmlife1(CMNAME,tempr,stress,flag)               
	  implicit real*8(a-h,m-z)                                          
      character*80 CMNAME                                              
c     Equation 2 used for material pointers 7,8,9 and 13
      eqn=0.0d0
c********************************************************************  
      IF (CMNAME.EQ.'CS_1-25_CR') Then
      ma0 = 41.44429d0
      ma1 = -0.0016608091d0
      ma2 = 0.0d0
      ma3 = -2.5842632d0
      ma4 = 0.0d0
      ma5 = 0.0d0
      ma6 = 0.0d0
      cmin = 20.0d0
      aa0 = 42.60525d0
      aa1 = 0.0d0
      aa2 = 0.0d0
      aa3 = -2.6236052d0
      aa4 = 0.0d0
      aa5 = 0.0d0
      aa6 = 0.0d0
      cavg = 20.0d0
      End if                     
c********************************************************************  
      IF (CMNAME.EQ.'CS_1-25_CR_WSRF') Then
      ma0 = 41.44429d0
      ma1 = -0.0016608091d0
      ma2 = 0.0d0
      ma3 = -2.5842632d0
      ma4 = 0.0d0
      ma5 = 0.0d0
      ma6 = 0.0d0
      cmin = 20.0d0
      aa0 = 42.60525d0
      aa1 = 0.0d0
      aa2 = 0.0d0
      aa3 = -2.6236052d0
      aa4 = 0.0d0
      aa5 = 0.0d0
      aa6 = 0.0d0
      cavg = 20.0d0
      End if                     
c********************************************************************  
	  xapi0=ma0+(ma1*stress)+(ma2*stress**2)+(ma3*log(stress))          
	  xapi1=aa0+(aa1*stress)+(aa2*stress**2)+(aa3*log(stress))          
	  xmin1=1000*xapi0/tempr-cmin                                        
	  xmin1=10.0d0**xmin1                                                 
	  xavg1=1000*xapi1/tempr-cavg                                        
	  xavg1=10.0d0**xavg1                                                
c********************************************************************  
	  xapi2n=(ma0+(ma2*stress**0.5)+(ma4*stress)+(ma6*stress**1.5))          
	  xapi2d=(1+(ma1*stress**0.5)+(ma3*stress)+(ma5*stress**1.5))          
	  xapi2=xapi2n/xapi2d          
	  xapi3n=(aa0+(aa2*stress**0.5)+(aa4*stress)+(aa6*stress**1.5))          
	  xapi3d=(1+(aa1*stress**0.5)+(aa3*stress)+(aa5*stress**1.5))          
	  xapi3=xapi2n/xapi2d          
	  xmin2=1000*xapi2/tempr-cmin                                        
	  xmin2=10.0d0**xmin2                                                 
	  xavg2=1000*xapi3/tempr-cavg                                        
	  xavg2=10.0d0**xavg2                                                 
c********************************************************************  
c min properties, flag=0                                               
	  IF (flag.EQ.0.0d0) Then                                           
	      IF (eqn.EQ.7.0d0) Then
	          xlmlife1=xmin2
	      ELSE IF (eqn.EQ.8.0d0) Then
	          xlmlife1=xmin2
	      ELSE IF (eqn.EQ.9.0d0) Then
	          xlmlife1=xmin2
	      ELSE IF (eqn.EQ.13.0d0) Then
	          xlmlife1=xmin2
	      ELSE IF (eqn.EQ.0.0d0) Then
	          xlmlife1=xmin1
	      End If
	  End if                                                            
c *********************************************************************
c avg properties, flag=1                                               
	  IF (flag.EQ.1.0d0) Then                                           
	      IF (eqn.EQ.7.0d0) Then
	      xlmlife1=xavg2
	      ELSE IF (eqn.EQ.8.0d0) Then
	      xlmlife1=xavg2
	      ELSE IF (eqn.EQ.9.0d0) Then
	      xlmlife1=xavg2
	      ELSE IF (eqn.EQ.13.0d0) Then
	      xlmlife1=xavg2
	      ELSE IF (eqn.EQ.0.0d0) Then
	      xlmlife1=xavg1
	      End If
	  End if                                                            
c                                                                      
90    RETURN                                                           
      End                                                              
c 	                                                                    
c *********************************************************************
c  Calculate life based on the MPC Omega Method                        
c *********************************************************************
	  real*8 function omegalife(CMNAME,tempr,stress,sig1,sig2,sig3)            
	  implicit real*8(a-h,m-z)                                            
	  character*80 CMNAME                                                 
c 	  MPC Omega parameter (beta_omega)	                                
	  beta_omega=0.33d0
c 	  Alpha parameter	                                                
	  alpha=2.0d0
c	MPC Omega strain rate data from API-579                             
      IF (CMNAME.EQ.'CS_1-25_CR') Then
      a0=-23.5d0
      a1=52610.0d0
      a2=-4500.0d0
      a3=-5190.0d0
      a4=825.0d0
      b0=-2.65d0
      b1=6110.0d0
      b2=3000.0d0
      b3=-4440.0d0
      b4=1375.0d0
c 	  MPC Omega adjustment factor	                                    
	  del_cd=0.0d0
	  del_sr=0.0d0
	  End if
***********************************************************************
      IF (CMNAME.EQ.'CS_1-25_CR_WSRF') Then
      a0=-23.5d0
      a1=52610.0d0
      a2=-4500.0d0
      a3=-5190.0d0
      a4=825.0d0
      b0=-2.65d0
      b1=6110.0d0
      b2=3000.0d0
      b3=-4440.0d0
      b4=1375.0d0
c 	  MPC Omega adjustment factor	                                    
	  del_cd=0.0d0
	  del_sr=0.0d0
	  End if
***********************************************************************
      s1=log(stress)/log(10.0d0)
      del_omega=beta_omega*((sig1+sig2+sig3)/stress-1.0d0)
      eta_bn=-((a2+2*a3*s1+3*a4*s1**2)/tempr)
      omega=10**((b0+del_cd)+((b1+b2*s1+b3*s1**2+b4*s1**3)/tempr))
      omega_n=max((omega-eta_bn),3.0d0)
      omega_m=omega_n**(del_omega+1.0d0)+alpha*eta_bn
      eps_dot=10**(-((a0+del_sr)+(a1+a2*s1+a3*s1**2+a4*s1**3)/tempr))
      omegalife=1/eps_dot/omega_m
50    return                                                           
      end                                                              
c 	                                                                    
c *********************************************************************
c  Subroutine to get Sinh model parameters based on MPC Omega Data
c *********************************************************************
	  subroutine matprop(CMNAME,acon,sig0,expn,wsrf)                     
	  implicit real*8(a-h,m-z)                                            
	  character*80 CMNAME                                                 	                                         
c	  Sinh model parameters calculated from MPC Omega model from API-579
c *********************************************************************
c     CS_1-25_CR at 20 ksi and 950 F
      IF (CMNAME.EQ.'CS_1-25_CR') Then
      acon=0.0030244057936187d0
      sig0=36.935856d0
      expn=8.958572d0
	  wsrf=1.0d0
	  End if
c                                                                      
***********************************************************************
c     CS_1-25_CR_WSRF at 20 ksi and . F
      IF (CMNAME.EQ.'CS_1-25_CR_WSRF') Then
      acon=0.0030244057936187d0
      sig0=36.935856d0
      expn=8.958572d0
	  wsrf=0.770d0
	  End if
c                                                                      
***********************************************************************
50    return                                                           
      end                                                              
