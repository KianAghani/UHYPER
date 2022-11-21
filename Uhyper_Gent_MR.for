************************************************************************
*  UHYPER sub. for a combination of a Mooney-Rivlin and Gent 
*  material model.
*
* U = Udeviatoric + Uvolumetric
*
* Udeviatoric = a1(I1-3) + a2(I1-3)^2 + a3(I1-3)^3 + a4(I1-3)^4
*             + b1(I2-3) + b2(I2-3)^2 + b3(I2-3)^3 + b4(I2-3)^4
*             - 0.5*mu*Jm*ln(1 - (I1-3/Jm))
*      
*
* Uvolumetric = 0.5*Bmod*(J-1)^2
*
************************************************************************
*
*     State Variables
*     --------------------------------------------------------------
*     statev(1) = I1 - 3
*     statev(2) = I2 - 3
*     statev(3) = sqrt(I1/3)
*
*
*     Material Properties 
*     --------------------------------------------------------------
*     a1    = props(1) --- I1 coefficient Mooney-Rivlin
*     a2    = props(2) --- I1 coefficient Mooney-Rivlin
*     a3    = props(3) --- I1 coefficient Mooney-Rivlin
*     a4    = props(4) --- I1 coefficient Mooney-Rivlin
*     b1    = props(5) --- I2 coefficient Mooney-Rivlin
*     b2    = props(6) --- I2 coefficient Mooney-Rivlin
*     b3    = props(7) --- I2 coefficient Mooney-Rivlin
*     b4    = props(8) --- I2 coefficient Mooney-Rivlin
*     mu    = props(9) --- Shear modulus Gent
*     Jm    = props(10) -- lock parameter Gent
*     Bmod  = props(11) -- Bulk modulus
*
************************************************************************
      SUBROUTINE uhyper(BI1,BI2,AJ,U,UI1,UI2,UI3,TEMP,NOEL,
     + CMNAME,INCMPFLAG,NUMSTATEV,STATEV,NUMFIELDV,FIELDV,
     + FIELDVINC,NUMPROPS,PROPS)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION U(2),UI1(3),UI2(6),UI3(6),STATEV(*),FIELDV(*),
     + FIELDVINC(*),PROPS(*)

      REAL*8 I1minus3,I2minus3,a1,a2,a3,a4,b1,b2,b3,b4,mu,Jm,Bmod,gentFac,effStr

      REAL*8 ZERO,ONE,TWO,HALF,THREE,FOUR,FIVE
      parameter(ZERO=0.d0,ONE=1.d0,TWO=2.d0,HALF=0.5d0,THREE=3.d0,FOUR=4.d0,FIVE=5.d0)


      ! Obtain material properties
      
      a1    = props(1)
      a2    = props(2)
      a3    = props(3)
      a4    = props(4)
      b1    = props(5)
      b2    = props(6)
      b3    = props(7)
      b4    = props(8)
      mu    = props(9)
      Jm    = props(10)
      Bmod  = props(11)


      U = ZERO
      UI1 = ZERO
      UI2 = ZERO
      UI3 = ZERO


      ! Avoid division by ZERO or infinite

      I1minus3 = BI1 - THREE
      I2minus3 = BI2 - THREE
	  
      if(Jm.le.1.d-6) then
         gentFac = ZERO
      else
         gentFac = I1minus3/Jm
      endif
	  
      if(gentFac .GT. 0.95d0) gentFac = 0.95d0


      U(1) = a1*I1minus3 + a2*(I1minus3**TWO) + a3*(I1minus3**THREE) +
     +     a4*(I1minus3**FOUR) + b1*I2minus3 + b2*(I2minus3**TWO) +
     +     b3*(I2minus3**THREE) + b4*(I2minus3**FOUR) -
     +     HALF*mu*Jm*dlog(ONE - gentFac) 
	 
      ! Check for compressible material
      
      if(incmpflag.eq.0) then
         U(1) = U(1) + HALF*Bmod*(AJ - ONE)**TWO
      endif

      ! First derivatives
      
      UI1(1) = a1 + TWO*a2*I1minus3 + THREE*a3*(I1minus3**TWO) +
     +     FOUR*a4*(I1minus3**THREE) + HALF*mu*((ONE - gentFac)**(-ONE)) 
	 
	 
      UI1(2) = b1 + TWO*b2*I2minus3 + THREE*b3*(I2minus3**TWO) +
     +     FOUR*b4*(I2minus3**THREE)
      
      ! Second derivatives
      
      if(Jm.le.1.d-6) then
         UI2(1) = TWO*a2 + 6.d0*a3*I1minus3 + 12.d0*a4*(I1minus3**TWO) 

      else
         UI2(1) = TWO*a2 + 6.d0*a3*I1minus3 + 12.d0*a4*(I1minus3**TWO) + 
     +        HALF*(mu/Jm)*((ONE - gentFac)**(-TWO)) 
	 
      endif
      UI2(2) = TWO*b2 + 6.d0*b3*I2minus3 + 12.d0*b4*(I2minus3**TWO)


      ! Derivatives for compressible material
      if(incmpflag.eq.0) then
         UI1(3) = Bmod*(AJ - ONE)
         UI2(3) = Bmod
      endif

      ! Compute the effective stretch
      
      effStr = dsqrt(dabs(BI1/THREE))


      ! Update the state variables
      
      statev(1) = I1minus3
      statev(2) = I2minus3
      statev(3) = effStr


      return
      end 

************************************************************************
