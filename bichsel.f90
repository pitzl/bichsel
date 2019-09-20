
! Hans Bichsel, Seattle, 1984
! Rev Mod Phys 60, 663 (1988)
! Straggling in thin silicon detectors

! Calculates the differential cross section for delta rays in silicon
! for incident particles of given mass, charge, and energy
! using tbles of dielectric 'constants', Generalised Oscillation Strengths,
! and 
! Calculates the "Landau" energy loss distribution for a silicon
! detector of given thickness by convolution

! Output:
! eesig.dat    differential cross section for delta rays
! CONV.OPA     Control debugging output 
! CONV.SPE     Folded E-loss spectrum      

! gfortran -O2 -o bichsel bichsel.f90
! bichsel (prompts for kinetic energy of e)

program LOSS

  ! Convolution program with realistic energy loss spectrum for silicon

  ! Input files:
  ! 14      heps.tab        in subr EPRED
  ! 15      MACOM.TAB       in subr AERED
  ! 16      EMERC.TAB       in subr EMRED

  ! Fortran is NOT case sensitive!
  ! nume and NUME are the same

  implicit none

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real saxk, Etop, bemx, FSG, zi, su(8)
  integer nels(2)
  common / NML / saxk, Etop, bemx, FSG, zi, su, nels

  real sig(6, 1252), stp(5), tsig(5), rM2(5), rim(1252)
  common / SPTT / sig, stp, tsig, rM2, rim

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry
  common / ABS / Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  real dec
  common / const / dec

  integer j

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  open( unit=3, file='CONV.OPA', status='replace' )

  call PREP
  call PREPE
  call EPRED
  call AERED
  call EMRED
  call SPECT

  open( unit=9, file='CONV.SPE', status='replace' )
  write( 9,*) betasq, tdedx
  write( 9,*) ' CONV.SPE , t=', exth

  print *
  print *, ptM, "  ", pkE
  do j = 1, nume
     print *, j, "  ", E(j), "  ", sig(5,j)
  enddo
  print *

  open( unit=1, file='eesig.dat',status='replace')

  write( 1, * ) "# mass ", ptM, " MeV"
  write( 1, * ) "# Ekin ", pkE, " MeV"
  write( 1, * ) "# path ", 1E4/FSG, " um"
  write( 1, * ) "# mean dE/dx ", 100*dEdx, " eV/um"
  write( 1, * ) "#         j      dE [eV]      dE*dE*dsig/ddE [eV*eV/cm]"
  write( 1, * ) "START"
  do j = 1, nume
     write( 1, * ) j, "  ", E(j), "  ", dec*sig(5,j)
  enddo

  call CONV

END program LOSS

!-----------------------------------------------------------------------
subroutine EPRED

  ! HEPS.TAB is the table of the dielectric constant for solid Si,
  ! epsilon = ep(1,j) + i*ep(2,j), as a function of energy loss E(j),
  ! section II.E in RMP, and rim is Im(-1/epsilon), Eq. (2.17), p.668.
  ! This is used for the cross section of small momentum transfer excitations.    

  implicit none

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real ep(2, 1252), dfdE(1252), BB(2)
  common / EP12 / ep, dfdE, BB

  real sig(6, 1252), stp(5), tsig(5), rM2(5), rim(1252)
  common / SPTT / sig, stp, tsig, rM2, rim

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer n2t, numt, j, jt
  real etbl

  open( unit = 14, file='bichdat/heps.tab', status='old' )

  read( 14, * ) n2t, numt
  print*, ' EPRED, n2t, numt=', n2t, numt
  if( nume .ne. numt ) print*, ' CAUTION: nume & numt differ'
  if( n2 .ne. n2t ) print*, ' CAUTION: n2 & n2t differ'
  do j = 1, numt
     read( 14, * ) jt, etbl, ep( 1, j ), ep( 2, j ), rim( j )
     dfdE( j ) = rim( j ) * 0.0092456 * E( j )
     !     if( ( j/20 )*20 .eq. j ) print 304, j, jt, E( j ), etbl, rim( j ), dfdE( j )
     !     304            format( ' EP:', 2i4, 2f11.2, 2f12.6 )
  enddo

END subroutine EPRED

!-----------------------------------------------------------------------
subroutine AERED

  ! MACOM.TAB is the table of the integrals over momentum transfer K of the
  ! generalized oscillator strength, summed for all shells, i.e. the A(E)
  ! of Eq. (2.11), p. 667 of RMP
  ! longitudinal excitation ( K+L shells ) with large momentum tansfers.   

  implicit none

  real sig(6, 1252), stp(5), tsig(5), rM2(5), rim(1252)
  common / SPTT / sig, stp, tsig, rM2, rim

  integer n2t, numt, j, jt
  real etbl

  open( unit = 15, file = 'bichdat/macom.tab', status='old' )

  read( 15, * ) n2t, numt
  print*, ' AERED, n2t, numt=', n2t, numt
  do j = 1, numt
     read( 15, * ) jt, etbl, sig( 6, j )
     !     if( ( j/20 )*20 .eq. j ) print 304, j, jt, E( j ), etbl, sig( 6, j )
     !     304            format( ' AE:', 2i4, 2f11.2, f12.6 )
  enddo

END subroutine AERED

!-----------------------------------------------------------------------
subroutine EMRED

  ! EMERC.TAB is the table of the integral over K of generalized oscillator
  ! strength for E < 11.9 eV with Im(-1/epsilon) from equations in the Appendix
  ! of Emerson et al., Phys Rev B7, 1798 (1973) (also see CCS-63)

  implicit none

  real sig(6, 1252), stp(5), tsig(5), rM2(5), rim(1252)
  common / SPTT / sig, stp, tsig, rM2, rim

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  character*8 tit( 9 )
  integer j, jt
  real etbl

  open( unit = 16, file = 'bichdat/emerc.tab', status='old' )

  do j=1, 4
     read( 16, 1615 ) tit
1615 format( 9a8 )
     print 1615, tit
  enddo
  do j=1, 175
     read( 16, * ) jt, etbl, sig( 6, j ), xkmn( j )
     if( jt .le. 31 ) print 304, j, jt, E( j ), etbl, sig( 6, j ), xkmn( j )
     ! if( ( j/20 )*20 .eq. j ) print 304, j, jt, E( j ), etbl, sig( 6, j )
304  format( ' EM:', 2i4, 2f11.2, 2f12.6 )
  enddo

END subroutine EMRED

!-----------------------------------------------------------------------
subroutine PREP

  ! Initialization routine for user parameter input & basic constants 
  ! Particle types NPM:
  ! = 1    proton
  ! = 2    pion
  ! = 3    alpha
  ! = 4    electron ( positron ) 
  ! = 5    kaon

  implicit none

  real Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry
  common / ABS / Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real saxk, Etop, bemx, FSG, zi, su(8)
  integer nels(2)
  common / NML / saxk, Etop, bemx, FSG, zi, su, nels

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real PMASS( 5 )
  Data PMASS / 938.256, 139.578, 3727.328, 0.511004, 497.034 /

  real emk

  pi = 3.14159265359
  Ry = 13.6058

  write( 3, 10 ) Ry
10 format( /, 1x, 'PREP:   Ry = ', F12.5, ' eV' )

  write( *, '( a, $ )' ) 'Particle type ( 1=P, 2=Pi, 3=Alpha, 4=e, 5=K ) : '
  !     Read*, npm
  npm = 4 ! DP
  print *, npm

  write( *, '( a, $ )' ) 'Silicon thickness ( microns ) : '
  !     Read*, exth
  exth = 150
  print *, exth

  ! convert microns into cm

  exth = exth / 1e4

  ! set particle mass values ( MeV ) according to input code number  

  PTM = PMASS( npm )

  ! charge of incident particle

  zi = 1.0
  if( npm .eq. 3 ) zi = 2.0

  ! properties of absober material ( Silicon ) 
  ! ZA = atomic number 
  ! AW = atomic weight
  ! rho = density  ( g/cm**3 ) 
  ! atnu = # of atoms/cm**3 

  ZA = 14.0
  AW = 28.086
  rho = 2.329
  atnu = 6.0222e23 * rho / Aw

  write( 3, 601 ) PTM, zi
601 format( 1x, 'PREP:  particle mass=', F10.3, ' MeV, charge=', f3.0 )

  ! Initialization kinematic parameters

  call EVANS

  Efin = Emax
  saxk = 153540. * zi**2 * rho / ( betasq*AW )
  ! Saxon Eq ( 3a ) for k
  Emk  = saxk * Emax

  write( 3, 602 ) Za, Aw, exth, saxk, Emk
602 format( /, 1x, 'Z=', F6.2, '  A=', f9.4, 3x, 't=', f9.5, 'cm', &
         '  k/Z=', f10.2, 3x, 'k*Emax/Z=', e12.5 )

END subroutine PREP

!-----------------------------------------------------------------------
subroutine EVANS

  ! Initialization of kinematic parameters 

  implicit none

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  integer jkm
  real xxx, w, pmom, telm, emx

  ! Evans, p 891, Uehling Eq ( 4a ).  Date : 15 June 1984

  print*, ' Evans: particle mass ( PTM ) = ', ptM, ' MeV'

  print*, ' Input option: 1=E_kin ( MeV ), 2=P_mom( MeV/c ), 3=bet*gam'
  !      write( *, '( a, $ )' ) 'Select option: '
  !      read*,    jkm
  jkm = 1 ! DP
  print *, jkm

  write( *, '( a, $ )' ) 'Input value  : '
  read*,    xxx

  go to( 11, 22, 33 ), jkm

  ! Calculate the following     W = Gamma  ( E_total/mass ) 
  ! pKe = Particle kinetic energy 
  ! bg = beta*gamma

  !---  input was kinetic energy

11 pkE = xxx
  W = xxx/ptM + 1.0
  bg = sqrt( W**2 - 1.0 )
  print *, pkE, "  ", bg
  go to 34

  !--   input was momentum  

22 pmom = xxx
  bg = xxx / ptM
  W = sqrt( bg**2 + 1.0 )
  pkE = ptM * ( W - 1.0 )
  go to 35 

  !--   input was beta*gamma 

33 bg = xxx                
  W = sqrt( bg**2 + 1.0 )
  pkE = ptM * ( W - 1.0 )

  !     further parameters:  pmom = particle momentum ( MeV )  
  !     ptE  = particle total energy ( MeV )
34 pmom = ptM * bg
35 betasq = bg**2 / ( 1 + bg**2 )
  ! beta = bg / W
  gam  = W
  ! ptE  = ptM * W

  ! Maximum energy transfer   Emax  ( MeV ) 
  ! Uehling, also Sternheimer & Peierls Eq.( 53 )

  telm = 2 * 0.511004
  Emax = ptM * ( W**2 - 1.0 ) / ( ptM/telm + telm/ptM + W )
  ! cannot distinguish i/o electrons
  if( npm .eq. 4 ) Emax = 0.5*pkE

  print*, 'particle type = ', npm, '   Emax = ', Emax, ' MeV'
  Emx  = telm * bg**2

  write( 3, 333 ) bg, pmom, pkE
333 format( /, 3x, 'beta*gamma=', f11.4, 3x, 'momentum=', f13.4, ' MeV/c', &
       3x, 'E kinetic of incident particle=', f15.2, ' MeV' )
  write( 3, 334 ) betasq, gam, Emax, Emx
334 format( 3x, 'beta**2=', f9.6, 3x, 'gamma=', f12.5, 3x, 'Emax=', &
       2e12.4, ' MeV'/ )

  Emax = 1e6 * Emax! [eV]

END subroutine EVANS

!-----------------------------------------------------------------------
subroutine PREPE

  ! Definitions of energy scale ( log ) bin size 

  implicit none

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  integer n1, nu
  real F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4
  common / ut / N1, NU, F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real saxk, Etop, bemx, FSG, zi, su(8)
  integer nels(2)
  common / NML / saxk, Etop, bemx, FSG, zi, su, nels

  integer ken, l
  real exs

  ! n2 = number of bins for each factor of 2 in energy 

  n2   = 64
  nume = 650
  if( n2 .eq. 64 ) nume = 1250
  u    = log( 2. ) / n2
  um   = exp( u )
  ken  = log( 1839. / 1.5 ) / u
  Emin = 1839. / 2**( 1.0*ken/n2 )
  E( 1 ) = Emin
  EXS = 1.

  print*, ' PREPE', n2, ken, E( 1 ), emin
  write( 3, 609 ) n2, Emin, u, um
609 format( /1x, 'PREPE:   N2=', I4, 3x, &
       3x, 'Emin=', f8.3, 3X, 'u=', F9.6, ' e**u=', f10.6 )

  lemx = nume + 450
  do L = 1, lemx
     EXS = EXS * um
     E( L+1 ) = E( L ) * um
     if( ( L/50 )*50 .eq. L ) print*, ' L, E=', L, E( L ), exs, um
     DI( L )  = -alog( 1.0 - 1.0/EXS ) / u
     DE( L )  = E( L+1 ) - E( L )
     if( L .le. nume ) H( L )   = 0.
     if( E( L ) .le. Emax ) leh = L
  enddo

  if( leh .gt. nume ) leh = nume
  Etop = E( nume ) * sqrt( um )
  if( Efin .gt. Etop ) Efin = Etop
  write( 3, * ) ' PREPE: Efin, Etop, Emax=', Efin, Etop, Emax
  print*, ' PREPE: Efin, Etop, Emax=', Efin, Etop, Emax

END subroutine PREPE

!-----------------------------------------------------------------------
subroutine SPECT

  ! generate collision spectrum from ep-1, 2 and ae

  implicit none

  real dec
  common / const / dec

  real Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry
  common / ABS / Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry

  real sig(6, 1252), stp(5), tsig(5), rM2(5), rim(1252)
  common / SPTT / sig, stp, tsig, rM2, rim

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real saxk, Etop, bemx, FSG, zi, su(8)
  integer nels(2)
  common / NML / saxk, Etop, bemx, FSG, zi, su, nels

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  real ep(2, 1252), dfdE(1252), BB(2)
  common / EP12 / ep, dfdE, BB

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  real elm, fac, blg, s0, s1, avi, avi1, pf
  real tmcb, uef, q1, qmin, epbe, thet, sgh, rmf, sgg
  integer jpr, jpd, l, j

  elm = 511004.
  fac = 8. * pi * Ry**2 * ( 0.529177e-8 )**2 / ( elm * betasq )
  DEC = zi**2 * atnu * fac

  blg = alog( ( 2.*elm ) * bg**2 ) - betasq
  write( 3, 307 ) betasq, atnu, blg
307 format( /4X, 'SPECT F.307:  beta**2=', F12.10, 4X, '# of', &
         ' atoms per cm**3=', e12.4, 3x, 'blg=', f9.4, / )
  !     write( 3, 308 )
  !     308    format( 3x, 'j', 5x, 'E/eV', 5x, 'df/dE ', 5x, 'sgg', 6x, 'sgh', 
  !     1  7x, 'S 1', 6x, 'S 3', 6x, 'S 4', 5x, 'sum S', 5x, 'S( 0 )', 3x, 'dE/dx'/ )

  jpr = 5
  jpd = 5
  do L=1, 5
     rM2( L )  = 0
     Tsig( L ) = 0
     STP( L )  = 0
  enddo
  S0 = 0
  S1 = 0
  avI = 0
  avI1 = 0
  bemx = betasq / Emax
  pf = pkE * 1e6
  tmcb = 2. * elm * betasq

  do 5 j = 1, nume

     if( E( j ) .gt. Emax ) go to 11
     if( npm .eq. 4 ) then
        uef = 1 + ( E( j )/( pf-E( j ) ) )**2 + ( ( ( gam-1 ) / gam ) &
             * E( j )/pf )**2 - ( 2*gam - 1 )*E( j )/( gam**2 * ( pf - E( j ) ) )
     else
        uef = 1 - E( j ) * bemx
     endif
     ! uef from Uehling Eqs. 9 & 2
     if( j .eq. 1 ) print*, ' uef=', uef
     S0   = S0   + dfdE( j ) * dE( j )
     avI  = avI  + dfdE( j ) * alog( E( j ) ) * dE( j )
     avI1 = avI1 + dfdE( j ) * E( j ) * alog( E( j ) ) * dE( j )
     S1   = S1   + dfdE( j ) * E( j ) * dE( j )
     Q1 = Ry
     !ee   red CCS-33, 39 & 47
     if( E( j ) .lt. 100. )  Q1 = 0.025**2 * Ry
     if( E( j ) .lt. 11.9 ) Q1 = xkmn( j )**2 * Ry
     Qmin = E( j )**2 / tmcb
     sig( 1, j ) = 0
     if( E( j ) .lt. 11.9 .and. Q1 .le. Qmin ) go to 14
     sig( 1, j ) = E( j ) * dfdE( j ) * alog( Q1 / Qmin ) 
14   epbe = 1 - betasq * ep( 1, j )
     ! Fano Eq 47
     if( epbe .eq. 0 ) epbe = 1e-20
     sgg = E( j ) * dfdE( j )*( -.5 )*alog( epbe**2+( betasq*ep( 2, j ) )**2 )
     thet = atan( ep( 2, j ) * betasq / epbe )
     if( thet .lt. 0 ) thet = thet + pi         
     ! plausible-otherwise I'd have a jump
     ! Fano says [p 21]: 'arctan approaches pi for betasq*eps1 > 1'
     sgh = E( j )**2 *( betasq-ep( 1, j ) / ( ep( 1, j )**2+ep( 2, j )**2 ) )*thet
     sgh = 0.0092456 * sgh
     sig( 3, j ) = sgg + sgh
     sig( 4, j ) = 2. * sig( 6, j ) * uef
     ! the integral was over  d lnK rather than  d lnQ
     !     if( ( j/10 )*10 .eq. j ) print 327, E( j ), sgg, sgh, ( sig( ii, j ), ii=1, 4 )
     !     327            format( 1x, f11.2, 1p6e11.3 )
     sig( 2, j ) = 0
     sig( 5, j ) = 0
     do  27 L=1, 4
        Tsig( L )  = Tsig( L )  + sig( L, j ) * dE( j ) / E( j )**2
        STP( L )   = STP( L )   + sig( L, j ) * dE( j ) / E( j )
        rM2( L )   = rM2( L )   + sig( L, j ) *  dE( j )
        sig( 5, j ) = sig( 5, j ) + sig( L, j )
27   enddo
     Tsig( 5 ) = Tsig( 5 ) + sig( 5, j ) * dE( j ) / E( j )**2
     STP( 5 )  = STP( 5 )  + sig( 5, j ) * dE( j ) / E( j )  
     rM2( 5 )  = rM2( 5 )  + sig( 5, j ) * dE( j )
     !     if( j .eq. 1 ) go to 28
     !     if( j .ge. 320 .and. j .le. 326 ) go to 28
     !     if( ( j/10 )*10 .ne. j ) go to 5
     !     28     write( 3, 608 ) j, E( j ), dfdE( j ), sgg, sgh, sig( 1, j ), ( sig( L, j ), L=3, 5 ), 
     !     1          S0, STP( 5 )
     !     608            format( 1x, i4, f9.1, 1pe11.3, 0p9f9.4 )

5 enddo

11 write( 3, * ) '  uef=', uef
  write( 3, 374 ) Tsig, STP, rm2
374 format( /9x, 'Integ. over sig =', 5F12.4 / 2( 28x, 5f12.3 / ) )
  write( 3, 375 ) S0, avI, S1, avI1
375 format( /9x, ' S( 0 )=', f9.5, 3x, 'ln( I )=', f10.5, 3x, &
       'S( 1 )=', f10.3, 3x, 'L( 1 )=', f10.3 / )
  write( 3, * ) '  following data without density effect'
  write( 3, * ) '  S( 0 )*blg=', S0*blg, '   2*L( 0 )=', 2*avI
  write( 3, * ) '  S( 1 )*blg=', S1*blg, '   2*L( 1 )=', 2*avI1
  print*, ' S( 0 )=', S0, '  L( 0 )=', avI
  FSG  = Tsig( 5 ) * DEC
  dEdx = STP( 5 ) * ( dec/1E6 )
  rmf  = rM2( 5 ) * ( dec/1E6 )
  write( 3, 388 ) S0, FSG, dEdx, rmf
388 format( /, 10X, 'Zeff=', F7.3, 4X, '# coll/cm=', f11.3, 4x, &
       'dE/dx=', F9.4, ' MeV/cm', 3x, 'M2=', f12.4, ' keV**2/cm' )

  write( 3, * ) ' DEC=', dec, '  # atoms/cm**3=', atnu, '  fac=', fac

  call SPTS

END subroutine SPECT

!-----------------------------------------------------------------------
subroutine SPTS

  implicit none

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  real ep(2, 1252), dfdE(1252), BB(2)
  common / EP12 / ep, dfdE, BB

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real dec
  common / const / dec

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real sig(6, 1252), stp(5), tsig(5), rM2(5), rim(1252)
  common / SPTT / sig, stp, tsig, rM2, rim

  real sgm, stpw, secm, eps, bbb, fft, he2, rm0, sbb
  integer jpr, ja, j, nlast

  !     write( 3, 333 )
  !333  format( /4X, 'SPTS, F.333:', /15X, 'E', 7x, 'sig*E**2', 6x, 'sig', 
  !     1       10X, 'M 0', 11X, 'M 1', 11X, 'M 2', 10X, '<E>', / )
  SGM  = 0
  Stpw = 0
  SECM = 0
  jpr = 20
  ja  = 20

  do j=1, nume
     if( E( j ) .gt. Emax ) go to 77
     nlast= j
     he2  = sig( 5, j ) * dec
     H( j ) = he2 / E( j )**2
     SGM  = SGM + H( j )*dE( j )
     STPW = STPW + H( j ) * E( j ) * dE( j )
     SECM = SECM + he2 * dE( j )
     eps  = STPW / SGM
     !     if( j .lt. 5 ) go to 11
     !     if( j .eq. nume ) go to 11
     !     if( jpr .ne. j ) go to 75
     !     jpr = jpr + ja
     !11   write( 3, 654 ) j, E( j ), he2, H( j ), SGM, STPW, SECM, eps
     !     654            format( 1x, i6, f12.2, 1p7e13.5 )
  enddo

77 write( 3, 610 )  nlast, SGM, STPW, SECM
610 format( /1X, 'SPTS: nlast  ', i5, ' total cross section=', E15.5, &
       3X, 'dE/dx=', E15.5, ', M2=', E15.5/2x, 'see FSR-99 and CCS-9' / )
  write( 3, * ) ' final E=', Efin, ' Emax=', Emax, ' he2=', he2

  !ee   FSR-99
  sbb = 720.
  bbb = 1. - sbb*betasq / Emax
  fft = 14. * dec / 1e6
  write( 3, * ) ' sbb=', sbb, ' eV,   bbb=', bbb, '  fft=', fft
  rm0 = bbb * ( ( 1/Efin - 1/Emax ) + 2 * ( 1/Efin**2 - 1/Emax**2 ) ) &
       - betasq * alog( Emax/Efin ) / Emax
  write( 3, * ) ' residual M0=', rm0, rm0*fft, '/cm'
  write( 3, * ) '  If residual M0 is large, look for error'

  call HPART( bbb, fft, sbb, stpw, secm )

END subroutine SPTS

!-----------------------------------------------------------------------
subroutine HPART( bbb, fft, sbb, stpw, secm )

  implicit none

  real bbb, fft, sbb, stpw, secm

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry
  common / ABS / Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry

  real te, rst, rm2p, secmv, del2

  if( npm .eq. 4 ) then ! e
     !ee   FSR-143 and Uehling Eq 9
     TE = pkE * 1e6
     print*, ' ele', efin, emax, TE, gam
     rst = alog( Emax/Efin ) + alog( Emax ) - alog( TE-Efin ) - 1.0 / ( 1.0 - Efin/TE ) &
          + 2 + ( ( gam-1 ) / gam )**2 * ( 1.0/8.0 - 0.5 * ( Efin/TE )**2 ) &
          + ( ( 2*gam - 1.0 ) / gam**2 ) * ( alog( Emax ) - alog( TE-Efin ) )
  else
     rst = bbb * alog( Emax/Efin ) + &
          sbb * ( 1.0/Efin - 1.0/Emax ) - betasq * ( 1.0 - Efin/Emax )
  endif

  write( 3, * ) ' residual dE/dx=', rst, rst*fft, ' MeV/cm'

  tdedx = stpw/1e6 + rst*fft
  write( 3, * ) ' dE/dx=', tdedx, ' MeV/cm ', tdedx/2.329, ' MeV cm**2/g'
  print*, ' dE/dx=', tdedx, tdedx/rho

  secmv = secm / 1e6
  rM2p  = Efin - 0.5 * betasq * Efin**2 / Emax
  del2  = secmv - fft * rM2p
  print *, ' M2=', secmv, '  M2"=', fft*rM2p, ' M2-M2"=', del2
  write( 3, * ) ' M2=', secmv, '  M2"=', fft*rM2p, ' M2-M2"=', del2

END subroutine HPART

! Hans Bichsel, 1984
! Convolution subroutines

subroutine CONV

  implicit none

  real Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry
  common / ABS /  Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry

  integer n1, nu
  real F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4
  common / ut /   N1, NU, F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND /  N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real saxk, Etop, bemx, FSG, zi, su(8)
  integer nels(2)
  common / NML / saxk, Etop, bemx, FSG, zi, su, nels

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real xmc, xx
  integer jxt, k, l, kl
  real s
  real eav, stpp, onrm
  real dk, dq
  real pb, dmmpl
  integer n

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  write( 3, 303 )
303 format( /1x, 50(' *'))
  write( 3, * ) '  CONV','  t=', exth,'  fsg=', fsg

  XMC = exth * FSG
  ! number of collisions in thickness exth
  jxt = log(xmc) / log(2.0) + 1
  print*, '  jxt=', jxt

  xx = xmc / 2.0**jxt
  CZ0 = xx / 2.**10
  NU  = jxt + 10 + 1
  write( 3,386 ) jxt, NU, exth, XMC, XX, CZ0
386 format( 3x, 'F.386: NU=', 2i3, 4x, 't=', f10.7, ' cm', &
         4X, '# coll=', F11.3, 4X, 'xx=', F9.5, 4X, 'CZ0=', 1pe12.5)
  write( 9, *) CZ0, Emax, betasq, pkE, bg, PTM, zi
  write( 9, *) lemx, NU, U, um, Emin
  do k=1, lemx, 6
     write( 9, 384) (E(kl), kl=k, k+5)
  enddo
384 format( 1p6e13.6)

  H0  = 0.0
  CM1 = 1.0
  CM2 = 1.0
  XN  = 1.0
  EX  = 1.e-15              
  MIE = 0                 
  MIF = 0
  MIH = 0

  call NORMAL

  xn = 1.0
  D1 = CM1
  D2 = CM2 + CM1**2
  D3 = CM3 + 3.0*CM2*CM1 + CM1**3
  D4 = CM4 + 4.0*CM3*CM1 + 6.0*CM1**2*CM2 + CM1**4
  S  = D2 / D1
  write( 3, 612) D1, S, um
612 format( 2x, 'conv  F.612:', 3X, 'Initial distribution', /,  &
         '  delta 1=', 1PE12.4, '  delta 2 =', E12.4, '  exp(u)=', E12.5)
  ! write( 3, 635) (H(J), J=1, nume, 5)

  ! set parameters from single collision spectrum

  EAV  = D1
  ONRM = CM0
  STPP = 2. * EAV * saxk * ONRM / ZA
  print *, ' Eav', d1, cm0, stpp, cma

  DK   = 2. * dEdx / STPP / ZA * CMA
  DQ   = DK / Emax - 1.
  ! write( 3, 615) ONRM, CMA, EAV, STPP, DK, DQ
  ! 615            format(/, ' Orig spec - M0, M2, EAV, STP, M2, D2', /1P6E15.6)

  call SHRINK

  H0 = 1. - CZ0
  do L=1, leh
     H(L) = H(L) * CZ0
  enddo
  ! first convolution: H0 + H(E)
  CN  = CZ0
  CM1 = CZ0 * D1
  CM2 = CZ0 * D2
  thi = exth / 2**(jxt+10)
  xi  = saxk * thi * Za
  rkap = xi / Emax
  write( 3, 401)  NU, leh, Emin, EAV, DQ, Emax
401 format(/2x, 'F 401: ', 2i5, 1p5e14.7)

  do N1 = 1, NU

     print*,  'convol: N1', N1
     thi  = 2. * thi
     xi   = 2. * xi
     rkap = 2. * rkap
     CN   = 2. * CN
     write( 3, 303)
     write( 3, 317) N1, leh, MIE, MIH, CN, CZ0
317  format( /, ' CONVOL NUMBER=', i3, ' leh=', i4, 2x, 'MIE, MIH=',  &
          2i5, 4x,  'mean collision number=', 1PE12.4, '  CZ0=', f9.6, /)
     ! write( 3, 635) (H(J), J=1, leh, 5)
     ! 635            format(' H= ', 1P10E12.5)
     N2P = N2

     call FOLD

     ! Landau-Vavilov parameters
     PB = -0.4227843351 - log(rkap) - betasq
     dmmpl = xi * ( pb + 0.225 )
     write( 3, 310) thi, xi, rkap, pb, dmmpl
310  format( /1x, 't =', 1pe12.5, 2x, 'xi=', e12.5, 2x, 'kappa=', e12.5,  &
          2x, '<lam>=', 0pf12.6/3x, 'Landau theory <del>-dmp=', f12.3)
     N = MIH - MIE
     ! write( 3, 408) N1, leh, CN, xi, rkap, PB, PZERO   ! PZERO undefined ?
     ! 408            format( ' CONV,  F.408', 2i5, f9.3, 1p4e11.4/)

     write( 3, 408) N1, leh, CN, xi, rkap, PB 
408  format( ' CONV,  F.408', 2i5, f9.3, 1p3e11.4/)

     call OUTPUT

  enddo

END subroutine CONV

!-------------------------------------------------------------------------------
subroutine OUTPUT               

  implicit none

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer n1, nu
  real F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4
  common / ut / N1, NU, F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  real Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry
  common / ABS / Aw, atnu, dEdx, exth, pi, rho, thi, Za, xi, rkap, Ry

  integer npm, nzch
  real ptM, bg, betasq
  common / EVA / npm, nzch, ptM, bg, betasq

  real ASP(1252), ASS(1252)

  real b, c, d, s1, s2, s3, s4
  real x, hmax
  integer n, k, l, j, kl, kk
  real bax, dmp, bax1, dmp1, hhun
  integer npp, nskip, j2, l2, llow, lup, lmax

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  B = CM2 / (CM1**2)
  C = CM3 / sqrt(CM2**3)
  D = CM4 / (CM2**2)
  S1 = CN * D1
  write( 3, *) ' OUT', b, c, d, s1, cn
  S2 = CN * D2 / S1**2
  S3 = D3 / sqrt(D2**3 * CN)
  S4 = D4 / (D2**2 * CN) + 3.
  write( 3, 350) H0, cm1, B, C, D, S1, S2, S3, S4
350 format( 2x, 'OUTP: zero component =', 1pe16.4/30x, 'mean',  &
         9x, 'variance/mean**2', 4x, 'central3/var**1.5',  &
         20H     central4/var**2/3x, 'actual values=', 4e20.4, / &
         3x, 'theoret values', 4e20.4, /)
  write( 3, 340) cm1, s1, cm1/s1
  print 340,  cm1, s1, cm1/s1
340 format( 1x, 'OUTPUT:', 1p2e13.6, '  ratio=', 0pf10.6)

  X = 1.
  N = MIH - MIE
  if( N1 .le. nu-1 ) RETURN ! ???
  write( 9, *) leh, N1, N2, N2P, N, H0, thi, xi, rkap
  if( N2 .eq. N2P ) go to 42
  write( 9, *) lemx
  do k=1, lemx, 8
     write( 9, 950) (E(kl), kl=k, k+7)
  enddo
950 format( 8f10.3)
42 if( N2 .ne. N2P ) N2P = N2
  ASP(1) = 0
  ASS(1) = 0
  hmax   = 0
  do L=1, leh
     if( h(L) .gt. hmax ) then
        hmax = h(L)
        lmax = L
     endif
     J = L + N
     ASP(L+1) = ASP(L) + H(L) * DE(J)
     ASS(L+1) = ASS(L) + H(L) * E(J) * DE(J)
  enddo
  print 346,   lmax+N, E(lmax+N), h(lmax)
346 format( 1x, 'lmax+N=', i4, '  E=', f11.2, '  h=', f12.6)
  write( 3, 346) lmax+N, E(lmax+N), h(lmax)
  bax = (h(lmax) - h(lmax-1)) / (h(lmax) - h(lmax+1) )
  dmp = E(lmax+N) + 0.5 * dE(lmax+N) * (bax - 1.) / (bax + 1.)
  bax1 = (h(lmax-1) - h(lmax-2)) / (h(lmax-1) - h(lmax) )
  dmp1 = E(lmax+N-1) + 0.5 * dE(lmax+N-1)*(bax1 - 1.)/(bax1 + 1.)
  write( 3, *) ' bax, dmp=', bax, dmp, '  lower:', bax1, dmp1
  hhun = 0.01 * h(lmax)
  do L=1, leh
     if( L .lt. lmax .and. h(L) .lt. hhun) llow = L
     if( L .gt. lmax .and. h(L) .gt. hhun) lup  = L
  enddo
  write( 3, 355 ) llow, E(llow+N), h(llow), lup, E(lup+N), h(lup)
  print 355,  llow, E(llow+N), h(llow), lup, E(lup+N), h(lup)
355 format( ' lower & upper cutoff:', 2(i6, 0pf9.1, 1pe12.4))

  write( 9, *) lmax+N, E(lmax+N), h(lmax), ' lmax+N,  E,  h'
  write( 9, *) leh, dmp, '  leh,  dmp'
  do kl=1, leh, 6
     write( 9, 944) (H(kk), kk=kl, kl+5)
  enddo
944 format( 1p6e12.5)
  NPP   = (lup + llow) / 2
  nskip = 1
  npm = (lup - llow) / 2
  if( npm .gt. 70) nskip = 2
  npm = npm + nskip
  print*,  ' N, npp, npm, nskip=', N, npp, npm, nskip
  write( 3, 344)
344 format(  2(7x, 'j', 4x, 'E/eV', 8x, 'phi(E)', 7x, 'dE/dx', 9x, 'M2', 3x))
  do L=llow, NPP, nskip
     J  = L + N
     J2 = J + npm
     L2 = L + npm
     write( 3, 356) J, E(J), H(L), ASP(L), ASS(L), J2, E(J2), H(L2),  ASP(L2), ASS(L2)
  enddo
356 format( 3x, i5, 0pF10.1, 1p3E13.4, 3x, i5, 0pF10.1, 1p3E13.4)

END subroutine OUTPUT

!-------------------------------------------------------------------------------
subroutine FOLD

  implicit none

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer n1, nu
  real F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4
  common / ut / N1, NU, F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  integer l, lh, jh, lf, jf, k, lff, le
  real flf, s

  do L = 1, 1250
     F(L) = H(L)
     H(L) = 0.
  enddo
  F0  = H0
  LEF = leh
  MIF = MIH
  if( LEF .lt. 80 ) call RESET
  H0  = F0**2
  MIH = MIF + N2
  LEH = LEF
  if( leh .gt. 1250) leh = 1250
  print*,      '  FOLD: MIH, lef, leh, nume=', MIH, lef, leh, nume
  write( 3, *) '  FOLD: MIH, lef, leh, nume=', MIH, lef, leh, nume

  do 61 LH = 1, leh

     JH = LH + MIH

     do 62 LF = 1, LH
        JF = LF + MIF
        K = JH - JF
        FLF = JH - MIF - DI(K) + 1.E-20
        LFF = FLF
        LE  = JF - MIE
        S   = FLF - LFF
        if( LFF .eq. 0 ) LFF = 1
        H(LH) = H(LH) + F(LF) * ((1.0-S)*F(LFF)+S*F(LFF+1)) * DE(LE)
62   enddo
  ! On 18 July 1984,  I have some questions whether ST 62 is correct

     H(LH) = H(LH) - F(LH)**2 * 0.5*DE(LE)

61 enddo

  do L = 1, leh
     H(L) = H(L) * 2.0
  enddo
  if( F0 .gt. EX ) call ZERO
  call SHRINK
  call NORMAL

END subroutine FOLD

!-------------------------------------------------------------------------------
subroutine RESET                

  implicit none

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real Emin, Efin, Emax, gam, pkE
  common / ENER / Emin, Efin, Emax, gam, pkE

  real u, s
  integer ll, l, n, j

  ! LEF = LEH initially

  if( N2 .ge. 128 ) go to 702

  N2 = N2 * 2
  U = log(2.) / float(N2)
  write( 3, 375) MIE, MIH, N2
375 format( 1x, 'RESET: ', 2i4, '  Coordinate change: doubling of ',  &
         'point grid', i4, ' points on a factor of 2', /)
  write( 3, *) ' LEH, MIH, MIE=', LEH, MIH, MIE
  do LL = 1, LEF                  
     ! this is from top down
     L = LEF + 1 - LL
     F(2*L) = F(L)
  enddo
  N = 2*LEF
  do L = 4, N, 2                   
     ! N is just some number
     F(L-1) = (F(L) + F(L-2)) / 2.
  enddo
  LEF = 2 * LEF + 1
  MIF = 2 * MIF
  MIE = MIF
  do J = 1, leh
     S     = float(J+MIE)
     E(J)  = exp(S*U) * Emin
     DE(J) = E(J)*U
     S = J
     DI(J) = - log(1. - exp(-S*U)) / U
  enddo

702 MIE = MIF                       

  do J = 1, leh
     S     = float(J+MIE)
     E(J)  = exp(S*U) * Emin
     DE(J) = E(J)*U
  enddo

  write( 3, *) ' LEH, MIH, MIE=', LEH, MIH, MIE

END subroutine RESET

!-------------------------------------------------------------------------------
subroutine ZERO

  implicit none

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  integer n1, nu
  real F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4
  common / ut / N1, NU, F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4

  real xs
  integer n, l, k, ll, la, lb, m, lex

  xs = 0.
  N = MIH - MIE
  write( 3, *) ' zero ', mih, mie, n, leh, F0
  do L = 1, leh
     K = L+N
     xs = xs + H(L)*DE(K)
  enddo
  xs = ( 1.0 - F0 )**2 / xs 
  write( 3, *) k, xs, H(l), dE(K)
  do  L = 1, leh
     H(L) = H(L) * xs
  enddo
  N   = MIH - MIF
  MIH = MIF
  LEH = LEH + N
  if( leh .gt. nume) leh = nume
  LEX = LEH + 1
  do LL=1, leh
     LA = LEX - LL
     K  = LA + N
     H(K) = H(LA)
  enddo
  do LB = 1, N
     H(LB) = 0.0
  enddo
  print 684,  N, LEF, leh, LEX, MIE, MIF, MIH
684 format( ' zero:', i3, 2(3x, 3i5))
  do M = 1, LEF
     H(M) = H(M) + 2*F0*F(M)
  enddo
  ! F(M) is the H of the previous convolution
  ! This seems to be correct (18 July 1984)

END subroutine ZERO

!-------------------------------------------------------------------------------
subroutine NORMAL

  implicit none

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  integer n1, nu
  real F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4
  common / ut / N1, NU, F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real cma, cmb, cmd, d2, d3, d4, tdedx, tDD(2,4), xkmn(200)
  common / MEAN / CMA, CMB, CMD, D2, D3, D4, tdedx, tDD, xkmn

  real y, z, s, ec, cmq, t
  integer n, l, le

  Y = CM1 * 2.0
  Z = CM2 * 2.0
  write( 3, 714) MIE, MIH, H0, xn, y, z
714 format(' NORM', 2i5, '  H0=', 1pe12.5, '  xn, y, z=', 0p2f12.3, 1pe12.5)
  CM0 = H0
  CM1 = 0.
  CMA = 0.0
  N = MIH - MIE
  do L = 1, leh
     LE  = L + N
     S   = H(L) * dE(LE)
     CM0 = CM0 + S
     CM1 = CM1 + S*E(LE)
     cma = cma + s*E(LE)**2
  enddo
  if( cm0 .ne. 0 ) cmq = cm1 / cm0
  write( 3, 635 ) CM0, CM1, CMQ
635 format( 7x, 'area=', 1pe12.5, 3x, 'straight mean=',  &
         1pe12.5, 3X, 'CM1/CM0=', 0pf14.4)
  if( cm0 - H0 .ne. 0 ) T = ( 1.0 - H0 ) / (cm0 - H0)
  CM1 = CM1*T
  CM2 = 0.0
  CM3 = 0.0
  CM4 = 0.0
  do L = 1, leh
     H(L) = H(L) * T
  enddo
  print*,      ' N, leh=', N, leh
  write( 3, *) ' N, leh=', N, leh
  do L = 1, leh
     LE = L + N
     EC = E(LE) - CM1
     S  = H(L) * DE(LE)
     CM2 = CM2 + S*EC**2
     CM3 = CM3 + S*EC**3
     CM4 = CM4 + S*EC**4
  enddo
  XN  = XN * CM0
  if( Y .ne. 0 ) Y = CM1 / Y
  if( z .ne. 0 ) Z = CM2 / Z
  write( 3, 332 ) XN, CM0, CM1, CM2, CM3, CM4, Y, Z
332 format( 7X, 'Precision control,  ', 1P6E13.5,  &
         /24x, 'mean=', e12.5, '   variance=', e12.5)

END subroutine NORMAL

!-------------------------------------------------------------------------------
subroutine SHRINK

  implicit none

  integer n1, nu
  real F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4
  common / ut / N1, NU, F0, H0, U, um, EX, CZ0, CN, CM0, CM1, CM2, CM3, CM4

  integer N2, N2P, MIE, MIF, MIH, LEF, LEH, nume, lemx
  real d1
  common / IND / N2, N2P, MIE, MIF, MIH, LEF, LEH, D1, nume, lemx

  real f(1405), h(1405), E(1705), DI(1705), dE(1705), xn
  common / barray / f, h, E, DI, dE, xn

  real s
  integer n, l, m, lla, la, k, kk

  S = 0.0
  N = MIH - MIE
  do L = 1, leh
     lla = L
     K = L + N
     S = S + H(L)*DE(K)
     if( S .gt. EX ) GOTO 42
  enddo
42 M = lla - 1
  MIH = MIH + M
  S = 0.0
  LA = LEH + 1
  do K=1, leh
     L = LA - K
     KK = L + N
     S = S + H(L) * DE(KK)
     if( S .gt. EX) goto 44
  enddo
44 LEH = L - M

  do L = 1, leh
     K = L + M
     H(L) = H(K)
  enddo
  K = LEH + 1

  do L = K, nume
     H(L) = 0.0
  enddo

END subroutine SHRINK
