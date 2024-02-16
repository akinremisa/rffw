        subroutine hrftn(rmodel, moddim, rayp, gaussalp, 
     1      delay, n, dt, x)
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: HRFTH96                                               c
c                                                                     c
c      COPYRIGHT 2001                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       20 DEC 2001 - create this program based on hspec96
c       21 NOV 2002 - modified KEVNM header as per 
c           Mejian An Dept Geophysics,
c            IAG, Sao Paulo, Br
c       11 APR 2008 - added more information in the help command
c       24 APR 2009 - added chftop to handle ray parameter in E format
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       19 APR 2011 - flipped sign of  Z R and RFTN for SV incident
c               to agree with expected particle motion
c       05 JUN 2013 - see note in excit for computation of Ur Uz
c               the forward synthetics now refer to a unit amplitude
c               incident P or SV wave. The change required multiplying
c               the previous amplitudes by a(mmax) or b(mmax), for
c               P or Sv incident, respectively. Also modified evgmat
c               so that the free surface effect for a near vertical 
c               incidence leads to a doubling of amplitude at the surface
c-----
c       This program creates a SAC file for the 
c           surface receiver function
c
c-----
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       earth model information
c-----
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer iunit, iiso, iflsph, idim, icnvel, ierr
        character mname*80, title*80
c-----
c       time series information
c-----
        integer NSAMP, NFREQ
        parameter (NSAMP=8192,NFREQ=4097)
        complex zdata(NSAMP)
        real x(NSAMP)
c
        complex Z
        complex ZFAC
c
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
c

        common/depref/refdep
c
        common/damp/alpha,ieqex
c-----
c       Adding the required variables to modify
c-----
         real*4 rmodel(moddim)
         character*40 kname
c-----
c       fref    R*4 - reference frequency for Causal Q
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c-----


c-----
c       command line arguments
c-----
        logical dop, dotwo
        real delay, dt, rayp, gaussalp
        integer n, iout

        logical ext

c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
c        call gcmdln(dop,dotwo,delay,dt,rayp,gaussalp,n,mname,iout)
c-----
c       fixed variables by IF
c-----
        IOUT = 1
        DOKJAR = .FALSE.
        ISHANK = .FALSE.
        DSTCOR = 0
        HNKARG = 3.0
        dop = .true.
        dotwo = .true.
c-----
c       get earth model information
c-----
c        if(mname.eq. ' ')call usage('Model not specified')  
c-----
c       get the earth model
c-----
c        inquire(file=mname,exist=ext)
c        if(.not. ext)call usage('Model file does not exist')
c        l = lgstr(mname)
c        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
c     1          idimen,icnvel,ierr,.false.)
c        if(ierr .lt. 0)stop
        mmax = (moddim+1)/4
        iunit =  0
        iiso = 0
        iflsph = 0
        idimen = 1
        icnvel = 0
c read the thickness
        do i=1,mmax-1
            d(i) = rmodel(mmax*3+i)
        end do
        do i=1,mmax
c read the VS
            b(i) = rmodel(i)
c read the VP
            a(i) = b(i) * rmodel(mmax+i)
            rho(i) = rmodel(mmax*2+i)
            qa(i) = 0.0
            qb(i) = 0.0
            etap(i) = 0.0 
            etas(i) = 0.0
            frefp(i) =  1.0
            frefs(i) = 1.0
c>>>>>>>>IF 20240216 -- Enable for debugging
c            write(*,"(A12,I6,4F10.5)") 'Layer #: ',i,d(i),
c     1      b(i),a(i),rho(i)
c<<<<<<<<IF 20240216 -- Enable for debugging
        end do
       
c-----
c       make sure that we use 1/Q
c-----
        do 3007 i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
 3007   continue
c-----
c       output control parameters
c-----
c>>>>>>>>IF 20240216 -- Enable for debugging
c        write(LOT,*)'dop          :',dop
c        write(LOT,*)'delay        :',delay
c        write(LOT,*)'dt           :',dt   
c        write(LOT,*)'rayp         :',rayp 
c        write(LOT,*)'gaussalp     :',gaussalp 
c        write(LOT,*)'n            :',n 
c        write(LOT,*)'Double Length:',dotwo
c        lm = lgstr(mname)
c        write(LOT,*)'Model        :',mname(1:lm)
c<<<<<<<<IF 20240216 -- Enable for debugging
c       ensure that the number of points is a power of 2
c-----
        call npow2(n)
        if(dotwo)then
            n = 2*n
        endif

        if(dokjar)then
        write(LOT,*)'Kjartansson Constant Q operator used'
        endif
c-----
c       generate a zero phase pulse, with a zero at the 
c           Nyquist Frequency
c-----
        delay = abs(delay)
c-----
c       define the alpha parameter for complex frequency
c-----
        alpha = 2.5/(n*dt)
c-----
c       note this will not work for symmetric pulses with 0 lag ???
c-----
        fac = 1.0
        dfac = exp(-alpha*dt)
        do 1000 i=1,n
            if(i.eq.2)then
                zdata(i) = cmplx(0.25/dt,0.0)*fac
            else if(i.eq.3)then
                zdata(i) = cmplx(0.50/dt,0.0)*fac
            else if(i.eq.4)then
                zdata(i) = cmplx(0.25/dt,0.0)*fac
            else
                zdata(i) = cmplx(0.0,0.0)
            endif
            fac = fac * dfac
 1000   continue
        call zfour(zdata,n,-1,dt,df)
c-----
c       now process in the frequency domain
c       include the desired time shift, but not that the 
c           pulse is already
c       shifted 2*dt units 
c-----
        n21 = n / 2 + 1
        do 4000 i=1,n21
                freq=(i-1)*df
            call excit(freq,dop,rayp,Z,iout)
            zdata(i) = zdata(i) * Z
            fac = - 6.2831853*freq*(delay - 2*dt)
            zdata(i) = zdata(i) * cmplx(cos(fac),sin(fac))
            if(i.gt.1)then
                zdata(n+2 -i ) = conjg(zdata(i))
            endif
 4000   continue
        call zfour(zdata,n,+1,dt,df)
c-----
c       undamp the time series
c-----
            fac = exp(-alpha*delay)
            dfac = exp(alpha*dt) 
            do 425 i = 1,n 
                zdata(i)= zdata(i) * fac
                fac = fac * dfac 
  425       continue 
c-----
c       now Gaussian filter
c-----
        call zfour(zdata,n,-1,dt,df)
        do 426 i=1,n21
            freq=(i-1)*df
            fac = (6.2831853*freq)/(2.0*gaussalp)
            if(fac.gt.25.0)then
                fac = 0.0
            else
                fac = exp( - fac * fac)
            endif
            zdata(i) = zdata(i) * fac
            if(i.gt.1)then
                zdata(n + 2 - i) = conjg(zdata(i))
            endif
  426   continue
c-----
c       The source pulse should have a zero at the Nyquist frequency
c-----
        zdata(n21) = cmplx(0.0,0.0)
        call zfour(zdata,n,+1,dt,df)
        do 427 i=1,n
            x(i) = real(zdata(i))
  427   continue
        return
        end

       subroutine zero(x,n)
c-----
c   zero a real array
c-----
c   x   R*4 - array to be zeroed
c   n   I*4 - number of points in array
c-----
       real x(n)
       integer i,n
       do i = 1,n
          x(i) = 0
       end do
c
       return
       end


        subroutine gcmdln(dop,dotwo,delay,dt,rayp,gaussalp,
     1      nsamp,mname,iout)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c       dop L   - .true. P-wave incident
c                 .false. SV-wave incident
c       dotwo   L   - .true. computations done with double length,
c                   but only single length is output
c       delay   R   - padding before the beginning of the 
c           receiver function
c       rayp    R   - ray parameter in sec/km
c       gaussalp R  - Gaussian filter parameter
c       nsamp   I   - number of samples (poer of 2)
c       mname   Ch* - earth model name - required
c       dokjar  L   - 
c       iout    I   - output 1 = rftn, 2 = z, 3 = r time series
c-----
        logical dop, dotwo
        real delay, dt, rayp, gaussalp
        integer nsamp, iout
        character mname*(*)

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        logical ishank,dokjar
        integer*4 dstcor
        real*4 hnkarg



        integer*4 mnmarg
        character*50 name

        IOUT = 1
        DOKJAR = .FALSE.
        ISHANK = .FALSE.
        DSTCOR = 0
        HNKARG = 3.0
        dop = .true.
        dotwo = .false.
        delay = 5.0
        dt = 1.0
        rayp = 0.05
        gaussalp = 1.0
        nsamp = 512
        mname = ' '
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-P' )then
                dop = .true.
            else if(name(1:2).eq.'-S' )then
                dop = .false.
            else if(name(1:2).eq.'-D' .and. name(1:3).ne.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')delay
            else if(name(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')dt
            else if(name(1:2).eq.'-N')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nsamp
            else if(name(1:5).eq.'-RAYP')then
                i = i + 1
                call mgtarg(i,name)
                call chtofp(name,rayp)
            else if(name(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')dt
            else if(name(1:4).eq.'-ALP')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')gaussalp
            else if(name(1:2).eq.'-M')then
                i = i + 1
                call mgtarg(i,mname)
            else if(name(1:2).eq.'-2')then
                dotwo = .true.
            else if(name(1:2).eq.'-z')then
                iout = 2
            else if(name(1:2).eq.'-r')then
                iout = 3
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'hrftn96 [-P] [-S] [-2] [-r] [-z] -RAYP p',
     2      ' -ALP alpha -DT dt -NSAMP nsamp',
     3      ' -M model'
        write(LER,*)
     1  '-P           (default true )    Incident P wave'
        write(LER,*)
     1  '-S           (default false)    Incident S wave'
        write(LER,*)
     1  '-RAYP p      (default 0.05 )    Ray parameter in sec/km'
        write(LER,*)
     1  '-DT dt       (default 1.0  )    Sample interval for synthetic'
        write(LER,*)
     1  '-NSAMP nsamp (default 512  )    Number samples for synthetic'
        write(LER,*)
     1  '-M   model   (default none )    Earth model name'
        write(LER,*)
     1  '-ALP alp     (default 1.0  )    Number samples for synthetic'
        write(LER,*)
     1  '     H(f) = exp( - (pi freq/alpha)**2) '
        write(LER,*)
     1  '     Filter corner ~ alpha/pi '
        write(LER,*)
     1  '-2           (default false)    Use 2x length internally'
        write(LER,*)
     1  '-r           (default false)    Output radial   time series'
        write(LER,*)
     1  '-z           (default false)    Output vertical time series'
        write(LER,*)
     1  '     -2  (default false) use double length FFT to'
        write(LER,*)
     1  '     avoid FFT wrap around in convolution '
        write(LER,*)
     1  '-D delay     (default 5 sec)    output delay sec before t=0'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        write(LER,*)
     1  ' SAC header values set by hrftn96'
        write(LER,*)
     1  '  B     :  delay'
        write(LER,*)
     1  '  USERO :  gwidth        KUSER0:  Rftn'
        write(LER,*)
     1  '  USER4 :  rayp (sec/km)'
        write(LER,*)
     1  '  USER5 :  fit in % (set at 100)'
        write(LER,*)
     1  '  KEVNM :  Rftn          KUSER1:  hrftn96'
        write(LER,*)
     1  'The program creates the file names hrftn96.sac'
        write(LER,*)
     1  'This is the receiver fucntion, Z or R trace according',
     2  ' to the command line flag'

        stop 
        end

        subroutine aten(om,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     2      frefp,frefs)
c-----
c       make velocities complex, using Futterman causality operator
c-----
        real*4 qa,qb,alpha,a,b
        complex*16 om,at,atna,atnb,xka,xkb
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        real*8 CDABS
        complex*16 CDLOG
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927d+00
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c           Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/pi
            gamb = atan(qb)/pi
            if(gama.le.0.0)then
                atna = cmplx(1.0,0.0)
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(b.gt.1.0e-04*a)then
                if(gamb.le.0.0)then
                    atnb = cmplx(1.0,0.0)
                else
                    fac = pi2*gamb
                    rfac = sin(fac)/cos(fac)
                    atnb = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1s)**dble(-gamb) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
                endif
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853d+00
            atna=dcmplx(1.0d+00,0.0d+00)
            atnb=dcmplx(1.0d+00,0.0d+00)
            if(qa.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1p)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1s)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        xka=om/(dble(a)*atna)
        if(b.le.1.0e-04*a)then
            iwat = 1
            xkb = dcmplx(0.0d+00,0.0d+00)
        else
            iwat = 0
            xkb=om/(dble(b)*atnb)
        endif
        return
        end

        subroutine cmult(e,ca,exa,exe)
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM EC where e(1x5) c(5x5)
        complex*16 ca(5,5)
        real*8 exa,exe,eval
        real *8 xnorm
        complex*16 e(5)
        complex*16 c, ee(5)

        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            c = dcmplx(0.0d+00,0.0d+00)
            do 1349 j=1,IUP
                c=c+ e(j) * ca(j,i)
 1349       continue
            ee(i)=c
 1350   continue
        exe = exe + exa
CRBH    call normc(ee,eval,xnorm)
CRBH    do 1351 i=1,IUP
CRBH        e(i) = ee(i)*xnorm
CRBH 1351   continue
CRBH    exe = exe + eval
        return
        end

        subroutine rcmult(y,c,exa,exe)
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM YC where y(5x5) c(5x5) RETURN Y
c-----
        complex*16 c(5,5)
        real*8 exa,exe,eval
        real *8 xnorm
        complex*16 y(5,5)
        complex*16 ztmp, ee(5,5)
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            do 1351 j=1,IUP
                ztmp = dcmplx(0.0d+00,0.0d+00)
                do 1349 k=1,IUP
                    ztmp=ztmp+ y(i,k) * c(k,j)
 1349           continue
                ee(i,j)=ztmp
 1351       continue
 1350   continue
        exe = exe + exa
CRBH    call rnormc(ee,eval,xnorm)
CRBH    do 1353 j=1,IUP
CRBH        do 1352 i=1,IUP
CRBH            y(i,j) = ee(i,j)*xnorm
CRBH 1352       continue
CRBH 1353   continue
CRBH    exe = exe + eval
        return
        end

        subroutine dmult(da,aa)
c-----
c       propagator up
c       FORM D = DA
c-----
        complex*16 aa(4,4)
        complex*16 sumd,ea(4,4),da(4,4)
        do 1360 i=1,4
            do 1361 j=1,4
                sumd = dcmplx(0.0d+00,0.0d+00)
                do 1362 jj=1,4
                    sumd=sumd+da(i,jj) * aa(jj,j)
 1362           continue
                ea(i,j)=sumd
 1361       continue
 1360   continue
        do 1363 j=1,4
            do 1364 i=1,4
                da(i,j)=ea(i,j)
 1364       continue
 1363   continue
        return
        end

        subroutine dnka(ca,wvno,wvno2,om2,gam,rho,iwat,w,x,cosp,ex)
        complex*16 ca(5,5)
        complex*16 wvno, wvno2, om2
        real*4 rho
        complex*16 gam
        complex*16 w,x,cosp
        real*8 ex
        common/ ovrflw / a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        real *8 a0
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz

        complex*16 gam2,gamm1,gamm2,a0c,xz2,wy2,temp
        real*8 dfac
        complex*16 cqww2, cqxw2, g1wy2, gxz2, g2wy2, g2xz2
        complex*16 gg1, a0cgg1
        complex*16 zrho, zrho2
c-----
c        A11     A12     A13    -A13     A15     A16
c        A21     A22     A23    -A23     A25     A15
c        A31     A32     A33    1-A33   -A23    -A13
c       -A31    -A32    1-A33    A33     A23     A13
c        A51     A52    -A32     A32     A22     A12
c        A61     A51    -A31     A31     A21     A11
c-----
c       this will be multipled on the left by the G matrix
c
c       [ G11   G12 G13 -G13    G15 G16 ]
c
c-----
c       or on the right by
c
c       [ H11   H21 H31 -H31    H51 H61  ] ^T
c-----
c       the number of multiplications can be reduced from 
c           36 to 25 if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c        2 A31   2 A32   2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the 
c           Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, 
c           we note that the old G14 = -old G14 = new G13
c-----
c-----
        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,5
                do 101 i=1,5
                    ca(i,j) = dcmplx(0.0d+00, 0.0d+00)
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            ca(3,3) = dfac
            ca(1,1) = cosp
            ca(5,5) = cosp
            ca(1,2) = - x/(zrho*om2)
            ca(2,1) = - zrho*w*om2
            ca(2,2) = cosp
            ca(4,4) = cosp
            ca(4,5) = ca(1,2)
            ca(5,4) = ca(2,1)
        else
c-----
c       elastic layer
c-----
            zrho2= dcmplx(dble(rho*rho),0.0d+00)
            gam2  = gam*gam
            gamm1 = gam-1.
            gamm2 = gamm1*gamm1
            cqww2 = cqw * wvno2
            cqxw2 = cqx / wvno2
            gg1 = gam*gamm1
            a0c  = dcmplx(2.0d+00,0.0d+00)*
     1          (dcmplx(a0,0.0d+00)-cpcq)
            xz2  = xz/wvno2
            gxz2 = gam*xz2
            g2xz2 = gam2 * xz2
            a0cgg1 = a0c*(gam+gamm1)
            wy2  = wy*wvno2
            g2wy2 = gamm2 * wy2
            g1wy2 = gamm1 * wy2

c-----
c       OK by symmetry
c----
            temp = a0c*gg1 + g2xz2 + g2wy2
            ca(3,3) = a0 + temp + temp
            ca(1,1) = cpcq-temp
            ca(1,2) = (-cqx + wvno2*cpy)/(zrho*om2)
            temp = dcmplx(0.5d+00,0.0d+00)*a0cgg1 + gxz2 + g1wy2
            ca(1,3) = wvno*temp/(zrho*om2)

            ca(1,4) = (-cqww2+cpz)/(zrho*om2)
            temp = wvno2*(a0c + wy2) + xz
            ca(1,5) = -temp/(zrho2*om2*om2)

            ca(2,1) = (-gamm2*cqw + gam2*cpz/wvno2)*zrho*om2
            ca(2,2) = cpcq
            ca(2,3) = (gamm1*cqww2 - gam*cpz)/wvno
            ca(2,4) = -wz
            ca(2,5)=ca(1,4)


            temp =dcmplx(0.5d+00,0.0d+00)*a0cgg1*gg1 
     1          + gam2*gxz2 + gamm2*g1wy2
            ca(3,1) = -dcmplx(2.0d+00,0.0d+00)*temp*zrho*om2/wvno
            ca(3,2) = -wvno*(gam*cqxw2 - gamm1*cpy)*
     1          dcmplx(2.0d+00,0.0d+00)

            ca(3,4)=-dcmplx(2.0d+00,00d+00)*ca(2,3)
            ca(3,5)=-dcmplx(2.0d+00,00d+00)*ca(1,3)

            ca(4,1) = (-gam2*cqxw2 + gamm2*cpy)*zrho*om2
            ca(4,2) = -xy
            ca(4,3)= -ca(3,2)/(dcmplx(2.0d+00,00d+00))
            ca(4,4)=ca(2,2)
            ca(4,5)=ca(1,2)

            temp = gamm2*(a0c*gam2 + g2wy2) + gam2*g2xz2
            ca(5,1) = -zrho2*om2*om2*temp/wvno2
            ca(5,2)=ca(4,1)
            ca(5,3)=-ca(3,1)/(dcmplx(2.0d+00,00d+00))
            ca(5,4)=ca(2,1)
            ca(5,5)=ca(1,1)
        endif
        return
        end

        subroutine hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,gam,gamm1,rho,
     1      iwat,ex,om2)
        complex*16 wvno,wvno2
        complex*16 aa(4,4),w,x,y,z,cosp,cosq,gam,gamm1
        complex*16 cpq, gcpq, zw2, gzw2, g1w, g1y, gx
        complex*16 zrho
        real*8 ex, dfac
        complex*16 om2
        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    aa(i,j) = dcmplx(0.0d+00,0.0d+00)
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            aa(1,1) = dfac
            aa(4,4) = dfac
            aa(2,2) = cosp
            aa(3,3) = cosp
            aa(2,3) = -x/(zrho*om2)
            aa(3,2) = - zrho*w*om2
        else
c-----
c       elastic layer
c-----
            cpq = cosp-cosq
            gcpq = gam*cpq
            zw2 = z/wvno2
            gzw2 = gam*zw2
            g1w = gamm1*w
            g1y = gamm1*y
            gx = gam*x
            aa(1,1)=   gcpq + cosq
                aa(1,3)= - wvno * cpq/(zrho*om2)
            aa(1,2)=   wvno*(-g1w+gzw2)
            aa(1,4)=   (wvno2*w-z)/(zrho*om2)
            aa(2,1)=   (gx - wvno2*g1y)/wvno
            aa(2,2)= - gcpq + cosp
            aa(2,3)=   (-x+wvno2*y)/(zrho*om2)
            aa(2,4)= - aa(1,3)
            aa(3,1)=   zrho*om2*gamm1*gcpq/wvno
            aa(3,2)=   zrho*om2*((-gamm1*g1w)+(gam*gzw2))
            aa(3,3)=   aa(2,2)
            aa(3,4)= - aa(1,2)
            aa(4,1)=   zrho*om2*(((gam*gx)/wvno2) - (gamm1*g1y))
            aa(4,2)= - aa(3,1)
            aa(4,3)= - aa(2,1)
            aa(4,4)=   aa(1,1)
        endif
        return
        end

        subroutine hskl(hl,cosql,yl,zl,h,iwat)
        complex*16 hl(2,2)
        complex*16 cosql,yl,zl,h
        integer iwat
c-----
c       cosql = cosh ( nu d )
c       zl    = nu sinh ( nu d )
c       yl    = sinh ( nu d ) / nu
c       h     = rho beta ^ 2
c-----
        if(iwat.eq.0)then   
            hl(1,1) = cosql
            hl(2,1) = zl*h
            hl(1,2) = yl/h
            hl(2,2) = cosql
        else
            hl(1,1) = dcmplx(1.0d+00,0.0d+00)
            hl(1,2) = dcmplx(0.0d+00,0.0d+00)
            hl(2,1) = dcmplx(0.0d+00,0.0d+00)
            hl(2,2) = dcmplx(1.0d+00,0.0d+00)
        endif
        return
        end 

        subroutine lmult(d11,d12,d21,d22,hl,iwat,exel,exb,icomp)
c-----
c       multiply SH matrix by a row vector on left
c-----
        complex*16 d11,d12,d21,d22,hl(2,2),e1,e2
        real*8 exel, exb
        logical icomp
c-----
c       fluid layer do nothing, just return, 
c           equivalent to multiplying by
c       identity matrix
c-----
        if(iwat.eq.0)then
c-----
c       elastic layer
c-----
            e1=d11
            e2=d12
c-----
c           a11 = cosql
c           a12 = yl
c           a21 = zl
c           a22 = cosql
c-----
            d11=e1*hl(1,1) + e2*hl(2,1)
            d12=e1*hl(1,2) + e2*hl(2,2)
            exel = exel + exb
            if(icomp)then
                e1=d21
                e2=d22
                d21=e1*hl(1,1) + e2*hl(2,1)
                d22=e1*hl(1,2) + e2*hl(2,2)
            endif
        endif
        return
        end

        subroutine normc(e,ex,xnorm)
        common/lwater/lfluid
        logical lfluid
        real*8 ex
        real *8 test,testt,x,y,fac,xnorm
        real*8 DREAL
        complex*16 e(*)
        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 2 i = 1,IUP
            if(dabs(dreal(e(i))).gt.testt) testt =dabs(dreal(e(i)))
            if(dabs(dimag(e(i))).gt.testt) testt =dabs(dimag(e(i)))
    2   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 1 i =1,IUP
            x=dreal(e(i))/testt
            y=dimag(e(i))/testt
            fac = x*x + y*y
            if(test.lt.fac) test = fac
    1   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine rnormc(e,ex,xnorm)
        common/lwater/lfluid
        logical lfluid
        real*8 ex
        real*8 DREAL
        real *8 test,testt,x,y,fac,xnorm
        complex*16 e(5,5)
        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 3 j=1,IUP
            do 2 i = 1,IUP
            if(dabs(dreal(e(i,j))).gt.testt) testt =dabs(dreal(e(i,j)))
            if(dabs(dimag(e(i,j))).gt.testt) testt =dabs(dimag(e(i,j)))
    2       continue
    3   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 4 j=1,IUP
            do 1 i =1,IUP
                x=dreal(e(i,j))/testt
                y=dimag(e(i,j))/testt
                fac = x*x + y*y
                if(test.lt.fac) test = fac
    1       continue
    4   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine evlmat(om,wvno,jbdrys,jbdryh,wvno2,om2)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/damp/alpha,ieqex
        complex*16 om,xka,xkb,ra,rb,gam,gamm1,p,q,w,x,y
        complex*16 cosq,z,cosp
        complex*16  ca(5,5), hl(2,2)
        complex*16 atna,atnb
        complex*16 yl,zl,cosql
        real *8 ex,exa,exb
        complex*16 wvno,wvno2,om2
        complex*16 h
        complex*16 aa(4,4)
        complex*16 zone
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 har(NL,4,4), dar(NL,5,5), hsr(2,5), gbr(2,5), 
     1      hal(NL,2,2), hsl(2,2), gbl(2,2)
        real*8 hex(NL), lex(NL), dex(NL), hexw(NL)
        common/hamat/har
        common/damat/dar
        common/hsrfr/hsr
        common/gbrfr/gbr
        common/hlmat/hal
        common/hsrfl/hsl
        common/gbrfl/gbl
        common/hexex/hex
        common/hexexw/hexw
        common/dexex/dex
        common/lexex/lex 
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        logical compute
        complex*16 CDSQRT
        zone  = dcmplx(1.0d+00,0.0d+00)
c-----
c       evaluate the G matrix 
c           gbr(1,x) is for normal stack
c           gbr(2,x) is for inverted stack
c-----
        call evalg(jbdryh,mmax,mmax-1,gbr,gbl,1,wvno,om,om2,wvno2)
        call evalg(jbdrys,1,   2,     gbr,gbl,2,wvno,om,om2,wvno2)
c-----
c       evaluate the H matrix
c           hsr(1,x) is for normal stack
c           hsr(2,x) is for inverted stack
c-----
        call evalh(jbdrys,1,   2,     hsr,hsl,1,wvno,om,om2,wvno2)
        call evalh(jbdryh,mmax,mmax-1,hsr,hsl,2,wvno,om,om2,wvno2)
c-----

c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 m = 1,mmax,1
c-----
c       first check to see if computations already done
c-----
            if(equald(m).and.m.lt.mmax)then
                compute = .false.
            else
                compute = .true.
            endif
            if(compute)then
                call aten(om,qa(m),qb(m),xka,xkb,
     1              alpha,a(m),b(m),atna,atnb,iwat,
     2              frefp(m),frefs(m))
                h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
                gam=dble(b(m))*(wvno/om)
                gam = gam * atnb
                gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
                gamm1 = gam - zone
                ra=CDSQRT(wvno2-xka*xka)
                rb=CDSQRT(wvno2-xkb*xkb)
                p=ra*dble(d(m))
                q=rb*dble(d(m))
                call var(p,q,ra,rb,w,x,y,z,cosp,cosq,
     1              ex,exa,exb,yl,zl,cosql,iwat)
                call dnka(ca,wvno, wvno2, om2,gam,rho(m),iwat,
     1              w,x,cosp,ex)
                call hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,gam,
     1              gamm1,rho(m),iwat,ex,om2)
                call hskl(hl,cosql,yl,zl,h,iwat)
            endif
            iwater(m) = iwat
            call copy5(ca,dar,m,0,dex,exa)
            call copy4(aa,har,m,0,hex,ex)
            call copy2(hl,hal,m,0,lex,exb)
 1340   continue
        return
        end

        subroutine evgmat(om, om2, wvno, wvno2, g, RA,RB)
        implicit none
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        real *8 ex,exa,exb
        complex*16 om,wvno,wvno2,om2
        complex*16 g(4,4)
        complex*16 zone
        integer i, j, m, iwat
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 CDSQRT
        zone  = dcmplx(1.0d+00,0.0d+00)
        m = mmax
            call aten(om,qa(m),qb(m),xka,xkb,
     1          alpha,a(m),b(m),atna,atnb,iwat,
     2          frefp(m),frefs(m))
            gam=dble(b(m))*(wvno/om)
            gam = gam * atnb
            gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
            gamm1 = gam - zone
            ra=CDSQRT(wvno2-xka*xka)
            rb=CDSQRT(wvno2-xkb*xkb)
            g(1,1) =   gam/wvno
            g(2,1) = - gamm1/rb
            g(3,1) =   g(1,1)
            g(4,1) = - g(2,1)
            g(1,2) = - gamm1 / ra
            g(2,2) =   g(1,1)
            g(3,2) =   gamm1 / ra
            g(4,2) =   g(1,2)
            g(1,3) = - 1.0d+00/(rho(m)*om2)
            g(2,3) =   wvno / ( rb * rho(m)*om2)
            g(3,3) =   g(1,3)
            g(4,3) = - g(2,3)
            g(1,4) =   wvno / ( ra * rho(m)*om2)
            g(2,4) =   g(1,3)
            g(3,4) = - wvno / ( ra * rho(m)*om2)
            g(4,4) =   g(2,4)
            do i=1,4
               do j=1,4
                  g(i,j) = g(i,j) * 0.5
               enddo
            enddo
        return
        end

        subroutine copy5(ca,dar,m,itofrm,dex,exa)
        parameter (NL=200)
        complex*16 dar(NL,5,5)
        complex*16 ca(5,5)
        integer itofrm
        real*8 dex(NL)
        real*8 exa
c-----
c       copy from ca to dar
c-----
        if(itofrm.eq.0)then
            do 100 j=1,5
                do 110 i=1,5
                    dar(m,i,j) = ca(i,j)
  110           continue
                dex(m) = exa
  100       continue
c-----
c       copy from dar to ca
c-----
        else
            do 200 j=1,5
                do 210 i=1,5
                    ca(i,j) = dar(m,i,j)
  210           continue
                exa = dex(m)
  200       continue
        endif
        return
        end

        subroutine copy4(aa,har,m,itofrm,hex,ex)
        parameter (NL=200)
        complex*16 har(NL,4,4)
        complex*16 aa(4,4)
        integer itofrm
        real*8 hex(NL)
        real*8 ex
c-----
c       copy from aa to har
c-----
        if(itofrm.eq.0)then
            do 100 j=1,4
                do 110 i=1,4
                    har(m,i,j) = aa(i,j)
  110           continue
                hex(m) = ex
  100       continue
c-----
c       copy from har to aa
c-----
        else
            do 200 j=1,4
                do 210 i=1,4
                    aa(i,j) = har(m,i,j)
  210           continue
                ex = hex(m)
  200       continue
        endif
        return
        end

        subroutine copy2(hl,hal,m,itofrm,lex,exb)
        parameter (NL=200)
        complex*16 hal(NL,2,2)
        complex*16 hl(2,2)
        integer itofrm
        real*8 lex(NL)
        real*8 exb
c-----
c       copy from hl to hal
c-----
        if(itofrm.eq.0)then
            do 100 j=1,2
                do 110 i=1,2
                    hal(m,i,j) = hl(i,j)
  110           continue
                lex(m) = exb
  100       continue
c-----
c       copy from hal to hl
c-----
        else
            do 200 j=1,2
                do 210 i=1,2
                    hl(i,j) = hal(m,i,j)
  210           continue
                exb = lex(m)
  200       continue
        endif
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,gbl,in,
     1      wvno,om,om2,wvno2)
        complex*16 gbr(2,5), gbl(2,2)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        complex*16 om,xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        complex*16 wvno,wvno2,om2
        complex*16 h
        common/lwater/lfluid
        logical lfluid
        complex*16 CDSQRT
c-----
c       set up halfspace conditions
c-----
        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
CRBH    WRITE(0,*)'RA :',RA
CRBH    WRITE(0,*)'RB :',RB
c-----
c       set up halfspace boundary conditions
c
c       jbdry   = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c
c-----
        if(jbdry.lt.0)then
c-----
c       RIGID - check properties of layer above
c-----
            if(b(m) .gt. 0.0)then
c-----
c               ELASTIC ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            else
c-----
c               FLUID ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(lfluid)then
                    gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,4) = dcmplx(1.0d+00,0.0d+00)
                endif
c-----
c               (pseudo SH)
c-----
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.0)then
c-----
c       HALFSPACE
c-----
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
c-----
c       multiply G of Herrmann 2001 by - rho^2 om^4 k^2 ra rb
c       should have no effect since it is in both numerator and
c       denominator -- however will not give the correct potential
c       coefficient -- so rethink?
c-----



                gbr(in,1)=dble(rho(m)*rho(m))*om2*om2*
     1              (-gam*gam*ra*rb+wvno2*gamm1*gamm1)
                gbr(in,2)=-dble(rho(m))*(wvno2*ra)*om2
                gbr(in,3)=-dble(rho(m))*(-gam*ra*rb+wvno2*gamm1)
     1              *om2*wvno
                gbr(in,4)=dble(rho(m))*(wvno2*rb)*om2
                gbr(in,5)=wvno2*(wvno2-ra*rb)
        gbr(in,1) = 0.25*gbr(in,1)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,2) = 0.25*gbr(in,2)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,3) = 0.25*gbr(in,3)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,4) = 0.25*gbr(in,4)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,5) = 0.25*gbr(in,5)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)

                gbl(in,1) =  dble(rho(m))*(dble(b(m)*b(m))
     1                  *atnb*atnb)*rb
                gbl(in,2) =  dcmplx(1.0d+00,0.0d+00)
            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(lfluid)then
                    gbr(in,1) = dble(0.5) / ra
                    gbr(in,2) = dcmplx(0.5d+00,0.0d+00)/
     1                  (-dble(rho(m))*om2)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dble(0.5*rho(m)*om2) / ra
                    gbr(in,5) = dcmplx(-0.5d+00,0.0d+00)
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.1)then
c-----
c       FREE - check properties of layer above
c-----
            if(b(m) .gt. 0.0)then
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
                
            else
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(lfluid)then
                    gbr(in,2) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        endif
CRBH        WRITE(0,*)'EVALG:'
CRBH        WRITE(0,*)'jbdry,in:',jbdry,in
CRBH        WRITE(0,*)'gbr:',gbr
CRBH        WRITE(0,*)'gbl:',gbl           
        return
        end

        subroutine evalh(jbdry,m,m1,hsr,hsl,in,wvno,om,om2,wvno2)
        complex*16 hsr(2,5), hsl(2,2)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        complex*16 om,xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        complex*16 wvno,wvno2,om2
        complex*16 h
        complex*16 CDSQRT
c-----
c       set up surface conditions
c-----
        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
c-----
c       do top surface condition
c
c           = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c-----
        if(iwat.eq.0)then
            if(jbdry.eq.0)then
c-----
c           ELASTIC ABOVE SOLID
c-----
                hsr(in,1) = (wvno2 - ra*rb)
                hsr(in,2) = dble(rho(m))*om2*rb
                hsr(in,3) = 2.0d+00*
     1              wvno*dble(rho(m))*om2*
     1              ( gamm1 - gam*ra*rb/wvno2)
                hsr(in,4) = -dble(rho(m))*om2*ra
                hsr(in,5) = dble(rho(m)*rho(m))*om2*om2
     1              *(gamm1*gamm1 - 
     1              gam*gam*ra*rb/wvno2)
                hsl(in,1) = dcmplx(0.5d+00,0.0d+00)
                hsl(in,2) = 0.5d+00*h*rb
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(1.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(0.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(1.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        else if(iwat.gt.0)then
            if(jbdry.eq.0)then
c-----
c           HALFSPACE ABOVE FLUID
c-----
                hsr(in,1) = ra
                hsr(in,2) = -dble(rho(m))*om2
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        endif
CRBH        WRITE(0,*)'EVALH:'
CRBH        WRITE(0,*)'jbdry,in:',jbdry,in
CRBH        WRITE(0,*)'hsr:',hsr
CRBH        WRITE(0,*)'hsl:',hsl            
        return
        end

        subroutine var(p,q,ra,rb,w,x,y,z,cosp,cosq,ex,
     1      exa,exl,yl,zl,cosql,iwat)
c     not modified for negative p,q
c     this assumes that real p and real q have same signs
        common/ovrflw/a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 p,q,ra,rb,w,x,y,z,cosp,cosq
        complex*16 yl,zl,cosql
        complex*16 eqp,eqm,epp,epm,sinp,sinq
        real *8 a0,pr,pi,qr,qi,fac,qmp,ex,exa,exl
c-----
c       form terms such as cos(p), sin(p), cos(q), sin(q)
c       and cos(p)*cos(q)
c
c       Introduce a factorization of exponentials to
c       make a pseudo floating point system
c
c       ex is the exponent in cosp
c       exl is the exponent in cosq for SH
c       exa is the exponent in cosp*cosq
c-----
        real*8 DREAL
      ex=0.0d+00
      exl = 0.0d+00
      a0=0.0d+00
      pr=dreal(p)
      pi=dimag(p)
      epp=dcmplx(dcos(pi),dsin(pi))/2.
      epm=dconjg(epp)
      ex=pr
      fac=0.0
      if(pr.lt.15.) fac=dexp(-2.*pr)
      cosp=epp + fac*epm
      sinp=epp - fac*epm
      w=sinp/ra
      x=ra*sinp
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            a0 = 1.0d+00
            exa = ex
            cosq = 1.0d+00
            y = 0.0d+00
            z = 0.0d+00
            cosql = 1.0d+00
            yl = 0.0d+00
            zl = 0.0d+00
            exl = 0.0d+00
        else
c-----
c       elastic layer
c-----
            qr=dreal(q)
            qi=dimag(q)
            eqp=dcmplx(dcos(qi),dsin(qi))/2.
            eqm=dconjg(eqp)
            exl=qr
            fac=0.0d+00
            if(qr.lt.15.) fac=dexp(-2.*qr)
            cosql=eqp + fac*eqm
            sinq=eqp - fac*eqm
            yl=sinq/rb
            zl=rb*sinq
c-----
c       form factors for compound P-SV matrix
c-----
            exa=pr + qr
            cpcq=cosp*cosql
            cpy=cosp*yl
            cpz=cosp*zl
            cqw=cosql*w
            cqx=cosql*x
            xy=x*yl
            xz=x*zl
            wy=w*yl
            wz=w*zl
            fac=0.0d+00
            qmp=qr-pr
            if(qmp.gt.-40.) fac=dexp(qmp)
            cosq=cosql*fac
            y=fac*yl
            z=fac*zl
            fac=0.0d+00
            if(exa.lt.60.) a0=dexp(-exa)
        endif
        return
        end

        subroutine dosph(mmax,v,h,ipsvsh)
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c           Fast surface wave and free
c       mode computations, in  
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c-----
c       mmax    I*4 number of layers
c       v() R*4 array of velocities
c       h() R*4 array of layer thicknesses
c       ipsvsh  I*4     1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c-----
        integer*4 mmax
        real*4 v(mmax), h(mmax)
        integer*4 ipsvsh

        double precision z0,z1,r0,r1,dr,ar,tmp

        if(ipsvsh.eq.3)then
            ifunc = 1
        else
            ifunc = 2
        endif

        ar=6370.0d0
        dr=0.0d0
        r0=ar
        h(mmax)=1.0
        do 10 i=1,mmax
C           dr=dr+dble(d(i))
            dr=dr+dble(h(i))
            r1=ar-dr
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
C           d(i)=z1-z0
            h(i)=z1-z0
            tmp=(ar/r1-ar/r0)*ar/(z1-z0)
            v(i)=v(i)*tmp
C           a(i)=a(i)*tmp
C           b(i)=b(i)*tmp
            if(ifunc.eq.1)then
                tmp=(r0*(r0/ar)**4-r1*(r1/ar)**4)
     1              /(dble(h(i))*5.0d0)
C     1             /(dble(d(i))*5.0d0)
            else if(ifunc.eq.2)then
                tmp=(r0*(r0/ar)-r1*(r1/ar))
     1              /(dble(h(i))*2.0d0)
C     1             /(dble(d(i))*2.0d0)
            endif
C           rho(i)=rho(i)*tmp
            r0 = r1
   10   continue
        h(mmax)=0.0
        return
        end

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        if(npts.ge.nsamp)return
        npts = 2*npts
        go to 1000
        end

        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end

        subroutine equlck()
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
c-----
c       To avoid repeated computation, 
c           check to see if neighboring layers
c       are identical, once for going up and another for going down
c-----
c       First check top down
c-----
        do 100 m=1,mmax
            if(m.eq.1)then
                equald(m) = .false.
            else if(m.gt.1
     1          .and. a(m).eq.a(m-1) 
     2          .and. b(m).eq.b(m-1)
     3          .and. d(m).eq.d(m-1) 
     4          .and. rho(m).eq.rho(m-1)
     5          .and. qa(m).eq.qa(m-1)
     6          .and. qb(m).eq.qb(m-1) )then
                equald(m) = .true.
            else
                equald(m) = .false.
            endif
  100   continue
c-----
c       check bottom up
c-----
        do 200 m=1,mmax
            if(m.eq.mmax)then
                equalu(m) = .false.
            else if(m.lt.mmax
     1          .and. a(m).eq.a(m+1) 
     2          .and. b(m).eq.b(m+1)
     3          .and. d(m).eq.d(m+1) 
     4          .and. rho(m).eq.rho(m+1)
     5          .and. qa(m).eq.qa(m+1)
     6          .and. qb(m).eq.qb(m+1) )then
                equalu(m) = .true.
            else
                equalu(m) = .false.
            endif
  200   continue
        return
        end

        subroutine copy(da,aa)
c-----
c       copy da array to aa
c-----
        implicit none
        complex*16 da(4,4), aa(4,4)
        integer i, j
        
        do 1000 j=1,4
            do 1100 i=1,4
                aa(i,j) = da(i,j)
 1100       continue
 1000   continue
        return
        end

        subroutine excit(freq,dop,rayp,Z,iout)
        implicit none
        real freq, omega, rayp
        logical dop
        complex Z
        integer iout

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
c-----
c       matrix components in layers and boundaries saved
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
            complex*16 har(NL,4,4)
        common/damat/dar
            complex*16 dar(NL,5,5)
        common/hsrfr/hsr
            complex*16 hsr(2,5)
        common/gbrfr/gbr
            complex*16 gbr(2,5)
        common/hlmat/hal
            complex*16 hal(NL,2,2)
        common/hsrfl/hsl
            complex*16 hsl(2,2)
        common/gbrfl/gbl
            complex*16 gbl(2,2)
        common/hexex/hex
            real*8 hex(NL)
        common/hexexw/hexw
            real*8 hexw(NL)
        common/dexex/dex
            real*8 dex(NL)
        common/lexex/lex 
            real*8 lex(NL)
        common/water/iwater, iwats, iwatb
            integer iwater(NL), iwats(2), iwatb(2)
        common/updnsm/equalu, equald
            logical equalu(NL), equald(NL)
         logical compute

c-----
c       internal variables
c-----
        complex*16 om, om2, wvno, wvno2
        complex*16 aa(4,4), da(4,4), g(4,4)
        real*8 ex, exe, exa, exl
        complex*16   cd(5),e(5),fr, ca(5,5)
        
        integer m, i, j
        integer jbdrys, jbdryh

        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        integer iwat

        omega = 6.2831853*freq
        om =  dcmplx(dble(omega), dble(-alpha))
        om2 = om * om
        wvno  = dcmplx(dble(rayp),0.0d+00)*om
        wvno2 = wvno * wvno

        exe = 0.0d+00
        exa = 0.0d+00
        exl = 0.0d+00
        ex  = 0.0d+00
c-----
c       First evaluate all the Haskell matrices
c       Surface is free  , jbdrys = 1
c       Base is halfspace, jbdryh = 0
c-----
        jbdrys = 1
        jbdryh = 0
c-----
c       Check the model to avoid excessive computation
c-----
        call equlck()
c-----
c       evaluate the matrices
c-----  
        call evlmat(om,wvno,jbdrys,jbdryh,wvno2,om2)
C       call evemat(om, om2, wvno, wvno2, g, RA,RB)
        call evgmat(om, om2, wvno, wvno2, g, RA,RB)
c       initialize the halfspace
c-----
        do 100 i=1,5
            e(i) = gbr(1,i)
            cd(i) = e(i)
  100   continue
c-----
c       now compute the Haskell Product from top down
c-----
        do 1000 j=1,4
            do 1100 i=1,4
                aa(i,j) = dcmplx(0.0d+00, 0.0d+00)
 1100       continue
            aa(j,j) = dcmplx(1.0d+00,0.0d+00)
 1000   continue
                
        do 1340 m=1,mmax-1
c-----
c           get da(m) matrix
c-----
            call copy4(da,har,m,1,hex,ex)
            exl = exl + ex
c-----
c           form da * aa
c-----
            call dmult(da,aa)
c-----
c           copy da to aa
c-----
            call copy(da,aa)
 1340   continue
        do 1341 m=mmax-1,1,-1
c-----
c       handle the compound matrices
c-----
            call copy5(ca,dar,m,1,dex,exa)
            call cmult(e,ca,exa,exe)
 1341   continue
c-----
c       put in the halfspace condition
c-----
        e(1) = e(1)*hsr(1,1) + e(2)*hsr(1,2) + e(3)*hsr(1,3)
     1      + e(4)*hsr(1,4) + e(5)*hsr(1,5)
            call dmult(g ,aa)
c-----
c       compute the receiver function
c-----
c       NOTES:
c         With respect to theory the vertical motions
c         have been multiplied by -1 to make Uz positive upward
c
c         To make the Ur and Uz for incident refer to a unit
c         incident P-wave amplitude, the value must be
c         multiplied by a(mmax)*atna, and
c         to make the Ur and Uz for incident refer to a unit
c         incident SV-wave amplitude, the value must be
c         multiplied by b(mmax)*atna
c         As of 06 JUN 2013 this has been introduced
c-----
            m = mmax
            call aten(om,qa(m),qb(m),xka,xkb,
     1          alpha,a(m),b(m),atna,atnb,iwat,
     2          frefp(m),frefs(m))
        if(dop)then
            if(iout.eq.1)then
c-----
c           Ur/Uz
c-----
            Z = - dcmplx(0.0d+00, 1.0d+00) * g(2,2) / g(2,1)
            else if(iout.eq.2)then
c-----
c           Uz
c-----
            Z = -  g(2,1)*exp(-exl)/(g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            Z = Z * a(mmax)*atna
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            Z =   dcmplx(0.0d+00, 1.0d+00) * g(2,2)*exp(-exl)/
     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            Z = Z * a(mmax)*atna
            endif
        else
            if(iout.eq.1)then
c-----
c           Uz/Ur
c-----
            Z =   dcmplx(0.0d+00, 1.0d+00) * g(1,1) / g(1,2)
            else if(iout.eq.2)then
c-----
c           Uz
c-----
            Z =   dcmplx(0.0d+00, 1.0d+00) * g(1,1)*exp(-exl)/
     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            Z = Z * b(mmax)*atnb
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            Z =                             g(1,2)*exp(-exl)/
     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            Z = Z * b(mmax)*atnb
            endif
        endif
        return
        end

        subroutine evemat(om, om2, wvno, wvno2, g, RA,RB)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        real *8 ex,exa,exb
        complex*16 om,wvno,wvno2,om2
        complex*16 g(4,4)
        complex*16 zone
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 CDSQRT
        zone  = dcmplx(1.0d+00,0.0d+00)
        m = mmax
            call aten(om,qa(m),qb(m),xka,xkb,
     1          alpha,a(m),b(m),atna,atnb,iwat,
     2          frefp(m),frefs(m))
            h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
            gam=dble(b(m))*(wvno/om)
            gam = gam * atnb
            gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
            gamm1 = gam - zone
            ra=CDSQRT(wvno2-xka*xka)
            rb=CDSQRT(wvno2-xkb*xkb)
            g(1,1) =   wvno
            g(2,1) =   ra
            g(3,1) =   rho(m)*om2*gamm1
            g(4,1) =   rho(m)*om2*gam*ra/wvno
            g(1,2) =   rb
            g(2,2) =   wvno
            g(3,2) =   rho(m)*om2*gam*rb/wvno
            g(4,2) =   rho(m)*om2*gamm1
            g(1,3) =   g(1,1)
            g(2,3) = - g(2,1)
            g(3,3) =   g(3,1)
            g(4,3) = - g(4,1)
            g(1,4) = - g(1,2)
            g(2,4) =   g(2,2)
            g(3,4) = - g(3,2)
            g(4,4) =   g(4,2)
        return
        end

        subroutine getmod(rlun,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,listmd)
c-----
c       HISTORY
c
c       09 08 2000  gave ierr an initial default value for g77
c       01 13 2001  put in close(lun) if file is not model file
c       03 MAY 2002     Modify to permit read from standard input
c       06 JUL 2005 moved inquire to permit use of STDIN
c
c-----
c       General purpose model input
c       This model specification is designed to be as 
c           general as possible
c
c       Input lines
c       Line 01: MODEL
c       Line 02: Model Name
c       Line 03: ISOTROPIC or ANISOTROPIC or 
c           TRANSVERSELY ANISOTROPIC
c       Line 04: Model Units, First character 
c           is length (k for kilometer
c           second is mass (g for gm/cc), third is time (s for time)
c       Line 05: FLAT EARTH or SPHERICAL EARTH
c       Line 06: 1-D, 2-D or 3-D
c       Line 07: CONSTANT VELOCITY
c       Line 08: open for future use
c       Line 09: open for future use
c       Line 10: open for future use
c       Line 11: open for future use
c       Lines 12-end:   These are specific to the model
c           For ISOTROPIC the entries are
c           Layer Thickness, P-velocity, S-velocity, Density, Qp, Qs,
c           Eta-P, Eta S (Eta is frequency dependence), 
c           FreqRefP, FreqRefP
c-----
cMODEL
cTEST MODEL.01
cISOTROPIC
cKGS
cFLAT EARTH
c1-D
cCONSTANT VELOCITY
cLINE08
cLINE09
cLINE10
cLINE11
c H  VP  VS   RHO   QP  QS   ETAP   ETAS REFP  REFS
c1.0    5.0 3.0 2.5 0.0 0.0 0.0 0.0 1.0 1.0
c2.0    5.1 3.1 2.6 0.0 0.0 0.0 0.0 1.0 1.0
c7.0    6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
c10.0   6.5 3.8 2.9 0.0 0.0 0.0 0.0 1.0 1.0
c20.0   7.0 4.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0
c40.0   8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0
c-----
c-----
c       rlun    I*4 - logical unit for reading model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name - if this is stdin or 
c           STDIN just read
c                 from standard input
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       title   C*(*)   - title of the model file
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       idimen  I*4 - 1 1-D
c               - 2 2-D
c               - 3 3-D
c       icnvel  I*4 - 0 constant velocity
c                 1 variable velocity
c       ierr    I*4 - 0 model file correctly read in
c               - -1 file does not exist
c               - -2 file is not a model file
c                 -3 error in the model file
c       listmd  L   - .true. list the model
c------

        implicit none
        character mname*(*), title*(*)
        integer rlun
        integer*4 mmax, iunit, iiso, iflsph, idimen, icnvel
        integer*4 ierr
        character string*80
        logical listmd
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        logical ext
        character ftype*80
        integer lun, j, i, irefdp

c-----
c       test to see if the file exists
c-----
        ierr = 0
c-----
c       test for input
c-----
        if(MNAME(1:5).eq.'stdin' .or. mname(1:5).eq.'STDIN')then
c-----
c           do not open anything, use standard output
c-----
            lun = LIN
        else
            lun = rlun
            inquire(file=mname,exist=ext)
            if(.not.ext)then
                ierr = -1
                write(LER,*)'Model file does not exist'
                return
            endif
c-----
c           open the file
c-----
            open(lun,file=mname,status='old',form='formatted',
     1          access='sequential')
            rewind lun
        endif
c-----
c       verify the file type
c-----
c-----
c       LINE 01
c-----
        read(lun,'(a)')ftype
        if(ftype(1:5).ne.'model' .and. ftype(1:5).ne.'MODEL')then
            ierr = -2
            write(LER,*)'Model file is not in model format'
            close(lun)
            return
        endif
c-----
c       LINE 02
c-----
        read(lun,'(a)')title
c-----
c       LINE 03
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'ISO' .or. string(1:3).eq.'iso')then
            iiso = 0
        else if(string(1:3).eq.'TRA' .or. string(1:3).eq.'tra')then
            iiso = 1
        else if(string(1:3).eq.'ANI' .or. string(1:3).eq.'ani')then
            iiso = 2
        endif
c-----
c       LINE 04
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'KGS' .or. string(1:3).eq.'kgs')then
            iunit = 0
        endif
c-----
c       LINE 05
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'FLA' .or. string(1:3).eq.'fla')then
            iflsph = 0
        else if(string(1:3).eq.'SPH' .or. string(1:3).eq.'sph')then
            iflsph = 1
        endif
c-----
c       LINE 06
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'1-d' .or. string(1:3).eq.'1-D')then
            idimen = 1
        else if(string(1:3).eq.'2-d' .or. string(1:3).eq.'2-D')then
            idimen = 2
        else if(string(1:3).eq.'3-d' .or. string(1:3).eq.'3-D')then
            idimen = 3
        endif
c-----
c       LINE 07
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'CON' .or. string(1:3).eq.'con')then
            icnvel = 0
        else if(string(1:3).eq.'VAR' .or. string(1:3).eq.'var')then
            icnvel = 1
        endif
c-----
c       get lines 8 through 11
c-----
        do 900 i=8,11
            read(lun,'(a)')string
  900   continue
c-----
c       get model specifically for 1-D flat isotropic
c-----
c-----
c       get comment line
c-----
        read(lun,'(a)')string
        mmax = 0
        refdep = 0.0
        irefdp = 0
        if(iiso.eq.0)then
 1000       continue
            j = mmax +1
                read(lun,*,err=9000,end=9000)d(j),a(j),b(j),
     1              rho(j),qa(j),qb(j),etap(j),etas(j),
     2              frefp(j),frefs(j)
                if(d(j).lt.0.0)then
                    d(j) = -d(j)
                    refdep = refdep + d(j)
                    irefdp = j
                endif
            mmax = j
            go to 1000
 9000       continue
        endif
    1   format(' LAYER             H      P-VEL     S-VEL   DENSITY  ')
    2   format(' ',i5,5x,4f10.3)
    3   format(' ','-SURFACE ','- - - - - ','- - - - - ',
     1      '- - - - - ','- - - - - -')
        if(mmax.gt.0)then
            if(listmd)then
            ierr = 0
            write(LOT,1)
            do 2000 i=1,mmax
                write(LOT,2)
     1              i,d(i),a(i),b(i),rho(i)
                if(i.eq.irefdp)write(LOT,3)
 2000       continue
            endif
        else 
            ierr = -3
            write(LER,*)'Error in model file'
        endif
        if(lun.ne.LIN)close (lun)
        return
        end

        subroutine chtofp(str,fout)
c------
c       routine to convert string to floating point
c       The E format is accepted as well as free form
c       input
c
c       This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c       If the string str contains an E or e, then
c       we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c       read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.0)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        function mnmarg()
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MNMARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        implicit none
        integer iargc
        integer mnmarg
            mnmarg = iargc()
        return
        end

        subroutine mgtarg(i,name)
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MGTARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c       return the i'th command line argument
c
c       This version works on SUN, IBM RS6000
c-----
        implicit none
        integer i
        character name*(*)
            call getarg(i,name)
        return
        end

      subroutine mchdep()
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: mchdep                                                c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        implicit none
      return
      end

        function lgstr(str)
c-----
c       function to find the length of a string
c       this will only be used with file system path names
c       thus the first blank 
c       indicates the end of the string
c-----
        implicit none
        character*(*) str
        integer*4 lgstr
        integer n, i
        n = len(str)
        lgstr = 1
        do 1000 i=n,1,-1
            lgstr = i
            if(str(i:i).ne.' ')goto 100
 1000   continue
  100   continue
        return
        end
