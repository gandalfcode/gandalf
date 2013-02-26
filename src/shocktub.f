      subroutine shocktub (RHOinL,
     *          RHOinR,
     *          UinL,
     *          UinR,
     *          PinL,
     *          PinR,
     *          xL,
     *          x0,
     *          xR,
     *          time,
     *          iMAX,
     *          gamma)
c     *********************************************************************
c     Exact Rieman solver (shock tube problem)
c     Reference: E.F. TORO, 1989, Intl. J. for Num. Methods in Fluid, 9
c
!f2py double precision intent(in) :: rhoinl, rhoinr, uinl, uinr
!f2py double precision intent(in) :: pinl, pinr, xl, x0, xr, time, gamma
!f2py integer intent(int) :: imax
      implicit double precision(a-h,l-z)
      dimension x(iMAX)
      dimension RHO(iMAX),U(iMAX),P(iMAX)
c
      common /gasconst/gm1,gp1,gmgp,gpgm,gfact1,gfact2

      gm1=gamma-1.d0
      gp1=gamma+1.d0
      gmgp=gm1/gp1
      gpgm=gp1/gm1
      gfact1=0.5*gm1/gamma
      gfact2=0.5*gp1/gamma
c     *********************************************************************
c
c..Initialize the grid
c
      call grid(x, iMAX, xL, xR)
c
c..Set up the initial conditions
c
      call icond(x,RHOl,RHOr,Ul,Ur,Pl,Pr,RHO,U,P, iMAX, RHOinL, RHOinR,
     *   UinL, UinR, PinL, PinR, xL, x0)
c
c..Calculate the star variables Pstar & Ustar
c   ******************
c   \ RHO*l | RHO*r /
c     \     Q*    /
c   Ql  \   |   /  Qr
c         \ | /
c---------------------------------
c
      call Qstar(RHOl,RHOr,Ul,Ur,Pl,Pr,Ustar,Pstar, gamma)
c
c..Solution of the Riemann problem
c
      Cl=sqrt(gamma*Pl/RHOl)
      Cr=sqrt(gamma*Pr/RHOr)
c
      do 100 i=1,iMAX
        xx=x(i)
c
        if(xx.lt.x0+Ustar*time) then
c     left of the contact discontinuty
c
c
          if(Pstar.le.Pl) then
c     left wave is a rarefaction
c
            PD=Pl/RHOl**gamma
            RHOs=(Pstar/PD)**(1./gamma)
            Cstar=sqrt(gamma*Pstar/RHOs)
c
          if(xx.lt.x0+(Ustar-Cstar)*time) then
c
            if(xx.lt.x0+(Ul-Cl)*time) then
c     left of the tail of the rarefaction
c
              RHO(i)=RHOl
              U(i)=Ul
              P(i)=Pl
            else
c     inside the rarefaction
c
             U(i)=2.0*((xx-x0)/time+Cl+0.5*gm1*Ul)/gp1
             Cs=Cl+0.5*gm1*(Ul-U(i))
             RHO(i)=(Cs*Cs/gamma/PD)**(1./gm1)
             P(i)=PD*RHO(i)**gamma
            endif
          else
c     right of the rarefaction
c
             RHO(i)=RHOs
             U(i)=Ustar
             P(i)=Pstar
            endif
          else
c     left wave is a shock
c
            sml2=1.0+gfact2*(Pstar/Pl-1.0)
            sml=-sqrt(sml2)
            Uls=Ul+Cl*sml
            if(xx.ge.x0+Uls*time) then
c     behind the shock
c
             RHO(i)=gp1*RHOl*sml2/(gm1*sml2+2.0)
             U(i)=Ustar
             P(i)=Pstar
            else
c     left of the shock
c
             RHO(i)=RHOl
             U(i)=Ul
             P(i)=Pl
            endif
          endif
c
        else
c     right of the contact discontinuty
c
          if(Pstar.le.Pr) then
c     right wave is a rarefaction
c
            if(xx.lt.x0+(Ur+Cr)*time) then
             PD=Pr/RHOr**gamma
             RHOs=(Pstar/PD)**(1./gamma)
             Cstar=sqrt(gamma*Pstar/RHOs)
c
            if(xx.lt.x0+(Ustar+Cstar)*time) then
c     left of the tail of the rarefaction
             RHO(i)=RHOs
             U(i)=Ustar
             P(i)=Pstar
           else
c     inside the rarefaction
c
             U(i)=2.0*((xx-x0)/time-Cr+0.5*gm1*Ur)/gp1
             Cs=Cr+0.5*gm1*(U(i)-Ur)
             RHO(i)=(Cs*Cs/gamma/PD)**(1./gm1)
             P(i)=PD*RHO(i)**gamma
            endif
          else
c     right of the rarefaction
c
            RHO(i)=RHOr
            U(i)=Ur
            P(i)=Pr
          endif
        else
c     right wave is a shock
c
            smr2=1.0+gfact2*(Pstar/Pr-1.0)
            smr=sqrt(smr2)
            Urs=Ur+Cr*smr
          if(xx.ge.x0+Urs*time) then
c     right of the shock
c
            RHO(i)=RHOr
            U(i)=Ur
            P(i)=Pr
          else
c     behind the shock
c
            RHO(i)=gp1*RHOr*smr2/(gm1*smr2+2.0)
            U(i)=Ustar
            P(i)=Pstar
          endif
        endif
c
                                                                               
      endif

c
  100 continue
c
c
      open(unit=50,file='sod.out')
      do i=1,iMAX
        write(50,1000) x(i),RHO(i),U(i),P(i),P(i)/RHO(i)/(gamma-1.)
c       write(50,1000) x(i),RHO(i),U(i)/dsqrt(gamma*P(i)/RHO(i)),
c    .                 P(i),P(i)/RHO(i)/(gamma-1.)
      enddo
      close(unit=50)
c
 1000 format(1x,5f18.8)
c

      end
c
      subroutine grid(x, iMAX, xR, xL)
c     *****************************************************************
c     Initialize the grid locations
      implicit double precision(a-h,l-z)
      dimension x(iMAX)
c     *****************************************************************
c
c..Equal spacing
c
      dx = (xR-xL)/(iMAX-1)
c
c..Load up the x-array
c
      do i=1,iMAX
        x(i) = xL+(i-1)*dx
      enddo
c
      return
      end
c
c
      subroutine icond(x,RHOl,RHOr,Ul,Ur,Pl,Pr,RHO,U,P, iMAX, RHOinL,
     *  RHOinR, UinL, UinR, PinL, PinR, xL, x0)
c     *********************************************************************
      implicit double precision(a-h,l-z)
      dimension x(iMAX)
      dimension RHO(iMAX),U(iMAX),P(iMAX)
c     *********************************************************************
c
      RHOl=RHOinL
      RHOr=RHOinR
      Ul=UinL
      Ur=UinR
      Pl=PinL
      Pr=PinR
      do i=1,iMAX
        xx=x(i)
      if(xx.ge.xL.and.xx.lt.x0) then
        RHO(i)=RHOl
        U(i)=Ul
        P(i)=Pl
      else
        RHO(i)=RHOr
        U(i)=Ur
        P(i)=Pr
      endif
c
      enddo
c
      return
      end
c
      subroutine Qstar(RHOl,RHOr,Ul,Ur,Pl,Pr,Ustar,Pstar, gamma)
c     *****************************************************************
c     Exact Riemann solver for the star region P & U
c     using the Newton-Raphson Method
      implicit double precision(a-h,l-z)
      common /gasconst/gm1,gp1,gmgp,gpgm,gfact1,gfact2
c     *****************************************************************
c
      delU=Ul-Ur
      Cl=sqrt(gamma*Pl/RHOl)
      Cr=sqrt(gamma*Pr/RHOr)
      tl=Cl/Pl**gfact1
      tr=Cr/Pr**gfact1
c..Initial guess for Pstar
      Pstar=((Cl+Cr+0.5*gm1*delU)/(tl+tr))**(1./gfact1)
c
c..start of iteration
      Ps0=Pstar
      do 10 i=0,100
        y=fpstar(Pl,RHOl,Pstar,gamma)+fpstar(Pr,RHOr,Pstar,gamma)+delU
        dy=dfpstar(Pl,RHOl,Pstar,gamma)+dfpstar(Pr,RHOr,Pstar,gamma)
        dPstar=-y/dy
        Pstar=Pstar+dPstar
        testps=abs(Pstar-Ps0)
        if(testps.le.1.d-6) go to 20
        if(Pstar.lt.1.d-6) Pstar=1.d-6
        Ps0=Pstar
   10 continue
c
   20 Ustar=0.5*(Ul+fpstar(Pl,RHOl,Pstar,gamma)
     >     +Ur-fpstar(Pr,RHOr,Pstar,gamma))
c
      return
      end
c
      function fpstar(P,RHO,Pstar, gamma)
c     *****************************************************************
      implicit double precision(a-h,l-z)
      common /gasconst/gm1,gp1,gmgp,gpgm,gfact1,gfact2
c     *****************************************************************
c
      H=Pstar/P
      C=sqrt(gamma*P/RHO)
      if(H.ge.1.) then
c..Shock relation
        fpstar=(1.d0-H)*sqrt(2.*P/gp1/RHO/(H+gmgp))
      else
c..Rarefaction fan relation
        fpstar=2.*C/gm1*(1.d0-H**gfact1)
      endif
c
      return
      end
c
      function dfpstar(P,RHO,Pstar, gamma)
c     *****************************************************************
      implicit double precision(a-h,l-z)
      common /gasconst/gm1,gp1,gmgp,gpgm,gfact1,gfact2
c     *****************************************************************
c
      H=Pstar/P
      C=sqrt(gamma*P/RHO)
      if(H.ge.1.) then
c..Shock relation
        f1=2./P/gp1/RHO
        f2=H+gmgp
        f3=1.d0+0.5*(1.d0-H)/f2
        dfpstar=-sqrt(f1/f2)*f3
c..Rarefaction fan relation
      else
        dfpstar=-C/gamma/P**gfact1/Pstar**gfact2
      endif
c
      return
      end
