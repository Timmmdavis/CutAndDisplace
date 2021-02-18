function [Stress]=LDstressLayeredElastic(Dn,Ds,xe,ye,Beta,a,x,y,mu1,mu2,nu1,nu2)
% LDstressLayeredElastic: LineDislocationInducedStressLayeredElastic,
%               Computes the influence (stress) of a single planar
%               line crack shearing and/or opening on the surronding
%               points. The medium is a stratified bonded elastic
%               (interface along y=0) where the upper layer is defined by
%               constants (nu1 and mu1) and the lower (nu2 and mu2).
%
%               Formulas from the following two publications:
%                 The tensile dislocation problem in a layered elastic
%                 medium. Maurizio Bonafede, Eleonora Rivalta. Geophysical
%                 Journal International, Volume 136, Issue 2, February 1999,
%                 Pages 341–356. 
%                 The edge dislocation problem in a layered elastic medium.
%                 Rivalta Eleonora, Mangiavillano Walter, Bonafede
%                 Maurizio. Geophysical Journal International, Volume 149, 
%                 Issue 2, May 2002, Pages 508–523
%
% Arguments: (input)
%      x & y  - The observation points locations in real Cartesian coords. 
%
%     xe & ye - The element midpoint location in real coords. 
%
%       a     - The elements half length
%
%       Beta  - The angle of the element away from the x-axis (radians).
%               When the normal in the -y axis down this is 0. In degrees
%               when the normal points east this is 90, west -90 and north
%               180.
%
%     Dn,Ds   - The defined displacement of each element.(normal and shear)
%               Dn+ is opening, Ds+ is left lateral shearing. 
%
%       nu    - The Poisson's ratio: For elastic 1 (upper) and elastic 2 (lower)
%
%       E     - The Young's modulus: For elastic 1 (upper) and elastic 2 (lower)
%
% Arguments: (output)
% Stress - Is the stress caused by the movement of the dislocation at the 
%               observataion points. [Sxx,Syy,Sxy].
%
% Dn=0;
% Ds=1;
% xe=0;
% ye=0;
% Beta=deg2rad(0);
% a=1;
% x=linspace(-2*a,2*a,20);
% y=linspace(-2*a,2*a,20);
% [x,y]=meshgrid(x,y);
% [dimx,dimy] = size(x); 
% nu1=0.25;
% nu2=0.25;
% mu1=1;
% mu2=1;
% %Layered case
% [Stress]=LDstressLayeredElastic(Dn,Ds,xe,ye,Beta,a,x,y,mu1,mu2,nu1,nu2);
% [Sxx,Syy,Sxy ] = ExtractCols( Stress );
% [SxxLayered,SyyLayered,SxyLayered]=ReshapeData2d( dimx,dimy,Sxx,Syy,Sxy);
% DrawContourFPlots2d(x,y,[],SxxLayered,SyyLayered,SxyLayered);
% %Homogenous case
% nu=nu1;
% mu=mu1;
% [ k,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );
% [Stress] = LDstressFS(x(:),y(:),xe,ye,a,Beta,Ds,Dn,nu,E);
% [Sxx,Syy,Sxy ] = ExtractCols( Stress );
% [SxxHomogenous,SyyHomogenous,SxyHomogenous]=ReshapeData2d( dimx,dimy,Sxx,Syy,Sxy);
% DrawContourFPlots2d(x,y,[],SxxHomogenous,SyyHomogenous,SxyHomogenous);
% 
%  Author: Tim Davis/Francesco Maccaferri
%  Copyright 2021, Tim Davis, Potsdam University
%  Modified from of Francesco's 2D fluid-filled fracture propagation scripts:
%  Modified so the notation is the same as in CutAndDisplace

a=a*2; %To make length
Beta=-Beta;%To rotate same as LD convention
delta=Beta+pi/2;


[~,~,~,~,c1,c2,d2]=SET_EPAR(mu1,mu2,nu1,nu2);

Ih1=zeros(4);  
Ih2=zeros(4);  
dI=zeros(4);  
dIfill=zeros(4);  

x=x(:);
y=y(:);
Sxxfill=zeros(size(x));
Syyfill=zeros(size(x));
Sxyfill=zeros(size(x));
for i=1:length(x)      
    Sxx = sxxUt(xe,ye,delta,a,x(i),y(i),mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*Dn + sxxUb(xe,ye,delta,a,x(i),y(i),mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*Ds;
    Syy = szzUt(xe,ye,delta,a,x(i),y(i),mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*Dn + szzUb(xe,ye,delta,a,x(i),y(i),mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*Ds;
    Sxy = sxzUt(xe,ye,delta,a,x(i),y(i),mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*Dn + sxzUb(xe,ye,delta,a,x(i),y(i),mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*Ds;
    Sxxfill(i)=Sxx;
    Syyfill(i)=Syy;
    Sxyfill(i)=Sxy;
end
Stress=[Sxxfill,Syyfill,Sxyfill];
      
end



function [Up,Us,Tp,Ts,c1,c2,d2]=SET_EPAR(mu1,mu2,nu1,nu2)
% Arguments: (input)
%      x & y  - The observation points locations in real Cartesian coords. 
%
%     xe & ye - The element midpoint location in real coords. 
%
%       a     - The elements half length
%
%       Beta  - The angle of the element away from the x-axis (radians).
%               When the normal in the -y axis down this is 0. In degrees
%               when the normal points east this is 90, west -90 and north
%               180.
%
%     Dn,Ds   - The defined displacement of each element.(normal and shear)
%               Dn+ is opening, Ds+ is left lateral shearing. 
%
%       nu    - The Poisson's ratio: For elastic 1 (upper) and elastic 2 (lower)
%
%       E     - The Young's modulus: For elastic 1 (upper) and elastic 2 (lower)
%
% Arguments: (output)
% Disp        - Is the displacement caused by the movement of the
%             dislocation at the observataion points. [Ux,Uy].

      pi = 4.D0*atan(1.D0);

      if mu2==0 %%Then % homogeneous half space
      
        nu2 = 0.D0;

        gamm1 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu1);
        gamm2 = 1.D0/(4.D0*pi);

        kap1  = 1.D0/(2.D0*pi*mu1);
        kap2  = 1.D0/(6.D0*pi*mu1);

        c1    = 0.D0;
        c2    = 1.D0/(2.D0*pi) * mu1/(1.D0-nu1);
        d2    = 0.D0;
      
      else
      
        if ((nu1~=nu2) || (mu1~=mu2)) %%Then % bounded medium

          a1    = (3.D0-4.D0*nu1)/(4.D0*mu1*mu1);
          a2    = (3.D0-4.D0*nu2)/(4.D0*mu2*mu2);

          delt  = 1.D0/(2.D0*pi) * (mu2 /(1.D0-nu2) - mu1 /(1.D0-nu1));

          gamm1 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu1);
          gamm2 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu2);
          gamm  = gamm2 - gamm1;

          kap1  = 1.D0/(2.D0*pi) * 1.D0/(mu1 + (3.D0 - 4.D0*nu1)*mu2);
          kap2  = 1.D0/(2.D0*pi) * 1.D0/(mu2 + (3.D0 - 4.D0*nu2)*mu1);

          cp    = (1.D0+(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2);
          cm    = (1.D0-(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2);

          d     = (1.D0-2.D0*nu2)/(2.D0*mu2)-(1.D0-2.D0*nu1)/(2.D0*mu1);
          e     = (1.D0-     nu2)/      mu2 +(1.D0-     nu1)/      mu1;

          c1    = 1.D0/(e*e-d*d) *( delt*(a1+cp) + gamm*d);
          c2    = 1.D0/(e*e-d*d) *(-delt*(a2+cp) + gamm*d);
          d2    = 1.D0/(e*e-d*d) *( delt*cm      - gamm*e);

        else %% homogeneous medium

          a1    = (3.D0-4.D0*nu1)/(4.D0*mu1*mu1);
          a2    = a1;

          delt  = 0.D0;

          gamm1 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu1);
          gamm2 = gamm1;
          gamm  = 0.D0;
          
          kap1  = 1.D0/(2.D0*pi) * 1.D0/(mu1 + (3.D0 - 4.D0*nu1)*mu2);
          kap2  = kap1;

          cp    = (1.D0+(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2);
          cm    = (1.D0-(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2);

          d     = 0.D0;
          e     = 2.D0*(1.D0-nu1)/mu1;

          c1    = 0.D0;
          c2    = 0.D0;
          d2    = 0.D0;

        end

      end

      Tp(1,1)=   mu1*kap1 - mu2*kap2;
      Tp(1,2)=   gamm1 - 2.D0*mu2*kap1;
      Tp(1,3)=  -gamm1 + 2.D0*mu1*kap1;
      Tp(1,4)=   2.D0*(gamm1 - 2.D0*mu2*kap1);

      Tp(2,1)=   mu1*kap1 - mu2*kap2;
      Tp(2,2)=   gamm2 - 2.D0*mu1*kap2;
      Tp(2,3)=  -gamm2 + 2.D0*mu1*kap1;
      Tp(2,4)=   0.D0;

      Tp(3,1)=   gamm1 - mu1*kap1 - mu2*kap2;
      Tp(3,2)=   gamm1 - 2.D0*mu2*kap1;
      Tp(3,3)=   gamm1 - 2.D0*mu1*kap1;
      Tp(3,4)=   2.D0*(gamm1 - 2.D0*mu2*kap1);

      Tp(4,1)=   gamm2 - mu1*kap1 - mu2*kap2;
      Tp(4,2)=  -gamm2 + 2.D0*mu1*kap2;
      Tp(4,3)=   gamm2 - 2.D0*mu1*kap1;
      Tp(4,4)=   0.D0;

      Ts(1,1)=   mu2*kap2 - mu1*kap1;
      Ts(1,2)=   gamm2 - 2.D0*mu1*kap2;
      Ts(1,3)=  -gamm2 + 2.D0*mu2*kap2;
      Ts(1,4)=   2.D0*(gamm2 - 2.D0*mu1*kap2);

      Ts(2,1)=   mu2*kap2 - mu1*kap1;
      Ts(2,2)=   gamm1 - 2.D0*mu2*kap1;
      Ts(2,3)=  -gamm1 + 2.D0*mu2*kap2;
      Ts(2,4)=   0.D0;

      Ts(3,1)=   gamm2 - mu2*kap2 - mu1*kap1;
      Ts(3,2)=   gamm2 - 2.D0*mu1*kap2;
      Ts(3,3)=   gamm2 - 2.D0*mu2*kap2;
      Ts(3,4)=   2.D0*(gamm2 - 2.D0*mu1*kap2);

      Ts(4,1)=   gamm1 - mu2*kap2 - mu1*kap1;
      Ts(4,2)=  -gamm1 + 2.D0*mu2*kap1;
      Ts(4,3)=   gamm1 - 2.D0*mu2*kap2;
      Ts(4,4)=   0.D0;

      Up(1,1)= -(1.D0-     nu1)/(     mu1) * c2 +...
               (1.D0-2.D0*nu1)/(2.D0*mu1) * d2;
      Up(1,2)=   1.D0/(2.D0*mu1) * (c2-d2);
      Up(1,3)=  (3.D0-4.D0*nu1)/(2.D0*mu1) * (c2-d2);
      Up(1,4)=  -1.D0/(     mu1) * (c2-d2);

      Up(2,1)=  (1.D0-     nu2)/(     mu2) * c1 +...
               (1.D0-2.D0*nu2)/(2.D0*mu2) * d2;
      Up(2,2)=   1.D0/(2.D0*mu2) * (c1+d2);
      Up(2,3)=  -1.D0/(2.D0*mu2) * (c1-d2);
      Up(2,4)=   0.D0;

      Up(3,1)=  (1.D0-     nu1)/(     mu1) * d2 -...
               (1.D0-2.D0*nu1)/(2.D0*mu1) * c2;
      Up(3,2)=  -1.D0/(2.D0*mu1) * (c2-d2);
      Up(3,3)=  (3.D0-4.D0*nu1)/(2.D0*mu1) * (c2-d2);
      Up(3,4)=   1.D0/(     mu1) * (c2-d2);

      Up(4,1)= -(1.D0-     nu2)/(     mu2) * d2 -...
               (1.D0-2.D0*nu2)/(2.D0*mu2) * c1;
      Up(4,2)=   1.D0/(2.D0*mu2) * (c1+d2);
      Up(4,3)=  -1.D0/(2.D0*mu2) * (c1-d2);
      Up(4,4)=   0.D0;

      Us(1,1)= -(1.D0-     nu2)/(     mu2) * c1 -...
               (1.D0-2.D0*nu2)/(2.D0*mu2) * d2;
      Us(1,2)=   1.D0/(2.D0*mu2) * (c1+d2);
      Us(1,3)=  (3.D0-4.D0*nu2)/(2.D0*mu2) * (c1+d2);
      Us(1,4)=  -1.D0/(     mu2) * (c1+d2);

      Us(2,1)=  (1.D0-     nu1)/(     mu1) * c2 -...
               (1.D0-2.D0*nu1)/(2.D0*mu1) * d2;
      Us(2,2)=   1.D0/(2.D0*mu1) * (c2-d2);
      Us(2,3)=  -1.D0/(2.D0*mu1) * (c2+d2);
      Us(2,4)=   0.D0;

      Us(3,1)= -(1.D0-     nu2)/(     mu2) * d2 -...
               (1.D0-2.D0*nu2)/(2.D0*mu2) * c1;
      Us(3,2)=  -1.D0/(2.D0*mu2) * (c1+d2);
      Us(3,3)=  (3.D0-4.D0*nu2)/(2.D0*mu2) * (c1+d2);
      Us(3,4)=   1.D0/(     mu2) * (c1+d2);

      Us(4,1)=  (1.D0-     nu1)/(     mu1) * d2 -...
               (1.D0-2.D0*nu1)/(2.D0*mu1) * c2;
      Us(4,2)=   1.D0/(2.D0*mu1) * (c2-d2);
      Us(4,3)=  -1.D0/(2.D0*mu1) * (c2+d2);
      Us(4,4)=   0.D0;

end %Subroutine SET_EPAR

%     xx comp of the stress field due to the tensile component of the burger vector
function [output]=sxxUt(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxxUt
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = sxxUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*cos(delta) -...
             sxxUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*sin(delta);

end %Function sxxUt

%     xz comp of the stress field due to the tensile component of the burger vector
function [output]=sxzUt(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxzUt
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = sxzUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*cos(delta) -...
             sxzUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*sin(delta);

end %Function sxzUt

%     zz comp of the stress field due to the tensile component of the burger vector
function [output]=szzUt(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: szzUt
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = szzUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*cos(delta) -...
             szzUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*sin(delta);

end %Function szzUt

%     xx comp of the stress field due to the shear component of the burger vector
function [output]=sxxUb(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxxUb
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = sxxUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*sin(delta) +...
             sxxUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*cos(delta);

end% Function sxxUb

%     xz comp of the stress field due to the shear component of the burger vector
function [output]=sxzUb(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxzUb
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = sxzUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*sin(delta) +...
             sxzUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*cos(delta);

end %Function sxzUb

%     zz comp of the stress field due to the shear component of the burger vector
function [output]=szzUb(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: szzUb
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = szzUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*sin(delta) +...
             szzUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,Ih1,Ih2,dI,dIfill)*cos(delta);

end% Function szzUb

%     xx component of the stress field due to x component of the burger vector
function [output]=sxxUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,...
    Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxxUX
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      %REAL(8) :: x1,z1,x2,z2

      
      z1 = z0 - (Hd*0.5D0)*cos(delta);
      z2 = z0 + (Hd*0.5D0)*cos(delta);

      x1 = x0 - (Hd*0.5D0)*sin(delta);
      x2 = x0 + (Hd*0.5D0)*sin(delta);

    if (abs(z1)<1.D-12)
        z1 = 0.D0;
    end
    if (abs(z2)<1.D-12)
        z2 = 0.D0;
    end

      if (z1>=0.D0) %Then

         [dI]=fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

         if (z<0.D0) %Then    % formula for CASE1 within half-space 2
            output =...
            -3.D0/(4.D0*(1.D0-nu2)) *...                                   
            2.D0*mu2/(3.D0*pi) *...                                       
            (  (z-z1)*(3.D0*(x-x1)^2+(z-z1)^2)/ ...                         
            ((x-x1)^2+(z-z1)^2)^2 - ...                                  
            (z-z2)*(3.D0*(x-x2)^2+(z-z2)^2)/ ...                         
            ((x-x2)^2+(z-z2)^2)^2) - ...                        
            (2.D0*c1+d2)*dI(2,1) - (c1+d2)*dI(2,2) + (c1-d2)*dI(2,3);
         else                   % formula for CASE1 within half-space 1
            output =...
            -3.D0/(4.D0*(1.D0-nu1)) *...                                   
            2.D0*mu1/(3.D0*pi) *...                                       
            (  (z-z1)*(3.D0*(x-x1)^2+(z-z1)^2)/ ...                         
            ((x-x1)^2+(z-z1)^2)^2 - ...                                  
            (z-z2)*(3.D0*(x-x2)^2+(z-z2)^2)/  ...                        
            ((x-x2)^2+(z-z2)^2)^2         ) + ...                        
            (2.D0*c2-d2)*dI(1,1) - (c2-d2)*dI(1,2) - 3.D0*(c2-d2)*dI(1,3) +...
            2.D0*(c2-d2)*dI(1,4);
         end
         
      else

         if (z2<=0.D0) %Then

            [dI]=fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

            if (z>0.D0) %Then    % formula for CASE2 within half-space 1
                output =...
                -3.D0/(4.D0*(1.D0-nu1)) *...                                   
                2.D0*mu1/(3.D0*pi) *      ...                                 
                (  (z-z1)*(3.D0*(x-x1)^2+(z-z1)^2)/ ...                         
                ((x-x1)^2+(z-z1)^2)^2 -         ...                          
                (z-z2)*(3.D0*(x-x2)^2+(z-z2)^2)/  ...                        
                ((x-x2)^2+(z-z2)^2)^2         ) - ...                        
                (2.D0*c2-d2)*dI(2,1) - (c2-d2)*dI(2,2) + (c2+d2)*dI(2,3);
            else                   % formula for CASE2 within half-space 2
                output =...
                -3.D0/(4.D0*(1.D0-nu2)) * ...                                  
                2.D0*mu2/(3.D0*pi) *...                                       
                ( (z-z1)*(3.D0*(x-x1)^2+(z-z1)^2)/  ...                        
                ((x-x1)^2+(z-z1)^2)^2 - ...                                  
                (z-z2)*(3.D0*(x-x2)^2+(z-z2)^2)/  ...                        
                ((x-x2)^2+(z-z2)^2)^2         ) + ...                        
                (2.D0*c1+d2)*dI(1,1) - (c1+d2)*dI(1,2) - 3.D0*(c1+d2)*dI(1,3) +...
                2.D0*(c1+d2)*dI(1,4);
            end
         else                      % formulas for mixed CASE (2-1)

         [Ih1]=fill_I(Ih1,x1,z1,x,z);
         [Ih2]=fill_I(Ih2,x2,z2,x,z);

              if (z<0.D0) %Then  % formula for mixed CASE within half-space 2
                output =...
                -3.D0/(4.D0*(1.D0-nu2)) *...                                   
                2.D0*mu2/(3.D0*pi) *...                                       
                ( (z-z1)*(3.D0*(x-x1)^2+(z-z1)^2)/...                          
                ((x-x1)^2+(z-z1)^2)^2 - ...                                  
                (z-z2)*(3.D0*(x-x2)^2+(z-z2)^2)/   ...                       
                ((x-x2)^2+(z-z2)^2)^2         ) + ...                        
                (2.D0*c1+d2)*(Ih1(1,1)+Ih2(2,1)) - ...                            
                (c1+d2)*(Ih1(1,2)+3.D0*Ih1(1,3)-2.D0*Ih1(1,4)-Ih2(2,2)) -  ...    
                (c1-d2)*Ih2(2,3);
              else                 % formula for mixed CASE within half-space 1
                output =...
                -3.D0/(4.D0*(1.D0-nu1)) *...                                   
                2.D0*mu1/(3.D0*pi) *...                                       
                ( (z-z1)*(3.D0*(x-x1)^2+(z-z1)^2)/ ...                         
                ((x-x1)^2+(z-z1)^2)^2 - ...                                  
                (z-z2)*(3.D0*(x-x2)^2+(z-z2)^2)/  ...                        
                ((x-x2)^2+(z-z2)^2)^2         ) -...                         
                (2.D0*c2-d2)*(Ih2(1,1)+Ih1(2,1)) + ...                            
                (c2-d2)*(Ih2(1,2)+3.D0*Ih2(1,3)-2.D0*Ih2(1,4)-Ih1(2,2)) +...      
                (c2+d2)*Ih1(2,3);
              end
              
         end

      end

end %Function sxxUX

%     xz component of the stress field due to x component of the burger vector
function [output]=sxzUX (x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,...
    Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxzUX
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      %REAL(8) :: x1,z1,x2,z2
      
      z1 = z0 - (Hd*0.5D0)*cos(delta);
      z2 = z0 + (Hd*0.5D0)*cos(delta);

      x1 = x0 - (Hd*0.5D0)*sin(delta);
      x2 = x0 + (Hd*0.5D0)*sin(delta);

    if (abs(z1)<1.D-12)
        z1 = 0.D0;
    end
    if (abs(z2)<1.D-12)
        z2 = 0.D0;
    end

      if (z1>=0.D0) %Then

           [dI]=fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

           if (z<0.D0) %Then        % formula for CASE1 within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) *...                               
                2.D0*mu2/(3.D0*pi) *...                                   
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/...                            
                ((x-x1)^2+(z-z1)^2)^2 -...                               
                (x-x2)*((z-z2)^2-(x-x2)^2)/...                           
                ((x-x2)^2+(z-z2)^2)^2   ) -...                           
                c1*dI(4,1) - (c1+d2)*dI(4,2) + (c1-d2)*dI(4,3);
           else                       % formula for CASE1 within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) *...                               
                2.D0*mu1/(3.D0*pi) * ...                                  
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/ ...                           
                ((x-x1)^2+(z-z1)^2)^2 -...                               
                (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
                ((x-x2)^2+(z-z2)^2)^2   ) -...                           
                c2*dI(3,1) + (c2-d2)*dI(3,2) + ...                           
                (c2-d2)*dI(3,3) - 2.D0*(c2-d2)*dI(3,4);
           end

      else

           if (z2<=0.D0) %Then

              [dI]=fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

              if (z>0.D0) %Then     % formula for CASE2 within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) *...                               
                2.D0*mu1/(3.D0*pi) * ...                                  
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/  ...                          
                ((x-x1)^2+(z-z1)^2)^2 -  ...                             
                (x-x2)*((z-z2)^2-(x-x2)^2)/...                            
                ((x-x2)^2+(z-z2)^2)^2   ) - ...                          
                c2*dI(4,1) - (c2-d2)*dI(4,2) + (c2+d2)*dI(4,3);
              else                    % formula for CASE2 within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) *  ...                             
                2.D0*mu2/(3.D0*pi) * ...                                  
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/  ...                          
                ((x-x1)^2+(z-z1)^2)^2 -  ...                             
                (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
                ((x-x2)^2+(z-z2)^2)^2   ) - ...                          
                c1*dI(3,1) + (c1+d2)*dI(3,2) + ...                           
                (c1+d2)*dI(3,3) - 2.D0*(c1+d2)*dI(3,4);
              end

           else                       % formulas for mixed CASE

              [Ih1]= fill_I (Ih1,x1,z1,x,z);
              [Ih2]= fill_I (Ih2,x2,z2,x,z);

              if (z<0.D0) %Then     % formula for mixed CASE within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) *  ...                             
                2.D0*mu2/(3.D0*pi) * ...                                  
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/ ...                          
                ((x-x1)^2+(z-z1)^2)^2 -  ...                             
                (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
                ((x-x2)^2+(z-z2)^2)^2   ) -  ...                        
                c1*(Ih1(3,1)-Ih2(4,1)) +  ...                                
                (c1+d2)*(Ih1(3,2)+Ih1(3,3)-2.D0*Ih1(3,4)+Ih2(4,2)) - ...      
                (c1-d2)*Ih2(4,3);
              else                    % formula for mixed CASE within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) *  ...                             
                2.D0*mu1/(3.D0*pi) * ...                                  
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/...                            
                ((x-x1)^2+(z-z1)^2)^2 - ...                              
                (x-x2)*((z-z2)^2-(x-x2)^2)/...                            
                ((x-x2)^2+(z-z2)^2)^2   ) + ...                          
                c2*(Ih2(3,1)-Ih1(4,1)) - ...                                 
                (c2-d2)*(Ih2(3,2)+Ih2(3,3)-2.D0*Ih2(3,4)+Ih1(4,2)) + ...      
                (c2+d2)*Ih1(4,3);
              end

           end

      end

end %Function sxzUX

%     zz component of the stress field due to x component of the burger vector
function [output]=szzUX(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,...
    Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: szzUX
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      %REAL(8) :: x1,z1,x2,z2
      
      z1 = z0 - (Hd*0.5D0)*cos(delta);
      z2 = z0 + (Hd*0.5D0)*cos(delta);

      x1 = x0 - (Hd*0.5D0)*sin(delta);
      x2 = x0 + (Hd*0.5D0)*sin(delta);

    if (abs(z1)<1.D-12)
        z1 = 0.D0;
    end
    if (abs(z2)<1.D-12)
        z2 = 0.D0;
    end

      if (z1>=0.D0) %Then

           dI= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

           if (z<0.D0) %Then        % formula for CASE1 within half-space 2
            output =-3.D0/(4.D0*(1.D0-nu2)) * ...                          
            2.D0*mu2/(3.D0*pi) * ...                              
            ( (z-z1)*((z-z1)^2-(x-x1)^2)/ ...                       
            ((x-x1)^2+(z-z1)^2)^2 -   ...                        
            (z-z2)*((z-z2)^2-(x-x2)^2)/ ...                      
            ((x-x2)^2+(z-z2)^2)^2   ) - ...                      
            d2*dI(2,1) + (c1+d2)*dI(2,2) - (c1-d2)*dI(2,3);
           else                       % formula for CASE1 within half-space 1
            output =-3.D0/(4.D0*(1.D0-nu1)) *...                           
            2.D0*mu1/(3.D0*pi) *...                               
            ( (z-z1)*((z-z1)^2-(x-x1)^2)/...                        
            ((x-x1)^2+(z-z1)^2)^2 -  ...                         
            (z-z2)*((z-z2)^2-(x-x2)^2)/ ...                      
            ((x-x2)^2+(z-z2)^2)^2   ) -...                       
            d2*dI(1,1) + (c2-d2)*dI(1,2) - ...                    
            (c2-d2)*dI(1,3) - 2*(c2-d2)*dI(1,4);
           end

      else
           if (z2<=0.D0) %Then

              [dI]= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

              if (z>0.D0) %Then     % formula for CASE2 within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) * ...                          
                2.D0*mu1/(3.D0*pi) * ...                              
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/ ...                       
                ((x-x1)^2+(z-z1)^2)^2 -  ...                         
                (z-z2)*((z-z2)^2-(x-x2)^2)/...                        
                ((x-x2)^2+(z-z2)^2)^2    ) + ...                     
                d2*dI(2,1) + (c2-d2)*dI(2,2) - (c2+d2)*dI(2,3);
              else                    % formula for CASE2 within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) * ...                          
                2.D0*mu2/(3.D0*pi) * ...                              
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/ ...                       
                ((x-x1)^2+(z-z1)^2)^2 - ...                          
                (z-z2)*((z-z2)^2-(x-x2)^2)/ ...                       
                ((x-x2)^2+(z-z2)^2)^2    ) + ...                     
                d2*dI(1,1) + (c1+d2)*dI(1,2) - ...                       
                (c1+d2)*dI(1,3) - 2*(c1+d2)*dI(1,4);
              end

           else                       % formulas for mixed CASE

              Ih1= fill_I (Ih1,x1,z1,x,z);
              Ih2= fill_I (Ih2,x2,z2,x,z);

              if (z<0.D0) %Then     % formula for mixed CASE within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) *...                           
                2.D0*mu2/(3.D0*pi) * ...                              
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/...                        
                ((x-x1)^2+(z-z1)^2)^2 - ...                          
                (z-z2)*((z-z2)^2-(x-x2)^2)/ ...                       
                ((x-x2)^2+(z-z2)^2)^2    ) +...                      
                d2*(Ih1(1,1)+Ih2(2,1)) + ...                             
                (c1+d2)*(Ih1(1,2)-Ih1(1,3)-2.D0*Ih1(1,4)-Ih2(2,2)) + ...  
                (c1-d2)*Ih2(2,3);
              else                    % formula for mixed CASE within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) *...                           
                2.D0*mu1/(3.D0*pi) *     ...                          
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/  ...                      
                ((x-x1)^2+(z-z1)^2)^2 -   ...                        
                (z-z2)*((z-z2)^2-(x-x2)^2)/ ...                       
                ((x-x2)^2+(z-z2)^2)^2    ) +  ...                    
                d2*(Ih2(1,1)+Ih1(2,1)) -  ...                            
                (c2-d2)*(Ih2(1,2)-Ih2(1,3)-2.D0*Ih2(1,4)-Ih1(2,2)) -  ... 
                (c2+d2)*Ih1(2,3);
              end

           end

      end

end% Function szzUX

%     xx component of the stress field due to z component of the burger vector
function [output]=sxxUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,...
    Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxxUZ
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      %REAL(8) :: x1,z1,x2,z2  

      z1 = z0 - (Hd*0.5D0)*cos(delta);
      z2 = z0 + (Hd*0.5D0)*cos(delta);

      x1 = x0 - (Hd*0.5D0)*sin(delta);
      x2 = x0 + (Hd*0.5D0)*sin(delta);

    if (abs(z1)<1.D-12)
        z1 = 0.D0;
    end
    if (abs(z2)<1.D-12)
        z2 = 0.D0;
    end

      if (z1>=0.D0) %Then

           [dI]= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

           if (z<0.D0) %Then        % formula for CASE1 within half-space 2
            output =-3.D0/(4.D0*(1.D0-nu2)) * ...                              
            2.D0*mu2/(3.D0*pi) *   ...                                
            ( (x-x1)*((z-z1)^2-(x-x1)^2)/ ...                           
            ((x-x1)^2+(z-z1)^2)^2 -     ...                          
            (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
            ((x-x2)^2+(z-z2)^2)^2    ) -  ...                        
            (c1+2.D0*d2)*dI(4,1) - (c1+d2)*dI(4,2) + (c1-d2)*dI(4,3);
           else                       % formula for CASE1 within half-space 1
            output =-3.D0/(4.D0*(1.D0-nu1)) *  ...                             
            2.D0*mu1/(3.D0*pi) *    ...                               
            ( (x-x1)*((z-z1)^2-(x-x1)^2)/ ...                           
            ((x-x1)^2+(z-z1)^2)^2 -    ...                           
            (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
            ((x-x2)^2+(z-z2)^2)^2    ) -   ...                       
            (c2-2.D0*d2)*dI(3,1) + (c2-d2)*dI(3,2) -   ...                
            3.D0*(c2-d2)*dI(3,3) + 2.D0*(c2-d2)*dI(3,4);
           end

      else

           if (z2<=0.D0) %Then

              [dI]= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

              if (z>0.D0) %Then     % formula for CASE2 within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) *...                               
                2.D0*mu1/(3.D0*pi) *...                                   
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/...                            
                ((x-x1)^2+(z-z1)^2)^2 - ...                              
                (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
                ((x-x2)^2+(z-z2)^2)^2    ) - ...                         
                (c2-2.D0*d2)*dI(4,1) - (c2-d2)*dI(4,2) + (c2+d2)*dI(4,3);
              else                    % formula for CASE2 within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) *...                               
                2.D0*mu2/(3.D0*pi) *...                                   
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/...                            
                ((x-x1)^2+(z-z1)^2)^2 - ...                              
                (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
                ((x-x2)^2+(z-z2)^2)^2    ) - ...                         
                (c1+2.D0*d2)*dI(3,1) + (c1+d2)*dI(3,2) - ...                  
                3.D0*(c1+d2)*dI(3,3) + 2.D0*(c1+d2)*dI(3,4);
              end

           else                       % formulas for mixed CASE

              Ih1= fill_I (Ih1,x1,z1,x,z);
              Ih2= fill_I (Ih2,x2,z2,x,z);

              if (z<0.D0) %Then     % formula for mixed CASE within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) * ...                              
                2.D0*mu2/(3.D0*pi) *...                                   
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/...                            
                ((x-x1)^2+(z-z1)^2)^2 - ...                              
                (x-x2)*((z-z2)^2-(x-x2)^2)/ ...                           
                ((x-x2)^2+(z-z2)^2)^2    ) - ...                         
                (c1+2.D0*d2)*(Ih1(3,1)-Ih2(4,1)) + ...                        
                (c1+d2)*(Ih1(3,2)-3.D0*Ih1(3,3)+2.D0*Ih1(3,4)+Ih2(4,2)) -  ...
                (c1-d2)*Ih2(4,3);
              else                    % formula for mixed CASE within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) * ...                              
                2.D0*mu1/(3.D0*pi) *         ...                          
                ( (x-x1)*((z-z1)^2-(x-x1)^2)/       ...                     
                ((x-x1)^2+(z-z1)^2)^2 -              ...                 
                (x-x2)*((z-z2)^2-(x-x2)^2)/         ...                   
                ((x-x2)^2+(z-z2)^2)^2    ) -           ...               
                (c2-2.D0*d2)*(Ih1(4,1)-Ih2(3,1)) -          ...               
                (c2-d2)*(Ih1(4,2)+Ih2(3,2)-3.D0*Ih2(3,3)+2.D0*Ih2(3,4)) +  ...
                (c2+d2)*Ih1(4,3);
              end

           end

      end

end% Function sxxUZ

%     xz component of the stress field due to z component of the burger vector
function [output]=sxzUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,...
    Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: sxzUZ
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      %REAL(8) :: x1,z1,x2,z2 

      z1 = z0 - (Hd*0.5D0)*cos(delta);
      z2 = z0 + (Hd*0.5D0)*cos(delta);

      x1 = x0 - (Hd*0.5D0)*sin(delta);
      x2 = x0 + (Hd*0.5D0)*sin(delta);

    if (abs(z1)<1.D-12)
        z1 = 0.D0;
    end
    if (abs(z2)<1.D-12)
        z2 = 0.D0;
    end

      if (z1>=0.D0) %Then

           [dI]= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

           if (z<0.D0) %Then        % formula for CASE1 within half-space 2
            output =-3.D0/(4.D0*(1.D0-nu2)) * ...                             
            2.D0*mu2/(3.D0*pi) * ...                                 
            ( (z-z1)*((z-z1)^2-(x-x1)^2)/ ...                          
            ((x-x1)^2+(z-z1)^2)^2 - ...                             
            (z-z2)*((z-z2)^2-(x-x2)^2)/...                           
            ((x-x2)^2+(z-z2)^2)^2    ) + ...                        
            d2*dI(2,1) + (c1+d2)*dI(2,2) - (c1-d2)*dI(2,3);
           else                       % formula for CASE1 within half-space 1
            output =-3.D0/(4.D0*(1.D0-nu1)) *...                              
            2.D0*mu1/(3.D0*pi) *    ...                              
            ( (z-z1)*((z-z1)^2-(x-x1)^2)/  ...                         
            ((x-x1)^2+(z-z1)^2)^2 -    ...                          
            (z-z2)*((z-z2)^2-(x-x2)^2)/   ...                        
            ((x-x2)^2+(z-z2)^2)^2    ) +     ...                    
            d2*dI(1,1) + (c2-d2)*dI(1,2) -        ...                    
            (c2-d2)*dI(1,3) + 2*(c2-d2)*dI(1,4);
           end

      else
           if (z2<=0.D0) %Then

              [dI]= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

              if (z>0.D0) %Then     % formula for CASE2 within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) * ...                             
                2.D0*mu1/(3.D0*pi) * ...                                 
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/...                           
                ((x-x1)^2+(z-z1)^2)^2 - ...                            
                (z-z2)*((z-z2)^2-(x-x2)^2)/ ...                          
                ((x-x2)^2+(z-z2)^2)^2    ) - ...                        
                d2*dI(2,1) + (c2-d2)*dI(2,2) - (c2+d2)*dI(2,3);
              else                    % formula for CASE2 within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) * ...                             
                2.D0*mu2/(3.D0*pi) *  ...                                
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/  ...                         
                ((x-x1)^2+(z-z1)^2)^2 -  ...                            
                (z-z2)*((z-z2)^2-(x-x2)^2)/...                           
                ((x-x2)^2+(z-z2)^2)^2    ) - ...                        
                d2*dI(1,1) + (c1+d2)*dI(1,2) - ...                           
                (c1+d2)*dI(1,3) + 2*(c1+d2)*dI(1,4);
              end

           else                       % formulas for mixed CASE

              Ih1= fill_I (Ih1,x1,z1,x,z);
              Ih2= fill_I (Ih2,x2,z2,x,z);

              if (z<0.D0) %Then     % formula for mixed CASE within half-space 2
                output =-3.D0/(4.D0*(1.D0-nu2)) *...                              
                2.D0*mu2/(3.D0*pi) * ...                                 
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/...                           
                ((x-x1)^2+(z-z1)^2)^2 -   ...                           
                (z-z2)*((z-z2)^2-(x-x2)^2)/...                           
                ((x-x2)^2+(z-z2)^2)^2   ) - ...                         
                d2*(Ih1(1,1)+Ih2(2,1)) + ...                                 
                (c1+d2)*(Ih1(1,2)-Ih1(1,3)+2.D0*Ih1(1,4)-Ih2(2,2)) +   ...   
                (c1-d2)*Ih2(2,3);
              else                    % formula for mixed CASE within half-space 1
                output =-3.D0/(4.D0*(1.D0-nu1)) *...                              
                2.D0*mu1/(3.D0*pi) * ...                                 
                ( (z-z1)*((z-z1)^2-(x-x1)^2)/ ...                          
                ((x-x1)^2+(z-z1)^2)^2 - ...                             
                (z-z2)*((z-z2)^2-(x-x2)^2)/ ...                          
                ((x-x2)^2+(z-z2)^2)^2    ) -...                         
                d2*(Ih1(2,1)+Ih2(1,1)) + ...                                 
                (c2-d2)*(Ih1(2,2)-Ih2(1,2)+Ih2(1,3)-2.D0*Ih2(1,4)) -   ...   
                (c2+d2)*Ih1(2,3);
              end

           end

      end

end %Function sxzUZ

%     zz component of the stress field due to z component of the burger vector
function [output]=szzUZ(x0,z0,delta,Hd,x,z,mu1,mu2,nu1,nu2,c1,c2,d2,...
    Ih1,Ih2,dI,dIfill)
      %IMPLICIT NONE

      %REAL(8) :: szzUZ
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      %REAL(8) :: x1,z1,x2,z2 

      z1 = z0 - (Hd*0.5D0)*cos(delta);
      z2 = z0 + (Hd*0.5D0)*cos(delta);

      x1 = x0 - (Hd*0.5D0)*sin(delta);
      x2 = x0 + (Hd*0.5D0)*sin(delta);

    if (abs(z1)<1.D-12)
        z1 = 0.D0;
    end
    if (abs(z2)<1.D-12)
        z2 = 0.D0;
    end

      if (z1>=0.D0) %Then

         [dI]= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

         if (z<0.D0) %Then        % formula for CASE1 within half-space 2
            output = 3.D0/(4.D0*(1.D0-nu2)) *...                                 
            (2.D0*mu2)/(3.D0*pi)* ...                                     
            ( (x-x1)*((x-x1)^2+3.D0*(z-z1)^2)/  ...                        
            ((x-x1)^2+(z-z1)^2)^2 - ...                                 
            (x-x2)*((x-x2)^2+3.D0*(z-z2)^2)/ ...                         
            ((x-x2)^2+(z-z2)^2)^2         ) - ...                       
            c1*dI(4,1) + (c1+d2)*dI(4,2) - (c1-d2)*dI(4,3)  ;              
         else                       % formula for CASE1 within half-space 1
            output = 3.D0/(4.D0*(1.D0-nu1)) *...                                 
            (2.D0*mu1)/(3.D0*pi)* ...                                     
            ( (x-x1)*((x-x1)^2+3.D0*(z-z1)^2)/  ...                        
            ((x-x1)^2+(z-z1)^2)^2 -  ...                                
            (x-x2)*((x-x2)^2+3.D0*(z-z2)^2)/  ...                        
            ((x-x2)^2+(z-z2)^2)^2         ) -  ...                      
            c2*dI(3,1) - (c2-d2)*dI(3,2) - (c2-d2)*dI(3,3) - ...             
            2.D0*(c2-d2)*dI(3,4);
         end
      else

         if (z2<=0.D0) %Then

            [dI]= fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z);

            if (z>0.D0) %Then     % formula for CASE2 within half-space 1
                output = 3.D0/(4.D0*(1.D0-nu1)) * ...                                
                (2.D0*mu1)/(3.D0*pi)*...                                      
                ( (x-x1)*((x-x1)^2+3.D0*(z-z1)^2)/  ...                        
                ((x-x1)^2+(z-z1)^2)^2 -...                                  
                (x-x2)*((x-x2)^2+3.D0*(z-z2)^2)/...                          
                ((x-x2)^2+(z-z2)^2)^2         ) -...                        
                c2*dI(4,1) + (c2-d2)*dI(4,2) - (c2+d2)*dI(4,3);
            else                    % formula for CASE2 within half-space 2
                output = 3.D0/(4.D0*(1.D0-nu2)) * ...                                
                (2.D0*mu2)/(3.D0*pi)*...                                      
                ( (x-x1)*((x-x1)^2+3.D0*(z-z1)^2)/ ...                          
                ((x-x1)^2+(z-z1)^2)^2 -  ...                                 
                (x-x2)*((x-x2)^2+3.D0*(z-z2)^2)/  ...                         
                ((x-x2)^2+(z-z2)^2)^2          ) - ...                       
                c1*dI(3,1) - (c1+d2)*dI(3,2) - (c1+d2)*dI(3,3) - ...             
                2.D0*(c1+d2)*dI(3,4);
            end
         else                       % formulas for mixed CASE

            Ih1= fill_I (Ih1,x1,z1,x,z);
            Ih2= fill_I (Ih2,x2,z2,x,z);

            if (z<0.D0) %Then     % formula for mixed CASE within half-space 2
                output = 3.D0/(4.D0*(1.D0-nu2)) * ...                                
                (2.D0*mu2)/(3.D0*pi)*...                                       
                ( (x-x1)*((x-x1)^2+3.D0*(z-z1)^2)/ ...                          
                ((x-x1)^2+(z-z1)^2)^2 - ...                                  
                (x-x2)*((x-x2)^2+3.D0*(z-z2)^2)/ ...                          
                ((x-x2)^2+(z-z2)^2)^2         ) - ...                        
                c1*(Ih1(3,1)-Ih2(4,1)) -...                                      
                (c1+d2)*(Ih1(3,2)+Ih1(3,3)+2.D0*Ih1(3,4)+Ih2(4,2)) +   ...        
                (c1-d2)*Ih2(4,3);
            else                    % formula for mixed CASE within half-space 1
                output = 3.D0/(4.D0*(1.D0-nu1)) * ...                                
                (2.D0*mu1)/(3.D0*pi)* ...                                      
                ( (x-x1)*((x-x1)^2+3.D0*(z-z1)^2)/...                           
                ((x-x1)^2+(z-z1)^2)^2 - ...                                  
                (x-x2)*((x-x2)^2+3.D0*(z-z2)^2)/   ...                        
                ((x-x2)^2+(z-z2)^2)^2        ) -...                          
                c2*(Ih1(4,1)-Ih2(3,1)) +  ...                                    
                (c2-d2)*(Ih1(4,2)+Ih2(3,2)+Ih2(3,3)+2.D0*Ih2(3,4)) -     ...      
                (c2+d2)*Ih1(4,3);
            end
            
         end

      end

end %Function szzUZ


function [dI]=fill_dI(dI,dIfill,x0,z0,delta,Hd,x,z)
      %IMPLICIT NONE
      
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dI(4,4),I1(4,4),I2(4,4)
      %REAL(8) :: x1,z1,x2,z2
      %INTEGER :: n,m
      %dI=zeros(4);    

    z1 = z0 - (Hd*0.5D0)*cos(delta);
    z2 = z0 + (Hd*0.5D0)*cos(delta);

    x1 = x0 - (Hd*0.5D0)*sin(delta);
    x2 = x0 + (Hd*0.5D0)*sin(delta);

    if (abs(z1)<1.D-12)
        z1 = 0.D0;
    end
    if (abs(z2)<1.D-12)
        z2 = 0.D0;
    end

    dI=fill_I(dI,x1,z1,x,z);
    dIfill=fill_I(dIfill,x2,z2,x,z);

    for n=1:4
       for m=1:4
          %dI(n,m) = 0.D0;
          dI(n,m) = dI(n,m) - dIfill(n,m);
       end
    end

end %Subroutine fill_dI


function [I]=fill_I(I,xd,zd,x,z)
     %IMPLICIT NONE

      %REAL(8) :: I(4,4)
      %REAL(8) :: xd,zd,x,z
      %I=zeros(4);   

      %Introducing some constants
      zdpz=zd+z;
      zdpz2=zdpz^2;
      xmxd=x-xd;
      xmxd2=xmxd^2;
      xmxd2pzdpz22=(xmxd2 + zdpz2)^2;
      zmzd=z-zd;
      zmzd2=zmzd^2;
      zdmz=zd-z;
      xmxd2pzmzd22=(xmxd2+zmzd2)^2;
      xmxd2pzmzd23=(xmxd2 + zmzd2)^3;
      xmxd2pzdpz23=(xmxd2 + zdpz2)^3;
      zdpz2mxmxd2=zdpz2-xmxd2;

      
     I(1,1)=(zdpz)/(xmxd2 + zdpz2);

     I(1,2)=z*(zdpz2mxmxd2)/...         
             xmxd2pzdpz22;

     I(1,3)=zd*(zdpz2mxmxd2)/...      
              xmxd2pzdpz22;

     I(1,4)=2.D0*z*zd*(zdpz)*...                  
           (zdpz2 - 3.D0*xmxd2)/...      
           xmxd2pzdpz23;

     I(2,1)=(zdmz)/(xmxd2 + zmzd2);

     I(2,2)=z*(zmzd2 - xmxd2)/...         
             xmxd2pzmzd22;

     I(2,3)=zd*(zmzd2 - xmxd2)/...        
              xmxd2pzmzd22;

     I(2,4)=2.D0*z*zd*(zdmz)*...                  
           (zmzd2 - 3.D0*xmxd2)/...      
           xmxd2pzmzd23;

     I(3,1)=(xmxd)/(xmxd2 + zdpz2);

     I(3,2)=2.D0*z*(xmxd)*(zdpz)/...              
              xmxd2pzdpz22;

     I(3,3)=2.D0*zd*(xmxd)*(zdpz)/...             
                xmxd2pzdpz22;

     I(3,4)=2.D0*z*zd*(xmxd)*...                  
           (3.D0*zdpz2mxmxd2)/...      
           xmxd2pzdpz23;

     I(4,1)=(xmxd)/(xmxd2 + zmzd2);

     I(4,2)=2.D0*z*(xmxd)*(zdmz)/...              
               xmxd2pzmzd22;

     I(4,3)=2.D0*zd*(xmxd)*(zdmz)/...             
                xmxd2pzmzd22;

     I(4,4)=2.D0*z*zd*(xmxd)*...                  
           (3.D0*(zdmz)^2 - xmxd2)/...      
           xmxd2pzmzd23;


end %Subroutine fill_I



