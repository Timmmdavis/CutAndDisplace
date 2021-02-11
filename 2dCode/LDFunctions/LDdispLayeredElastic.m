
function [Disp]=LDdispLayeredElastic(Dn,Ds,xe,ye,Beta,a,x,y,mu1,mu2,nu1,nu2)
% LDdispLayeredElastic: LineDislocationInducedDisplacementLayeredElastic,
%               Computes the influence (displacement) of a single planar
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
% Disp        - Is the displacement caused by the movement of the
%             dislocation at the observataion points. [Ux,Uy].
%
% Example usage (compare to homogenous case):
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
% [Disp]=LDdispLayeredElastic(Dn,Ds,xe,ye,Beta,a,x,y,mu1,mu2,nu1,nu2);
% [ Ux,Uy ] = ExtractCols( Disp );
% [Ux,Uy]=ReshapeData2d( dimx,dimy,Ux,Uy );
% figure,quiver(x(:),y(:),Ux(:),Uy(:));
% title('Layered elastic')
% 
% %Homogenous case
% nu=nu1;
% mu=mu1;
% [ k,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );
% [Disp] = LDdispFS(x(:),y(:),xe,ye,a,Beta,Ds,Dn,nu);
% [ Ux,Uy ] = ExtractCols( Disp );
% [Ux,Uy]=ReshapeData2d( dimx,dimy,Ux,Uy );
% figure,quiver(x(:),y(:),Ux(:),Uy(:));
% title('Homogenous elastic')
% 
%  Author: Tim Davis/Francesco Maccaferri
%  Copyright 2021, Tim Davis, Potsdam University
%  Modified from of Francesco's 2D fluid-filled fracture propagation scripts:
%  Modified so the notation is the same as in CutAndDisplace


a=a*2; %To make length same as LD convention
Beta=-Beta;%To rotate same as LD convention
delta=Beta+pi/2;

[Up,Us,Tp,Ts,~,~,~]=SET_EPAR(mu1,mu2,nu1,nu2);

x=x(:);
y=y(:);
Uxfill=zeros(size(x));
Uyfill=zeros(size(x));
for i=1:length(x)
    Ux = uxUt(xe,ye,delta,a,x(i),y(i),Up,Us,Tp,Ts,nu1,nu2)*Dn + uxUb(xe,ye,delta,a,x(i),y(i),Up,Us,Tp,Ts,nu1,nu2)*Ds;
    Uy = uzUt(xe,ye,delta,a,x(i),y(i),Up,Us,Tp,Ts,nu1,nu2)*Dn + uzUb(xe,ye,delta,a,x(i),y(i),Up,Us,Tp,Ts,nu1,nu2)*Ds;
    Uxfill(i)=Ux;
    Uyfill(i)=Uy;
end
Disp=[Uxfill,Uyfill];
      
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

      mu1=mu1;
      mu2=mu2;
      nu1=nu1;
      nu2=nu2;
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

%     x comp of the displacement field due to the tensile component of the burger vector
function [output]=uxUt(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uxUt
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = uxUX(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*cos(delta) -...
            uxUZ(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*sin(delta);

end %Function uxUt

%     x comp of the displacement field due to the shear component of the burger vector
function [output]=uxUb(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uxUb
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = uxUX(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*sin(delta) +...
            uxUZ(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*cos(delta);

end %Function uxUb

%     z comp of the displacement field due to the tensile component of the burger vector
function [output]=uzUt(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uzUt
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = uzUX(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*cos(delta) -...
            uzUZ(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*sin(delta);

end %Function uzUt

%     z comp of the displacement field due to the shear component of the burger vector
function [output]=uzUb(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uzUb
      %REAL(8) :: x0,z0,delta,Hd,x,z

      output = uzUX(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*sin(delta) +...
            uzUZ(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)*cos(delta);

end %Function uzUb

%     x component of displacement due to x component of the burger vector
function [output]=uxUX(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uxUX
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      %REAL(8) :: x1,z1,x2,z2,dGx
      %INTEGER :: j
 

      output = 0.D0;

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

      [dGx]= fill_dGx(x0,z0,delta,Hd,x,z,nu1,nu2);

      output = output + dGx;

      if (z1>=0.D0) %Then

        [dY]=fill_dY(x0,z0,delta,Hd,x,z);

        if (z>=0.D0) %Then    % formula for CASE1 within half-space 1
          for j=1:4
            output = output + (-Up(1,j)*dY(1,j));
          end
        else                   % formula for CASE1 within half-space 2
          for j=1:4
            output = output + (-Up(2,j)*dY(2,j));
          end
        end
      else
        if (z2<=0.D0) %Then

          [dY]=fill_dY (x0,z0,delta,Hd,x,z);

          if (z<=0.D0) %Then  % formula for CASE2 within half-space 2
            for j=1:4
                output = output + (-Us(1,j)*dY(1,j));
            end
          else                 % formula for CASE2 within half-space 1
            for j=1:4
                output = output + (-Us(2,j)*dY(2,j));
            end
          end
        else                   % formulas for mixed CASE (2-1)

          [Y1]=fill_Y(x1,z1,Hd,x,z);
          [Y2]=fill_Y(x2,z2,Hd,x,z);

          if (z>=0.D0) %Then  % formula for mixed CASE within half-space 1
            for j=1:4
                output = output + (-Us(2,j)*Y1(2,j) + Up(1,j)*Y2(1,j));
            end
          else                 % formula for mixed CASE within half-space 2
            for j=1:4
                output = output + (-Us(1,j)*Y1(1,j) + Up(2,j)*Y2(2,j));
            end
          end
        end
      end

end %Function uxUX

%     x component of displacement due to z component of the burger vector
function [output]=uxUZ(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uxUZ
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      %REAL(8) :: x1,z1,x2,z2,dHx
      %INTEGER :: j

      output = 0.D0;

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

      [dHx]= fill_dHx (x0,z0,delta,Hd,x,z,nu1,nu2);

      output = output - dHx;

      if (z1>=0.D0) %Then

        [dY]=fill_dY (x0,z0,delta,Hd,x,z);

        if (z>=0.D0) %Then    % formula for CASE1 within half-space 1
          for j=1:4
          output = output - (-Tp(3,j)*dY(3,j));
          end
        else                   % formula for CASE1 within half-space 2
          for j=1:4
          output = output - (-Tp(4,j)*dY(4,j)); %+
          end
        end
      else
        if (z2<=0.D0) %Then

          [dY]=fill_dY (x0,z0,delta,Hd,x,z);

          if (z<=0.D0) %Then  % formula for CASE2 within half-space 2
            for j=1:4
            output = output - (-Ts(3,j)*dY(3,j));
            end
          else                 % formula for CASE2 within half-space 1
            for j=1:4
            output = output - (-Ts(4,j)*dY(4,j)); %+
            end
          end
        else                   % formulas for mixed CASE (2-1)

          [Y1]=fill_Y(x1,z1,Hd,x,z);
          [Y2]=fill_Y(x2,z2,Hd,x,z);

          if (z>=0.D0) %Then  % formula for mixed CASE within half-space 1
            for j=1:4
            output = output - (-Ts(4,j)*Y1(4,j) + Tp(3,j)*Y2(3,j));
            end
          else                 % formula for mixed CASE within half-space 2
            for j=1:4
            output = output - (-Ts(3,j)*Y1(3,j) + Tp(4,j)*Y2(4,j));
            end
          end
        end
      end

end %Function uxUZ

%     z component of displacement due to x component of the burger vector
function [output]=uzUX(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uzUX
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      %REAL(8) :: x1,z1,x2,z2,dGz
      %INTEGER :: j
      
      output = 0.D0;

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

      [dGz]= fill_dGz (x0,z0,delta,Hd,x,z,nu1,nu2);

      output = output + dGz;

      if (z1>=0.D0) %Then

        [dY]= fill_dY (x0,z0,delta,Hd,x,z);

        if (z>=0.D0) %Then    % formula for CASE1 within half-space 1
          for j=1:4
            output = output + (-Up(3,j)*dY(3,j));
          end
        else                   % formula for CASE1 within half-space 2
          for j=1:4
            output = output + (-Up(4,j)*dY(4,j));
          end
        end
      else
        if (z2<=0.D0) %Then

          [dY]=fill_dY (x0,z0,delta,Hd,x,z);

          if (z<=0.D0) %Then  % formula for CASE2 within half-space 2
            for j=1:4
                output = output + (-Us(3,j)*dY(3,j));
            end
          else                 % formula for CASE2 within half-space 1
            for j=1:4
                output = output + (-Us(4,j)*dY(4,j));
            end
          end
        else                   % formulas for mixed CASE (2-1)

          Y1= fill_Y(x1,z1,Hd,x,z);
          Y2= fill_Y(x2,z2,Hd,x,z);

          if (z>=0.D0) %Then  % formula for mixed CASE within half-space 1
            for j=1:4
                output = output + (-Us(4,j)*Y1(4,j) + Up(3,j)*Y2(3,j));
            end
          else                 % formula for mixed CASE within half-space 2
            for j=1:4
                output = output + (-Us(3,j)*Y1(3,j) + Up(4,j)*Y2(4,j));
            end
          end
        end
      end

end %Function uzUX

%     z component of displacement due to z component of the burger vector
function [output]=uzUZ(x0,z0,delta,Hd,x,z,Up,Us,Tp,Ts,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: uzUZ
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      %REAL(8) :: x1,z1,x2,z2,dHz
      %INTEGER :: j
      
      output = 0.D0;

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

      [dHz]= fill_dHz (x0,z0,delta,Hd,x,z,nu1,nu2);

      output = output - dHz;

      if (z1>=0.D0) %Then

        [dY]=fill_dY (x0,z0,delta,Hd,x,z);

        if (z>=0.D0) %Then    % formula for CASE1 within half-space 1
          for j=1:4
            output = output - (-Tp(1,j)*dY(1,j));
          end
        else                   % formula for CASE1 within half-space 2
          for j=1:4
          output = output - (-Tp(2,j)*dY(2,j));
          end
        end
      else
        if (z2<=0.D0) %Then

          [dY]=fill_dY (x0,z0,delta,Hd,x,z);

          if (z<=0.D0) %Then  % formula for CASE2 within half-space 2
            for j=1:4
                output = output - (-Ts(1,j)*dY(1,j));
            end
          else                 % formula for CASE2 within half-space 1
            for j=1:4
                output = output - (-Ts(2,j)*dY(2,j));
            end
          end
        else                   % formulas for mixed CASE (2-1)

          [Y1]=fill_Y(x1,z1,Hd,x,z);
          [Y2]=fill_Y(x2,z2,Hd,x,z);

          if (z>=0.D0) %Then  % formula for mixed CASE within half-space 1
            for j=1:4
                output = output - (-Ts(2,j)*Y1(2,j) + Tp(1,j)*Y2(1,j));
            end
          else                 % formula for mixed CASE within half-space 2
            for j=1:4
                output = output - (-Ts(1,j)*Y1(1,j) + Tp(2,j)*Y2(2,j));
            end
          end
        end
      end

end %Function uzUZ



function [dY]=fill_dY(x0,z0,delta,Hd,x,z)
      %IMPLICIT NONE
      
      %REAL(8) :: x0,z0,delta,Hd,x,z

      %REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      %REAL(8) :: x1,z1,x2,z2
      %INTEGER :: n,m
      dY=zeros(4);    

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

    Y1=fill_Y(x1,z1,Hd,x,z);
    Y2=fill_Y(x2,z2,Hd,x,z);

    for n=1:4
       for m=1:4
          dY(n,m) = 0.D0;
          dY(n,m) = Y1(n,m) - Y2(n,m);
       end
    end

end %Subroutine fill_dY


function [I]=fill_I(xd,zd,x,z)
      %IMPLICIT NONE

      %REAL(8) :: I(4,4)
      %REAL(8) :: xd,zd,x,z
      I=zeros(4);   
      
     I(1,1)=(zd+z)/((x-xd)^2 + (z+zd)^2);

     I(1,2)=z*((z+zd)^2 - (x-xd)^2)/...         
             ((x-xd)^2 + (z+zd)^2)^2;

     I(1,3)=zd*((z+zd)^2 - (x-xd)^2)/...      
              ((x-xd)^2 + (z+zd)^2)^2;

     I(1,4)=2.D0*z*zd*(zd+z)*...                  
           ((z+zd)^2 - 3.D0*(x-xd)^2)/...      
           ((x-xd)^2 + (z+zd)^2)^3;

     I(2,1)=(zd-z)/((x-xd)^2 + (z-zd)^2);

     I(2,2)=z*((z-zd)^2 - (x-xd)^2)/...         
             ((x-xd)^2 + (z-zd)^2)^2;

     I(2,3)=zd*((z-zd)^2 - (x-xd)^2)/...        
              ((x-xd)^2 + (z-zd)^2)^2;

     I(2,4)=2.D0*z*zd*(zd-z)*...                  
           ((z-zd)^2 - 3.D0*(x-xd)^2)/...      
           ((x-xd)^2 + (z-zd)^2)^3;

     I(3,1)=(x-xd)/((x-xd)^2 + (z+zd)^2);

     I(3,2)=2.D0*z*(x-xd)*(zd+z)/...              
              ((x-xd)^2 + (z+zd)^2)^2;

     I(3,3)=2.D0*zd*(x-xd)*(zd+z)/...             
                ((x-xd)^2 + (z+zd)^2)^2;

     I(3,4)=2.D0*z*zd*(x-xd)*...                  
           (3.D0*(zd+z)^2 - (x-xd)^2)/...      
           ((x-xd)^2 + (z+zd)^2)^3;

     I(4,1)=(x-xd)/((x-xd)^2 + (z-zd)^2);

     I(4,2)=2.D0*z*(x-xd)*(zd-z)/...              
               ((x-xd)^2 + (z-zd)^2)^2;

     I(4,3)=2.D0*zd*(x-xd)*(zd-z)/...             
                ((x-xd)^2 + (z-zd)^2)^2;

     I(4,4)=2.D0*z*zd*(x-xd)*...                  
           (3.D0*(zd-z)^2 - (x-xd)^2)/...      
           ((x-xd)^2 + (z-zd)^2)^3;

end %Subroutine fill_I


function [Y]=fill_Y(xd,zd,Hd,x,z)
      %IMPLICIT NONE

      %REAL(8) :: Y(4,4)
      %REAL(8) :: xd,zd,x,z,Hd
      Y=zeros(4);

      %REAL(8) :: arctanP,arctanM,segnoP,segnoM

      if (abs(x-xd)<1.D-12) %Then % IMPORTANT: this will happen every time fill_Y
                                    % is called for computing the displacement at the
                                    % CENTER of a VERTICAL elementary disl. (that is
                                    % needed for computing the asimmetry in the opening
                                    % and shear displacement)
        if (zd+z>0.D0) %Then
          arctanP =  pi*0.5D0;
        else
          arctanP = -pi*0.5D0;
        end

        if (zd-z>0.D0) %Then
          arctanM =  pi*0.5D0;
        else
          arctanM = -pi*0.5D0;
        end

        if (abs(zd)>1.D-12) %Then

          if (zd>1.D-12) %Then
            segnoP = +1.D0;
            segnoM = +1.D0;
          else
            segnoP = -1.D0;
            segnoM = -1.D0;
          end

        else
          segnoP = +signFORTRAN(1.D0,z);
          segnoM = -signFORTRAN(1.D0,z);
        end

      else % (when abs(x-xd) > 1.D-12)

        arctanP = atan((zd+z)/(x-xd));
        arctanM = atan((zd-z)/(x-xd));

        if (abs(zd)<1.D-12) %Then
          segnoP = +signFORTRAN(1.D0,(x-xd)*z);%+(x-xd)/abs(x-xd)*signFORTRAN(1.D0,z)
          segnoM = -signFORTRAN(1.D0,(x-xd)*z);%-(x-xd)/abs(x-xd)*signFORTRAN(1.D0,z)
        else
          segnoP =  signFORTRAN(1.D0,(x-xd)*zd);%(x-xd)*zd/(abs((x-xd)*zd))
          segnoM =  signFORTRAN(1.D0,(x-xd)*zd);%(x-xd)*zd/(abs((x-xd)*zd))
        end

      end

%       NUMERICAL LIMIT INSTEAD OF ANALYTICAL (SEE ABOVE)
%       zdp = zd
%       zdm = zd
% 
%       If (abs(x-xd)<1.D-14) %Then %(disl. verticale)
%         x = x  + 1.D-12
%       EndIf
%       If (abs(z-zd)<1.D-14) %Then %(disl. orizzontale)
%         If (abs(zd)<1.D-14) %Then
%           zdp = zd + 1.D-12
%           zdm = zd + 1.D-12
%         EndIf
%       Else
%         If (abs(zd)<1.D-14) %Then
%           zdp = zd + 1.D-12*z/abs(z)
%           zdm = zd - 1.D-12*z/abs(z)
%         EndIf
%       EndIf
% 
%                  Y(1,1)= -(-pi*0.5D0*(x-xd)*zdp/abs((x-xd)*zdp)
%      &                   + atan((z+zd)/(x-xd)))

                 Y(1,1)= pi*0.5D0*segnoP - arctanP;

                 Y(1,2)=  z*(x-xd)/...                        
                           ((x-xd)^2 + (z+zd)^2);

                 Y(1,3)= zd*(x-xd)/...                        
                           ((x-xd)^2 + (z+zd)^2);

                 Y(1,4)=2.D0*z*zd*(z+zd)*(x-xd)/...           
                          (((x-xd)^2 + (z+zd)^2)^2);

%                  Y(2,1)= -(-pi*0.5D0*(x-xd)*zdm/abs((x-xd)*zdm)
%      &                   + atan((zd-z)/(x-xd)))

                 Y(2,1)= pi*0.5D0*segnoM - arctanM;

                 Y(2,2)= z*(x-xd)/...                         
                        ((x-xd)^2 + (z-zd)^2);

                 Y(2,3)= zd*(x-xd)/...                        
                        ((x-xd)^2 + (z-zd)^2);          

                 Y(2,4)= 2.D0*z*zd*(zd-z)*(x-xd)/...          
                        (((x-xd)^2 + (z-zd)^2)^2);

                 Y(3,1)= 0.5D0*log( ((x-xd)^2 + (z+zd)^2)/...    
                                   (Hd^2) );

                 Y(3,2)= -z*(z+zd)/...                        
                        ((x-xd)^2 + (z+zd)^2);

                 Y(3,3)= -zd*(z+zd)/...                       
                        ((x-xd)^2 + (z+zd)^2);

                 Y(3,4)= z*zd*((x-xd)^2 - (z+zd)^2)/...     
                               (((x-xd)^2 + (z+zd)^2)^2);

                 Y(4,1)= 0.5D0*log( ((x-xd)^2 + (z-zd)^2)/...    
                                   (Hd^2) );

                 Y(4,2)= -z*(zd-z)/...                        
                        ((x-xd)^2 + (z-zd)^2);

                 Y(4,3)= -zd*(zd-z)/...                       
                        ((x-xd)^2 + (z-zd)^2);

                 Y(4,4)= z*zd*((x-xd)^2 - (z-zd)^2)/...     
                               (((x-xd)^2 + (z-zd)^2)^2);

end %Subroutine fill_Y


function [dGx]=fill_dGx(x0,z0,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: dGx,x0,z0,x,z,Hd,delta
      dGx=zeros(1);
      %REAL(8) :: z1,z2,x1,x2
      %REAL(8) :: Gx1,Gx2

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

      [Gx1]=fill_Gx (x1,z1,delta,Hd,x,z,nu1,nu2);
      [Gx2]=fill_Gx (x2,z2,delta,Hd,x,z,nu1,nu2);

      dGx = Gx1 - Gx2;

end  %Subroutine fill_dGx


function [dGz]=fill_dGz(x0,z0,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: dGz,x0,z0,x,z,Hd,delta

      %REAL(8) :: z1,z2,x1,x2
      %REAL(8) :: Gz1,Gz2
      dGz=zeros(1);

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

      [Gz1]=fill_Gz(x1,z1,delta,Hd,x,z,nu1,nu2);
      [Gz2]=fill_Gz(x2,z2,delta,Hd,x,z,nu1,nu2);

      dGz = Gz1 - Gz2;

end %Subroutine fill_dGz


function [dHx]=fill_dHx(x0,z0,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: dHx,x0,z0,x,z,Hd,delta
      dHx=zeros(1);
      %REAL(8) :: z1,z2,x1,x2
      %REAL(8) :: Hx1,Hx2

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

      [Hx1]=fill_Hx(x1,z1,delta,Hd,x,z,nu1,nu2);
      [Hx2]=fill_Hx(x2,z2,delta,Hd,x,z,nu1,nu2);

      dHx = Hx1 - Hx2;

end %Subroutine fill_dHx


function [dHz]=fill_dHz(x0,z0,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: dHz,x0,z0,x,z,Hd,delta

      %REAL(8) :: z1,z2,x1,x2
      %REAL(8) :: Hz1,Hz2
      dHz=zeros(1);

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

      [Hz1]=fill_Hz (x1,z1,delta,Hd,x,z,nu1,nu2);
      [Hz2]=fill_Hz (x2,z2,delta,Hd,x,z,nu1,nu2);

      dHz = Hz1 - Hz2;

end %Subroutine fill_dHz


function [Gx]=fill_Gx(xd,zd,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: Gx,xd,zd,x,z,Hd,delta

      %REAL(8) :: nu,Arctan
        Gx=zeros(1);
      
      if (z>=0.D0) %Then
        nu = nu1;
      else
        nu = nu2;
      end

      if ((x-xd)*cos(delta)-(z-zd)*sin(delta)>1.D-12) %Then
          Arctan =  pi*0.5D0 +...                                           
                   atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/...         
                         ((x-xd)*cos(delta)-(z-zd)*sin(delta)) );
          else
          if ((x-xd)*cos(delta)-(z-zd)*sin(delta)<-1.D-12) %Then
              Arctan = -pi*0.5D0 +...                                           
                       atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/...         
                             ((x-xd)*cos(delta)-(z-zd)*sin(delta)) );
              else
              if ((x-xd)*sin(delta)+(z-zd)*cos(delta)>1.D-12) %Then
                  Arctan = pi;
                  else
                  if ((x-xd)*sin(delta)+(z-zd)*cos(delta)<-1.D-12) %Then
                    Arctan = 0.D0;
                  else
                    Arctan = 3.D0*pi/4.D0;
                  end
              end
          end
      end

      Gx = 1.D0/(2.D0*pi)*...                                           
          (Arctan + 1.D0/(2.D0*(1.D0-nu))*...                          
                     (((x-xd)*(z-zd))/((x-xd)^2+(z-zd)^2)) );

end %Subroutine fill_Gx


function [Gz]=fill_Gz(xd,zd,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: Gz,xd,zd,x,z,Hd,delta

      %REAL(8) :: nu
        Gz=zeros(1);
      
      if (z>=0.D0) %Then
        nu = nu1;
      else
        nu = nu2;
      end

      Gz = -1.D0/(4.D0*pi*(1.D0-nu))*...                                
           ( (1.D0-2.D0*nu)*0.5D0*log( ((x-xd)^2+(z-zd)^2)/...       
                                        (Hd^2) ) -...                 
             (z-zd)^2/((x-xd)^2+(z-zd)^2) );

end %Subroutine fill_Gz


function [Hx]=fill_Hx(xd,zd,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: Hx,xd,zd,x,z,Hd,delta

      %REAL(8) :: nu
      Hx=zeros(1);

      if (z>=0.D0) %Then
        nu = nu1;
      else
        nu = nu2;
      end
                                                                     
      Hx = -1.D0/(4.D0*pi*(1.D0-nu))*...                                
           ( (1.D0-2.D0*nu)*0.5D0*log( ((x-xd)^2+(z-zd)^2)/...       
                                        (Hd^2) ) + ...                
             (z-zd)^2/((x-xd)^2+(z-zd)^2) );

end %Subroutine fill_Hx


function [Hz]=fill_Hz(xd,zd,delta,Hd,x,z,nu1,nu2)
      %IMPLICIT NONE

      %REAL(8) :: Hz,xd,zd,x,z,Hd,delta
       Hz=zeros(1);
      %REAL(8) :: nu,Arctan

      if (z>=0.D0) %%Then
        nu = nu1;
      else
        nu = nu2;
      end

      if ((x-xd)*cos(delta)-(z-zd)*sin(delta)>1.D-12) 
          Arctan =  pi*0.5D0 +...                                            
                   atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/...          
                         ((x-xd)*cos(delta)-(z-zd)*sin(delta)) );
          else
          if ((x-xd)*cos(delta)-(z-zd)*sin(delta)<-1.D-12) %Then
              Arctan = -pi*0.5D0 +...                                            
                       atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/...          
                             ((x-xd)*cos(delta)-(z-zd)*sin(delta)) );
              else
              if ((x-xd)*sin(delta)+(z-zd)*cos(delta)>1.D-12) %Then
                  Arctan = pi;
                  else
                  if ((x-xd)*sin(delta)+(z-zd)*cos(delta)<-1.D-12) %Then
                      Arctan = 0.D0;
                      else
                      Arctan = 3.D0*pi/4.D0;
                  end
              end
          end
      end

      Hz = -1.D0/(2.D0*pi)*...                                           
          (Arctan - 1.D0/(2.D0*(1.D0-nu))*...                           
                      (((x-xd)*(z-zd))/((x-xd)^2+(z-zd)^2)) );
     
end

function [A]=signFORTRAN(A,B)
%Fortran signFORTRAN(A,B) returns the value of A with the signFORTRAN of B.
if (B<0) && (0<A)
    A=-A;
end
if (A<0) && (0<B)
    A=-A;
end

end

      %End Module DISL2D
