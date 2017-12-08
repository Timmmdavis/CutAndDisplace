function [X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,Ux,Uy,Uz]=MengEshelbyFunc(a,b,c,x,y,z,sxx,syy,szz,sxy,sxz,syz,E,nu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the data structure needed for an ellipsoid inhomogeniety problem;
% solve the problem with an equevilent Eshelby inclusion problem;
% demonstrate the outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sharma's ellpdoidal heterogeneity with observation line
%clear all;
%poisson ratio of matrix
incl.vm= nu;
%elastic modulus of matrix
incl.Em=E;
%hetergeneity poisson ratio
incl.vh=0;
%hetergeneity elastic modulus
incl.Eh=0;

% dimensiona of the ellipsoid.
incl.dim=[a b c]; %incl.dim=[1 1 1];
% ortation angles in radian applied to the ellipsoid in sequence of [x y z];.
incl.ang = [0 0 0];

% remote stress ordered by sigma11, sigma12, sigma13, sigma22,sigma23,sigma33
%incl.stressvec=[0;0;0;0;0;0];
incl.stressvec=[sxx;sxy;sxz;syy;syz;szz];
incl.eigp=[0;0;0;0;0;0];
% first observaton grid
x = x+eps; 
y = y;
z = z;
incl.grid{1} = {x,y,z};


% call Eshelby solver,arg: 'disp','stress','strain' only output displacements by default
 incl.sol = Esh_sol_NOPARA(incl,'disp','stress','strain');
%incl.sol = Esh_sol(incl,'disp','stress','strain');

% % demonstrate the third grid's results for the first and second components of
% % the output displacement and stress.
X = squeeze(incl.sol.grid{1,1});
Y = squeeze(incl.sol.grid{1,2});
Z = squeeze(incl.sol.grid{1,3});
u = squeeze(incl.sol.u{:});

stress = squeeze(incl.sol.stress{:});

if z==0; %2d plane, this check could be improved so it catches any 2d plane not just flat ones
Sxx=squeeze(stress(:,:,1));
Syy=squeeze(stress(:,:,4));
Szz=squeeze(stress(:,:,6));
Sxy=squeeze(stress(:,:,2));
Sxz=squeeze(stress(:,:,3));
Syz=squeeze(stress(:,:,5));
Ux=squeeze(u(:,:,1));
Uy=squeeze(u(:,:,2));
Uz=squeeze(u(:,:,3));    
else    
Sxx=squeeze(stress(:,:,:,1));
Syy=squeeze(stress(:,:,:,4));
Szz=squeeze(stress(:,:,:,6));
Sxy=squeeze(stress(:,:,:,2));
Sxz=squeeze(stress(:,:,:,3));
Syz=squeeze(stress(:,:,:,5));
Ux=squeeze(u(:,:,:,1));
Uy=squeeze(u(:,:,:,2));
Uz=squeeze(u(:,:,:,3));
end
% Sxx(R<1.5)=nan;
% Sxy(R<1.5)=nan;
% Syy(R<1.5)=nan;
% Ux(R<1.5)=nan;
% Uy(R<1.5)=nan;
% 
% figure;quiver(X,Y,Ux,Uy);
% figure;scatter(X(:),Y(:),30,Ux(:));title('Ux');colorbar;
% figure;scatter(X(:),Y(:),30,Uy(:));title('Uy');colorbar;
% 
% 
% Sxx=reshape(Sxx,(size(X)));
% Syy=reshape(Syy,(size(X)));
% Sxy=reshape(Sxy,(size(X)));
% 
% figure;subplot(2,2,1),contourf(X,Y,Sxx);%caxis([-0.55 0.55])
% xlabel('x'); ylabel('z'); axis('equal');  title('Sxx');colorbar;
% subplot(2,2,2),contourf(X,Y,Syy);%caxis([-0.55 0.55])
% xlabel('x'); ylabel('z'); axis('equal');  title('Syy');colorbar;
% subplot(2,2,3),contourf(X,Y,Sxy);%caxis([0.65 1.15])
% xlabel('x'); ylabel('z'); axis('equal');  title('Sxy');colorbar;
% 
% figure;subplot(2,2,1),scatter(X(:),Y(:),15,Sxx(:));%caxis([-0.4 0.4])
% xlabel('x'); ylabel('z'); axis('equal');  title('Sxx');colorbar;
% subplot(2,2,2),scatter(X(:),Y(:),15,Syy(:));%caxis([-0.4 0.4])
% xlabel('x'); ylabel('z'); axis('equal');  title('Syy');colorbar;
% subplot(2,2,3),scatter(X(:),Y(:),15,Sxy(:));%caxis([0.7 1.1])
% xlabel('x'); ylabel('z'); axis('equal');  title('Sxy');colorbar;
% 
% 
% S1 = (Sxx+Syy)/2 + sqrt( ((Sxx-Syy)/2).^2 + Sxy.^2);
% S2 = (Sxx+Syy)/2 - sqrt( ((Sxx-Syy)/2).^2 + Sxy.^2);
% 
% figure;subplot(1,2,1),scatter(X(:),Y(:),15,S1(:));%caxis([-2 2])
% xlabel('x'); ylabel('z'); axis('equal'); title('S1');colorbar;
% subplot(1,2,2),scatter(X(:),Y(:),15,S2(:));%caxis([-2 2])
% xlabel('x'); ylabel('z'); axis('equal'); title('S2');colorbar;
% 
% 
