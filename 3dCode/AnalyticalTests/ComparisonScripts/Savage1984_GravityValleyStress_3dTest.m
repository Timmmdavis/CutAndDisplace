% Test: Gravitational unloading stresses after erosion of a valley into an
% elastic medium
% 
% Solution: Coded a MATLAB version of the FORTRAN code for the analytical
% solution. This code is from: Savage, William Z., Philip S. Powers, and
% Henri S. Swolfs. In situ geomechanics of crystalline and sedimentary
% rocks; Part V, RVT, a Fortran program for an exact elastic solution for
% tectonics and gravity stresses in isolated symmetric ridges and valleys.
% No. 84-827. US Geological Survey,, 1984. This runs through the Savage
% code calculating the stresses at points below the valley. A surface is
% loaded that reflects this 2d valley profile at Y=0 and extends far in
% either Y direction with the same profile. The weight of the material is
% set to 1. A plot is drawn of for the analytical and BEM solution along
% Y=0. The test checks the mean residual stress is below 0.01mpa.
% 
% Proof: This tests the implementation of the elastic stresses due to
% topography as described in Martel and Muller 2000.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%Clearing old figures and arrays. 
clear;close all

%Loading a colourmap to be used for figures.
cmap = colormap_cpt('Ccool-warm');
cmap2 = colormap_cpt('Ccool-warm2');
		
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading Gocad ascii data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 string='ValleySurface_ProperEq_1000Faces.ts';  
 %string='ValleySurface_ProperEq_5700Faces.ts'; 
 [ Points,Triangles ] = GoCadAsciiReader( string );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define full/halfspace and elastic constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Defining halfspace. 1 is on, 0 is off. The code will run much faster with
% this off as the calculations are simpler in a full space.
halfspace = 0; 

% Defining the height of the freesurface relative to the value of 0 of
% your imported surfaces. 
freesurface_height = 0;

% Defining elastic constants
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

%Locked Els
Fdisp= zeros(size(Triangles(:,1))); 

%TEST locking edges
TR = triangulation(Triangles,Points(:,2:4));
[TotalConnectionDistance,SortedTriangles,Noconnections] = ConnectedTrianglesFinder(TR,MidPoint);
Flag= Noconnections(:,2)<3;
Fdisp(Flag)=1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%
%If Cartesian stress tensor components are not known but principal stresses
%are use the function 'TensorsFromPrincipal' to calculate the tensors in a
%Cartesian reference frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;


P=1; 	%Density of elastic material 
g=1;		%Acceleration due to gravity
Tectsx=0; %Tectonic stress (xx)
Tectsy=0; %Tectonic stress (yy)
Sxx = ((nu/(1-nu))*P*g*+MidPoint(:,3))+Tectsx;	%Equation 1b, Martel and Muller
Szz = P*g*+MidPoint(:,3);			%Equation 1a, Martel and Muller
Syy = nu*(Sxx+Szz);	%Plane strain conditions
Sxy = 0;          						%Positive is right lateral movement
Sxz = 0;
Syz = 0;
Option='B'; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(Option,'A')    
    [Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
     Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
     FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);
end

Sxx=Tectsx; %Tectonic stress (xx)
Syy=Tectsy; %Tectonic stress (xx)

    % Removing any fixed elements
if any(Fdisp)==1
    [Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP A: Analytical solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%
P=1; %Density
ts=Tectsx; %Can only do X as its a 2d solution
rg= P*g; %Weight (Dens*AccelGrav)
pr=0.25;
%A and B describing slope, See paper for diagram
a=2;
b=-1;
%Setting up grid points U and V, these are mapped to a different location
%later
u = linspace(0,4,50);
v = linspace(-0.3265,-4,46); %Points slightly below 0 as we do not want observation points on the elements in the TDE code 
[u,v] = meshgrid(u,v);
[sigx,sigy,sigxy,xobs,yobs] = Savage1984_GravityValleyStress(ts,rg,pr,a,b,u,v);
%%%%%%%%%%%%%%%%%%%%%

%Getting points from this solution to calc numerically. 
X=xobs;
X=X+50; %Imported surfaces valley mid points are at 50. 
Z=yobs;
Y=zeros(size(X));
ArrayReshapeSize=size(X);
X=X(:);
Y=Y(:);
Z=Z(:);

%Drawing figure of surface and the observation points
figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceColor', [0.5 0 0.9 ],'FaceAlpha',(.2));
hold on
scatter3(X(:),Y(:),Z(:),5,'k')  %Showing surface and obs points
xlabel('x'); ylabel('z'); axis('equal'); title('Surface and Obs Points');




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
 [StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]=...
     CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d( Dss,Dds,Dn, nu, X,Y,Z, P1, P2, P3,halfspace);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

Sxx=StressTChg(:,1);
Szz=StressTChg(:,3);
Sxz=StressTChg(:,5);

X=reshape(X,ArrayReshapeSize);
Z=reshape(Z,ArrayReshapeSize);
Sxx=reshape(Sxx,ArrayReshapeSize);
Szz=reshape(Szz,ArrayReshapeSize);
Sxz=reshape(Sxz,ArrayReshapeSize);


%Adding vertical Stress profile: 
SxxGrav = (nu/(1-nu))*P*g*Z;	%Equation 1b, Martel and Muller
SzzGrav = P*g*Z;		%Equation 1a, Martel and Muller
SxzGrav =zeros(size(SxxGrav));

Sxx=SxxGrav+Sxx;
Szz=SzzGrav+Szz;
Sxz=SxzGrav+Sxz;


%Setting contour steps, the steps set below are the same as the figures in the Martel
%2000 paper.
LvlstpX=0.2;
LvlstpY=0.5;
LvlstpXY=0.02;

%Setting axis limits so supurious values do not change the overall trend. 
caxissxx=[min(sigx(:)),max(sigx(:))];
caxissyy=[min(sigy(:)),max(sigy(:))];
caxissxy=[min(sigxy(:)),max(sigxy(:))];
%Rounding values, rounds the axis values to the set step.
Xrnd=1/LvlstpX;
Yrnd=1/LvlstpY;
XYrnd=1/LvlstpXY;
caxissxx = (round(caxissxx*Xrnd))/Xrnd; %Rounding to nearest 0.2
caxissyy = (round(caxissyy*Yrnd))/Yrnd; %Rounding yy to nearest 0.5 (same as levelstep)
caxissxy = (round(caxissxy*XYrnd))/XYrnd; %Rounding to nearest 0.02

%Creating specific contour levels, See MATLABS ''highlight-specific-contour-levels'' Doc. 
xmin = floor(min(sigx(:)));
xmax = ceil(max(sigx(:)));
xindex = xmin:LvlstpX:xmax;

ymin = floor(min(sigy(:)));
ymax = ceil(max(sigy(:)));
yindex = ymin:LvlstpY:ymax;

xymin = floor(min(sigxy(:)));
xymax = ceil(max(sigxy(:)));
xyindex = xymin:LvlstpXY:xymax;

%Plotting figures for Savage 2d solution
figure,figure,subplot(2,3,4),[C,h] =contourf(xobs,yobs,sigx);
caxis(caxissxx);colormap(cmap2),colorbar;
hold on
contour(xobs,yobs,sigx,xindex,'k');
set (h, 'levelstep', LvlstpX);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('z'); axis('equal'); title('SxxSavage');
hold off
%Syy
subplot(2,3,5),[C,h] =contourf(xobs,yobs,sigy);
caxis(caxissyy);colormap(cmap2),colorbar;
hold on
contour(xobs,yobs,sigy,yindex,'k');
set (h, 'levelstep', LvlstpY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('z'); axis('equal'); title('SzzSavage');
hold off
%Sxy
subplot(2,3,6),[C,h] =contourf(xobs,yobs,sigxy);
caxis(caxissxy);colormap(cmap2),colorbar;
hold on
contour(xobs,yobs,sigxy,xyindex,'k');
set (h, 'levelstep', LvlstpXY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('z'); axis('equal'); title('SxzSavage');
hold off

%Plotting figures of TDE result
%Sxx
subplot(2,3,1),[C,h] =contourf(X,Z,Sxx);
caxis(caxissxx);colormap(cmap2),colorbar;
hold on
contour(X,Z,Sxx,xindex,'k');
set (h, 'levelstep', LvlstpX);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('z'); axis('equal'); title('SxxScript');
hold off

%Szz
subplot(2,3,2),[C,h] =contourf(X,Z,Szz);
caxis(caxissyy);colormap(cmap2),colorbar;
hold on
contour(X,Z,Szz,yindex,'k');
set (h, 'levelstep', LvlstpY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('z'); axis('equal'); title('SzzScript');
hold off

%Sxz
subplot(2,3,3),[C,h] =contourf(X,Z,Sxz);
caxis(caxissxy);colormap(cmap2),colorbar;
hold on
contour(X,Z,Sxz,xyindex,'k');
set (h, 'levelstep', LvlstpXY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('z'); axis('equal'); title('SxzScript');
hold off

SxxRes=Sxx-sigx;
bad=isnan(SxxRes);SxxRes=SxxRes(~bad);%removing nans to get proper mean
%disp('MeanResidualSxx');disp(mean(SxxRes(:)))
fprintf('MeanResidualSxx %i.\n',mean(SxxRes(:))) %prints on one line unlike 'disp
SzzRes=Szz-sigy;
bad=isnan(SzzRes);SzzRes=SzzRes(~bad);%removing nans to get proper mean
fprintf('MeanResidualSzz %i.\n',mean(SzzRes(:))) %prints on one line unlike 'disp
SxzRes=Sxz-sigxy;
bad=isnan(SxzRes);SxzRes=SxzRes(~bad);%removing nans to get proper mean
fprintf('MeanResidualSxz %i.\n',mean(SxzRes(:))) %prints on one line unlike 'disp


P=[mean(SxxRes(:)),mean(SzzRes(:)),mean(SxzRes(:))];
if any(P>0.15)
    error('Error in surrounding grid is too large. This no longer matches analytical solutions, check the comparative images')
else
    disp('Everything looks good, tolerance checks mean stress residuals and flags errors above 0.15')
end


%Percentages
% SlipResPercx=mean(abs((100./sigx(:)).*Sxx(:))); 
% SlipResPercy=mean(abs((100./sigy(:)).*Syy(:))); 
% SlipResPercxy=mean(abs((100./sigxy(:)).*Sxy(:))); 
%Mean Percentage Error
SlipResPercx=mean(100-(abs((100./sigx(:)).*Sxx(:))));
SlipResPercy=mean(100-(abs((100./sigy(:)).*Syy(:)))); 
bad=abs((100./sigxy(:)))==inf;
SlipResPercxy=100-(abs((100./sigxy(:)).*Sxy(:))); 
SlipResPercxy(bad)=[];
SlipResPercxy=mean(SlipResPercxy);


