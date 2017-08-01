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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

clear;close all

% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/'); end 
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir


	%Script to run dislocations on a fault surface

	% This follows the same unit convention as Poly3d so ''assumes dimensionally consistent units for physical quantities with
	% dimensions of length or stress. Thus if you use kilometers for any quantity with
	% dimensions of length (e.g. a coordinate position), you must use kilometers for all
	% quantities with dimensions of length (e.g. displacement discontinuity).''
	% Andrew Lyell Thomas, Masters thesis, Stanford University. 

cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

	%Conventions: 
		%See figure 6.14 in David Pollards book for tensors. 
		%Geological convention is that a confining pressure is a positive stress. Geological convention is used throughout this script. 
		%For a normal stress: positive if it produces compression in the material.
		%For a shear stress: negative if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is - if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear stress Sxy will cause RightLateral movement.
		
		%Strain uses the same convention as engineering stress convention where extension is a positve value. When entering strain values remember this. 
		%For a normal stress: positive (negative) if it produces tension (compression) in the material.
		%For a shear stress: positive if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is + if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear strain Exy will cause Leftlateral movement.
		
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %The file must list the traingle vertex's then triangles. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 string='ValleySurface_ProperEq_1000Faces.ts';  
 %string='ValleySurface_ProperEq_5700Faces.ts'; 
 [ Points,Triangles ] = GoCadAsciiReader( string,pathstring );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
freesurface_height = 0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=1;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
lambda=1;%4000;  %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;

	%Equations you can use to calculate elastic parameters that you need above. 
	
% E = 10000;                        %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
% lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard


	%Function to extract the three vertex's for each triangle and 
	%triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp= zeros(size(Triangles(:,1))); 

%TEST locking edges
TR = triangulation(Triangles,Points(:,2:4));
[TotalConnectionDistance,SortedTriangles,Noconnections] = ConnectedTrianglesFinder(TR,MidPoint);
Flag= Noconnections(:,2)<3;
Fdisp(Flag)=1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface
	%Runs a constant slip on every triangle. A mixture of two slips will give an inclined slip.
	%Output 'stresses' are blank arrays of 0's as this option is not driven by boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option E = Finding stresses beneath topographic surface. Adapted from 2d script from Martel and Muller (2000) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);
Sxx=Tectsx; %Tectonic stress (xx)
Syy=Tectsy; %Tectonic stress (xx)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
if any(Fdisp)==1
[Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

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
[sigx,sigy,sigxy,xobs,yobs] = SavageSlope1984FortranCodeMATLABFunc(ts,rg,pr,a,b,u,v);
%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses  on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%TotalStress,Stress and StressChange are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
	%Strain is the first 6 coloumns and stress the last 6. Stress is the regional stress 
	%StressChange is the stress on the points from the event.
	%TotalStress is the driving stress and stress change from the event added. 
	
 [TotalStrainStress,StrainStressChange,StrainStressRemote,Syy,Szz,Sxy,Sxz,Syz]=CalculateStressOnSurroundingPoints(StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,mu,lambda,X,Y,Z,Sxx,Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints( StrikeSlipDisp,DipSlipDisp,TensileSlipDisp, nu, X,Y,Z, P1, P2, P3,halfspace);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%XX YY ZZ XY XZ YZ
Sxx=TotalStrainStress(:,7);
Szz=TotalStrainStress(:,9);
Sxz=TotalStrainStress(:,11);

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


