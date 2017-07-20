% Test: Gravitational unloading stresses after erosion of a valley into an
% elastic medium: Checks to see if the script can reproduce solutions for
% stresses caused by a valley profile in an elastic medium.
% 
% Solution: Coded a MATLAB version of the FORTRAN code for the analytical
% solution. This code is from: Savage, William Z., Philip S. Powers, and
% Henri S. Swolfs. In situ geomechanics of crystalline and sedimentary
% rocks; Part V, RVT, a Fortran program for an exact elastic solution for
% tectonics and gravity stresses in isolated symmetric ridges and valleys.
% No. 84-827. US Geological Survey,, 1984. This runs through the Savage
% code calculating the stresses at points below the valley. These points
% are then passed to the gravitational stress code using the parameters set
% up in the savage part. The weight of the material is set to 1. A plot is
% drawn of for the analytical and BEM solution. The test checks the mean
% residual stress is below 0.01mpa.
%  
% Proof: This tests the implementation of the elastic stresses due to
% topography as described in Martel and Muller 2000.

%   Copyright 2017, Tim Davis, The University of Aberdeen
%===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir

close all; clear 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution, Savage 1984.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for Savage code, change below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tectonic Stress
ts=0; 

%Weight and elastic constants
P=1; 	%Density of elastic material 
g=1;	%Acceleration due to gravity
rg= P*g; %Weight (Dens*AccelGrav)
pr=0.25;

%A and B describing slope, See paper for diagram
a=2;
b=-1;

%Setting up grid points U and V, these are mapped to a different location
%later
u = linspace(0,4,50);
v = linspace(-0.3265,-4,46); %Points slightly below 0 as we do not want observation points on the elements for the twodd code 
% u = linspace(-0,4,50);
% v = linspace(0,-4,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sigx,sigy,sigxy,xobs,yobs] = SavageSlope1984FortranCodeMATLABFunc(ts,rg,pr,a,b,u,v);

%Equations 7a and b of
%Martel, S.J. and Muller, J.R., 2000. A two-dimensional boundary element
%method for calculating elastic gravitational stresses in slopes. pure and
%applied geophysics, 157(6-8), pp.989-1007.
%Calculating slope profile
u2 = linspace(-40,40,800);
v2 = zeros(size(u2)); 
x=u2+((a.*b.*u2)./((u2.^2)+((v2-a).^2)));
y=v2-((a.*b.*(v2-a))./((u2.^2)+((v2-a).^2)));


		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY defined above
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1))));
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms

    %Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
Fdisp=zeros(size(Pointsxy(:,1)));
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
    %Use this if the algorithm is throwing errors for vertex's above 0m.
    
freesurface_height = 0;
Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=1;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
%lambda=4000;%4000;  %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.
nu = pr;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;

	%Equations you can use to calculate elastic parameters that you need above. 
	
E = mu*(2*(1+nu)) ;                     %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard



    %Reshaping the list of XY into usable variables
[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng,NUM,XBEG,XEND,YBEG,YEND]...
 = CreateElements2d( Pointsxy,mystruct );


line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
title('fractures'), xlabel('x'), ylabel('y');axis equal;

   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option E = Finding stresses beneath topographic surface. Adapted from 2d script from Martel and Muller (2000) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=P; 	%Density of elastic material 
g=g;%9.8;		%Acceleration due to gravity
tects=0; %Tectonic stress (xx)
% % Creating gravitational stress on each elements centre presuming the
% % ground surface was originally at 0m. 
Sxx = ((nu/(1-nu))*P*g*ye)+tects;	%Equation 1b, Martel and Muller
Syy = P*g*ye;		%Equation 1a, Martel and Muller
Sxy = zeros(NUM,1);  %Positive is right lateral movement
Option='B'; %Valley surface deformation after erosion

[ShearDisp,TensileDisp,Sxx,Syy,Sxy]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);
%The resultant stresses in the body from these dislocations is taken away
%from a calculated far field stress that pervades the body, this far field stress is calculated using the constants P&g. 
Sxx=tects; %Tectonic stress (xx)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X=xobs;
Y=yobs;
dimx = length(X(:,1));
dimy = length(X(1,:));

    %Replacing XY of any obs points for Nan if these lie on top of the end points of a
    %dislocation, without this it can result in spurious values. 
[ X,Y  ] = NullPointsLyingOnElement( X,Y,XBEG,YBEG,XEND,YEND,NUM );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses and displacements on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Displacements on Observation points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%StressChange is the stress on the points from the event.                       Sxx Syy Sxy
	%TotalStress is the remote stress and stress change from the event added.       Sxx Syy Sxy
    %StressReg is the remote stress                                                 Sxx Syy Sxy
    %InducedDisplacements is the displacements Ux Uy from the dislocations          Ux Uy
[StressChange,TotalStress,StressReg,InducedDisplacements]=StressDispOnSurroundingPoints(X,Y,xe,ye,HalfLength,Sxx,Syy,Sxy,nu,E,halfspace,ShearDisp,TensileDisp,Beta);





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
Pxx=Sxx(1,1);
Pyy=Syy(1,1);
Pxy=Sxy(1,1);
%Drawing Figs
Sxx=TotalStress(:,1);Syy=TotalStress(:,2);Sxy=TotalStress(:,3);


Sxx = reshape(Sxx,dimx,dimy);
Syy = reshape(Syy,dimx,dimy);
Sxy = reshape(Sxy,dimx,dimy);

%Adding vertical Stress profile: 
SxxGrav = (nu/(1-nu))*P*g*Y;	%Equation 1b, Martel and Muller
SyyGrav = P*g*Y;		%Equation 1a, Martel and Muller
SxyGrav =zeros(size(SxxGrav));

Sxx=SxxGrav+Sxx;
Syy=SyyGrav+Syy;
Sxy=SxyGrav+Sxy;

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
xmin = floor(min(sigx(isfinite(sigx(:)))));
xmax = ceil(max(sigx(isfinite(sigx(:)))));
xindex = xmin:LvlstpX:xmax;

ymin = floor(min(sigxy(isfinite(sigy(:)))));
ymax = ceil(max(sigxy(isfinite(sigy(:)))));
yindex = ymin:LvlstpY:ymax;

xymin = floor(min(sigxy(isfinite(sigxy(:)))));
xymax = ceil(max(sigxy(isfinite(sigxy(:)))));
xyindex = xymin:LvlstpXY:xymax;

%Plotting figures from the analytical solution
figure,subplot(2,3,4),[C,h] =contourf(X,Y,sigx);
caxis(caxissxx);colormap(cmap2),colorbar;
hold on
contour(X,Y,sigx,xindex,'k');
set (h, 'levelstep', LvlstpX);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x Analytical'); 
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','k','LineStyle','--')
ax = gca;
set( ax,'YLim', [min(Y(:)) 0]);
set( ax,'XLim', [0 max(X(:))]);
%ax.YLim = [min(Y(:)) 0];
%ax.XLim = [0 max(X(:))];
hold off

%Syy
subplot(2,3,5),[C,h] =contourf(X,Y,sigy);
caxis(caxissyy);colormap(cmap2),colorbar;
hold on
contour(X,Y,sigy,yindex,'k');
set (h, 'levelstep', LvlstpY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y Analytical');
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','k','LineStyle','--')
ax = gca;
set( ax,'YLim', [min(Y(:)) 0]);
set( ax,'XLim', [0 max(X(:))]);
%ax.YLim = [min(Y(:)) 0];
%ax.XLim = [0 max(X(:))];
hold off

%Sxy
subplot(2,3,6),[C,h] =contourf(X,Y,sigxy);
caxis(caxissxy);colormap(cmap2),colorbar;
hold on
contour(X,Y,sigxy,xyindex,'k');
set (h, 'levelstep', LvlstpXY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y Analytical');
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','k','LineStyle','--')
ax = gca;
set( ax,'YLim', [min(Y(:)) 0]);
set( ax,'XLim', [0 max(X(:))]);
%ax.YLim = [min(Y(:)) 0];
%ax.XLim = [0 max(X(:))];
hold off



%Plotting figures from the dislocation solution
%Sxx
subplot(2,3,1),[C,h] =contourf(X,Y,Sxx);
caxis(caxissxx);colormap(cmap2),colorbar;
hold on
contour(X,Y,Sxx,xindex,'k');
set (h, 'levelstep', LvlstpX);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x Numerical');
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','k','LineStyle','--')
ax = gca;
set( ax,'YLim', [min(Y(:)) 0]);
set( ax,'XLim', [0 max(X(:))]);
%ax.YLim = [min(Y(:)) 0];
%ax.XLim = [0 max(X(:))];
hold off

%Syy
subplot(2,3,2),[C,h] =contourf(X,Y,Syy);
caxis(caxissyy);colormap(cmap2),colorbar;
hold on
contour(X,Y,Syy,yindex,'k');
set (h, 'levelstep', LvlstpY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y Numerical');
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','k','LineStyle','--')
ax = gca;
set( ax,'YLim', [min(Y(:)) 0]);
set( ax,'XLim', [0 max(X(:))]);
%ax.YLim = [min(Y(:)) 0];
%ax.XLim = [0 max(X(:))];
hold off

%Sxy
subplot(2,3,3),[C,h] =contourf(X,Y,Sxy);
caxis(caxissxy);colormap(cmap2),colorbar;
hold on
contour(X,Y,Sxy,xyindex,'k');
set (h, 'levelstep', LvlstpXY);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y Numerical');
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','k','LineStyle','--')
ax = gca;
set( ax,'YLim', [min(Y(:)) 0]);
set( ax,'XLim', [0 max(X(:))]);
%ax.YLim = [min(Y(:)) 0];
%ax.XLim = [0 max(X(:))];
hold off
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Printing results (residual)
SxxRes=Sxx-sigx;
bad=isnan(SxxRes);SxxRes=SxxRes(~bad);%removing nans to get proper mean
%disp('MeanResidualSxx');disp(mean(SxxRes(:)))
fprintf('MeanResidualSxx %i.\n',mean(SxxRes(:))) %prints on one line unlike 'disp
SyyRes=Syy-sigy;
bad=isnan(SyyRes);SyyRes=SyyRes(~bad);%removing nans to get proper mean
%disp('MeanResidualSyy');disp(mean(SyyRes(:)))
fprintf('MeanResidualSyy %i.\n',mean(SyyRes(:))) %prints on one line unlike 'disp
SxyRes=Sxy-sigxy;
bad=isnan(SxyRes);SxyRes=SxyRes(~bad);%removing nans to get proper mean
%disp('MeanResidualSxy');disp(mean(SxyRes(:)))
fprintf('MeanResidualSxy %i.\n',mean(SxyRes(:))) %prints on one line unlike 'disp

P=[mean(SxxRes(:)),mean(SyyRes(:)),mean(SxyRes(:))];
if any(P>0.01)
    error('Error in surrounding grid is too large. This no longer matches analytical solutions, check the comparative images')
else
    disp('Everything looks good, checks the mean Sxx Syy Sxy residual of all points is below 0.01')
end
%disp('Compare the two plots of Sxx+Syy induced by the valley without graviational loading, ignore features above valley surface')
