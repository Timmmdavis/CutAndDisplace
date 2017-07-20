% Test: Constant shear dislocation: Comparison of stresses and
% displacements in an observation grid surrounding a flat lying dislocation
% with a constant shearing discontinuity (Glide dislocation).
% 
% Solution: Stresses induced in the surrounding medium are given by the
% formulas for a glide dislocation from Barber, J., 2010. Elasticity (ed.,
% Vol. 172).(G. Gladwell, Ed.) Michigan. Displacements are found using
% equations 8.63-8.37 from Pollard, D.D. and Fletcher, R.C.,
% 2005.Fundamentals of structural geology. Cambridge University Press.
% (these displacement equations are originally from Weertman and Weertman
% 1964) . This solution uses the superposition of two analytical solutions
% for dislocations extending to infinity from a point in a whole space. The
% residual stress and displacement at each observation point is calculated
% and this throws an error if this is too high.
% 
% Proof: Each element in the 2d TWODD code is essentially this solution
% underneath but with coordinate transformations to get Cartesian stress
% components that are dependant on the elements orientation. Comparing the
% code to the analytical solution gives a good idea that the basic TWODD
% maths is correctly coded for a constant dislocation. its a basic test and
% only really shows the calculations are correct when not having to use the
% complex TWODD maths for the coordinate transformations.

%   Copyright 2017, Tim Davis, The University of Aberdeen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir


clear;close all
cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

		
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


%Single flat user defined fracture
x=linspace(-1,1,2);                        %List of x points
% x = fliplr(x);
y= zeros(1,numel(x));%x*0.5;
Pointsxy=[x;y]';
% mystruct.line1=(1:(length(Pointsxy(:,1)))-1);
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
	
mu=10000;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
%lambda=4000;%4000;  %Lame's constant  0 perfectly compressible like cork, Infinate for incompressible material like rubber.
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;

	%Equations you can use to calculate elastic parameters that you need above. 
	
E = mu*(2*(1+nu)) ;                     %Young's Modulus
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard



    %Reshaping the list of XY into usable variables
[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng,NUM,XBEG,XEND,YBEG,YEND]...
 = CreateElements2d( Pointsxy,mystruct );

    %Plotting figure of the the loaded fracture/fractures
hold on    
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
nx=cos(NormAng);
ny=cos((pi/2)-NormAng);
quiver(xe,ye,nx,ny)
title('fractures'), xlabel('x'), ylabel('y');axis equal
%Fixed elements shown as blue
line([Points(logical(Fdisp),1)';Points(logical(Fdisp),2)'],[Points(logical(Fdisp),3)';Points(logical(Fdisp),4)'],'color','b')
hold off

   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface
	%Runs a constant slip on every element. A mixture of two slips will give an inclined slip.
	%Output 'stresses' are blank arrays of 0's as this option is not driven by boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
cc=zeros(NUM,1); 
ShearDisp  = 0.0001;  ShearDisp	= cc+(ShearDisp*-1);   %For fractures lying along Xaxis positive stress drives right lateral movement, along Y is left lat. Pollard fig6.13c
TensileDisp =0;       TensileDisp	= cc+(TensileDisp*-1);  %Positive = extensional movement
Sxx = 0;  Syy = 0;  Sxy = 0;   	
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

spacing=0.1;
%minx=-4; maxx=4;
minx=-4+0.05; maxx=4+0.05; %Moving off dislocation slightly, MATLAB handles it on the dislocation but octave pulls a hissy fit
[X,Y] = meshgrid(minx:spacing:maxx); %large grid
dimx = length(X(:,1));
dimy = length(Y(:,1));

    %Replacing XY of any obs points for Nan if these lie on top of the end points of a
    %dislocation, without this it can result in spurious values. 
[ X,Y  ] = NullPointsLyingOnElement( X,Y,XBEG,YBEG,XEND,YEND,NUM );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses  on dispersed XYZ
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
%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Drawing figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Drawing Figs
Sxx = reshape(StressChange(:,1),dimx,dimy);
Syy = reshape(StressChange(:,2),dimx,dimy);
Sxy = reshape(StressChange(:,3),dimx,dimy);

Ux=InducedDisplacements(:,1);Uy=InducedDisplacements(:,2);
Ux = reshape(Ux,dimx,dimy);Uy = reshape(Uy,dimx,dimy);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating Barber Dislocation
%%

PR = 0.25;
%k = (3-PR)/(1+PR); %Kolosov's constant for plane stress
k = (3-(4*PR)); %Kolosov's constant for plane strain
U = mu; %100; %ShearMod
%spacing=0.05;     %Set above
%minx=-4; maxx=4;  %Set above
[x,y] = meshgrid(minx:spacing:maxx); %large grid
a = 1;  % Unit half-length displacement discontinuity
%B=0.0001; % Length of the Burger's vector.
B=-ShearDisp; %Convention

[SxxBarb,SyyBarb,SxyBarb,uxBarb,uyBarb]=Mode2_Timoshenko_Barber(k,U,minx,maxx,spacing,a,B,PR);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cleaning data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Removing spurios values for the code that are right next to the crack. The Pollard analytical
% %solution is far more accurate at this at the tips. 
% bad=abs(X)<1+(spacing*1);
% bad2=abs(Y)<1*spacing;
% bad3=bad+bad2;
% bad4=bad3==2;
% Sxx(bad4)=nan;Syy(bad4)=nan;Sxy(bad4)=nan;
% %Sxx_An(bad4)=nan;Syy_An(bad4)=nan;Sxy_An(bad4)=nan;
% Ux(bad4)=nan;Uy(bad4)=nan;
% %u(bad4)=nan;v(bad4)=nan;

%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
figure; quiver(X,Y,uxBarb,uyBarb);
xlabel('x'); ylabel('y'); axis('equal'); title('Displacement dislocation solution');
figure;quiver(X,Y,Ux,Uy);
xlabel('x'); ylabel('y'); axis('equal'); title('Displacement code solution');

%Octave ContourF doesn't like 'nans', replacing with '-infs'
%MATLAB unsuprisingly handles this fine
Sxx(isnan(Sxx)) = -inf;    SxxBarb(isnan(SxxBarb)) = -inf;    
Syy(isnan(Syy)) = -inf;    SyyBarb(isnan(SyyBarb)) = -inf;    
Sxy(isnan(Sxy)) = -inf;    SxyBarb(isnan(SxyBarb)) = -inf;  

%Setting contour steps, the steps set below are the same as the figures in the Martel
%2000 paper.
LvlstpX =0.1;
LvlstpY =0.1;
LvlstpXY=0.1;

%Setting axis limits so supurious values do not change the overall trend. 
%Only looking at finite values, No axis scaling to infs or nans. 
caxissxx=[min(SxxBarb(isfinite(SxxBarb(:)))),max(SxxBarb(isfinite(SxxBarb(:))))];caxissxx=caxissxx/4;
caxissyy=[min(SyyBarb(isfinite(SyyBarb(:)))),max(SyyBarb(isfinite(SyyBarb(:))))];caxissyy=caxissyy/4;
caxissxy=[min(SxyBarb(isfinite(SxyBarb(:)))),max(SxyBarb(isfinite(SxyBarb(:))))];caxissxy=caxissxy/4;
%Rounding values, rounds the axis values to the set step.
Xrnd=1/LvlstpX;
Yrnd=1/LvlstpY;
XYrnd=1/LvlstpXY;
caxissxx = [-0.5 0.5];% (round(caxissxx*Xrnd))/Xrnd; %Rounding to nearest 0.2
caxissyy = [-0.5 0.5];%(round(caxissyy*Yrnd))/Yrnd; %Rounding yy to nearest 0.5 (same as levelstep)
caxissxy = [-0.5 0.5];%(round(caxissxy*XYrnd))/XYrnd; %Rounding to nearest 0.02


%Creating specific contour levels, See MATLABS ''highlight-specific-contour-levels'' Doc. 
% xmin = floor(min(SxxBarb(isfinite(SxxBarb(:)))));
% xmax = ceil(max(SxxBarb(isfinite(SxxBarb(:)))));
% xindex = xmin:LvlstpX:xmax;
% 
% ymin = floor(min(SyyBarb(isfinite(SyyBarb(:)))));
% ymax = ceil(max(SyyBarb(isfinite(SyyBarb(:)))));
% yindex = ymin:LvlstpY:ymax;
% 
% xymin = floor(min(SxyBarb(isfinite(SxyBarb(:)))));
% xymax = ceil(max(SxyBarb(isfinite(SxyBarb(:)))));
% xyindex = xymin:LvlstpXY:xymax;

fntsz=25;
titlesz=23;
%Sxx
figure,
subplot(2,3,1),[C,h]= contourf(X,Y,Sxx); caxis(caxissxx)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x numerical'),colorbar; 
ChangeFontSizes(titlesz,fntsz);set (h, 'levelstep', LvlstpX);
set(colorbar,'YTick',-0.5:0.25:0.50)

pause(3); %Octave has issues with timing drawing this, this slows it down enough so it doesn't fail

subplot(2,3,4),[C,h]= contourf(X,Y,SxxBarb);
caxis(caxissxx)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x analytical'),colorbar; 
ChangeFontSizes(titlesz,fntsz);set (h, 'levelstep', LvlstpX);
set(colorbar,'YTick',-0.5:0.25:0.50)

%Syy
subplot(2,3,2),[C,h]= contourf(X,Y,Syy);
caxis(caxissyy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y numerical'),colorbar; 
ChangeFontSizes(titlesz,fntsz);set (h, 'levelstep', LvlstpY);
set(colorbar,'YTick',-0.5:0.25:0.50)

subplot(2,3,5),[C,h]= contourf(X,Y,SyyBarb);
caxis(caxissyy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y analytical'),colorbar; 
ChangeFontSizes(titlesz,fntsz);set (h, 'levelstep', LvlstpY);
set(colorbar,'YTick',-0.5:0.25:0.50)

%Sxy
subplot(2,3,3),[C,h]= contourf(X,Y,Sxy);
caxis(caxissxy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y numerical'),colorbar; 
ChangeFontSizes(titlesz,fntsz);set (h, 'levelstep', LvlstpXY);
set(colorbar,'YTick',-0.5:0.25:0.50)

subplot(2,3,6),colormap(cmap2),[C,h]= contourf(X,Y,SxyBarb);%
caxis(caxissxy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y analytical'),colorbar; 
ChangeFontSizes(titlesz,fntsz);set (h, 'levelstep', LvlstpXY);
set(colorbar,'YTick',-0.5:0.25:0.50)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Printing results (residual)
SxxRes=Sxx-SxxBarb;
%disp('MaxResidualSxx');disp(max(SxxRes(:)))
%fprintf('MaxResidualSxx %i.\n',max(abs(SxxRes(:)))) %prints on one line unlike 'disp
SyyRes=Syy-SyyBarb;
%disp('MaxResidualSyy');disp(max(SyyRes(:)))
%fprintf('MaxResidualSyy %i.\n',max(abs(SyyRes(:)))) %prints on one line unlike 'disp
SxyRes=Sxy-SxyBarb;
%disp('MaxResidualSxy');disp(max(SxyRes(:)))
%fprintf('MaxResidualSxy %i.\n',max(abs(SxyRes(:)))) %prints on one line unlike 'disp

DispRes=(uxBarb-Ux)+(Uy-uyBarb);
%disp('MaxResidualDisp');disp(max(DispRes(:)))
fprintf('MaxResidualDisp %i.\n',max(abs(DispRes(:)))) %prints on one line unlike 'disp

%Error Percent
Error_percSxy=((SxyBarb-Sxy)./SxyBarb)*100; MaxError_percSxy=max(max(Error_percSxy));
Error_percSxx=((SxxBarb-Sxx)./SxxBarb)*100; MaxError_percSxx=max(max(Error_percSxx));
Error_percSyy=((SyyBarb-Syy)./SyyBarb)*100; MaxError_percSyy=max(max(Error_percSyy));
fprintf('BiggestPercentErrorSxx %% %i.\n',MaxError_percSxx); %prints on one line unlike 'disp
fprintf('BiggestPercentErrorSyy %% %i.\n',MaxError_percSyy);
fprintf('BiggestPercentErrorSxy %% %i.\n',MaxError_percSxy);

P=[MaxError_percSxx,MaxError_percSyy,MaxError_percSxy];

if any(P>1e-9)
    error('Error in surrounding grid for stress is too large. Something has changed in Twodd Calculation. This no longer matches analytical solutions')
else
    disp('Everything looks good, checks the max % error between analytical and TWODD Sxx Syy Sxy of all points. Values above 1e-9% flagged')
end

