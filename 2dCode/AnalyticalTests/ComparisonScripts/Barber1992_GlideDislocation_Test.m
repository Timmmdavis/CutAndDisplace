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
% only really shows the calculations are correct for the TWODD maths
% without any the coordinate transformations.


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
%STEP 1: import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Single flat user defined fracture
x=linspace(-1,1,2);                        %List of x points
% x = fliplr(x);
y= zeros(1,numel(x));%x*0.5;
Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1))));
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms

%Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
    Fdisp=[];
end


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
Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;

% Defining elastic constants
mu=10000;           %Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

%Reshaping the list of XY into usable variables
[MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp  ]...
 = CreateElements2d( Pointsxy,mystruct,Fdisp );

%Drawing fracture
%Normal els as red.
PlotFracture( P1,P2,'r' )
%Fixed elements shown as blue
id=find(Fdisp); 
PlotFracture(P1(id,:),P2(id,:),'b' )
%Titles
title('fractures'), xlabel('x'), ylabel('y');axis equal
%Drawing Normals
quiver(MidPoint(:,1),MidPoint(:,2),LineNormalVector(:,1),LineNormalVector(:,2))



   
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface.
	%Runs a constant slip on every element. Output 'stresses' are blank
	%arrays of 0's as this option uses displacement boundary conditions.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
cc=zeros(numel(MidPoint(:,1)),1); 
Ds  = -0.0001;  Ds	= cc+(Ds);   
Dn =0;       Dn	= cc+(Dn);
Sxx = 0;  Syy = 0;  Sxy = 0;   	
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

spacing=0.1;
%Moving off dislocation slightly, MATLAB handles it on the dislocation but octave pulls a hissy fit
minx=-4+0.05; maxx=4+0.05; 
[X,Y] = meshgrid(minx:spacing:maxx); %large grid
[dimx,dimy] = size(X);  




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
= CalculateStressesOnSurroundingPoints2d(X,Y,MidPoint,HalfLength...
,Sxx,Syy,Sxy,nu,E,halfspace,Ds,Dn,LineNormalVector);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d...
(X,Y,MidPoint,HalfLength,nu,E,halfspace,Ds,Dn,LineNormalVector);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%Reshape data
[Sxx,Syy,Sxy,Ux,Uy]=ReshapeData2d(dimx,dimy,StressTChg(:,1),StressTChg(:,2),StressTChg(:,3),Ux,Uy);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating Barber Dislocation
%%

nu = 0.25;
%k = (3-PR)/(1+PR); %Kolosov's constant for plane stress
k = (3-(4*nu));     %Kolosov's constant for plane strain
U = mu;             %ShearMod
[x,y] = meshgrid(minx:spacing:maxx); 
a = 1;              %Unit half-length displacement discontinuity
b=-Ds;              %Convention


[x,y,SxxBarb,SyyBarb,SxyBarb,uxBarb,uyBarb]=Barber1992_GlideDislocation(k,U,x,y,a,b,nu);


%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
figure; quiver(x,y,uxBarb,uyBarb);
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
caxissxx = [-0.5 0.5];
caxissyy = [-0.5 0.5];
caxissxy = [-0.5 0.5];


%Sxx
figure,

subplot(2,3,1),[~,h]= contourf(X,Y,Sxx); caxis(caxissxx)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x numerical'),colorbar; 
set (h, 'levelstep', LvlstpX);
set(colorbar,'YTick',-0.5:0.25:0.50)

pause(3); %Octave has issues with timing drawing this, this slows it down enough so it doesn't fail

subplot(2,3,4),[~,h]= contourf(X,Y,SxxBarb);
caxis(caxissxx)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x analytical'),colorbar; 
set (h, 'levelstep', LvlstpX);
set(colorbar,'YTick',-0.5:0.25:0.50)

%Syy
subplot(2,3,2),[~,h]= contourf(X,Y,Syy);
caxis(caxissyy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y numerical'),colorbar; 
set (h, 'levelstep', LvlstpY);
set(colorbar,'YTick',-0.5:0.25:0.50)

subplot(2,3,5),[~,h]= contourf(X,Y,SyyBarb);
caxis(caxissyy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y analytical'),colorbar; 
set (h, 'levelstep', LvlstpY);
set(colorbar,'YTick',-0.5:0.25:0.50)

%Sxy
subplot(2,3,3),[~,h]= contourf(X,Y,Sxy);
caxis(caxissxy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y numerical'),colorbar; 
set (h, 'levelstep', LvlstpXY);
set(colorbar,'YTick',-0.5:0.25:0.50)

subplot(2,3,6),colormap(cmap2),[C,h]= contourf(X,Y,SxyBarb);%
caxis(caxissxy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y analytical'),colorbar; 
set (h, 'levelstep', LvlstpXY);
set(colorbar,'YTick',-0.5:0.25:0.50)


titlesz=25;
fntsz=21;
ChangeFontSizes(titlesz,fntsz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Printing results (residual)
SxxRes=Sxx-SxxBarb;
SyyRes=Syy-SyyBarb;
SxyRes=Sxy-SxyBarb;


DispRes=(uxBarb-Ux)+(Uy-uyBarb);
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

