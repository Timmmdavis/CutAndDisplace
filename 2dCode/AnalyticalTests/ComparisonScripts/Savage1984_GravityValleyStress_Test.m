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
%STEP AA: Analytical solution, Savage 1984.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for Savage code, change below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tectonic Stress
ts=0; 

%Weight and elastic constants
P=1; 	%Density of elastic material 
g=1;	%Acceleration due to gravity
rg= P*g; %Weight (Dens*AccelGrav)
nu=0.25;

%A and B describing slope, See paper for diagram
a=2;
b=-1;

%Setting up grid points U and V, these are mapped to a different location
%later
u = linspace(0,4,50);
v = linspace(-0.3265,-4,46); %Points slightly below 0 as we do not want observation points on the elements for the twodd code 
[u,v] = meshgrid(u,v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sigx,sigy,sigxy,xobs,yobs] = Savage1984_GravityValleyStress(ts,rg,nu,a,b,u,v);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
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
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
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

[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;

P=P; 	 %Density of elastic material 
g=g;	 %Acceleration due to gravity
tects=0; %Tectonic stress (xx)
% % Creating gravitational stress on each elements centre, the
% % ground surface was originally at 0m. 
Sxx = ((nu/(1-nu))*P*g*MidPoint(:,2))+tects;	%Equation 1b, Martel and Muller
Syy = P*g*MidPoint(:,2);		%Equation 1a, Martel and Muller
Sxy = zeros(numel(MidPoint(:,1)),1);  %Positive is right lateral movement
Option='B'; %Valley surface deformation after erosion

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ds,Dn,~,Syy,Sxy]=SlipCalculator2d(...
    MidPoint,HalfLength,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);

%The resultant stresses in the body from these dislocations is taken away
%from a calculated far field stress that pervades the body, this far field stress is calculated using the constants P&g. 
Sxx=tects; %Tectonic stress (xx)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

X=xobs;
Y=yobs;
dimx = length(X(:,1));
dimy = length(X(1,:));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
    = CalculateStressesOnSurroundingPoints2d(X,Y,MidPoint,HalfLength,Sxx,Syy,Sxy,nu,E,halfspace,Ds,Dn,LineNormalVector);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d(X,Y,MidPoint,HalfLength,nu,E,halfspace,Ds,Dn,LineNormalVector);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Sxy,Ux,Uy]=ReshapeData2d(dimx,dimy,StressTTotal(:,1),StressTTotal(:,2),StressTTotal(:,3),Ux,Uy);


%Adding vertical Stress profile: 
SxxGrav = (nu/(1-nu))*P*g*Y;	%Equation 1b, Martel and Muller
SyyGrav = P*g*Y;                %Equation 1a, Martel and Muller
SxyGrav =zeros(size(SxxGrav));

Sxx=SxxGrav+Sxx;
Syy=SyyGrav+Syy;
Sxy=SxyGrav+Sxy;

%Setting contour steps, the steps set below are the same as the figures in
%the Martel (2000) paper.
LvlstpX=0.2;
LvlstpY=0.5;
LvlstpXY=0.02;

%Axis limits=
yl=[min(Y(:)) 0];
xl=[0 max(X(:))];


%Plotting figures from the analytical solution
figure;
%Sxx
subplot(2,3,4); hold on;
PlotFracture( P1,P2,'k' )
DrawContourFPlot2dWithContours( X,Y,sigx,LvlstpX,cmap2 );
title('\sigma_x_x Analytical');
xlim(xl);ylim(yl)
hold off

%Syy
subplot(2,3,5); hold on;
PlotFracture( P1,P2,'k' )
DrawContourFPlot2dWithContours( X,Y,sigy,LvlstpY,cmap2 );
title('\sigma_y_y Analytical');
xlim(xl);ylim(yl)
hold off

%Sxy
subplot(2,3,6); hold on;
PlotFracture( P1,P2,'k' )
DrawContourFPlot2dWithContours( X,Y,sigxy,LvlstpXY,cmap2 );
title('\sigma_x_y Analytical');
xlim(xl);ylim(yl)
hold off

%Plotting figures from the dislocation solution
%Sxx
subplot(2,3,1); hold on;
PlotFracture( P1,P2,'k' )
DrawContourFPlot2dWithContours( X,Y,Sxx,LvlstpX,cmap2 );
title('\sigma_x_x Numerical');
xlim(xl);ylim(yl)
hold off

%Syy
subplot(2,3,2); hold on;
PlotFracture( P1,P2,'k' )
DrawContourFPlot2dWithContours( X,Y,Syy,LvlstpY,cmap2 );
title('\sigma_y_y Numerical');
xlim(xl);ylim(yl)
hold off

%Sxy
subplot(2,3,3); hold on;
PlotFracture( P1,P2,'k' )
DrawContourFPlot2dWithContours( X,Y,Sxy,LvlstpXY,cmap2 );
title('\sigma_x_y Numerical');
xlim(xl);ylim(yl)
hold off

titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


SlipResPercx=mean(100-(abs((100./sigx(:)).*Sxx(:))));
SlipResPercy=mean(100-(abs((100./sigy(:)).*Syy(:)))); 
bad=abs((100./sigxy(:)))==inf;
SlipResPercxy=100-(abs((100./sigxy(:)).*Sxy(:))); 
SlipResPercxy(bad)=[];
SlipResPercxy=mean(SlipResPercxy);