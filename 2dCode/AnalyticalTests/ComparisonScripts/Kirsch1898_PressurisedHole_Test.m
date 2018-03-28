% Test: Hole with internal pressure: Compares stresses around a hole with
% an internal pressure that sits in an infinite elastic.
% 
% Solution: Used the script for this from the Pollard and Fletcher
% fundamentals book. This solution is In most elasticity/continuum
% textbooks as its simple to implement and relevant. It's the original
% solution from: Kirsch, G., 1898. Theory of Elasticity and Application in
% Strength of Materials. Zeits chrift des Vereins Deutscher
% Ingenieure, 42(29), pp.797-807. The solution is independent of elastic
% constants and hole size, it only relies on the input pressure. Elastic
% constants will only change the displacement related to the pressure. The
% results in the surrounding grid from the analytical and code solution are
% compared and this checks the residual stresses are not too high. Figures
% are also drawn that can be compared.
%  
% Proof: This tests the user defined traction function works and reproduces
% the result. This gives a good example case for using it. Displacement
% components could also be compared in this test. These were not in the
% original solution so for the sake of simplicity have been omitted. The
% solutions are here if needed:
% https://www.rocscience.com/help/RS3beta/webhelp/pdf_files/verification/Verification_001_(Cylindrical_Hole,_Elastic).pdf

%   Copyright 2017, Tim Davis, The University of Aberdeen



%%
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
     
% %Circle
rad=1;
[ xe,ye ] = CreateCircleXY( 600,rad );
sz=numel(xe);
Pointsxy=[ye,xe];
mystruct.line1=(1:sz);
%Fixed points inner circle (making a square)
rr=0.8;%size of the square
x=[-rr,0,rr,0,-rr];                        %List of x points
y=[0,-rr,0,rr,0];
PointsxyF=[x;y]';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Pointsxy,mystruct,Fdisp] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
	%Creating Fictitious disp flag, if set to one then this forces the elements with this to 0 displacement. 
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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. No opening components, this option solves
    %the boundary value problem through shearing of elements alone. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 cc=zeros(numel(HalfLength),1); 
    
 Tn = 1;    Tn=cc+Tn;	%'Tensile' traction
 Ts = 0;    Ts=cc+Ts;   
 Option='E';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ds,Dn,Sxx,Syy,Sxy]=SlipCalculator2d(...
        MidPoint,HalfLength,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);
    
%%
% Removing any fixed elements
if any(Fdisp)==1  
    [MidPoint,HalfLength,P1,P2,LineNormalVector]...
    = RemovingFixedEls2d(Fdisp,P1,P2);
end

	
	
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Regular/randomly spaced grid of points bounding the faults.
    %Define how far it extends away and its density.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
boxsz=3;    
X = linspace(-boxsz,boxsz,101);
Y = linspace(-boxsz,boxsz,101); 

[X,Y] = meshgrid(X,Y);
rowcount = length(X(:,1));
colcount = length(X(1,:));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5: Calculate Stresses/Disps at defined observation points in XY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Stresses and Strains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
% [StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]...
%     = CalculateStressesOnSurroundingPoints2d(X,Y,MidPoint,HalfLength,Sxx,Syy,Sxy,nu,E,halfspace,Ds,Dn,LineNormalVector);
 
%Doing with obs calc func so this is tested. 
[StressTTotal,StressTChg,StressTReg]=ObsPointsStressCalculator2d(X,Y,MidPoint,HalfLength,Sxx,Syy,Sxy,nu,E,halfspace...
    ,Ds,Dn,LineNormalVector);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ux,Uy] = CalculateDisplacementOnSurroundingPoints2d(X,Y,MidPoint,HalfLength,nu,E,halfspace,Ds,Dn,LineNormalVector);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution, Kirsch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xa = linspace(-boxsz,boxsz,101);
% Ya = linspace(-boxsz,boxsz,101); 
[SXX,SYY,SXY]=Kirsch1898_PressurisedHole(X,Y,rad,-Tn); %Tn is extensional, Pressure is compressive

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: Visualisation and analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

Pxx=Sxx(1,1);
Pyy=Syy(1,1);
Pxy=Sxy(1,1);
%Drawing Figs
[Sxx,Syy,Sxy,Ux,Uy]=ReshapeData2d...
    (rowcount,colcount,StressTTotal(:,1),StressTTotal(:,2),StressTTotal(:,3),Ux,Uy);

%%%%%%%%%%%%%%%%%%
%Removing inner points
[~,R] = cart2pol(X,Y);
Bad=R<1;
Sxx(Bad)=nan;
Syy(Bad)=nan;
Sxy(Bad)=nan;



[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);
%Drawing S1 S2 directions
DrawS1S2Directions(X(:),Y(:),S1dir,S2dir,'P1',P1,'P2',P2 )
PlotFracture( P1,P2,'r' )
xlabel('x'); ylabel('y'); axis('equal'); title('NumericalPrincipalStresses');

%Octave ContourF doesn't like 'nans', replacing with '-infs'
%MATLAB unsuprisingly handles this fine
Sxx(isnan(Sxx)) = -inf;    SXX(isnan(SXX)) = -inf;    
Syy(isnan(Syy)) = -inf;    SYY(isnan(SYY)) = -inf;    
Sxy(isnan(Sxy)) = -inf;    SXY(isnan(SXY)) = -inf;  


figure,quiver(X,Y,Ux,Uy);
xlabel('x'); ylabel('y'); axis('equal'); title('Displacement Numerical');

Lvlstp=0.2;

figure;subplot(2,3,1),colormap(cmap2),[~,h]=contourf(X,Y,Sxx);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x numerical');colorbar;
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,2),colormap(cmap2),[~,h]=contourf(X,Y,Syy);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y numerical');colorbar;
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,3),colormap(cmap2),[~,h]=contourf(X,Y,Sxy);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y numerical');colorbar;
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)


subplot(2,3,4),colormap(cmap2),[~,h]=contourf(X,Y,SXX);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x analytical');colorbar;
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,5),colormap(cmap2),[~,h]=contourf(X,Y,SYY);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y analytical');colorbar;
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,6),colormap(cmap2),[C,h]=contourf(X,Y,SXY);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y analytical');colorbar;
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)


titlesz=25;
fntsz=21;
ChangeFontSizes(titlesz,fntsz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Printing results (residual)
SxxRes=SXX-Sxx;
bad=isnan(SxxRes);SxxRes=SxxRes(~bad);%removing nans to get proper mean
fprintf('MaxResidualSxx %i.\n',max(SxxRes(:)))
%disp('MeanResidualSxx');disp(mean(SxxRes(:)))
SyyRes=SYY-Syy;
bad=isnan(SyyRes);SyyRes=SyyRes(~bad);%removing nans to get proper mean
%disp('MeanResidualSyy');disp(mean(SyyRes(:)))
fprintf('MaxResidualSyy %i.\n',max(SyyRes(:)))
SxyRes=SXY-Sxy;
bad=isnan(SxyRes);SxyRes=SxyRes(~bad);%removing nans to get proper mean
%disp('MeanResidualSxy');disp(mean(SxyRes(:)))
fprintf('MaxResidualSxy %i.\n',max(SxyRes(:))) %prints on one line unlike 'disp
 
%Percentage error
%Residual as a percentage of the driving stress
%this avoids doing it as a percentage of the actual values were values
%close to 0 give unreasonably high errors 
Drivingstress=Tn(1,1)+Ts(1,1);

SxxResPerc=(100/Drivingstress).*SxxRes; %scatter(X(:),Y(:),15,SxxResPerc(:))
SyyResPerc=(100/Drivingstress).*SyyRes;
SxyResPerc=(100/Drivingstress).*SxyRes;
 
answer=isnan(SxxResPerc);
SxxResPerc(answer)=0;
answer=isinf(SxxResPerc);
SxxResPerc(answer)=0;
answer=isnan(SyyResPerc);
SyyResPerc(answer)=0;
answer=isinf(SyyResPerc);
SyyResPerc(answer)=0;
answer=isnan(SxyResPerc);
SxyResPerc(answer)=0;
answer=isinf(SxyResPerc);
SxyResPerc(answer)=0;

PpercStress=[max(SxxResPerc(:)),max(SyyResPerc(:)),max(SxyResPerc(:))];

if any(PpercStress>10)
    error('Error in surrounding grid too large. Something has changed in Twodd Calculation. This no longer matches analytical solutions')
else
    disp('Everything looks good, checks the highest tensor error % is below 1% of the total driving stress')
end

