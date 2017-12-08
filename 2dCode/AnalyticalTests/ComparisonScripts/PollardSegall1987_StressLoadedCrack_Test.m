% Test: Fracture loaded with remote stress: Comparison of slip
% distributions along a planar fracture loaded remotely. Also compares the
% induced stresses in a surrounding grid of observation points.
% 
% Solution: Pollard, D.D. and Segall, P., 1987. Theoretical displacements
% and stresses near fractures in rock: with applications to faults, joints,
% veins, dikes, and solution surfaces.'Fracture mechanics of rock,277(349),
% pp.277-349. This solution is an analytical solution using complex stress
% functions. It loads a planar 2d fracture again lying along the X axis
% with far field stress Sxx, Syy or Sxy. It calculates the slip
% distribution and induced Cartesian stress components for points
% surrounding the fracture. To quote Steve Martel on this complex stress
% method ''They also have a key advantage over numerical methods such as
% constant displacement-discontinuity boundary element techniques (Crouch
% and Starfield, 1974) in that the stresses and displacements can be
% calculated at points arbitrarily close to the patch walls.''. This can be
% seen in the results where the DDM methods match is less accurate closer
% to the fault. Especially near the tips where the slip distribution should
% increases exponentially and the DDM is blocky. This test compares both
% the residual slip and residual stress in the points for the maximum
% deviation from the analytical solution. This test drives displacement
% with a tensile stress but this is set up so the user could change to a
% driving shear stress is they need. Comparative figures are also plotted.
% 
% Proof: This gives a good idea that the influence matrix collation and
% linear equations are working correctly. It also allows for an estimation
% of sampling needed before the solution is accurate enough. There are
% plenty of papers on TWODDs accuracy already. As the workflow given by
% Crouch and Starfield has changed in my code and it's  a new code base it's
% a good proof of the basics.

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

%Single flat user defined fracture
%x=linspace(-1,1,400);
x=-1:0.001:1;                        %List of x points
y=zeros(1,numel(x));
Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1)))-1);
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



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. 
    %Choose to define in stress or strain, the stresses can be defined as
    %single values that are repeated for every element or defined as a
    %vector that is the stress at every element: i.e.  
    %Syy = linspace(0,1,NUM); %(linearly increasing for each element).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%

%Put to 1 to define the stresses defined in 'stress input' as strain values    
strain=0;   

    %%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
%Input Stress
Sxx = 0; 					
Syy = 0.01;
Sxy = 0; 
Option='B';


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ds,Dn,Sxx,Syy,Sxy]=SlipCalculator2d(...
    MidPoint,HalfLength,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,LineNormalVector,strain,Fdisp,Mu,Sf,Option);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Function that calls eq 8.35a and 8.35b from the Pollard/Segall paper
[OpeningOrSlip]=PollardSegall1987_FractureSlipProfile(mu,nu,Syy,Sxy,MidPoint(:,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and drawing figures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SlipDispRes=abs(OpeningOrSlip)-abs(Dn+Ds);
SlipDispResPerc=((abs(OpeningOrSlip)-abs(Dn+Ds))./abs(OpeningOrSlip))*100;

%disp('MaxResidualSlip');disp(max(SlipDispRes(:)));
fprintf('MaxResidualSlip %i.\n',max(SlipDispRes(:))) %prints on one line unlike 'disp
fprintf('MaxResidualSlipAsAPercentage %i.\n',max(SlipDispResPerc(:))) %prints on one line unlike 'disp
lnwidth=2.5;
plot(MidPoint(:,1),OpeningOrSlip,'g',MidPoint(:,1),SlipDispRes,'r','LineWidth',lnwidth);title('OpeningDistributionAlongCrack','LineWidth',lnwidth);
xe=MidPoint(:,1);
hold on
step=xe(1,1)-xe(2,1); stairs([xe(1,1)+step;xe;xe(end,end)-step],[0;(Dn+Ds);0],'m','LineWidth',lnwidth);
xlabel('Distance from crack centre')
ylabel({'Opening, analytical (green)','numerical (pink) residual(red)'})
xlim ([-1 1])

fntsz=18;
titlesz=22;
ChangeFontSizes(titlesz,fntsz);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Drawing a different way:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subsetting data - would be better to interpolate at points I want to draw
gd=xe>0;
DnP=Dn(gd);
xP=xe(gd);
N=8;
DnN = DnP(1:N:end,:);
xN = xP(1:N:end,:);

%trying interpolation
xN = linspace(0,0.9,10);
DnN = interp1(xe,Dn,xN,'spline');

lnwidth=2.5;
figure;
subplot(1,2,1), %create subplot with 2 parts (horz alligned)
hold on
title('Opening of crack walls')
xlabel('Distance from crack centre')
ylabel({'Crack wall displacements'})

%an- 
plot(abs(xe),OpeningOrSlip,'g','LineWidth',2);
%num-
scatter(abs(xN),abs(DnN),24,'k','filled');
%res-
plot(x(gd),SlipDispRes(gd),'r','LineWidth',lnwidth)
rectangle('Position',[0.9 0 0.1 0.0015]) %zoom rectangle
grid on
hold off

%Zoomed in plot
subplot(1,2,2), 
hold on
%an- 
plot(abs(xe),OpeningOrSlip,'g','LineWidth',2);
%num-
scatter(xP,abs(DnP),24,'k','filled');
%res-
plot(x(gd),SlipDispRes(gd),'r','LineWidth',lnwidth)
xlim ([0.9 1]);
ylim ([0 0.0015]);
grid on
hold off
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);

%getting % error for these inp pnts
[OpeningOrSlip]=PollardSegall1987_FractureSlipProfile(mu,nu,Syy,Sxy,xP,1);
PercentErrorOpeningInterpPnts=((100./-OpeningOrSlip).*DnP)-100;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
spacing=0.1; 
minx=-4; maxx=4;
[X,Y] = meshgrid(minx:spacing:maxx); 
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
%STEP BB: Analytical solution stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sxx_An,Syy_An,Sxy_An,u,v]=PollardSegall1987_FractureStressDispInBody(mu,nu,Syy,Sxy,X,Y);
[Sxx_An,Syy_An,Sxy_An,u,v]=ReshapeData2d(dimx,dimy,Sxx_An,Syy_An,Sxy_An,u,v);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cleaning data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Removing spurios values for the code that are right next to the crack. The Pollard analytical
%solution is far more accurate at this at the tips. 
bad=abs(X)<1+(spacing*1);
bad2=abs(Y)<1*spacing;
bad3=bad+bad2;
bad4=bad3==2;
SxxTwodd=StressTChg(:,1); 
SyyTwodd=StressTChg(:,2);
SxyTwodd=StressTChg(:,3);
%Reshape
[SxxTwodd,SyyTwodd,SxyTwodd,Ux,Uy]=ReshapeData2d(dimx,dimy,SxxTwodd,SyyTwodd,SxyTwodd,Ux,Uy);
SxxTwodd(bad4)=nan;SyyTwodd(bad4)=nan;SxyTwodd(bad4)=nan;
Ux(bad4)=nan;Uy(bad4)=nan;

figure,quiver(X,Y,Ux,Uy);
xlabel('x'); ylabel('y'); axis('equal'); title('displacement from code, central patch removed');
%Drawing Figs
Sxxd=Sxx; Sxx = SxxTwodd;
Syyd=Syy; Syy = SyyTwodd;
Sxyd=Sxy; Sxy = SxyTwodd;

%Setting contour steps, the steps set below are the same as the figures in the Martel
%2000 paper.
LvlstpX=0.001;
LvlstpY=0.001;
LvlstpXY=0.001;

%Setting axis limit's  so supurious values do not change the overall trend. 
%Only looking at finite values, No axis scaling to infs or nans. 
caxissxx=[min(Sxx(isfinite(Sxx(:)))),max(Sxx(isfinite(Sxx(:))))];
caxissyy=[min(Syy(isfinite(Syy(:)))),max(Syy(isfinite(Syy(:))))];
caxissxy=[min(Sxy(isfinite(Sxy(:)))),max(Sxy(isfinite(Sxy(:))))];
%Rounding values, rounds the axis values to the set step.
Xrnd=1/LvlstpX;
Yrnd=1/LvlstpY;
XYrnd=1/LvlstpXY;
caxissxx = (round(caxissxx*Xrnd))/Xrnd; %Rounding to nearest 0.2
caxissyy = (round(caxissyy*Yrnd))/Yrnd; %Rounding yy to nearest 0.5 (same as levelstep)
caxissxy = (round(caxissxy*XYrnd))/XYrnd; %Rounding to nearest 0.02
%Getting round linearly increasing numbers for the cbar
length=caxissxx(1)-caxissxx(2); %the distance between the max and min limit's 
vectoroflocsX=linspace(caxissxx(1),caxissxx(2),abs(length/LvlstpX)+1); %now using this length to create a spaced vector for labeling the colourbar
length=caxissyy(1)-caxissyy(2);
vectoroflocsY=linspace(caxissyy(1),caxissyy(2),abs(length/LvlstpY)+1);
length=caxissxy(1)-caxissxy(2);
vectoroflocsXY=linspace(caxissxy(1),caxissxy(2),abs(length/LvlstpXY)+1);

%Sxx
figure,        
subplot(2,3,1), colormap(cmap2);
[~,h]=contourf(X,Y,Sxx); 
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x numerical');
set (h, 'levelstep', LvlstpX);
colorbar(); cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', vectoroflocsX)
caxis([-0.005 0.005]);set(colorbar,'YTick',-0.005:0.0025:0.0050);

pause(3)%octave issue

subplot(2,3,4),
[~,h]=contourf(X,Y,Sxx_An); 
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x analytical');
set (h, 'levelstep', LvlstpX);
colorbar(); cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', vectoroflocsX)
caxis([-0.005 0.005]);set(colorbar,'YTick',-0.005:0.0025:0.0050);

%Syy
subplot(2,3,2),
[~,h]=contourf(X,Y,Syy); 
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y numerical');
set (h, 'levelstep', LvlstpX);
colorbar(); cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', vectoroflocsY)
caxis([-0.01 0.01]);

subplot(2,3,5),
[~,h]=contourf(X,Y,Syy_An); 
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y analytical');
set (h, 'levelstep', LvlstpY);
colorbar(); cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', vectoroflocsY)
caxis([-0.01 0.01]);

%Sxy
subplot(2,3,3),
[~,h]=contourf(X,Y,Sxy); 
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y numerical');
set (h, 'levelstep', LvlstpXY);
colorbar(); cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', vectoroflocsXY)
caxis(caxissxy); set(colorbar,'YTick',-0.005:0.0025:0.0050);


subplot(2,3,6);
[~,h]=contourf(X,Y,Sxy_An); 
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y analytical');
set (h, 'levelstep', LvlstpXY);
colorbar(); cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', vectoroflocsXY)
caxis(caxissxy); set(colorbar,'YTick',-0.005:0.0025:0.0050);

titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Printing results (residual)
SxxRes=SxxTwodd-Sxx_An;
%disp('MaxResidualSxx');disp(max(SxxRes(:)))
fprintf('MaxResidualSxx %i.\n',max(abs(SxxRes(:))))
SyyRes=SyyTwodd-Syy_An;
%disp('MaxResidualSyy');disp(max(abs(SyyRes(:)))
fprintf('MaxResidualSyy %i.\n',max(abs(SyyRes(:))))
SxyRes=SxyTwodd-Sxy_An;
%disp('MaxResidualSxy');disp(max(abs(SxyRes(:)))
fprintf('MaxResidualSxy %i.\n',max(abs(SxyRes(:))))
UxRes=Ux-u;
%disp('MaxResidualUx');disp(max(abs(UxRes(:)))
fprintf('MaxResidualUx %i.\n',max(abs(UxRes(:))))
UyRes=Uy-v;
%disp('MaxResidualUy');disp(max(abs(UyRes(:)))
fprintf('MaxResidualUy %i.\n',max(abs(UyRes(:))))


%Percentage error
%Residual as a percentage of the driving stress
%this avoids doing it as a percentage of the actual values were values
%close to 0 give unreasonably high errors 
Drivingstress=Sxxd+Syyd+Sxyd;

SxxResPerc=((100./Sxx_An(:)).*SxxTwodd(:))-100; %scatter(X(:),Y(:),15,SxxResPerc(:))
SyyResPerc=((100./Syy_An(:)).*SyyTwodd(:))-100;
SxyResPerc=((100./Sxy_An(:)).*SxyTwodd(:))-100;
 
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

PpercStress=[max(abs(SxxResPerc(:))),max(abs(SyyResPerc(:))),max(abs(SxyResPerc(:)))];

% if any(PpercStress>1);
%     error('Error in surronding grid too large. Something has changed in Twodd Calculation. This no longer matches analytical solutions')
% else
%     disp('Everything looks good, checks the highestn tensor error % is below 1% of the total driving stress')
% end
