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

%   Copyright 2017, Tim Davis, The University of Aberdeen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bit's  you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/');  end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    %Reshaping the list of XY into usable variables

    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
    %Use this if the algorithm is throwing errors for vertex's above 0m.
    
freesurface_height = 0;
Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;

%Elastic parameters
mu=1;            	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
E = mu*(2*(1+nu)) ;                     %Young's Modulus
lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard

    %Reshaping the list of XY into usable variables
[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng,NUM,XBEG,XEND,YBEG,YEND]...
 = CreateElements2d( Pointsxy,mystruct );



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%
    
strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    
    %%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
%Input Stress
Sxx = 0; 					%Positive is tension
Syy = 0.01;
Sxy = 0; 
 Option='B';

%Loading fault that lies in XY and calculating slip
[ShearDisp,TensileDisp,Sxx,Syy,Sxy]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Function that calls eq 8.35a and 8.35b from the Pollard/Segall paper
[OpeningOrSlip]=Func_SlipPollardSegall(mu,nu,Syy,Sxy,xe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and drawing figures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SlipDispRes=abs(OpeningOrSlip)-abs(TensileDisp+ShearDisp);
SlipDispResPerc=((abs(OpeningOrSlip)-abs(TensileDisp+ShearDisp))./abs(OpeningOrSlip))*100;

%disp('MaxResidualSlip');disp(max(SlipDispRes(:)));
fprintf('MaxResidualSlip %i.\n',max(SlipDispRes(:))) %prints on one line unlike 'disp
fprintf('MaxResidualSlipAsAPercentage %i.\n',max(SlipDispResPerc(:))) %prints on one line unlike 'disp
lnwidth=2.5;
plot(xe,OpeningOrSlip,'g',x,SlipDispRes,'r','LineWidth',lnwidth);title('OpeningDistributionAlongCrack','LineWidth',lnwidth);
hold on
step=x(1,1)-x(2,1); stairs([x(1,1)+step;x;x(end,end)-step],[0;(TensileDisp+ShearDisp);0],'m','LineWidth',lnwidth);
xlabel('Distance from crack centre')
ylabel({'Opening, analytical (green)','numerical (pink) residual(red)'})
xlim ([-1 1])

fntsz=18;
titlesz=22;
ChangeFontSizes(titlesz,fntsz);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Drawing a different way:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subsetting data - would be better to interpolate at points I want to draw
gd=x>0;
TensileDispP=TensileDisp(gd);
xP=x(gd);
N=8;
TensileDispN = TensileDispP(1:N:end,:);
xN = xP(1:N:end,:);

%trying interpolation
xN = linspace(0,0.9,10);
TensileDispN = interp1(x,TensileDisp,xN,'spline');

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
scatter(abs(xN),abs(TensileDispN),24,'k','filled');
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
scatter(xP,abs(TensileDispP),24,'k','filled');
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
[OpeningOrSlip]=Func_SlipPollardSegall(mu,nu,Syy,Sxy,xP);
PercentErrorOpeningInterpPnts=((100./-OpeningOrSlip).*TensileDispP)-100;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
spacing=0.1; %0.01
minx=-4; maxx=4;
[X,Y] = meshgrid(minx:spacing:maxx); 
dimx = length(X(:,1));
dimy = length(X(1,:));

%Twodd grid solution
SinB=sin(Beta);	

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

Ux=InducedDisplacements(:,1);Uy=InducedDisplacements(:,2);
Ux = reshape(Ux,dimx,dimy);Uy = reshape(Uy,dimx,dimy);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP BB: Analytical solution stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sxx_An,Syy_An,Sxy_An,u,v]=Func_StressDisplacementGridPollardPaperTrigFunctions(mu,nu,Syy,Sxy,spacing,minx,maxx,cmap);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cleaning data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Removing spurios values for the code that are right next to the crack. The Pollard analytical
%solution is far more accurate at this at the tips. 
bad=abs(X)<1+(spacing*1);
bad2=abs(Y)<1*spacing;
bad3=bad+bad2;
bad4=bad3==2;
SxxTwodd=StressChange(:,1); 
SyyTwodd=StressChange(:,2);
SxyTwodd=StressChange(:,3);
SxxTwodd=reshape(SxxTwodd,dimx,dimy);
SyyTwodd=reshape(SyyTwodd,dimx,dimy);
SxyTwodd=reshape(SxyTwodd,dimx,dimy);
SxxTwodd(bad4)=nan;SyyTwodd(bad4)=nan;SxyTwodd(bad4)=nan;
%Sxx_An(bad4)=nan;Syy_An(bad4)=nan;Sxy_An(bad4)=nan;
Ux(bad4)=nan;Uy(bad4)=nan;
%u(bad4)=nan;v(bad4)=nan;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,quiver(X,Y,Ux,Uy);
xlabel('x'); ylabel('y'); axis('equal'); title('displacement from code, central patch removed');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Drawing Figs
Sxxd=Sxx; Sxx = SxxTwodd;
Syyd=Syy; Syy = SyyTwodd;
Sxyd=Sxy; Sxy = SxyTwodd;
%Octave ContourF doesn't like 'nans', replacing with '-infs'
%Matlab unsuprisingly handles this fine
% Sxx(isnan(Sxx)) = -inf;    Sxx_An(isnan(Sxx_An)) = -inf;    
% Syy(isnan(Syy)) = -inf;    Syy_An(isnan(Syy_An)) = -inf;    
% Sxy(isnan(Sxy)) = -inf;    Sxy_An(isnan(Sxy_An)) = -inf;  
% Sxx(isnan(Sxx)) = 0;    Sxx_An(isnan(Sxx_An)) = 0;    
% Syy(isnan(Syy)) = 0;    Syy_An(isnan(Syy_An)) = 0;    
% Sxy(isnan(Sxy)) = 0;    Sxy_An(isnan(Sxy_An)) = 0; 

%Setting contour steps, the steps set below are the same as the figures in the Martel
%2000 paper.
LvlstpX=0.001;
LvlstpY=0.001;
LvlstpXY=0.001;

%Setting axis limit's  so supurious values do not change the overall trend. 
%Only looking at finite values, No axis scaling to infs or nans. 
%caxissxx=[min(Sxx_An(isfinite(Sxx_An(:)))),max(Sxx_An(isfinite(Sxx_An(:))))];
%caxissyy=[min(Syy_An(isfinite(Syy_An(:)))),max(Syy_An(isfinite(Syy_An(:))))];
%caxissxy=[min(Sxy_An(isfinite(Sxy_An(:)))),max(Sxy_An(isfinite(Sxy_An(:))))];
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


fntsz=25;
titlesz=23;
%Sxx
figure,         %create figure
subplot(2,3,1), %create subplot with 6 parts
colormap(cmap), %import the blue red colourmap that was created at the beggining of the file
[C,h]= contourf(X,Y,Sxx); %now draw the data in the figure
caxis(caxissxx);          %set the min and max contour limit's    
xlabel('x'); ylabel('y'); %create the X label and Y label on the axes
axis('equal');            %make the X and Y axis equal
title('\sigma_x_x numerical');   %create a title over the top of the plot
set (h, 'levelstep', LvlstpX);  %set the levels for the contours based on the vector created above
length=caxissxx(1)-caxissxx(2); %the distance between the max and min limit's 
vectoroflocs=linspace(caxissxx(1),caxissxx(2),abs(length/LvlstpX)+1); %now using this length to create a spaced vector for labeling the colourbar
colorbar();               %drawing the colourbar
cbh = findobj( gcf(), 'tag', 'colorbar');   %Finding this object we just created
set( cbh,'ytick', [vectoroflocs])%Setting the ticks to nice round numbers
ChangeFontSizes(titlesz,fntsz);caxis([-0.005 0.005]);set(colorbar,'YTick',-0.005:0.0025:0.0050);
pause(3)%octave issue

subplot(2,3,4),colormap(cmap),[C,h]=contourf(X,Y,Sxx_An);
caxis(caxissxx)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x analytical');
set (h, 'levelstep', LvlstpX);
length=caxissxx(1)-caxissxx(2);
vectoroflocs=linspace(caxissxx(1),caxissxx(2),abs(length/LvlstpX)+1);
colorbar();
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', [vectoroflocs])
ChangeFontSizes(titlesz,fntsz);caxis([-0.005 0.005]);set(colorbar,'YTick',-0.005:0.0025:0.0050);

%Syy
subplot(2,3,2),colormap(cmap),[C,h]=contourf(X,Y,Syy);
caxis(caxissyy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y numerical');
set (h, 'levelstep', LvlstpY);
length=caxissyy(1)-caxissyy(2);
vectoroflocs=linspace(caxissyy(1),caxissyy(2),abs(length/LvlstpY)+1);
colorbar();
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', [vectoroflocs])
ChangeFontSizes(titlesz,fntsz);caxis([-0.01 0.01]);

subplot(2,3,5),colormap(cmap),[C,h]=contourf(X,Y,Syy_An);
caxis(caxissyy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y analytical');
set (h, 'levelstep', LvlstpY);
length=caxissyy(1)-caxissyy(2);
vectoroflocs=linspace(caxissyy(1),caxissyy(2),abs(length/LvlstpY)+1);
colorbar();
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', [vectoroflocs])
ChangeFontSizes(titlesz,fntsz);caxis([-0.01 0.01]);

%Sxy
subplot(2,3,3),colormap(cmap),[C,h]=contourf(X,Y,Sxy);
caxis(caxissxy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y numerical'),
set (h, 'levelstep', LvlstpXY);
length=caxissxy(1)-caxissxy(2);
vectoroflocs=linspace(caxissxy(1),caxissxy(2),abs(length/LvlstpXY)+1);
colorbar();
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', [vectoroflocs])
ChangeFontSizes(titlesz,fntsz);set(colorbar,'YTick',-0.005:0.0025:0.0050);

subplot(2,3,6),colormap(cmap2),[C,h]=contourf(X,Y,Sxy_An);
caxis(caxissxy)
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y analytical');
set (h, 'levelstep', LvlstpXY);
vectoroflocs=linspace(caxissxy(1),caxissxy(2),abs((caxissxy(1)-caxissxy(2))/LvlstpXY)+1);
colorbar();
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,'ytick', [vectoroflocs])
ChangeFontSizes(titlesz,fntsz);set(colorbar,'YTick',-0.005:0.0025:0.0050);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% SxxTwodd=SxxTwodd(~bad4);Sxx_An=Sxx_An(~bad4);
% SyyTwodd=SyyTwodd(~bad4);Syy_An=Syy_An(~bad4);
% SxyTwodd=SxyTwodd(~bad4);Sxy_An=Sxy_An(~bad4);
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
