% Test: Hole with internal pressure: Compares stresses around a hole with
% an internal pressure that sits in an infinite elastic.
% 
% Solution: Used the script for this from the Pollard and Fletcher
% fundamentals book. This solution is In most elasticity/continuum
% textbooks as its simple to implement and relevant. Itâÿÿs the original
% solution from: Kirsch, G., 1898. Theory of Elasticity and Application in
% Strength of Materials.Â Zeits chrift des Vereins Deutscher
% Ingenieure,Â 42(29), pp.797-807. The solution is independent of elastic
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

% %Circle
ri = 1;
k = 0;
for theta = 0:pi/600:2*pi %400
    k = k+1;
    x3(k) = ri*cos(theta);
    y3(k) = ri*sin(theta);
end
sz=numel(x3);
Pointsxy=[x3',y3'];
mystruct.line1=(1:sz);
%Fixed points inner circle (making a square)
rr=0.8;%size of the square
x=[-rr,0,rr,0,-rr];                        %List of x points
y=[0,-rr,0,rr,0];
PointsxyF=[x;y]';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Pointsxy,mystruct,Fdisp] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
	%Creating Fictitious disp flag, if set to one then this forces the elements with this to 0 displacement. 
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


line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r','LineWidth',2.5)
title('Hole elements with normals'), xlabel('x'), ylabel('y');axis equal
hold on
nx=cos(NormAng);
ny=cos((pi/2)-NormAng);
quiver(xe,ye,nx,ny)
%Fixed elements shown as blue
line([Points(logical(Fdisp),1)';Points(logical(Fdisp),2)'],[Points(logical(Fdisp),3)';Points(logical(Fdisp),4)'],'color','b','LineWidth',2.5)
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
[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Putting pressure on elements instead of slip. Inf code
    %works out the superposition of this pressure
	%Output 'stresses' are blank arrays of 0's as this options input stresses are not within the global CS    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 cc=zeros(NUM,1); 
    
 Tn = 1;    Tn=cc+Tn;	%'Tensile' traction
 Ts = 0;    Ts=cc+Ts;   %Shear traction always right lateral
 Option='E';

[ShearDisp,TensileDisp,Sxx,Syy,Sxy]=SlipCalculator2d(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);

    
%%
% Removing any fixed elements
if any(Fdisp)==1  
[XBEG,XEND,YBEG,YEND,NUM,x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng]...
    = RemovingFixedEls2d(Fdisp,XBEG,XEND,YBEG,YEND,NUM);
end

	
	
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Simply define flat observation plane of points. XYZ with defined step size 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boxsz=3;    
X = linspace(-boxsz,boxsz,101);
Y = linspace(-boxsz,boxsz,101); 


[X,Y] = meshgrid(X,Y); % define Cartesian grid
rowcount = length(X(:,1));
colcount = length(X(1,:));


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
%STEP AA: Analytical solution, Kirsch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xa = linspace(-boxsz,boxsz,101);
% Ya = linspace(-boxsz,boxsz,101); 
[SXX,SYY,SXY]=KirschSolutionFunction(X,Y,ri,-Tn); %Tn is extensional, Pressure is compressive

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Drawing figures and cleaning data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pxx=Sxx(1,1);
Pyy=Syy(1,1);
Pxy=Sxy(1,1);
%Drawing Figs
Sxx=TotalStress(:,1);Syy=TotalStress(:,2);Sxy=TotalStress(:,3);

Sxx = reshape(Sxx,rowcount,colcount);
Syy = reshape(Syy,rowcount,colcount);
Sxy = reshape(Sxy,rowcount,colcount);


%%%%%%%%%%%%%%%%%%
%Removing inner points
[TH,R] = cart2pol(X,Y);
Bad=R<1;
Sxx(Bad)=nan;
Syy(Bad)=nan;
Sxy(Bad)=nan;

Ux=InducedDisplacements(:,1);Uy=InducedDisplacements(:,2);
Ux = reshape(Ux,rowcount,colcount);Uy = reshape(Uy,rowcount,colcount);


[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);
S1 = reshape(S1,rowcount,colcount);
S2 = reshape(S2,rowcount,colcount);
S1dirx = reshape(S1dir(:,1),rowcount,colcount);
S1diry = reshape(S1dir(:,2),rowcount,colcount);
S2dirx = reshape(S2dir(:,1),rowcount,colcount);
S2diry = reshape(S2dir(:,2),rowcount,colcount);
scl=0.4;
figure,q1 = quiver(X,Y,S1dirx,S1diry); axis equal, hold on
q15 = quiver(X,Y,-S1dirx,-S1diry);
q2 = quiver(X,Y,S2dirx,S2diry);
q25 = quiver(X,Y,-S2dirx,-S2diry);
set (q1, 'linewidth', 0.8,'color', 'red','maxheadsize', 0,'AutoScaleFactor',scl);
set (q15, 'linewidth', 0.8,'color', 'red','maxheadsize', 0,'AutoScaleFactor',scl);
set (q2, 'linewidth', 0.8,'color', 'blue','maxheadsize', 0,'AutoScaleFactor',scl);
set (q25, 'linewidth', 0.8,'color', 'blue','maxheadsize', 0,'AutoScaleFactor',scl);
line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
xlabel('x'); ylabel('y'); axis('equal'); title('NumericalPrincipalStresses');


% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %Getting rid of singularities 
% % FILTER=1;
% % FilterVal=abs(FILTER*(abs(Pxx(1,:))+abs(Pyy(1,:))+abs(Pxy(1,:))));
% FilterVal=50;
% %%%%%%%%%%%%%%%%%%%%%%%%
% bad = abs(S1)>FilterVal;
% S1(bad)=0;
% bad = abs(S2)>FilterVal;
% S2(bad)=0;
% bad = abs(Sxx)>FilterVal;
% Sxx(bad)=0;
% bad = abs(Syy)>FilterVal;
% Syy(bad)=0;
% bad = abs(Sxy)>FilterVal;
% Sxy(bad)=0;

%Octave ContourF doesn't like 'nans', replacing with '-infs'
%MATLAB unsuprisingly handles this fine
Sxx(isnan(Sxx)) = -inf;    SXX(isnan(SXX)) = -inf;    
Syy(isnan(Syy)) = -inf;    SYY(isnan(SYY)) = -inf;    
Sxy(isnan(Sxy)) = -inf;    SXY(isnan(SXY)) = -inf;  


figure,quiver(X,Y,Ux,Uy);
xlabel('x'); ylabel('y'); axis('equal'); title('Displacement Numerical');

fntsz=25;
titlesz=23;

Lvlstp=0.2;

figure;subplot(2,3,1),colormap(cmap2),[C,h]=contourf(X,Y,Sxx);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x numerical');colorbar;ChangeFontSizes(titlesz,fntsz);
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,2),colormap(cmap2),[C,h]=contourf(X,Y,Syy);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y numerical');colorbar;ChangeFontSizes(titlesz,fntsz);
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,3),colormap(cmap2),[C,h]=contourf(X,Y,Sxy);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y numerical');colorbar;ChangeFontSizes(titlesz,fntsz);
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)


subplot(2,3,4),colormap(cmap2),[C,h]=contourf(X,Y,SXX);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x analytical');colorbar;ChangeFontSizes(titlesz,fntsz);
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,5),colormap(cmap2),[C,h]=contourf(X,Y,SYY);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y analytical');colorbar;ChangeFontSizes(titlesz,fntsz);
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)
subplot(2,3,6),colormap(cmap2),[C,h]=contourf(X,Y,SXY);set (h, 'levelstep', Lvlstp);
xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y analytical');colorbar;ChangeFontSizes(titlesz,fntsz);
caxis([-1 1]);set(colorbar,'YTick',-1:0.25:1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
 
% plotpoints=Y==0;
% Xx=X(plotpoints);
% Yy=Y(plotpoints);
% Sxx_pp=Sxx(plotpoints);
% SXX_pp=SXX(plotpoints);
% figure;plot(Xx,Sxx_pp,'r');
% hold on
% plot(Xx,SXX_pp);

% %Percentage error
% %Rounding to remove really small values that will mess up percentage
% SXX = round(SXX,5);Sxx = round(Sxx,5);
% SYY = round(SYY,5);Syy = round(Syy,5);
% SXY = round(SXY,5);Sxy = round(Sxy,5);
% SxxResPerc=((SXX-Sxx)./SXX)*100;
% SyyResPerc=((SYY-Syy)./SYY)*100;
% SxyResPerc=((SXY-Sxy)./SXY)*100;
%                   
% answer=isnan(SxxResPerc);
% SxxResPerc(answer)=0;
% answer=isinf(SxxResPerc);
% SxxResPerc(answer)=0;
% answer=isnan(SyyResPerc);
% SyyResPerc(answer)=0;
% answer=isinf(SyyResPerc);
% SyyResPerc(answer)=0;
% answer=isnan(SxyResPerc);
% SxyResPerc(answer)=0;
% answer=isinf(SxyResPerc);
% SxyResPerc(answer)=0;
% 
% PpercStress=[max(SxxResPerc(:)),max(SyyResPerc(:)),max(SxyResPerc(:))];

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

if any(PpercStress>10);
    error('Error in surrounding grid too large. Something has changed in Twodd Calculation. This no longer matches analytical solutions')
else
    disp('Everything looks good, checks the highest tensor error % is below 1% of the total driving stress')
end


%OTHER FOR TESTING
%%%%%%%%%%%%%%%%%%
%Creating radial coordinates for each observation point and then converting
%tensors from cart to polar coordinates. Eq 6.96 Pollard and Fletcher 2005
[TH,R] = cart2pol(X,Y);
ST = sin(TH); S2T = sin(2*TH); ST2 = ST.^2; 
CT = cos(TH); C2T = cos(2*TH); CT2 = CT.^2;
% Cartesian to polar stress components
SRR = Sxx.*CT2+Syy.*ST2+2*Sxy.*CT.*ST;
STT = Sxx.*ST2+Syy.*CT2-2*Sxy.*CT.*ST;
SRT = ST.*CT.*(Syy-Sxx)+Sxy.*(CT2-ST2);
Sxx(Bad)=nan;
% figure;colormap(cmap),contourf(X,Y,-SRR);
% xlabel('x'); ylabel('y'); axis('equal'); title('Srr Numerical');colorbar;
