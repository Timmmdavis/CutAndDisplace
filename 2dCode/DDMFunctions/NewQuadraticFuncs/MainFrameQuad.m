%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you don’t need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%   Copyright 2017, Tim Davis, The University of Aberdeen

% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
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

clear;close all
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
		
		%Strain uses the same convention as engineering stress convention where extension is a positive value. When entering strain values remember this. 
		%For a normal stress: positive (negative) if it produces tension (compression) in the material.
		%For a shear stress: positive if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is + if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear strain Exy will cause Leftlateral movement.
		
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading XY ascii data or manually creating fractures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% %Single flat user defined fracture
% x=linspace(-1,1,27);    
% y=zeros(1,numel(x)); 
% Pointsxy=[x;y]';
% mystruct.line1=(1:(length(Pointsxy(:,1))));
% Fnms=fieldnames(mystruct);
% nf=numel(Fnms);clear Fnms


% % % %Circle
xscl=1;
yscl=1;
% %Circle
ri = 1;
k = 0;
for theta = 0:pi/5:2*pi %600
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
[Pointsxy,mystruct,Fdisp] = DataAppender2d( Pointsxy,PointsxyF,mystruct );



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Shp
    %Loads Points that represent the file. 
    %If there is no Y component the objects Y values have been exported as Z!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Pointsxy,mystruct]=m_shaperead('CurvedLine');
% Fnms=fieldnames(mystruct);
% nf=numel(Fnms);clear Fnms

%Shifting this line to the right place
%Pointsxy(:,1)=Pointsxy(:,1)-50; 
%Pointsxy(:,2)=Pointsxy(:,2)-2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Appending and flagging the fixed data points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[Pointsxy,mystruct,nf,Fdisp] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
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


    %Defining elastic constants,  See equations below if you need to calculate these 
	%from other constants. 
	
mu=1;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
E = mu*(2*(1+nu)) ;                     %Young's Modulus

	%Equations you can use to calculate elastic parameters that you need above. 
% K = E/3*(1-2*nu);                 %Bulk Modulus.    Equation 8.25 Pollard
% mu = E/(2*(1+nu));                %Shear Modulus.   Equation 8.26 Pollard  which is Geo convention 
lambda= E*nu/((1+nu)*(1-2*nu));   %Lamé's  constant,  Equation 8.27 Pollard
% nu =lambda/(2*(mu+lamda);            %Poisson's ratio, Equation 8.28 Pollard

    %Reshaping the list of XY into usable variables
[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng,NUM,XBEG,XEND,YBEG,YEND]...
 = CreateElements2dQuad( Pointsxy,mystruct );



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
%STEP 3: Define how you want to calculate slip - Constant, Remote (Friction), Remote
%no opening, tractions on the elements (internal pressure) or gravitational induced stresses. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Sxy,Tn,Ts,Mu,Sf,strain ] = CreateBlankVars;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Constant Slip across the fault surface
	%Runs a constant slip on every element. A mixture of two slips will give an inclined slip.
	%Output 'stresses' are blank arrays of 0's as this option is not driven by boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Option='A';     
% cc=zeros(NUM*3,1); 
% ShearDisp  = 0;       ShearDisp	  = cc+(ShearDisp*-1);   %Always right lateral any orientation for this option
% TensileDisp =1;       TensileDisp = cc+(TensileDisp*-1); %Positive = extensional movement
% Sxx = 0;  Syy = 0;  Sxy = 0;   	
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %the stresses can be defined as single values that are repeated for
    %every element or defined as a vector that is the stress at every
    %element: i.e.  Syy = linspace(0,1,NUM); monatomic increasing Syy to 1 at the last element 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%
    
 strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
     
 	%%%%%%%%%%%%%% 
     %StressInput
 	%%%%%%%%%%%%%%
 	
  Sxx = 0; Sxxd=Sxx;					%Positive is tension
  Syy = -1; Syyd=Syy;
  Sxy = 0; Sxyd=Sxy;                  %For fractures lying along Xaxis positive stress drives right lateral movement, along Y is left lat. Pollard fig6.13c
  Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%
    %StrainInput       
	%%%%%%%%%%%%%%
   
% strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
%     
% 	%%%%%%%%%%%%%% 
%     %StressInput
% 	%%%%%%%%%%%%%%
% 	
%  Sxx = 0; 					%Positive is tension
%  Syy = -2;
%  Sxy = -1;                   %For fractures lying along Xaxis positive stress drives right lateral movement, along Y is left lat. Pollard fig6.13c
% 
% 
% % % %Frictional parameters, these can be changed to variable friction
%  Mu  = 0.2;     Mu=repmat(Mu,NUM*3,1); %Mu=[Mu1;zeros(freesz,1);Mu1]; %Coefficient of friction
%  Sf  = 0;       Sf=repmat(Sf,NUM*3,1); %Sf=[Sf1;zeros(freesz,1);Sf1]; %Frictional strength
% 
%  
% Option='C'; %slip from uniform remote stress with friction


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option D = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. This option has no opening components.  
	%Choose to define in stress or strain.  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
  	
% strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
%     
% 	%%%%%%%%%%%%%
%     %StressInput
% 	%%%%%%%%%%%%%
%     
% Sxx = 0;         			%Positive is tension 
% Syy = 0; 
% Sxy = 1;                    %Lying along Xaxis positive is right lateral movement, along Y is left lat. Pollard fig6.13c
% Option='D'; %slip no opening


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option E = Putting pressure on elements instead of slip. Inf code
    %works out the superposition of this pressure
	%Output 'stresses' are blank arrays of 0's as this options input stresses are not within the global CS    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%  cc=zeros(NUM,1); 
%     
%  Tn = 0;    Tn=cc+Tn;	%Tensile traction
%  Ts = 1;    Ts=cc+Ts;   %Shear traction always right lateral
%  Option='E'; %Traction defined on elements (pressure)



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(Option,'A')
%if ~Option=='A'
[ShearDisp,TensileDisp,Sxx,Syy,Sxy]=SlipCalculator2dQuad(x,y,xe,ye,HalfLength,Beta,Sxx,Syy,Sxy,Tn,Ts,nu,E,halfspace,NUM,NormAng,strain,Fdisp,Mu,Sf,Option);
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3 Drawing Figures: Drawing Figures of slip distribution on the crack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Removing any fixed elements
if any(Fdisp)==1  
[XBEG,XEND,YBEG,YEND,NUM,x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng]...
    = RemovingFixedEls2dQuad(Fdisp,XBEG,XEND,YBEG,YEND,NUM);
end


%Actually needs to use func to be an accurate plot. 
GraphSlipVsElNumQuad( NUM,TensileDisp,ShearDisp )


%Creating directions and magnitudes of slip for plotting
PlotOpeningVsShearOnElsQuad( NormAng,TensileDisp,ShearDisp,Points,x,y,HalfLength )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Animate Fault movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To do this just call 'Animate2d' in the cmd window. To change the
%parameters, save movie or the camera angle go into the 'Animate' function


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Regular Grid that sits around the fault surfaces. Define how far it
    %extends away and its density
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cells=5;   %Define the number of observations points on the square grid
% padding=10;  %How much extra bumph you want to add away from the fault sticks, reduce to 0 if having half space issues 
% [maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents(Points,cells,padding);
% 
% %Uniform grid spacing
% [X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);
% rowcount = size(X,1);
% colcount = size(X,2);

% boxsz=3;    
% X = linspace(-boxsz,boxsz,101);
% Y = linspace(-boxsz,boxsz,101); 
% 
% [X,Y] = meshgrid(X,Y); % define Cartesian grid
% rowcount = length(X(:,1));
% colcount = length(X(1,:));

% %For comp to An
% spacing=0.5; %0.01
% minx=-4; maxx=4;
% [X,Y] = meshgrid(minx:spacing:maxx); 
% rowcount = length(X(:,1));
% colcount = length(X(1,:));

% %maxgriX=2;mingriX=-2;maxgriY=2;mingriY=-2;
% %%Random points within the calculated bounds
% NumPnts=15000; %number of points, don’t have cells 
% lengthx=maxgriX-mingriX;
% lengthy=maxgriY-mingriY;
% xmv=(maxgriX+mingriX)/2; ymv=(maxgriY+mingriY)/2;
% x=rand(1,NumPnts)*(lengthx*2);
% y=rand(1,NumPnts)*(lengthy*2);
% X=x-lengthx+xmv;
% Y=y-lengthy+ymv;

%At el coeff pnts, moving along normal slightly
nx=cos(NormAng);
ny=cos((pi/2)-NormAng);
[ nxq,nyq ] = RepeatQuadValues( nx,ny );
%Grabbing col pnts
X=x+(nxq*1e-12);
Y=y+(nyq*1e-12);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = User defined XY grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Creating my own centred grid with how many points
% X = linspace(-10,15,25)+eps; % add eps to avoid singularity at origin
% Y = linspace(-10,15,25); 
% [X,Y] = meshgrid(X,Y); % define Cartesian grid
% dimx = length(X(:,1));
% dimy = length(X(1,:));   

% %Creating a flat line
%  X = linspace(0,15,100)+eps; % add eps to avoid singularity at origin
%  Y = zeros(size(X)); 

 
% % Random data at on n*n grid centred at xmv,ymv
% n=4;
% xmv=0; ymv=0;
% x=rand(1,10201)*n;
% y=rand(1,10201)*n;
% X=x-n+xmv;
% Y=y-n+ymv;

% %points at midpoints
% X=x;
% Y=y;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = Import XY from a .text file within the folder you are
    %running this script from
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Loads the dat file, This contains 3 columns XY
	%XY are the points the stresses/strains and displacements will be calculated on.

% S = load('Points.dat');	
% X = S(:,1);                   
% Y = S(:,2);                         


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing and fixing Obs Point data just created
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Replacing XY of any obs points for Nan if these lie on top of the end points of a
    %dislocation, without this it can result in spurious values. 
[ X,Y  ] = NullPointsLyingOnElement( X,Y,XBEG,YBEG,XEND,YEND,NUM );


    %Drawing these and the location compared to the Obs Points
figure;line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
title('fractures and obs points'), xlabel('x'), ylabel('y');axis equal
hold on
scatter(X(:),Y(:))
hold off

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
[StressChange,TotalStress,StressReg,InducedDisplacements]=StressDispOnSurroundingPointsQuad(X,Y,xe,ye,HalfLength,Sxx,Syy,Sxy,nu,E,halfspace,ShearDisp,TensileDisp,Beta);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 6: Finite strain computation, error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Grabbing Displacements
%[ Ux,Uy ] = ExtractCols( InducedDisplacements );
   

%[exx,eyy,exy,ExxInfErrorPerc,EyyInfErrorPerc,ExyInfErrorPerc]=FiniteStrainLagrangian2d(Ux,Uy,X,Y);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Drawing figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Grabbing Displacements
[ Ux,Uy ] = ExtractCols( InducedDisplacements );
   
%Remote stress
Pxx=Sxx(1,1);
Pyy=Syy(1,1);
Pxy=Sxy(1,1);


%Extracting stress
[ Sxx,Syy,Sxy ] = ExtractCols( TotalStress );


%Stress2strain
[ Exx,Eyy,Exy ] = Hooke'sLaw2dStress2Strain( Sxx,Syy,Sxy,lambda,mu );

%Calclating 2d EigenValues
[S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);

%Calclating 2d EigenValues
[E1,E2,E1dir,E2dir]=EigCalc2d(Exx,Eyy,Exy);
Dilatation2d=E1+E2;

% %If the stresses don’t draw well use this function
% FilterValue=2.5;
% [S1,S2,Sxx,Syy,Sxy] = NanOrRemoveBadPoints( FilterValue,5,1,S1,S2,Sxx,Syy,Sxy );


%Draw data with quiver, doesn't matter if its gridded or not
figure,quiver(X(:),Y(:),Ux,Uy);
xlabel('x'); ylabel('y'); axis('equal'); title('Disp'); 

%Drawing S1 S2 directions
DrawS1S2Directions(X(:),Y(:),S1dir,S2dir )


%Check if data is on a uniform or non uniform grid
[ Uniform ] = UniformGridCheck2d( X,Y );


%If data is not uniform plot as non grid data. 
if Uniform==0
    
DrawScatterPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );

%else reshapefrom col vectors to meshes then draw with nice grid functions
else

%Reshaping stresses to grid dimensions
[Sxy,Syy,Sxx,S1,S2]=ReshapeData2d( rowcount,colcount,Sxy,Syy,Sxx,S1,S2 );

DrawContourFPlots2d( X,Y,cmap2,Sxy,Syy,Sxx,S1,S2 );

Scl=1;
DrawDeformedGrid2d( X,Y,Ux,Uy,Scl,cmap2,(E1+E2)); %Dilatation2d %(abs(Ux)+abs(Uy)

end


%KirschComp    
pm=0;sH=0;sh=-1;    ri=1;
%CalcAnResult    
[SxxAn,SyyAn,SxyAn,SRRan,STTan,SRTan]=KirschSolutionFunction(nxq,sin(acos(nxq)),ri,pm,sH,sh);
figure
scatter(X(:),Y(:),15,STTan(:))

%polar coords
[ SRR,STT,SRT ] = StressTensorTransformation2d(Sxx,Syy,Sxy,nyq,nxq);
[TH,R] = cart2pol(X./xscl,Y./yscl);
figure;
plot(rad2deg(TH),STTan,'b','LineWidth',2)
hold on
N=7;
scatter(rad2deg(TH),SRR,12,'k')
hold off
title('Max hoop stress 2d numerical');
xlabel('Theta'); ylabel('HoopStress'); 

Res=STTan-SRR;
disp('max');max(abs(Res))
disp('mean');nanmean(abs(Res))


% [Sxy,Syy,Sxx,S1,S2]=ReshapeData2d( rowcount,colcount,Sxy,Syy,Sxx,S1,S2 );
% 
% %%comparing to an solution
% [OpeningOrSlip]=Func_SlipPollardSegall(mu,nu,Syyd,Sxyd,x);
% Rel=OpeningOrSlip./TensileDisp;
% scatter(x,TensileDisp);
% hold on
% scatter(x,OpeningOrSlip);
% 
% 
% [Sxx_An,Syy_An,Sxy_An,Ux_An,Uy_An]=Func_StressDisplacementGridPollardPaperTrigFunctions(mu,nu,Syyd,Sxyd,spacing,minx,maxx,cmap);
% %Draw data with quiver, doesn't matter if its gridded or not
% figure,quiver(X(:),Y(:),Ux_An(:),Uy_An(:));
% xlabel('x'); ylabel('y'); axis('equal'); title('AnDisp'); 
% 
% % Ux_An=Ux_An(:);
% % Uy_An=Uy_An(:);
% 
% %Drawing Figs
% %Setting contour steps, the steps set below are the same as the figures in the Martel
% %2000 paper.
% LvlstpX=0.001;
% LvlstpY=0.001;
% LvlstpXY=0.001;
% 
% %Setting axis limits so supurious values don’t change the overall trend. 
% %Only looking at finite values, No axis scaling to infs or nans. 
% %caxissxx=[min(Sxx_An(isfinite(Sxx_An(:)))),max(Sxx_An(isfinite(Sxx_An(:))))];
% %caxissyy=[min(Syy_An(isfinite(Syy_An(:)))),max(Syy_An(isfinite(Syy_An(:))))];
% %caxissxy=[min(Sxy_An(isfinite(Sxy_An(:)))),max(Sxy_An(isfinite(Sxy_An(:))))];
% caxissxx=[min(Sxx(isfinite(Sxx(:)))),max(Sxx(isfinite(Sxx(:))))];
% caxissyy=[min(Syy(isfinite(Syy(:)))),max(Syy(isfinite(Syy(:))))];
% caxissxy=[min(Sxy(isfinite(Sxy(:)))),max(Sxy(isfinite(Sxy(:))))];
% %Rounding values, rounds the axis values to the set step.
% Xrnd=1/LvlstpX;
% Yrnd=1/LvlstpY;
% XYrnd=1/LvlstpXY;
% caxissxx = (round(caxissxx*Xrnd))/Xrnd; %Rounding to nearest 0.2
% caxissyy = (round(caxissyy*Yrnd))/Yrnd; %Rounding yy to nearest 0.5 (same as levelstep)
% caxissxy = (round(caxissxy*XYrnd))/XYrnd; %Rounding to nearest 0.02
% 
% %%
% %Removing spurios values for the code that are right next to the crack. The Pollard analytical
% %solution is far more accurate at this at the tips. 
% % bad=abs(X)<1+(spacing*1);
% % bad2=abs(Y)<1*spacing;
% % bad3=bad+bad2;
% % bad4=bad3==2;
% % 
% % Sxx=reshape(Sxx,rowcount,colcount);
% % Syy=reshape(Syy,rowcount,colcount);
% % Sxy=reshape(Sxy,rowcount,colcount);
% % Sxx(bad4)=nan;Syy(bad4)=nan;Sxy(bad4)=nan;
% 
% %%%
% fntsz=25;
% titlesz=23;
% %Sxx
% figure,         %create figure
% subplot(2,3,1), %create subplot with 6 parts
% colormap(cmap), %import the blue red colourmap that was created at the beggining of the file
% [C,h]= contourf(X,Y,Sxx); %now draw the data in the figure
% caxis(caxissxx);          %set the min and max contour limits   
% xlabel('x'); ylabel('y'); %create the X label and Y label on the axes
% axis('equal');            %make the X and Y axis equal
% title('\sigma_x_x numerical');   %create a title over the top of the plot
% set (h, 'levelstep', LvlstpX);  %set the levels for the contours based on the vector created above
% length=caxissxx(1)-caxissxx(2); %the distance between the max and min limits 
% vectoroflocs=linspace(caxissxx(1),caxissxx(2),abs(length/LvlstpX)+1); %now using this length to create a spaced vector for labeling the colourbar
% colorbar();               %drawing the colourbar
% cbh = findobj( gcf(), 'tag', 'colorbar');   %Finding this object we just created
% set( cbh,'ytick', [vectoroflocs])%Setting the ticks to nice round numbers
% ChangeFontSizes(titlesz,fntsz);caxis([-0.005 0.005]);set(colorbar,'YTick',-0.005:0.0025:0.0050);
% pause(3)%octave issue
% 
% subplot(2,3,4),colormap(cmap),[C,h]=contourf(X,Y,Sxx_An);
% caxis(caxissxx)
% xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_x analytical');
% set (h, 'levelstep', LvlstpX);
% length=caxissxx(1)-caxissxx(2);
% vectoroflocs=linspace(caxissxx(1),caxissxx(2),abs(length/LvlstpX)+1);
% colorbar();
% cbh = findobj( gcf(), 'tag', 'colorbar');
% set( cbh,'ytick', [vectoroflocs])
% ChangeFontSizes(titlesz,fntsz);caxis([-0.005 0.005]);set(colorbar,'YTick',-0.005:0.0025:0.0050);
% 
% %Syy
% subplot(2,3,2),colormap(cmap),[C,h]=contourf(X,Y,Syy);
% caxis(caxissyy)
% xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y numerical');
% set (h, 'levelstep', LvlstpY);
% length=caxissyy(1)-caxissyy(2);
% vectoroflocs=linspace(caxissyy(1),caxissyy(2),abs(length/LvlstpY)+1);
% colorbar();
% cbh = findobj( gcf(), 'tag', 'colorbar');
% set( cbh,'ytick', [vectoroflocs])
% ChangeFontSizes(titlesz,fntsz);caxis([-0.01 0.01]);
% 
% subplot(2,3,5),colormap(cmap),[C,h]=contourf(X,Y,Syy_An);
% caxis(caxissyy)
% xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_y_y analytical');
% set (h, 'levelstep', LvlstpY);
% length=caxissyy(1)-caxissyy(2);
% vectoroflocs=linspace(caxissyy(1),caxissyy(2),abs(length/LvlstpY)+1);
% colorbar();
% cbh = findobj( gcf(), 'tag', 'colorbar');
% set( cbh,'ytick', [vectoroflocs])
% ChangeFontSizes(titlesz,fntsz);caxis([-0.01 0.01]);
% 
% %Sxy
% subplot(2,3,3),colormap(cmap),[C,h]=contourf(X,Y,Sxy);
% caxis(caxissxy)
% xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y numerical'),
% set (h, 'levelstep', LvlstpXY);
% length=caxissxy(1)-caxissxy(2);
% vectoroflocs=linspace(caxissxy(1),caxissxy(2),abs(length/LvlstpXY)+1);
% colorbar();
% cbh = findobj( gcf(), 'tag', 'colorbar');
% set( cbh,'ytick', [vectoroflocs])
% ChangeFontSizes(titlesz,fntsz);set(colorbar,'YTick',-0.005:0.0025:0.0050);
% 
% subplot(2,3,6),colormap(cmap2),[C,h]=contourf(X,Y,Sxy_An);
% caxis(caxissxy)
% xlabel('x'); ylabel('y'); axis('equal'); title('\sigma_x_y analytical');
% set (h, 'levelstep', LvlstpXY);
% vectoroflocs=linspace(caxissxy(1),caxissxy(2),abs((caxissxy(1)-caxissxy(2))/LvlstpXY)+1);
% colorbar();
% cbh = findobj( gcf(), 'tag', 'colorbar');
% set( cbh,'ytick', [vectoroflocs])
% ChangeFontSizes(titlesz,fntsz);set(colorbar,'YTick',-0.005:0.0025:0.0050);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Calculating residual and checking for errors. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Printing results (residual)
% SxxRes=Sxx-Sxx_An;
% %disp('MaxResidualSxx');disp(max(SxxRes(:)))
% fprintf('MaxResidualSxx %i.\n',max(abs(SxxRes(:))))
% SyyRes=Syy-Syy_An;
% %disp('MaxResidualSyy');disp(max(abs(SyyRes(:)))
% fprintf('MaxResidualSyy %i.\n',max(abs(SyyRes(:))))
% SxyRes=Sxy-Sxy_An;
% %disp('MaxResidualSxy');disp(max(abs(SxyRes(:)))
% fprintf('MaxResidualSxy %i.\n',max(abs(SxyRes(:))))
% 
% %Printing results (residual)
% UxRes=Ux-Ux_An(:);
% %disp('MaxResidualSyy');disp(max(abs(SyyRes(:)))
% fprintf('MaxResidualUy %i.\n',max(abs(UxRes(:))))
% UyRes=Uy-Uy_An(:);
% %disp('MaxResidualSxy');disp(max(abs(SxyRes(:)))
% fprintf('MaxResidualUx %i.\n',max(abs(UyRes(:))))
% 
