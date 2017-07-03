% Test: Check to see if the results from a flat dislocation lying in the half space is the same coming out of the 3d and 2d code. Runs on a planar fracture lying in X -1 to 1  and 25m below ground surface and checks the displacements along Y=0 X=0-100 
% 
% Solution: The 3d code for a single triangle is the analytical solution of Nikkhoo, M. and Walter, T.R., 2015. Triangular dislocation: an analytical, artefact-free solution.Geophysical Journal International,'201(2), pp.1117-1139
% I am using 4 triangles for the fault plane but these are meshed correctly so the results will not change
% Comparing this to the 2d code requires that we extend the fault along the Y axis far enough that end affects of non plain strain components don't mess with the '2d observations' result. I have extended the fault 500m in either direction along Y. The displacement is 0.01m to avoid the infinitesimal strain breaking down. 
% 
% Proof: This test that the HS solution is correct for the crack in 2d as this code is from disparate parts of the Crouch and Starfield book and is not explicitly written. Steve Martel provided an excellent base code for the 2d hs solution to start from. 

%   Copyright 2017, Tim Davis, The University of Aberdeen


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2D Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you don’t need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   
clear;close all
%===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir

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
x=-1:0.004:1;                        %List of x points

y= zeros(1,numel(x))-25;%x*0.5;
Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1)))-1);
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

    %Reshaping the list of XY into usable variables

XBEG=[];
XEND=[];
YBEG=[];
YEND=[];
for i=1:nf
    stringname{i} = strcat('line', num2str(i));
    tin=Pointsxy(mystruct.(stringname{i})(1,1:1:end-1),1);
    tan=Pointsxy(mystruct.(stringname{i})(1,2:1:end),1);
    ten=Pointsxy(mystruct.(stringname{i})(1:1:end-1),2); 
    tun=Pointsxy(mystruct.(stringname{i})(2:1:end),2);                          
    XBEG = [XBEG ; tin];
    XEND = [XEND ; tan];
    YBEG = [YBEG ; ten];
    YEND = [YEND ; tun];
    clear tin tan ten tun
end
NUM=numel(XBEG);      %Number of individual elements (segments)

    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
    
halfspace = 1; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
    %Use this if the algorithm is throwing errors for vertex's above 0m.
    
freesurface_height = 0;

YBEG=YBEG-freesurface_height;
YEND=YEND-freesurface_height;
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


[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng ] = MidPoint_Orientation( XBEG,XEND,YBEG,YEND,NUM );

line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
title('fractures'), xlabel('x'), ylabel('y')


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
ShearDisp  = 0;      ShearDisp	= cc+(ShearDisp*-1);   %For fractures lying along Xaxis positive stress drives right lateral movement, along Y is left lat. Pollard fig6.13c
TensileDisp =0.01;       TensileDisp	= cc+(TensileDisp*-1);  %Positive = extensional movement
Sxx = 0;  Syy = 0;  Sxy = 0;   	
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XYZ observation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

X=linspace(0,100,100);
Y=ones(1,numel(X))*-0.01;
dimx = length(X(:,1));
dimy = length(Y(:,1));



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

Ux2d=InducedDisplacements(:,1);Uy2d=InducedDisplacements(:,2);

	
	
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3D Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %The file must list the traingle vertex's then triangles. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
 string='FlatBlade_1a_4Faces.ts' ;
 [ Points,Triangles ] = GoCadAsciiReader( string,pathstring );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simplier in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 1; 

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
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);

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
	%Runs a constant slip on every triangle. A mixture of two slips will give an inclined slip.
	%Output 'stresses' are blank arrays of 0's as this option is not driven by boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
a=size(P1);
b=a(1,1);
c=zeros(b,1); 
    
StrikeSlipDisp  =0;      StrikeSlipDisp	= c+StrikeSlipDisp;   %Positive = left lat movement V
DipSlipDisp     =0;       DipSlipDisp		= c+DipSlipDisp;      %Positive = reverse movement
TensileSlipDisp =0.01;       TensileSlipDisp	= c+TensileSlipDisp;  %Positive = extensional movement
clear a b c             	%a b and c are used to create an array of slips the size of the number of triangles on the fault surface.

Sxx = 0; 					%Positive is compression 
Syy = 0; 
Szz = 0; 
Sxy = 0; 					%Positive is right lateral movement
Sxz = 0; 					
Syz = 0;   


X=linspace(0,100,100);
Y=zeros(1,numel(X));
Z=ones(1,numel(X))*-0.01;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 6: Calculate Displacements on dispersed XYZ
%observation points. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Calculate Displacements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints( StrikeSlipDisp,DipSlipDisp,TensileSlipDisp, nu, X,Y,Z, P1, P2, P3,halfspace);

 Ux3d=Ux;
 Uz3d=Uz;
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 7: VISULISATION AND ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

plot(X,Ux3d,'LineWidth',2,'color','r')
hold on
plot(X,Uz3d,'LineWidth',2,'color','g')
scatter(X(:),Ux2d(:),24,'k','filled')
scatter(X(:),Uy2d(:),24,'k','filled')
ylabel({'Displacement, 3d Ux = red, 2d Ux =dots';'3d Uz = pink, 2d Uy/z = dots'})
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Residual=(Ux3d-Ux2d)+(Uz3d-Uy2d);

Error_percUx=((Ux3d-Ux2d)/Ux3d)*100;
MaxError_percUx=max(max(Error_percUx));
Error_percUy=((Uz3d-Uy2d)/Uz3d)*100;
MaxError_percUy=max(max(Error_percUy));

fprintf('MaxResidualDispUx %i.\n',MaxError_percUx) %prints on one line unlike 'disp
fprintf('MaxResidualDispUy %i.\n',MaxError_percUy) %prints on one line unlike 'disp


if  MaxError_percUx > 0.25 && MaxError_percUy > 0.25  %making sure the results are really comparible. 
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, tolerance checks sum of UxUy displacement residuals and will throw error if any value is above 0.25%')
end

