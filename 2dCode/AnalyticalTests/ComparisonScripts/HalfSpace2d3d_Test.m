% Test: Check to see if the results from a flat dislocation lying in the half space is the same coming out of the 3d and 2d code. Runs on a planar fracture lying in X -1 to 1  and 25m below ground surface and checks the displacements along Y=0 X=0-100 
% 
% Solution: The 3d code for a single triangle is the analytical solution of Nikkhoo, M. and Walter, T.R., 2015. Triangular dislocation: an analytical, artefact-free solution.Geophysical Journal International,'201(2), pp.1117-1139
% I am using 4 triangles for the fault plane but these are meshed correctly so the results will not change
% Comparing this to the 2d code requires that we extend the fault along the Y axis far enough that end affects of non plain strain components don't mess with the '2d observations' result. I have extended the fault 500m in either direction along Y. The displacement is 0.01m to avoid the infinitesimal strain breaking down. 
% 
% Proof: This test that the HS solution is correct for the crack in 2d as this code is from disparate parts of the Crouch and Starfield book and is not explicitly written. Steve Martel provided an excellent base code for the 2d hs solution to start from. 



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2D Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
x=-1:0.004:1;                        %List of x points

y= zeros(1,numel(x))-25;%x*0.5;
Pointsxy=[x;y]';
mystruct.line1=(1:(length(Pointsxy(:,1)))-1);
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define full/halfspace and elastic constants.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Defining halfspace. 1 is on, 0 is off. The code will run much faster with
% this off as the calculations are simpler in a full space.
halfspace = 1; 

% Defining the height of the freesurface relative to the value of 0 of
% your imported surfaces. 
freesurface_height = 0;
Pointsxy(:,2)=Pointsxy(:,2)-freesurface_height;


% Defining elastic constants
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,~,nu,~ ] = ElasticConstantsCheck( mu,nu );

%Reshaping the list of XY into usable variables
[MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,~,Fdisp  ]...
 = CreateElements2d( Pointsxy,mystruct,[] );

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
    
Option='A';     
cc=zeros(numel(HalfLength),1); 
Ds  = 0;       Ds	  = cc+(Ds);     
Dn =0.01;       Dn = cc+(Dn); 
Sxx = 0;  Syy = 0;  Sxy = 0;   
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 4: Define dispersed XY observation points to calculate stress strain
%and displacement on. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

X=linspace(0,100,50);
Y=ones(1,numel(X))*-0.01;
dimx = length(X(:,1));
dimy = length(Y(:,1));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate Displacements.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ux2d,Uy2d] = CalculateDisplacementOnSurroundingPoints2d...
    (X,Y,MidPoint,HalfLength,nu,E,halfspace,Ds,Dn,LineNormalVector);

	
	
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3D Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %The file must list the traingle vertex's then triangles. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
 string='FlatBlade_1a_4Faces.ts' ;
 [ Points,Triangles ] = GoCadAsciiReader( string );

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
    
StrikeSlipDisp  =1;      StrikeSlipDisp     = c+StrikeSlipDisp;     %Positive = RightLatMovement (any orientation)
DipSlipDisp     =0;      DipSlipDisp		= c+DipSlipDisp;        %Positive = Reverse movement
TensileSlipDisp =0.01;   TensileSlipDisp	= c+TensileSlipDisp;    %Positive = Opening movement

clear a b c             	%a b and c are used to create an array of slips the size of the number of triangles on the fault surface.

Sxx = 0; 					%Positive is compression 
Syy = 0; 
Szz = 0; 
Sxy = 0; 					%Positive is right lateral movement
Sxz = 0; 					
Syz = 0;   


X=linspace(0,100,50);
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
    
 [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d( StrikeSlipDisp,DipSlipDisp,TensileSlipDisp, nu, X,Y,Z, P1, P2, P3,halfspace);

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
ylabel({'Ground Displacement'})
xlabel({'X distance from the problem centre'})
title({'Half space ground displacement comparison, flat lying sill with constant opening'})
legend('show')
legend('Ux 3D solution','Uz 3D solution','2D solution results')
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

