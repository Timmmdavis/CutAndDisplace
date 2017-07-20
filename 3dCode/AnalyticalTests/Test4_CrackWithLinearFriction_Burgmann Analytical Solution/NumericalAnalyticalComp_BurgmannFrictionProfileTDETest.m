
%   Copyright 2017, Tim Davis, The University of Aberdeen%

%Test: 2d Planar fracture with increasing frictional profile: Comparison of
% slip distributions along a planar fracture with linearly increasing
% friction (sliding friction/frictional strength) from 0 at the centre to a
% set value at the tips.
% Also testing slip reduces with normal stress (mu). Havent yet found a
% analytical solution for this.
% 
% Solution: From Burgmann, R., Pollard, D.D. and Martel, S.J., 1994. Slip
% distributions on faults: effects of stress gradients, inelastic
% deformation, heterogeneous host-rock stiffness, and fault interaction.
% Journal of Structural Geology, 16(12), pp.1675-1690. This test calculates
% the analytical slip profile for a planar fracture with its friction
% increasing linearly to its tips. This fracture is subjected to a set
% shear stress. In the code a planar fracture lying along the X axis is
% loaded with a shear stress (Sxy) with 100 elements where each has its own
% set frictional strength. The shearing slip distribution is calculated.
% Checks the residual value is smaller than a threshold.  Res=Slip from
% solver ' Slip from solution.
% 
% Proof: This shows the frictional properties for each element in the code
% is working correctly. It shows with enough elements this results in a
% accurate slip profile for a fracture with friction. This does not need to
% test the resultant stress in the surrounding material as this is
% dependant on each slip contribution and is accurate as shown in the other
% tests.


% Setup: I have put this into 2 parts each with two subparts. Firstly I run
% a tall vertical fracture under Sxy and check it matches the analytical
% slip profile for its strike slip shearing. Secondly I check this fractures
% slip is reduced with a normal stress.
% Third and fourth I run the same two tests but for a flat lying fracture
% long in YY that is loaded under shearing Sxz. For this I check the dip
% slip shearing not strike slip. 








%%
%PART 1.1 Checking a vertical fracture matches the analytical frictional
%profile when loaded under Sxy.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
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
	
	%Script to run dislocations on a fault surface

	% This follows the same unit convention as Poly3d so ''assumes dimensionally consistent units for physical quantities with
	% dimensions of length or stress. Thus if you use kilometers for any quantity with
	% dimensions of length (e.g. a coordinate position), you must use kilometers for all
	% quantities with dimensions of length (e.g. displacement discontinuity).''
	% Andrew Lyell Thomas, Masters thesis, Stanford University. 

    
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
		
		%Strain uses the same convention as engineering stress convention where extension is a positve value. When entering strain values remember this. 
		%For a normal stress: positive (negative) if it produces tension (compression) in the material.
		%For a shear stress: positive if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is + if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear strain Exy will cause Leftlateral movement.
		
        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 string='BladeFaultShort_500Tris_Vertical_NorthDip_0.5.ts'; %low sampling, still passes
 %string='BladeFaultShort_500Tris_Vertical_SouthDip_0.5.ts'; %low sampling, still passes
 %string='BladeFaultShort_2000Tris_Vertical_NorthDip_0.5.ts'; %low sampling, still passes
 %string='BladeFaultShort_2000Tris_Vertical_SouthDip_0.5.ts'; %low sampling, still passes
 [ Points,Triangles ] = GoCadAsciiReader( string,pathstring );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 0; 

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
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp= zeros(size(Triangles(:,1))); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%% 
    %StressInput - Friction
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is tension 
Syy = -0.05; 
Szz = 0; %Syy*nu - plane strain, not too important as this test is main shr stress
Sxy = 1;        			%Positive is right lateral movement  
Sxz = 0;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
X=(MidPoint(:,1));
ne=size(MidPoint(:,1));
Mu  = 0.0;  Mu=repmat(Mu,ne);     %Coefficient of friction
%Creating a linear friction gradient across fractures X axis that increases
%with its X distance. As the fracture has a half length of 1 the max value
%will be 1 at the tips. 
Sf  = abs(X); 
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);
Option='C'; %slip from uniform remote stress with friction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);

%redefine driving stress
Sxy = -1;        			%Positive is right lateral movement  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDispNegDr,DipSlipDispNegDr,TensileSlipDispNegDr,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);

    %%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
    
strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    

    
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
Sxx = 0;       			%Positive is tension
Syy = -0.05; 
Szz = 0;
Sxy = 1;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;
Syz = 0; 
Option='B'; %slip No Fric
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDispNoFr,DipSlipDispNoFr,TensileSlipDispNoFr,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)


%Filtering values from the calculation as  fractures tips along Y have less
%slip so will not match the 2d analytical 'plane strain/stress'  profile 
FracMid=(MidPoint(:,3))<1 & (MidPoint(:,3))>-1; %Only finding slips in the central Y strip of fault
Profile=-StrikeSlipDisp(FracMid); %Grabbing values of slip with friction
Profile2=-StrikeSlipDispNoFr(FracMid); %Grabbing values of slip for no friction
Profile3=StrikeSlipDispNegDr(FracMid); %Grabbing values of slip for no friction
Srt=(MidPoint((FracMid),1)); %Grabbing X points within Y distance to sort on
C = [Srt,Profile,Profile2,Profile3];
jnk = sortrows(C); %Sorted along X axis for slips. 

%Drawing figure of calculates slips
figure;colormap(cmap),scatter(jnk(:,1),abs(jnk(:,2)),24,'k');
hold on
scatter(jnk(:,1),abs(jnk(:,3)),24,'k'); 
scatter(jnk(:,1),abs(jnk(:,4)),18,'b','.'); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
%Linear increasing fricstion model of Burgmann Pollard and Martel 1994.
%Using Points from calculation as observation points. 
%%%%%%
Sg=1;
a = 1;  %Unit half length        
Sr =  abs(Sxy); %Driving stress (positive is extension)
G=mu; %ShearMod
pr=0.25;%Poissions ratio
%Points along crack
x=jnk(:,1);
%Pollard Burgmann_1994_eq5_UniformRemoteStress
Slip_UniformRemote=(2*Sr).*((1-pr)/G).*(sqrt(a.^2-x.^2));
%Pollard Burgmann_1994_eq14_LinearIncreasingFriction
Slip_IncreasingFriction=((1-pr)/G).*(2.*(((Sr.*(sqrt(a.^2-x.^2))))-((Sg.*(1./pi)).*(((sqrt(a.^2-x.^2))+(x.^2./a).*acosh(a./x))))));
Slip_IncreasingFriction=real(Slip_IncreasingFriction);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual, plotting and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Calculating residual
FrictionResidual=Slip_IncreasingFriction-abs(jnk(:,2));
FrictionResidualNegDirStress=Slip_IncreasingFriction-abs(jnk(:,4));

%Drawing figures of analytical profiles
plot(x,Slip_IncreasingFriction,'g','LineWidth',2.5)
plot(x,Slip_UniformRemote,'b','LineWidth',2.5)
plot(x,(FrictionResidual),'--r')

title('Frictional Slip Distribution')
xlabel('Distance from crack centre')
ylabel('Slip')
xlim ([-1 1]);

title({'SlipDistribution';'The numerical result shown here is the strike slip displacement for a vertical lying fracture'});
xlabel('Distance from crack centre')
ylabel({'Magnitude of slip';'No Friction Analytical = green';'Linear Fric Analytical = blue, Residual = dashed red';'ComplimentarySolver&NoFricResult = black + cyan dots'})
xlim ([-1 1]);
hold off

fprintf('MaxSlipDifference %i.\n',max(abs(FrictionResidual(:))))

%Printing error if the match is too far from an expected tolerance. 
if  max(abs(FrictionResidual)) > 0.12 || max(abs(FrictionResidualNegDirStress)) > 0.12 
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, tolerance checks frictional slip residuals and flags differences to analytical profile above 0.11')
end

disp('Now checking mu changes strike slip disp when there is a compressional normal stress on the vertical fracture')   












%%
%PART 1.2 Checking a vertical fractures slip is reduced when I load this
%with a compressional normal stress along with shearing and give it a
%coefficient of friction.

%%


    %%%%%%%%%%%%%% 
    %StressInput - Friction 2, checking normal stress induces changes in
    %crack wall disp
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is tension 
Syy = -2; 
Szz = 0;
Sxy = 1;        			%Positive is right lateral movement  
Sxz = 0;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
ne=size(MidPoint(:,1));
Mu  = 0.2;  Mu=repmat(Mu,ne);     %Coefficient of friction
Sf  = zeros(size(Mu));  
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);
Option='C'; %slip from uniform remote stress with friction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)

MaxDisp=max(abs(StrikeSlipDisp));

if  MaxDisp < 0.9 || MaxDisp > 1
    error('Not a good match, the normal stress does not reduce the friction corrrectly')
else
    disp('Everything looks good. Mu and Normal stress relate')
end

clear a AdRs0 AdRs05 AdRs1 C DipSlipDisp DipSlipDispNegDr DipSlipDispNoFr DirPart FaceNormalVector Fdisp FracMid 
clear freesurface_height FrictionResidualNegDirStress G halfspace jnk lambda mu MidPoint Mu ne nu P1 P2 P3 parts Points pr Profile 
clear Profile2 Profile3 Sf Sg Slip_IncreasingFriction Slip_UniformRemote Sr Srt strain StrikeSlipDisp StrikeSlipDispNegDr 
clear StrikeSlipDispNoFr string Sxx Sxy Sxz Syy Syz Szz TensileSlipDisp TensileSlipDispNegDr TensileSlipDispNoFr Triangles x X














%%
%PART 2.1 Checking a flat fracture matches the analytical frictional
%profile when loaded under Sxz.
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %string='BladeFaultShort_500Tris_nodip.ts';
 %string='BladeFaultShort_500Tris_Flat_WestDip_0.1.ts'; %low sampling, still passes
 string='BladeFaultShort_500Tris_Flat_EastDip_0.1.ts'; %low sampling, still passes
 %string='BladeFaultShort_2000Tris_Flat_WestDip_0.1.ts'; %low sampling, still passes
 %string='BladeFaultShort_2000Tris_Flat_EastDip_0.1.ts'; %low sampling, still passes

 [ Points,Triangles ] = GoCadAsciiReader( string,pathstring );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 0; 

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
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp= zeros(size(Triangles(:,1))); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
    
strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    

    
	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
Sxx = 0;       			%Positive is tension
Syy = 0; 
Szz = -0.05;
Sxy = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 1;
Syz = 0; 
Option='B'; %slip from uniform remote stress with no friction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDispNoFr,DipSlipDispNoFr,TensileSlipDispNoFr,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);

    %%%%%%%%%%%%%% 
    %StressInput - Friction
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is compression 
Syy = 0; 
Szz = -0.05;
Sxy = 0;        			%Positive is right lateral movement  
Sxz = 1;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
X=(MidPoint(:,1));
ne=size(MidPoint(:,1));
Mu  = 0.0;  Mu=repmat(Mu,ne);     %Coefficient of friction
%Creating a linear friction gradient across fractures X axis that increases
%with its X distance. As the fracture has a half length of 1 the max value
%will be 1 at the tips. 
Sf  = abs(X); 
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);
Option='C'; %slip from uniform remote stress with no friction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);

%redefine driving stress
Sxz = -1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDispNegDr,DipSlipDispNegDr,TensileSlipDispNegDr,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)


%Filtering values from the calculation as  fractures tips along Y have less
%slip so will not match the 2d analytical 'plane strain/stress'  profile 
FracMid=(MidPoint(:,2))<1 & (MidPoint(:,2))>-1; %Only finding slips in the central Y strip of fault
Profile=DipSlipDisp(FracMid); %Grabbing values of slip with friction
Profile2=DipSlipDispNoFr(FracMid); %Grabbing values of slip for no friction
Profile3=DipSlipDispNegDr(FracMid); %Grabbing values of slip for no friction
Srt=(MidPoint((FracMid),1)); %Grabbing X points within Y distance to sort on
C = [Srt,Profile,Profile2,Profile3];
jnk = sortrows(C); %Sorted along X axis for slips. 

%Drawing figure of calculates slips
figure;colormap(cmap),scatter(jnk(:,1),abs(jnk(:,2)),24,'k');
hold on
scatter(jnk(:,1),abs(jnk(:,3)),24,'k'); 
scatter(jnk(:,1),abs(jnk(:,4)),18,'b','.'); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
%Linear increasing fricstion model of Burgmann Pollard and Martel 1994.
%Using Points from calculation as observation points. 
%%%%%%
Sg=1;
a = 1;  %Unit half length        
Sr =  abs(Sxz); %Driving stress (positive is tension)
G=mu; %ShearMod
pr=0.25;%Poissions ratio
%Points along crack
x=jnk(:,1);
%Pollard Burgmann_1994_eq5_UniformRemoteStress
Slip_UniformRemote=(2*Sr).*((1-pr)/G).*(sqrt(a.^2-x.^2));
%Pollard Burgmann_1994_eq14_LinearIncreasingFriction
Slip_IncreasingFriction=((1-pr)/G).*(2.*(((Sr.*(sqrt(a.^2-x.^2))))-((Sg.*(1./pi)).*(((sqrt(a.^2-x.^2))+(x.^2./a).*acosh(a./x))))));
Slip_IncreasingFriction=real(Slip_IncreasingFriction);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual, plotting and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Calculating residual
FrictionResidual=Slip_IncreasingFriction-abs(jnk(:,2));
FrictionResidualNegDirStress=Slip_IncreasingFriction-abs(jnk(:,4));

%Drawing figures of analytical profiles
plot(x,Slip_IncreasingFriction,'g','LineWidth',2.5)
plot(x,Slip_UniformRemote,'b','LineWidth',2.5)
plot(x,(FrictionResidual),'--r')

title({'SlipDistribution';'The numerical result shown here is the dip slip displacement for a flat lying fracture'});
xlabel('Distance from crack centre')
ylabel({'Magnitude of slip';'No Friction Analytical = green';'Linear Fric Analytical = blue, Residual = dashed red';'ComplimentarySolver&NoFricResult = black dots'})
xlim ([-1 1]);
hold off
titlesz=25;
fntsz=21;
grid on
ChangeFontSizes(fntsz,titlesz);

fprintf('MaxSlipDifference %i.\n',max(abs(FrictionResidual(:))))

%Printing error if the match is too far from an expected tolerance. 
if  max(abs(FrictionResidual)) > 0.12 || max(abs(FrictionResidualNegDirStress)) > 0.12 
    error('Not a good match to analytical solution')
else
    disp('Everything looks good, tolerance checks frictional slip residuals and flags differences to analytical profile above 0.11')
end

disp('Now checking mu changes dip slip disp when there is a compressional normal stress on the flat lying fracture') 









%%
%PART 2.2 Checking a flat fractures slip is reduced when I load this
%with a compressional normal stress along with shearing and give it a
%coefficient of friction.


    %%%%%%%%%%%%%% 
    %StressInput - Friction 2, checking normal stress induces changes in
    %crack wall disp
	%%%%%%%%%%%%%%
    
Sxx = 0;       			%Positive is tension 
Syy = 0; 
Szz = -2;
Sxy = 0;        			%Positive is right lateral movement  
Sxz = 1;
Syz = 0; 
% %Finding distance of points along X from 0,0. 
ne=size(MidPoint(:,1));
Mu  = 0.2;  Mu=repmat(Mu,ne);     %Coefficient of friction
Sf  = zeros(size(Mu));  
%trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Sf);

Option='C'; %slip from uniform remote stress with no friction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)

MaxDisp=max(abs(DipSlipDisp));

if  MaxDisp < 0.9 || MaxDisp > 1
    error('Not a good match, the normal stress does not reduce the friction corrrectly')
else
    disp('Everything looks good. Mu and Normal stress relate')
end

