% Test: Penny cracks remote stress induced displacement profile
%   
% Solution: Eshelby, J.D., 1957, August. The determination of the elastic
% field of an ellipsoidal inclusion, and related problems. In'Proceedings
% of the Royal Society of London A: Mathematical, Physical and Engineering
% Sciences'(Vol. 241, No. 1226, pp. 376-396). The Royal Society. Rewritten
% in Segall, P., 2010.'Earthquake and volcano deformation. Princeton
% University Press. eq 4.74 3d Analytical solution for the relative
% displacement/slip across a flat circular crack under a remote load. Both
% the shear and normal displacements are calculated. For the TDE solution
% the calculated displacements are converted to radial coordinates so 'r'
% is the distance form the crack centre. These are then plotted against the
% analytical solution.
% 
% Proof: This tests that the remote stress function is calculating the
% correct displacements across the cracks face. its quite slow even for a
% fairly low sampling ~14 mid points in each radial 16th of the circle
% (pi/8 rad). The match gets better with increased sampling much like in 2d
% BEM models.

%   Copyright 2017, Tim Davis, The University of Aberdeen
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
 
 string='CircleMesh_1a_500Faces.ts'; 
 %string='CircleMesh_1a_8000Faces.ts'; 
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
	%%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
    
strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    

	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
Sxx = 0;       			%Positive is tension
Syy = 0; 
Szz = 0;
Sxy = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;
Syz = 1; 
Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);

%Code solution for tensile Disp
Sxx2 = 0;       			%Positive is extension
Syy2 = 0; 
Szz2 = 1;
Sxy2 = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz2 = 0;
Syz2 = 0; 

Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StrikeSlipDisp2,DipSlipDisp2,TensileSlipDisp2,Sxx2,Syy2,Szz2,...
Sxy2,Sxz2,Syz2]=SlipCalculator3d(MidPoint,Sxx2,Syy2,Szz2,Sxy2,Sxz2,Syz2,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);	

TotalShearing = abs(StrikeSlipDisp)+abs(DipSlipDisp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting the calculated slips.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Finding distance from 0,0 for code midpoints for plotting. 
X=(MidPoint(:,1));
Y=(MidPoint(:,2));
[TH,R] = cart2pol(X,Y);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Eshelby Displacement solution for a penny shaped crack in a whole space
%subject to a shearing load
%Equations from Segall, P., 2010. Earthquake and volcano deformation. Princeton University Press.
Ts=Syz(1,1); %Driving shear, could use Sxz for flat penny also
Tn=Szz2(1,1); %Driving normal traction
a=1; %Radius of penny
r=linspace(0,1,250); %Obs points

[us]=Func_SlipPennyCrackSneddon(mu,nu,0,Ts,r);

[ut]=Func_SlipPennyCrackSneddon(mu,nu,Tn,0,r);

%%For plots
%Spacing and points we will interpolate on
RN = linspace(0,0.95,20);
%Making sure we have unique points or the interpolation is unhappy. 
[R2, ~] = uniquify(R);
[TotalShearing, okFlag] = uniquify(TotalShearing);
%Interpolating the numerical results for plotting (Makes it easier to see
%errors)
TensileDispN = interp1(R2,TensileSlipDisp2,RN,'pchip');
TotalShearingN = interp1(R2,TotalShearing,RN,'pchip');
%Drawing the figure
figure;
hold on
%Plotting analytical profiles
plot(r,us,'b','LineWidth',2.5);
plot(r,ut,'g','LineWidth',2.5);
%Plotting interpolated numerical disps every 20th halflength
scatter(RN,TensileDispN,24,'k','filled');
scatter(RN,TotalShearingN,24,'k','filled');
%Adding titles etc
title({'Penny crack slip distributions'})
xlabel('Distance from crack centre')
ylabel('Crack Wall Displacement')
grid on
legend('show')
legend('Shear displacement','Tensile displacement','Numerical results')
%Making the plot nicer
titlesz=25;
fntsz=21;
ChangeFontSizes(fntsz,titlesz);



%Max error of the interpolated points (%)
%Eshelby Displacement solution for a penny shaped crack in a whole space
%subject to a shearing load
%Equations from Segall, P., 2010. Earthquake and volcano deformation. Princeton University Press.
Ts=Syz(1,1); %Driving shear, could use Sxz for flat penny also
Tn=Szz2(1,1); %Driving normal traction
a=1; %Radius of penny
r=linspace(0,1,250); %Obs points
%One side of crack (eq 4.74 segall)
[us]=Func_SlipPennyCrackSneddon(mu,nu,0,Ts,RN);
[ut]=Func_SlipPennyCrackSneddon(mu,nu,Tn,0,RN);

PercentErrorOpeningInterpPnts=((100./ut).*TensileDispN)-100;
PercentErrorShrInterpPnts=((100./us).*TotalShearingN)-100;


fprintf('MaxPercentFromAnalyticalSolutionTs %i.\n',max(abs(PercentErrorOpeningInterpPnts(:)))) %prints on one line unlike 'disp
fprintf('MaxPercentFromAnalyticalSolutionSS %i.\n',max(abs(PercentErrorShrInterpPnts(:)))) %prints on one line unlike 'disp
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Returns error if distance is too high
P=[max(abs(PercentErrorOpeningInterpPnts(:))),max(abs(PercentErrorShrInterpPnts(:)))];
if any(P>20)
    error('Calculated displacements are a long way from the analytical solutions > 20%')
else
    disp('Everything looks good, tolerance checks slip/opening residuals and flags errors above 0.1')
end



