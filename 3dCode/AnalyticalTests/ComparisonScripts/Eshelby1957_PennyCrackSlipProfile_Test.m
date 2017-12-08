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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

%Clearing old figures and arrays. 
clear;close all

%Loading a colourmap to be used for figures.
cmap = colormap_cpt('Ccool-warm');
cmap2 = colormap_cpt('Ccool-warm2');
		
		
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Import the fault surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option A = Loading Gocad ascii data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
 string='CircleMesh_1a_500Faces.ts'; 
 %string='CircleMesh_1a_8000Faces.ts'; 
 [ Points,Triangles ] = GoCadAsciiReader( string );


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

% Defining elastic constants
mu=1;           	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

% Getting surface ready for script. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 

%Locked Els
Fdisp= zeros(size(Triangles(:,1))); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, remote (with
%friction), remote with no opening components, tractions on the elements
%(internal pressure) etc.
%
%If Cartesian stress tensor components are not known but principal stresses
%are use the function 'TensorsFromPrincipal' to calculate the tensors in a
%Cartesian reference frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;

	%%%%%%%%%%%%%%
    %StrainInput        
	%%%%%%%%%%%%%%
    
%Put to 1 to define the stresses defined in 'stress input' as strain values    
strain=0;   

	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
Sxx = 0;       			
Syy = 0; 
Szz = 0;
Sxy = 0;        		
Sxz = 0;
Syz = 1; 
Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);

TotalShearing = sqrt((Dss).^2+(Dds).^2);


%Code solution for tensile Disp
Sxx2 = 0;       			%Positive is extension
Syy2 = 0; 
Szz2 = 1;
Sxy2 = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz2 = 0;
Syz2 = 0; 

Option='B'; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculating slip due to boundary conditions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dss2,Dds2,Dn2,~,~,~,~,~,~]=SlipCalculator3d(MidPoint,Sxx2,...
Syy2,Szz2,Sxy2,Sxz2,Syz2,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);


%Finding distance from 0,0 for code midpoints for plotting. 
X=(MidPoint(:,1));
Y=(MidPoint(:,2));
[TH,R] = cart2pol(X,Y);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Eshelby Displacement solution for a penny shaped crack in a whole space
%subject to a shearing load
%Equations from Segall, P., 2010. Earthquake and volcano deformation. Princeton University Press.
Ts=Syz(1,1); %Driving shear, could use Sxz for flat penny also
Tn=Szz2(1,1); %Driving normal traction
a=1; %Radius of penny
r=linspace(0,1,250); %Obs points

[us]=Eshelby1957_PennyCrackSlipProfile(mu,nu,0,Ts,r);
[ut]=Eshelby1957_PennyCrackSlipProfile(mu,nu,Tn,0,r);

%%For plots
%Spacing and points we will interpolate on
RN = linspace(0,0.95,20);
%Making sure we have unique points or the interpolation is unhappy. 
[R2, ~] = uniquify(R);
[TotalShearing, okFlag] = uniquify(TotalShearing);
%Interpolating the numerical results for plotting (Makes it easier to see
%errors)
TensileDispN = interp1(R2,Dn2,RN,'pchip');
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
[us]=Eshelby1957_PennyCrackSlipProfile(mu,nu,0,Ts,RN);
[ut]=Eshelby1957_PennyCrackSlipProfile(mu,nu,Tn,0,RN);

PercentErrorOpeningInterpPnts=((100./ut).*TensileDispN)-100;
PercentErrorShrInterpPnts=((100./us).*TotalShearingN)-100;


fprintf('MaxPercentFromAnalyticalSolutionTs %i.\n',max(abs(PercentErrorOpeningInterpPnts(:)))) %prints on one line unlike 'disp
fprintf('MaxPercentFromAnalyticalSolutionSS %i.\n',max(abs(PercentErrorShrInterpPnts(:)))) %prints on one line unlike 'disp
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating residual and checking for errors. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Returns error if distance is too high
P=[max(abs(PercentErrorOpeningInterpPnts(:))),max(abs(PercentErrorShrInterpPnts(:)))];
if any(P>20)
    error('Calculated displacements are a long way from the analytical solutions > 20%')
else
    disp('Everything looks good, tolerance checks slip/opening residuals and flags errors above 0.1')
end



