%Runs all the analytical 3d comparison tests and prints results to cmd window

%===== SUPRESSING FILE PATH WARNING==========================
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    warning('off', 'Octave:data-file-in-path');
end
%===== SUPRESSING FILE PATH WARNING==========================

%starting timer
tic

disp('Starting 3d ground surface displacement comparison, Okada solution vs TDE')
run Okada1985_RectangularDislocation_Test
disp('End, TDE and Okada solutions match')
disp(' ') %line break

TakeScreenShotOfFigure(1)
	
disp(' ') %line break	
disp('Starting Mogi ground surface displacement comparison')
run Mogi1958_SphericalCavity_Test
disp('End, TDE Ground surface displacement matches Mogi point source analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(2)

disp(' ') %line break	
disp('Starting Penny shaped crack crack wall displacements. Segall Analytical solution')
run Eshelby1957_PennyCrackSlipProfile_Test
disp('End, Crack wall displacements fit closely to the analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(3)
	
disp(' ') %line break	
disp('Starting friction analytical solution comparison (3d against 2d solution)')
run Burgmann1994_LinearFrictionCrack_3dTest
disp('End, 3d Fracture wall displacement profile matches frictional Burgmann/Pollard analytical solution')
disp(' ') %line break

close;close;close; %getting to the actual comp to an solution
TakeScreenShotOfFigure(4)
	
disp(' ') %line break
disp('Starting Savage Gravitational stress in valley solution (3d against 2d solution)')
run Savage1984_GravityValleyStress_3dTest
disp('End, Stress due to gravational loading in valley matches 2d Savage analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(5)
	
disp(' ') %line break	
disp('Starting Inhomogeneous elastic test, (annulus)')
run CrouchStar1983_InhomogeneousAnnulus_3dTest
disp('End, compares well to the Crouch analytical solution of a dual elastic annulus (2d)')
disp(' ') %line break

TakeScreenShotOfFigure(6)


disp('Starting stress intensity tests')
run StressIntensity3dTest
disp('End, compares well to analyical stress intensity for inclined 3d/curved 2d crack')
disp(' ') %line break

TakeScreenShotOfFigure(7)


    
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    disp('eshelby solution doesnt work in octave');
return
end

    %%%Will need to put meng/pollard eshelby functions in folder before you
    %%%turn this on
% disp(' ') %line break	
% disp('Eshelby Void test')
% run NumericalAnalyticalComp_EshelbyVoid
% disp('End, compares well to the Meng analytical code for the stress in the matrix surronging a void (Eshelby)')
% disp(' ') %line break
% TakeScreenShotOfFigure(~)

disp(' ') %line break	
disp('Starting Traction free surface test')
run TractionFreeSurfaceTest
disp('End, The surfaces are all "free" of tractions')
disp(' ') %line break

clear
%ending timer
elapsedTime = toc;
