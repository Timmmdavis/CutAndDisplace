%Runs all the analytical 2d comparison tests and prints results to cmd window

%===== SUPRESSING FILE PATH WARNING==========================
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    warning('off', 'Octave:data-file-in-path');
end
%===== SUPRESSING FILE PATH WARNING==========================


%starting timer
tic

disp('Starting 2d and 3d halfspace solution comparison')
run HalfSpace2d3d_Test
disp('End, 2d and 3d halfspace solutions match')
disp(' ') %line break

TakeScreenShotOfFigure(1)

disp(' ') %line break
disp('Starting comparison to analytical solution for a glide dislocation')
run Barber1992_GlideDislocation_Test
disp('End, 2d codes basic TWODD functions match analytical solution for a glide dislocation')
disp(' ') %line break

TakeScreenShotOfFigure(2)
	
disp(' ') %line break	

disp('Starting test of Pollard and Segall stress loaded crack')
run PollardSegall1987_StressLoadedCrack_Test
disp('End, Displacement and stress pattern for stress loaded fracture match P&S analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(3)
	
disp(' ') %line break	

disp('Starting Kirsch pressurised hole solution')
run Kirsch1898_PressurisedHole_Test
disp('End, Displacement and stress pattern for pressurised hole Kirsch analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(4)
	
disp(' ') %line break	

disp('Starting friction analytical solution comparison')
run Burgmann1994_LinearFrictionCrack_Test
disp('End, Fracture wall displacement profile matches frictional Burgmann/Pollard analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(5)
	
disp(' ') %line break

disp('Starting Savage Gravitational stress in valley solution')
run Savage1984_GravityValleyStress_Test
disp('End, Stress due to gravational loading in valley matches Savage analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(6)
	
disp(' ') %line break

disp('Starting Inhomogeneous elastic test, (annulus)')
run CrouchStar1983_InhomogeneousAnnulus_Test
disp('End, compares well to the Crouch analytical solution of a dual elastic annulus')
disp(' ') %line break

TakeScreenShotOfFigure(7)

%Sets the maxmimised figure size back to default
set(gcf, 'Position', getfield(get(groot,'default'), 'defaultFigurePosition')); 

clear
%ending timer
elapsedTime = toc;


