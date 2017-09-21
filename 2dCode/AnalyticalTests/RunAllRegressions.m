%Runs all the analytical 2d comparison tests and prints results to cmd window
%I recommend you take a look at these individually to see what the
%residual values compare to in the figures. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
%===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/');  end
[~,n] = find(~cellfun(@isempty,strfind(parts,'BEM_DDM_MATLAB'))); %finding the scripts root directory and its location (n)
if ispc; addpath(genpath(strjoin(parts(1,1:n),'\'))); %adding all folders to the path that are in this dir 
else; addpath(genpath(strjoin(parts(1,1:n),'/'))); end;%mac/linux
cd(pathstring)                                       %jumping back to the current dir
%===== add file paths ==========================


%===== SUPRESSING FILE PATH WARNING==========================
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
warning('off', 'Octave:data-file-in-path');
end
%===== SUPRESSING FILE PATH WARNING==========================


%starting timer
tic

disp('Starting 2d and 3d halfspace solution comparison')
run NumericalAnalyticalComp_HalfSpace2d3dTest
disp('End, 2d and 3d halfspace solutions match')
disp(' ') %line break

TakeScreenShotOfFigure(1)

disp(' ') %line break
disp('Starting comparison to analytical solution for a glide dislocation')
run NumericalAnalyticalComp_BarberGlideDislocation_Test
disp('End, 2d codes basic TWODD functions match analytical solution for a glide dislocation')
disp(' ') %line break

TakeScreenShotOfFigure(2)
	
disp(' ') %line break	

disp('Starting test of Pollard and Segall stress loaded crack')
run NumericalAnalyticalComp_XaxisCrack_Opening_StressAndDispInObs
disp('End, Displacement and stress pattern for stress loaded fracture match P&S analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(3)
	
disp(' ') %line break	

disp('Starting Kirsch pressurised hole solution')
run NumericalAnalyticalComp_KirschSolution
disp('End, Displacement and stress pattern for pressurised hole Kirsch analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(4)
	
disp(' ') %line break	

disp('Starting friction analytical solution comparison')
run NumericalAnalyticalComp_LinearFrictionTest
disp('End, Fracture wall displacement profile matches frictional Burgmann/Pollard analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(5)
	
disp(' ') %line break

disp('Starting Savage Gravitational stress in valley solution')
run NumericalAnalyticalComp_GravityComparisonWithSavage_Fortran
disp('End, Stress due to gravational loading in valley matches Savage analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(6)
	
disp(' ') %line break

disp('Starting Inhomogeneous elastic test, (annulus)')
run NumericalAnalyticalComp_InHomogeneousAnnulusCrouch
disp('End, compares well to the Crouch analytical solution of a dual elastic annulus')
disp(' ') %line break

TakeScreenShotOfFigure(7)

%Sets the maxmimised figure size back to default
set(gcf, 'Position', getfield(get(groot,'default'), 'defaultFigurePosition')); 

clear
%ending timer
elapsedTime = toc;

function TakeScreenShotOfFigure(num)
%num is the figure number you want

%grabbing and printing at fullscreensize
set(gcf, 'Position', get(0,'Screensize')); 
print(strcat('Test',num2str(num)),'-dpng')
end
