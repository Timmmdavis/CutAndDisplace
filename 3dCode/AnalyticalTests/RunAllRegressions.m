%Runs all the analytical 2d comparison test and prints results to cmd window
%I recommend you take a look at these individually to see what the
%residual values compare to in the figures. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
% %===== add file paths ==========================
pathstring = pwd;                                   %Get the address of the current working directory
if ispc; parts = strsplit(mfilename('fullpath'), '\');       %Getting the address of the script working directory and splitting into cell array
else; parts = strsplit(mfilename('fullpath'), '/'); end 
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

disp('Starting 3d ground surface displacement comparison, Okada solution vs TDE')
run NumericalAnalyticalComp_OkadaGroundDispTDETest
disp('End, TDE and Okada solutions match')
disp(' ') %line break

TakeScreenShotOfFigure(1)
	
disp(' ') %line break	
disp('Starting Mogi ground surface displacement comparison')
run NumericalAnalyticalComp_MogiGroundDispTDETest
disp('End, TDE Ground surface displacement matches Mogi point source analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(2)

disp(' ') %line break	
disp('Starting Penny shaped crack crack wall displacements. Segall Analytical solution')
run NumericalAnalyticalComp_Mode1_2_PennyTDETest
disp('End, Crack wall displacements fit closely to the analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(3)
	
disp(' ') %line break	
disp('Starting friction analytical solution comparison (3d against 2d solution)')
run NumericalAnalyticalComp_BurgmannFrictionProfileTDETest
disp('End, 3d Fracture wall displacement profile matches frictional Burgmann/Pollard analytical solution')
disp(' ') %line break

close;close;close;close;close; %getting to the actual comp to an solution
TakeScreenShotOfFigure(4)
	
disp(' ') %line break
disp('Starting Savage Gravitational stress in valley solution (3d against 2d solution)')
run NumericalAnalyticalComp_SavageTDETest
disp('End, Stress due to gravational loading in valley matches 2d Savage analytical solution')
disp(' ') %line break

TakeScreenShotOfFigure(5)
	
disp(' ') %line break	
disp('Starting Inhomogeneous elastic test, (annulus)')
run NumericalAnalyticalComp_InHomogeneousAnnulusCrouch3d
disp('End, compares well to the Crouch analytical solution of a dual elastic annulus (2d)')
disp(' ') %line break

TakeScreenShotOfFigure(6)
    
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
    if  isOctave==1;
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
% TakeScreenShotOfFigure(7)

disp(' ') %line break	
disp('Starting Traction free surface test')
run TractionFreeSurfaceTest
disp('End, The surfaces are all free of tractions')
disp(' ') %line break

clear
%ending timer
elapsedTime = toc;

function TakeScreenShotOfFigure(num)
%num is the figure number you want

%grabbing and printing at fullscreensize
set(gcf, 'Position', get(0,'Screensize')); 
print(strcat('Test',num2str(num)),'-dpng')
end
