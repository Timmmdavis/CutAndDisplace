function [ CartSxx,CartSyy,CartSzz,CartSxy,CartSxz,CartSyz]...
    =TensorsFromPrincipal3d( S1,S2,S3,DirS1 )
% TensorsFromPrincipal3d: Function that computes the Cartesian stress
%                   tensor from the magnitude and directions of defined
%                   principal stresses. Use to compute tensors for the BEM
%                   boundary conditions. This prints these to your console.
%
%                   The 3rd stress direction (most compressive is presumed
%                   to be vertical). 
%
%               
% usage #1: Prints results to cmd window
% TensorsFromPrincipal3d( S1,S2,S3,DirS1 )
%
% usage #2: Outputs these for futher use
% [ CartSxx,CartSyy,CartSzz,CartSxy,CartSxz,CartSyz]...
%    =TensorsFromPrincipal3d( S1,S2,S3,DirS1 )
%
% Arguments: (input)
%     S1,S2,S3       - The magnitudes of the principal stresses. 
%
%     DirS1       - The azimuth of the S1 component (degrees). 
%                  (angle away from N clockwise)
%
% Arguments: (output)
% CartSxx,CartSyy,CartSzz
% CartSxy,CartSxz,CartSyz- The magnitudes of the calculated Cartesian
%                          stress tensor.
%
% Example usage:
%
%  S1=0;
%  S2=-2;
%  S3=0;
%  DirS1=45;
%  TensorFromPrincipal2d( S1,S2,S3,DirS1 );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%% DefineMagnitude
%Horizontal stresses mags
Sz1= S1;    %Magnitude of component 1
Sz2= S2;	%Magnitude of component 2
%Vertical stress mag
Sz3= S3; 	%Magnitude of component 3

%% DefineDirection
Dir1=90-DirS1; 	%direction (az) of component 1 (degrees)
Dir2=90+Dir1; 	%direction (az) of component 2 (degrees)
%Dir 3 we presume is vertical so do not define this

%% Checking Dir 1&2 are 90 degrees apart
Test=round((asind(abs(sind(Dir1)))+asind(abs(sind(Dir2)))),8)==90; %rounding to 8 digits as precision can be an issue with this test
if Test==0
    error('Input Azimuth angles are not 90 degrees apart')
end

%% Calculating direction cosines
Ax1=abs(sind(Dir1)); %Angle away from X axis for Dir 1
Ay1=abs(cosd(Dir1)); %Angle away from Y axis for Dir 1
Az1=0;
if Dir1>180
	Ax1=-Ax1;
end
%Correcting sign for Y
if 	Dir1>90 && Dir1<270
	Ay1=-Ay1;
end

Ax2=abs(sind(Dir2)); %Angle away from X axis for Dir 1
Ay2=abs(cosd(Dir2)); %Angle away from Y axis for Dir 1
Az2=0;
if Dir2>180
	Ax2=-Ax2;
end
%Correcting sign for Y
if 	Dir2>90 && Dir2<270
	Ay2=-Ay2;
end

%Dipping directly down cosines
Ax3=0;
Ay3=0;
Az3=1;

%accumulating new axes to direction cosines
Xaxis=[Ax1,Ax2,Ax3];
Yaxis=[Ay1,Ay2,Ay3];
Zaxis=[Az1,Az2,Az3];

[ CartSxx,CartSyy,CartSzz,CartSxy,CartSxz,CartSyz ]...
    = StressTensorTransformation3d(Sz1,Sz2,Sz3,0,0,0,Xaxis,Yaxis,Zaxis);


if nargout==0
    %Printing to command window, you can just cut and paste this into the
    %script to use as the driving stress
    disp('Calculated cartesian tensors:')
    fprintf('Sxx= %i;\n',CartSxx)
    fprintf('Syy= %i;\n',CartSyy) 
    fprintf('Szz= %i;\n',CartSzz) 
    fprintf('Sxy= %i;\n',CartSxy) 
    fprintf('Sxz= %i;\n',CartSxz) 
    fprintf('Syz= %i;\n',CartSyz) 
else
    %do nothing, you get these out the script and probably dont need these
    %displaying
end
end
