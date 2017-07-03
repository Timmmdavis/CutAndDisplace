function [StressChange,TotalStress,StressReg,InducedDisplacements]=StressDispOnSurroundingPoints(x,y,xe,ye,a,Pxx,Pyy,Pxy,nu,E,halfspace,Sd,Nd,Beta)
%CalculateStressOnSurroundingPoints calculates the stress by superposition
%   This runs the calculated/defined slip for each element and works out
%   the induced stress in the surrounding points. The stresses for each
%   element are summed to a the total stress for each point within the
%   loop. 
%
%   INPUTS
%   x & y are observation points
%   xe & ye are also element midpoints, the coefficients function looks at
%   how each xe affects each x etc. 
%   a is an array of each each elements half length
%   Pxx, Pyy & Pxy are the remote stress inputs
%   nu is the Poisson's ratio
%   E is the Young's modulus
%   halfspace defines if we work out the coefficientsin a half or whole
%   space
%   Nd&Sd are the Shear and Normal displacements for each element
%   Beta is the orientation of each element relative to the X axis in the
%   C&S convention (See C&S section 2.8 p22 for Beta definition)
%
%   OUTPUTS
%   StressChange is the stress on the points from the event.                       Sxx Syy Sxy
%   TotalStress is the remote stress and stress change from the event added.       Sxx Syy Sxy
%   StressReg is the remote stress                                                 Sxx Syy Sxy
%   InducedDisplacements is the displacements Ux Uy from the dislocations          Ux Uy

%   Copyright 2017, Tim Davis, The University of Aberdeen
% Getting size of observation point array
 n = numel(x);
 
%Setting up loop where the stresses from each element are appended to a master array
 DDM_Obs = zeros(n,5);
 if halfspace==1
    progressbar('Calculating Change in Stress on Observation XYZ HS') % Create figure and set starting time
    for i = 1:size(xe,1)
    DDM_Obs =DDM_Obs + coeffhs_func(x,y,xe(i),ye(i),a(i),Beta(i),Sd(i),Nd(i),nu,E);    
    progressbar(i/size(xe,1)) % Update figure
    end
 else
    progressbar('Calculating Change in Stress on Observation XYZ FS') % Create figure and set starting time
    for i = 1:size(xe,1)
    DDM_Obs =DDM_Obs + coeff_func(x,y,xe(i),ye(i),a(i),Beta(i),Sd(i),Nd(i),nu,E); 
    progressbar(i/size(xe,1)) % Update figure
    end
 end   
 
%Extracting two arrays from the master array of induced stress/displacements 
StressChange=DDM_Obs(:,1:3); 
InducedDisplacements=DDM_Obs(:,4:5); 
 
%Using the regional/driving stress and adding this on top of the stresses
%induced by the elements. This is stored in the array 'Total Stress' 
StressReg = [Pxx(1,:),Pyy(1,:),Pxy(1,:)];
TotalStress=DDM_Obs(:,1:3)+repmat(StressReg,n,1);


% %  Turn on if you want to write an ascii table that can be used elsewhere (i.e. excel) 
% %  Grid nodes in column vector form
%  xg = x(:);			
%  yg = y(:);	
%  StressTensorsTable = table(xg,yg,TotalStress,InducedDisplacements); 
%  writetable(StressTensorsTable);                  %Writes the table


end
