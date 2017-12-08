function [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d(Ss,Ds,Ts,nu,X,Y,Z,P1,P2,P3,halfspace)
% CalculateDisplacementOnSurroundingPoints3d: calculates the displacement by
%				superposition.
%               This runs the calculated/defined slip for each element and
%               works out the induced disp in the surrounding points. The
%               disps for each element are summed to a the total disp
%               for each point within the loop.
%
%               Commented lines at the base of the script allow writing of
%               a table of the stress tensors. 
%
% usage #1:
% [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d(Ss,Ds,Ts,nu,X,Y,Z...
% ,P1,P2,P3,halfspace)
%
% Arguments: (input)
%    X,Y & Z  - The observation points
%
% P1,P2,P3    - n*3 Column vectors where each 'P' represents the
%               different corner points of one of the triangles (XYZ).
%
%  Dss,Dds,Dn - Vectors that describe how much the elements displace in the
%               normal (Dn) and strike slip (Dss) and dipslip (Dds)
%               directions on the elements.
%
%       nu    - The materials Poisson's ratio
%
%
%  halfspace  - Defines if we work out the coefficientsin a half or whole
%              space
%
% Arguments: (output)
%  Ux,Uy,Uz   - Displacement in X,Y and Z at each observation point (XYZ
%  input).
%
%
% Example usage:
% %Assuming the inputs are already defined.
% [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d(Ss,Ds,Ts,nu,X,Y,Z...
% ,P1,P2,P3,halfspace)
%
% figure,quiver3(X(:),Y(:),Z(:),Ux(:),Uy(:),Uz(:));
% xlabel('x'); ylabel('y'); zlabel('z'); axis('equal');
% title('Vectors showing displacement of points') ;
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

n = numel(X);
DispObs=[zeros(n,1),zeros(n,1),zeros(n,1)];

%Setting up progress bar
if halfspace==1
    progressbar('Calculating Displacement U HS') % Create figure and set starting time
else 
    progressbar('Calculating Displacement U FS') % Create figure and set starting time
end

NumPnts=size(P1,1);
%Runs the script to create displacement across everypoint in the XYZ grid
%defined. Runs for each triangle for each point and adds the tensors together (superposition). 
if halfspace==1
    for i = 1:NumPnts
        DispObs = DispObs + TDdispHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss(i,:),Ds(i,:),Ts(i,:),nu);
        progressbar(i/NumPnts) % Update figure 
    end
else
    for i = 1:size(P1,1)
        DispObs = DispObs + TDdispFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss(i,:),Ds(i,:),Ts(i,:),nu);
        progressbar(i/NumPnts) % Update figure 
    end
end   

 %Spliiting the output from the calcuation into 3 column vectors. 
Ux = DispObs(:,1);
Uy = DispObs(:,2);
Uz = DispObs(:,3);


% %Reshapes the matrix's into columns. Only really needed when the XYZ is a square array. If its vectors already this function does nothing.
% X2 = X(:);            
% Y2 = Y(:);
% Z2 = Z(:);

%Splitting the output tensors into a table with column headers for XYZ and
%the tensors
%DisplacementTable = table(X2,Y2,Z2,Ux,Uy,Uz); 
%writetable(DisplacementTable);                  %Writes the table

end

