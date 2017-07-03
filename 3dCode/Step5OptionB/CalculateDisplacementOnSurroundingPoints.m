function [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints( Ss,Ds,Ts, nu, X,Y,Z, P1, P2, P3,halfspace)
%CalculateDisplacementOnSurroundingPoints This function imports the fault surfaces with the defined/calculated slip and the observation points. These are run through the displacement calculations
%from Mehdi Nikkhoo and Thomas R. Walter published in Geophys. J. Int. (2015) 201, 1119âÿÿ1141. 
%   The calculation output is an array of Xx3 columns where column1 is X column2 is Y and columns is Z. These are then reshaped and drawn with the MATLAB function quiver3. A text file is also
%	created called 'DisplacementTable' that can be used by other software.  

%   Copyright 2017, Tim Davis, The University of Aberdeen

n = numel(X);
DisplacementXYZ=[zeros(n,1),zeros(n,1),zeros(n,1)];

%Runs the script to create displacement across everypoint in the XYZ grid
%defined. Runs for each triangle for each point and adds the tensors together (superposition). 
 if halfspace==1
    progressbar ('Calculating Displacement U HS') % Create figure and set starting time
    for i = 1:size(P1,1)
    DisplacementXYZ = DisplacementXYZ + TDdispHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss(i,:),Ds(i,:),Ts(i,:),nu);
    progressbar(i/size(P1,1)) % Update figure 
    end
 else
    progressbar ('Calculating Displacement U FS') % Create figure and set starting time
    for i = 1:size(P1,1)
    DisplacementXYZ = DisplacementXYZ + TDdispFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Ss(i,:),Ds(i,:),Ts(i,:),nu);
    progressbar(i/size(P1,1)) % Update figure 
    end
 end   

 %Spliiting the output from the calcuation into 3 column vectors. 
Ux = DisplacementXYZ(:,1);
Uy = DisplacementXYZ(:,2);
Uz = DisplacementXYZ(:,3);


%Reshapes the matrix's into columns. Only really needed when the XYZ is a square array. If its vectors already this function does nothing.
X2 = X(:);            
Y2 = Y(:);
Z2 = Z(:);

%Draws a quiver plot showing the displacement around the fault surface. These vectors are scaled 
%automatically by MATLAB and do not correspond to the size of the actual displacement.
% figure('name','Step 5, option B'); quiver3(X2,Y2,Z2,Ux,Uy,Uz);xlabel('x'); ylabel('y'); axis('equal'); title ('U');
% title('Vectors showing displacement of points') 

%Splitting the output tensors into a table with column headers for XYZ and
%the tensors
%DisplacementTable = table(X2,Y2,Z2,Ux,Uy,Uz); 
%writetable(DisplacementTable);                  %Writes the table

end

