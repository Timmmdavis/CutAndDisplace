function [P1,P2,P3] = CreateP1P2P3( Triangles,Points )
%CreateP1P2P3 From the imported surface this creates an array for each point of each triangle (P1P2P3) where P1 has three columns XYZ. 
% Each row in P1P2P3 correspond to the vertex's of one triangle that was originally defined by indexing in the array 'triangles'. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
MaxZ = max(Points(:,4));
%disp 'Maximum Z value of fault surface points', disp (MaxZ);
fprintf('Maximum Z value of fault surface points %i.\n',MaxZ) %prints on one line unlike 'disp

TrirowCount =  size(Triangles,1);
P1P2P3 = cell(TrirowCount,3);

%  For function that creates a cell array where each row is a triangle and
%  with each cell holding the XYZ of a vertex
for i = 1:size(Triangles,1) %Assigns 'i' number of rows in triangles array
        P1 = (Triangles(i,1));  %Assigns P1 as all rows in the first column of triangles
        for j = 1:size(Points,1) %Finds number of rows in imported table 'S'
            if Points(j,1) == P1 
                val = Points(j,2:4); %Checks that if j = i for any rows XYZ in the tables
                P1P2P3{i,1} = val; %Fills the Grid ID values that correspond in the interpolated table
                break; %Stops so everytime it finds a duplicate XYZ in the tables it fills GRID id and searches again
            end
        end
end
for i = 1:size(Triangles,1) %Assigns 'i' number of rows in triangles array
        P1 = (Triangles(i,2));  %Assigns P1 as all rows in the first column of triangles
        for j = 1:size(Points,1) %Finds number of rows in imported table 'S'
            if Points(j,1) == P1 
                val = Points(j,2:4); %Checks that if j = i for any rows XYZ in the tables
                P1P2P3{i,2} = val; %Fills the Grid ID values that correspond in the interpolated table
                break; %Stops so everytime it finds a duplicate XYZ in the tables it fills GRID id and searches again
            end
        end
end
for i = 1:size(Triangles,1) %Assigns 'i' number of rows in triangles array
        P1 = (Triangles(i,3));  %Assigns P1 as all rows in the first column of triangles
        for j = 1:size(Points,1) %Finds number of rows in imported table 'S'
            if Points(j,1) == P1 
                val = Points(j,2:4); %Checks that if j = i for any rows XYZ in the tables
                P1P2P3{i,3} = val; %Fills the Grid ID values that correspond in the interpolated table
                break; %Stops so everytime it finds a duplicate XYZ in the tables it fills GRID id and searches again
            end
        end
end

%Creates an array for each point defining the triangle (p1p2p3) with 3
%columns that are the point on the traingles 3 vertex's in XYZ.
numericVector = cell2mat(P1P2P3); %extracts the variables P1-P3 into a array
P1=numericVector(1:i,1:3);        %extracts the variables P1-P3 into a array
P2=numericVector(1:i,4:6);
P3=numericVector(1:i,7:9);

end

