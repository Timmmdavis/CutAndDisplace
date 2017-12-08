function [P1,P2,P3] = CreateP1P2P3( Triangles,Points )
% STLReader: From the imported surface each row of the outputs is one of
%                   the corner points of one of the triangles on the mesh.
%                   Each triangle is defined by a seperate row. Unlike the
%                   inputs this is an expanded form where points are
%                   duplicated.
%                   
%               
% usage #1:
% [P1,P2,P3] = CreateP1P2P3( Triangles,Points )
%
% Arguments: (input)
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index.
%                    Not needed unless you want to draw a surface too. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                    locations in "Points" which contains the XYZ location
%                    of each corner of the triangle.
%                    Not needed unless you want to draw a surface too. 
%
% Arguments: (output)
% P1,P2,P3          - n*3 Column vector where each 'P' represents the
%                    different corner points of one of the triangles (XYZ).
%                    Not as efficient in terms of storage but easier to
%                    understand.
%
%
% Example usage:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% [P1,P2,P3] = CreateP1P2P3( Triangles,Points );
% trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),z);
% figure; %Drawing with new indexing:
% sz=numel(P1)/3;
% Mono=(1:1:sz)';
% NewTriangulation=[Mono,Mono+sz,Mono+sz.*2];
% XPnts=[P1(:,1);P2(:,1);P3(:,1)];
% YPnts=[P1(:,2);P2(:,2);P3(:,2)];
% ZPnts=[P1(:,3);P2(:,3);P3(:,3)];
% %Now drawing with a clean monotomic indexing.
% trisurf(NewTriangulation,XPnts,YPnts,ZPnts);
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

TrirowCount =  size(Triangles,1);
P1P2P3 = cell(TrirowCount,3);

%Internal function
[P1P2P3] = CreateP1P2P3_Int( P1P2P3,Triangles,Points,1 );
[P1P2P3] = CreateP1P2P3_Int( P1P2P3,Triangles,Points,2 );
[P1P2P3] = CreateP1P2P3_Int( P1P2P3,Triangles,Points,3 );


%Creates an array for each point defining the triangle (p1p2p3) with 3
%columns that are the point on the traingles 3 vertex's in XYZ.
numericVector = cell2mat(P1P2P3); %extracts the variables P1-P3 into a array
P1=numericVector(1:end,1:3);        %extracts the variables P1-P3 into a array
P2=numericVector(1:end,4:6);
P3=numericVector(1:end,7:9);

%Internal function
function [P1P2P3] = CreateP1P2P3_Int( P1P2P3,Triangles,Points,Val )
%  For function that creates a cell array where each row is a triangle and
%  with each cell holding the XYZ of a vertex
for i = 1:size(Triangles,1) %Assigns 'i' number of rows in triangles array
        P = (Triangles(i,Val));  %Assigns P1 as all rows in the first column of triangles
        for j = 1:size(Points,1) %Finds number of rows in imported table 'S'
            if Points(j,1) == P 
                val = Points(j,2:4); %Checks that if j = i for any rows XYZ in the tables
                P1P2P3{i,Val} = val; %Fills the Grid ID values that correspond in the interpolated table
                break; %Stops so everytime it finds a duplicate XYZ in the tables it fills GRID id and searches again
            end
        end
end
end %end internal func

end

