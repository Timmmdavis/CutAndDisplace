function [AppendedPoints,AppendedTriangles,Flag] = DataAppender3d( PntsInput1,PntsInput2,TrisInput1,TrisInput2,Flag,Value )
% DataAppender2d: Appending lists of elements from different fractures, for
%                 getting data ready for the function 'MidPointCreate'
%                 and later functions in the script. New surface variables
%                 will have the correct indexing.
%
%                 Appending inputs 2 to the back of inputs 1. Output flag says
%                 where the data in two lies.
%               
% usage #1: First use, appending some points to some others.
%[AppendedPoints,AppendedTriangles,Flag] = DataAppender3d(...
%PntsInput1,PntsInput2,TrisInput1,TrisInput2)
%
% usage #2: 2nd use, appending for a second time and adding to flag. 
%[AppendedPoints,AppendedTriangles,Flag] = DataAppender3d(...
%PntsInput1,PntsInput2,TrisInput1,TrisInput2,Flag)
%
% usage #3: nth use, you care about the value of the flag.
%[AppendedPoints,AppendedTriangles,Flag] = DataAppender3d(...
%PntsInput1,PntsInput2,TrisInput1,TrisInput2,Flag,Value )
%
% Arguments: (input)
% PntsInput1       - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index.
%                    This is the inputs for the first surface.
%
% PntsInput2       - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index.
%                    This is the inputs for the second surface.
%
% TrisInput1/TrisInput2  -  TrisInput1 and TrisInput2 is a list where each row
%                   contains 3 index locations in "Points" which contains
%                   the XYZ location of each corner of the triangle. These
%                   correspond to the PntsInput1 & 2 respectivly. 
%
% Flag             -  Flag values increase with each
%                   use of this function if a flag from previous calls is
%                   brought in. If not the places containing the data in
%                   two are flagged as 1's.
%
% Value            - If you want the value of the 'Flag' for the new points
%                   to be a certain number set this with this variable. 
%
%
% Arguments: (output)
%  AppendedPoints  - The two input points vectors appended and indexed
%                   correctly for the new triangles array. 
%
% AppendedTriangles- The new array of triangles for both surfaces merged. 
%
%  Flag            - The output flag with numbering chosen by the user. For
%                   the first use of this function with this not called on
%                   the input this will be a flag of ones at the index
%                   locations of input 'PntsInput2'
%
% Example usage:
% 
% % Appending two seamount surfaces and drawing, colour based on 'Flag'
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Traingles2=Triangles;
% Points=[[1:numel(x)]',x(:),y(:),z(:)+1];
% Points2=[[1:numel(x)]',x(:),y(:),-z(:)];
% [AppendedPoints,AppendedTriangles,Flag] = DataAppender3d( Points,Points2,Triangles,Traingles2 );
% PlotSlipDistribution3d(AppendedTriangles,AppendedPoints,[],Flag);
% 
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Number of original points that make the tris for Input 1
Sz=numel(PntsInput1(:,1));    

%Now adding this to the 2nd set so the 2nd set have the right numbers
PntsInput2(:,1)= PntsInput2(:,1)+Sz;
TrisInput2=TrisInput2+Sz; 

%Appending for output
AppendedPoints=[PntsInput1;PntsInput2];
AppendedTriangles=[TrisInput1;TrisInput2];

%Getting sizes
Sz1=numel(TrisInput1(:,1));
Sz2=numel(TrisInput2(:,1));

if nargin<5
    %Zero CV
    Flag=zeros(Sz1,1);
    Value=1;    
end
FlagAppend=zeros(Sz2,1);


%Vaues of 2nd input are now n. 
Flag=[Flag;FlagAppend+Value];


end

