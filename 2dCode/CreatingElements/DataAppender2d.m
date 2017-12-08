function [AppendedPoints,mystruct,Flag] = DataAppender2d( PntsInput1,PntsInput2,mystruct,Flag,Value )
% DataAppender2d: Appending lists of elements from different fractures, for
%                   getting data ready for the function 'CreateElements2D'.
%               
% usage #1: First use, appending some points to some others.
%[AppendedPoints,mystruct,Flag] = DataAppender2d(...
%PntsInput1,PntsInput2,mystruct)
%
% usage #2: 2nd use, appending for a second time and adding to flag. 
%[AppendedPoints,mystruct,Flag] = DataAppender2d(...
%PntsInput1,PntsInput2,mystruct,Flag )
%
% usage #3: nth use, you care about the value of the flag.
%[AppendedPoints,mystruct,Flag] = DataAppender2d(...
%PntsInput1,PntsInput2,mystruct,Flag,Value )
%
% Arguments: (input)
% PntsInput1       - An n*2 vector of the XY locations of each elements end
%                   points, each set of consecutive rows are one element.
%                   The total list is for a single fracture. 
%
% PntsInput2       - An n*2 vector of the XY locations of each elements end
%                   points, each set of consecutive rows are one element.
%                   The total list is for a single fracture. 
%
% mystruct         - Index's of each seperate fracture, fields such as
%                   line1 = 1:601 say that the first 601 points in Pointsxy
%                   are elements, i.e. connected to the rows above and
%                   below so a single 'connected' fracture. We assume this
%                   has only been created for the points of input 'PntsInput1'
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
%  AppendedPoints  - The two input points vectors appended
%
%  mystruct        - The mystruc brought in the but with additional indexes
%                   of points two in a new field in this structure. 
%
%  Flag            - The output flag with numbering chosen by the user. For
%                   the first use of this function with this not called on
%                   the input this will be a flag of ones at the index
%                   locations of input 'PntsInput2'
%
% Example usage:
%
% % Creating a circle with a square of fixed inner points. 
% % Circle
% ri = 1;k = 0;
% for theta = 0:pi/300:2*pi;k = k+1;x(k) = ri*cos(theta);y(k) = ri*sin(theta);end
% sz=numel(x);Pointsxy=[x',y'];
% mystruct.line1=(1:sz);
% % Square of points inside circle
% rr=0.8;
% x=[-rr,0,rr,0,-rr];y=[0,-rr,0,rr,0];
% PointsxyF=[x;y]';
%
% [Pointsxy,mystruct,Fdisp] = DataAppender2d( Pointsxy,PointsxyF,mystruct );
%
% [xe,ye,HalfLength,Points,LineNormalVector] = CreateElements2d( Pointsxy,mystruct );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

if nargin>=3
%Calculating the number of fractures now
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms
else 
nf=0;
end

%Number of original points that make the tris for Input 1
Sz1=numel(PntsInput1(:,1));    
Sz2=numel(PntsInput2(:,1));

%Appending for output
AppendedPoints=[PntsInput1;PntsInput2];

stringname = strcat('line', num2str(nf+1));
mystruct.(stringname)=(Sz1+1:Sz1+Sz2);


if nargin<4
    %Zero CV
    Flag=zeros(Sz1-1,1);
    Value=1;
end

FlagAppend=zeros(Sz2-1,1);
% %Finding the current max of the flag, the 2nd input has to be 1 higher. 
% n=max(Flag)+1;
%Vaues of 2nd input are now n. 
Flag=[Flag;FlagAppend+Value];




end

