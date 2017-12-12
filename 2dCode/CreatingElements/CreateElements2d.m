function [ MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp ] ...
    = CreateElements2d( Pointsxy,mystruct,Fdisp,FlipNormalsFlag )
% CreateElements2d: Extracts lists of elements from
%               the input points and calculates variables
%               needed for calculations
%               
% usage #1:
% [ MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp]...
% = CreateElements2d( Pointsxy,mystruct,Fdisp,FlipNormalsFlag )
%
% [MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp...
% ] = CreateElements2d( Pointsxy,mystruct,Fdisp )
%
% Arguments: (input)
% Pointsxy         - An n*2 vector of the XY locations of the elements end
%                   points, each set of consecutive rows are one element
%                   when these are indexed by a single field in
%                   'mystruct'. 
%
% mystruct         - Index's of each seperate fracture, fields such as
%                   line1 = 1:601 say that the first 601 points in Pointsxy
%                   are elements, i.e. connected to the rows above and
%                   below so a single 'connected' fracture. 
%
% FlipNormalsFlag  - Flag the length of fields's in mystruct. If this is one
%                   the normals of this fracture is flipped.
%
% Fdisp           - Flag to say if elements will be locked. If this comes
%                   in as an empty array it is created. 
%
%
% Arguments: (output)
%  MidPoint        - The x and y locations of each elements midpoint [x,y]. 
%
%  HalfLength      - The half length of each element
%
%  P1,P2           - The start (P1) and end (P2) points of each separate 
%                   element in [x,y] coordiantes. First col is x. 
%
% LineNormalVector - The direction cosines, CosAx (Nx) and CosAy in a list
%                   for each element. 
%
% Points            - Columns 2 3 are the XYZ locations of one the
%                    corner points of a line. Column 1 is the index.
%
% Lines             -  Lines is a list where each row contains 2 index
%                    locations in "Points" which contains the XY location
%                    of each corner of the line.
%
% Fdisp           - Flag to say if elements will be locked. If this comes
%                   in as an empty array it is created. 
%
%
% Example usage:
%
% x=linspace(-0.5,0.5,15);    
% y=zeros(1,numel(x)); 
% Pointsxy=[x;y]';
% mystruct.line1=(1:(length(Pointsxy(:,1))));
%
% [ xe,ye,HalfLength,Points,LineNormalVector...
% ] = CreateElements2d( Pointsxy,mystruct )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

if nargin < 4
    names = fieldnames(mystruct);
    numstructs=numel(names);
    FlipNormalsFlag=zeros(numstructs,1); %nothing is flipped
end


%Calculating the number of fractures now
Fnms=fieldnames(mystruct);
nf=numel(Fnms);clear Fnms

%initilise vars
XBEG=[];
XEND=[];
YBEG=[];
YEND=[];
Lines=[];
for i=1:nf
    %Grab the data out the structure
    stringname = strcat('line', num2str(i));
    tin=Pointsxy(mystruct.(stringname)(1,1:1:end-1),1);
    tan=Pointsxy(mystruct.(stringname)(1,2:1:end),1);
    ten=Pointsxy(mystruct.(stringname)(1:1:end-1),2); 
    tun=Pointsxy(mystruct.(stringname)(2:1:end),2);
    
    if FlipNormalsFlag(i)==0
        XBEG = [XBEG ; tin]; XEND = [XEND ; tan];
        YBEG = [YBEG ; ten]; YEND = [YEND ; tun];
    else %Put change Xbeg for Xend etc..
        XBEG = [XBEG ; tan]; XEND = [XEND ; tin];
        YBEG = [YBEG ; tun]; YEND = [YEND ; ten];
    end
    %Build Index for PointsXY. 
    if i==1
        Lines=[Lines;[(1:numel(tin))',(2:numel(tin)+1)']];
    else
        Lines=[Lines;[(1:numel(tin))',(2:numel(tin)+1)']+1+numel(Lines)/2];
    end
    
    
    clear tin tan ten tun
end

%Prep for output
Lines=[(1:1:numel(Lines)/2)',Lines];

%Reshape PointsXY
bad=isnan(Pointsxy(:,1));
Pointsxy(bad,:)=[];
Points=reshape(Pointsxy,[],2);

%Call Func
[ MidPoint,HalfLength,P1,P2,LineNormalVector ]...
    = MidPoint_Orientation( XBEG,XEND,YBEG,YEND );

%Check to see if any line segments are invalid. 
if any(HalfLength==0)
    error('Some line segments have 0 length')
end  


if isempty(Fdisp)
    Fdisp=zeros(numel(P1(:,1)),1);
end

end

