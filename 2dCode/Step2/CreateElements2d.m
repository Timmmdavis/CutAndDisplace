function [ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng,NUM,XBEG,XEND,YBEG,YEND ] = CreateElements2d( Pointsxy,mystruct,FlipNormalsFlag )
%CreateElements2d Extracts lists of elements from the input points and calculates variables
%needed for calculations 

%   Copyright 2017, Tim Davis, The University of Aberdeen
if nargin < 3
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
for i=1:nf
    stringname{i} = strcat('line', num2str(i));
    tin=Pointsxy(mystruct.(stringname{i})(1,1:1:end-1),1);
    tan=Pointsxy(mystruct.(stringname{i})(1,2:1:end),1);
    ten=Pointsxy(mystruct.(stringname{i})(1:1:end-1),2); 
    tun=Pointsxy(mystruct.(stringname{i})(2:1:end),2);
    
    if FlipNormalsFlag(i)==0
    XBEG = [XBEG ; tin];
    XEND = [XEND ; tan];
    YBEG = [YBEG ; ten];
    YEND = [YEND ; tun];
    else
    XBEG = [XBEG ; tan];
    XEND = [XEND ; tin];
    YBEG = [YBEG ; tun];
    YEND = [YEND ; ten];
    end
    clear tin tan ten tun
end
NUM=numel(XBEG);      %Number of individual elements (segments)

[ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng ] = MidPoint_Orientation( XBEG,XEND,YBEG,YEND,NUM );


end

