function [AppendedPoints,mystruct,Flag] = DataAppender2d( PntsInput1,PntsInput2,mystruct,Flag,Val )
%DataAppender3d Appending 2 to the back of 1. Flag says where the data in two lies.  
%   PntsInput1,PntsInput2,TrisInput1,TrisInput2 are the point and triangle
%   inputs. 1 is the data in 1. 2 is two. Flag values increase with each
%   use of this function if a flag from previous calls is brought in. If
%   not the places containing the data in two are flagged as 1's. 

%   Copyright 2017, Tim Davis, The University of Aberdeen


%We append 2 to the back of 1. 

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
    Val=1;
end

FlagAppend=zeros(Sz2-1,1);
% %Finding the current max of the flag, the 2nd input has to be 1 higher. 
% n=max(Flag)+1;
%Vaues of 2nd input are now n. 
Flag=[Flag;FlagAppend+Val];




end

