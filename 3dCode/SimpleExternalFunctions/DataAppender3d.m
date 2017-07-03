function [AppendedPoints,AppendedTriangles,Flag] = DataAppender3d( PntsInput1,PntsInput2,TrisInput1,TrisInput2,Flag,Val )
%DataAppender3d Appending 2 to the back of 1. Flag says where the data in two lies.  
%   PntsInput1,PntsInput2,TrisInput1,TrisInput2 are the point and triangle
%   inputs. 1 is the data in 1. 2 is two. Flag values increase with each
%   use of this function if a flag from previous calls is brought in. If
%   not the places containing the data in two are flagged as 1's. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
%We append 2 to the back of 1. 

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
    Val=1;    
end
FlagAppend=zeros(Sz2,1);

%Finding the current max of the flag, the 2nd input has to be 1 higher. 
%n=max(Flag)+1;

%Vaues of 2nd input are now n. 
Flag=[Flag;FlagAppend+Val];


end

