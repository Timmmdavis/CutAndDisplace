function [TotalConnectionDistance,SortedTriangles,Noconnections] = ConnectedTrianglesFinder(TR,MidPoint)

%   Copyright 2017, Tim Davis, The University of Aberdeen
%THIS SCRIPT IS NOT OPTIMISED IN ANY SENSE OF THE WORD

%SortedTriangles First Third and 5th col - the tri number, other cols, the
%connected tris. 
%Noconnections 1st col the tri number, other col the number of connected
%tris
%TotalConnectionDistance distance between this and the other tris midpoints


%if you need inputs as structured in the BEM script use: 
%TR = triangulation(Triangles,Points(:,2:4)); 
%MidPoint= incenter(TR);

%Say you wanted to average variable StrikeSlipDisp at every tri face:
% % for i=1:numel(StrikeSlipDisp)
% % Mid=StrikeSlipDisp(i); %value at tri 1
% % one=StrikeSlipDisp(SortedTriangles(i,2)); %value at 1st connected tri
% % two=StrikeSlipDisp(SortedTriangles(i,4)); %value at 2nd connected tri
% % three=StrikeSlipDisp(SortedTriangles(i,6)); %etc
% % avg(i)=mean([Mid,one,two,three]);
% % end
% % StrikeSlipDisp=avg;
% %[ StrikeSlipDisp] = RowVecToCol( StrikeSlipDisp);




E = edges(TR); %Gives all connected edges. 
NoEdges=numel(E)/2;
AllTriangleNos=zeros(2,NoEdges);
%Will need loop
%Going to create a list the size of E that says which two triangles are connected to each other. 
for i=1:NoEdges;
    Lia = (sum((ismember(TR,E(i,:)))'))'; 
    Lia=Lia==2 ; %Creating indentity matrix where triangles have two adjacent edges at E1
    TriangleNos = find(Lia);%Finding the two connected triangles. 
    AllTriangleNos(1,i)=TriangleNos(1); 
    if numel(TriangleNos) == 2
    AllTriangleNos(2,i)=TriangleNos(2);   %Gives a matrix where each column is an edge connection with a triangle or TWO. 
    end
        
end

length=size(TR(:,1));length=length(1,1);

SortedTriangles=zeros(length,6);
SortTri=zeros(3,length); %Just associated vertex's

%To find number of connections for each triangle
for i=1:length;
    indices = find(any(AllTriangleNos==i));     
    
    Flag=AllTriangleNos(:,indices(1))==i;
    TriNo=AllTriangleNos(:,indices(1));
    TriNo(Flag,:)=[];
    SortTri(1,i)=TriNo;
    
    b=numel(indices);    
    if b==1             %Stops unconnected triangles breaking loop
        continue
    end
    Flag=AllTriangleNos(:,indices(2))==i;
    TriNo=AllTriangleNos(:,indices(2));
    TriNo(Flag,:)=[];
    SortTri(2,i)=TriNo;


    if b==2           %Stops triangles with two connections breaking loop
        continue
    end
    Flag=AllTriangleNos(:,indices(3))==i;
    TriNo=AllTriangleNos(:,indices(3));
    TriNo(Flag,:)=[];
    SortTri(3,i)=TriNo;
    %Gives list of each triangles connection in each set of two rows. 0 is
    %no connection/an edge. 
end

SortTri=SortTri';
SortTri = sort(SortTri,2,'descend'); %Sorting so 0's are always last row. Useful for loops later. 
%Ie if only 2 connections only loop through size twice not three times. 
SortedTriangles(:,2)=SortTri(:,1);SortedTriangles(:,4)=SortTri(:,2);SortedTriangles(:,6)=SortTri(:,3);
Ascend=(1:length)';
SortedTriangles(:,1)=Ascend;SortedTriangles(:,3)=Ascend;SortedTriangles(:,5)=Ascend;
%List of every triangle and connections. 

%Creating list of how many connections each tri has (can be used in loops
%later) 
Noconnections = SortTri<=0;
Noconnections=(zeros(1,length)+3)'- sum(Noconnections,2);
Noconnections=[Ascend,Noconnections];

ConnectionDist1=zeros(length,3);
ConnectionDist2=zeros(length,3);
ConnectionDist3=zeros(length,3);

for i=1:length
    ConnectionDist1=MidPoint(SortedTriangles(:,1),:)-MidPoint(SortedTriangles(:,2),:);
    if Noconnections(i,2) == 2
        ConnectionDist2(i,:)=MidPoint(SortedTriangles(i,3),:)-MidPoint(SortedTriangles(i,4),:);
    end
    if Noconnections(i,2) == 3
        ConnectionDist2(i,:)=MidPoint(SortedTriangles(i,3),:)-MidPoint(SortedTriangles(i,4),:);
        ConnectionDist3(i,:)=MidPoint(SortedTriangles(i,4),:)-MidPoint(SortedTriangles(i,5),:);
    end
end
ConnectionDistListXYZ=[ConnectionDist1;ConnectionDist2;ConnectionDist3];
TotalConnectionDistance=sqrt((ConnectionDistListXYZ(:,1).^2)+(ConnectionDistListXYZ(:,2).^2)+(ConnectionDistListXYZ(:,3).^2));
TotalConnectionDistance=reshape(TotalConnectionDistance,[],3);



end

