function [TotalConnectionDistance,SortedTriangles,Noconnections] = ConnectedTrianglesFinder(TR,MidPoint)
% ConnectedTrianglesFinder: Finds a list of index's of connected
%                   triangles for each triangle on the surface, the
%                   distance between these and the total number of edges.
%
%                   It would be nice to optimise this as currently its a
%                   little slow. 
%               
% usage #1:
%[TotalConnectionDistance,SortedTriangles,Noconnections]...
% = ConnectedTrianglesFinder(TR,MidPoint)
%
% Arguments: (input)
%   Tr              - Output of MATLAB's triangulation function. 
%
% MidPoint          - Midpoints of the triangles of your mesh. 
%
%
% Arguments: (output)
% TotalConnectionDistance  - Each row is the distance between the triangle
%                           on this row and the other tris midpoints
%
%  SortedTriangles         - First Third and 5th column - the triangles index,
%                           other cols, the are the connected triangle
%                           index's.
%
%  Noconnections           - 1st column is the the triangles indexr, the
%                           other col the number of connected triangles
%
% Example usage (1):
%
% %If you need inputs as structured as in this BEM script use: 
%TR = triangulation(Triangles,Points(:,2:4)); 
%MidPoint= incenter(TR);
% %Call func:
%[TotalConnectionDistance,SortedTriangles,Noconnections]...
% = ConnectedTrianglesFinder(TR,MidPoint)
% % Then say you wanted to average variable StrikeSlipDisp (Dss) at every
% % tri face:
% for i=1:numel(Dss)
% Mid=Dss(i); %value at tri 1
% one=Dss(SortedTriangles(i,2)); %value at 1st connected tri
% two=Dss(SortedTriangles(i,4)); %value at 2nd connected tri
% three=Dss(SortedTriangles(i,6)); %etc
% avg(i)=mean([Mid,one,two,three]);
% end
% Dss=avg;
%[ Dss] = RowVecToCol( Dss);
%
%
% Example usage (2):
% %Meshing surface to plotting the average distance
% %between tris:
%
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[[1:numel(x)]',x(:),y(:),z(:)];
% TR = triangulation(Triangles,Points(:,2:4)); 
% MidPoint= incenter(TR);
% [TotalConnectionDistance,SortedTriangles,Noconnections]...
% = ConnectedTrianglesFinder(TR,MidPoint);
% AvgConnectionDist=(sum(TotalConnectionDistance,2))./Noconnections(:,2);
% PlotSlipDistribution3d(Triangles,Points,[],AvgConnectionDist);
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen




E = edges(TR); %Gives all connected edges. 
NoEdges=numel(E)/2;
AllTriangleNos=zeros(2,NoEdges);
%Will need loop
%Going to create a list the size of E that says which two triangles are connected to each other. 
for i=1:NoEdges
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
for i=1:length
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
    if Noconnections(i,2)==0
        continue
    end
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

