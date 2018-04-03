function [TriNo,P1P2FreeFlg,P2P3FreeFlg,P3P1FreeFlg]=EdgeCons(P1,P2,P3)
% EdgeCons: 
%               
% usage:
% [TriNo]=EdgeCons(Pa,Pb);
%
% Arguments: (input)
% Pa,Pb          - The corner point of each triangle in 'Triangles'.
%                    Arranged so the row index's correspond exactly to
%                    'Triangles' and 'MidPoint'. 
%
%
% Arguments: (output)
% TriNo           - The triangle index that is the connected triangle to
%                   that edge. Location of index is the triangle in question.
%
% P1P2FreeFlg     
% P2P3FreeFlg
% P3P1FreeFlg     - Flags of free edges between Pa and Pb (location
%                   specific)
%                             
% Example usage:
%
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University


[SixPntsP1P2]=CreateSortedEdgeVec(P1,P2);
[SixPntsP2P3]=CreateSortedEdgeVec(P2,P3);
[SixPntsP3P1]=CreateSortedEdgeVec(P3,P1);

%Flags for free edges:
P1P2FreeFlg=zeros(numel(P1(:,1)),1);
P2P3FreeFlg=zeros(numel(P1(:,1)),1);
P3P1FreeFlg=zeros(numel(P1(:,1)),1);

%Assuming there is max 3 connections to each tri:
TriNo=zeros(numel(P1(:,1)),15);
%Do Pa Pb connections
for j=1:numel(P1(:,1))
	for i=1:numel(P1(:,1))
     
    %find index of first 0
    [L]=find(TriNo(j,:)==0);
    if ~isempty(L)
        L=L(1);
    end
    [J]=find(TriNo(i,:)==0);
    if ~isempty(L)
        J=J(1);
    end    

   
    %Diff edge cons:
    if isequal(SixPntsP1P2(i,:),SixPntsP2P3(j,:))
        TriNo(j,L)=(i);  
        P1P2FreeFlg(i,:)=P1P2FreeFlg(i,:)+1; 
        [L]=find(TriNo(j,:)==0);if ~isempty(L);L=L(1);end 
        TriNo(i,J)=(j);  
        P2P3FreeFlg(j,:)=P2P3FreeFlg(j,:)+1; 
        [J]=find(TriNo(i,:)==0);if ~isempty(J);J=J(1);end 
    end
    
    if isequal(SixPntsP1P2(i,:),SixPntsP3P1(j,:))
        TriNo(j,L)=(i);  
        P1P2FreeFlg(i,:)=P1P2FreeFlg(i,:)+1; 
        [L]=find(TriNo(j,:)==0);if ~isempty(L);L=L(1);end 
        TriNo(i,J)=(j); 
        P3P1FreeFlg(j,:)=P3P1FreeFlg(j,:)+1; 
        [J]=find(TriNo(i,:)==0);if ~isempty(J);J=J(1);end 
    end
    
    if isequal(SixPntsP2P3(i,:),SixPntsP3P1(j,:))
        TriNo(j,L)=(i);  
        P2P3FreeFlg(i,:)=P2P3FreeFlg(i,:)+1; 
        [L]=find(TriNo(j,:)==0);if ~isempty(L);L=L(1);end 
        TriNo(i,J)=(j);  
        P3P1FreeFlg(j,:)=P3P1FreeFlg(j,:)+1; 
        [J]=find(TriNo(i,:)==0);if ~isempty(J);J=J(1);end      
    end
        
    if i==j
        continue
    end
    
    %Self cons: (e.g. only P1P2 to P1P2 edges)
    if isequal(SixPntsP1P2(i,:),SixPntsP1P2(j,:))
        TriNo(j,L)=(i);  
        P1P2FreeFlg(i,:)=P1P2FreeFlg(i,:)+1; 
        [L]=find(TriNo(j,:)==0);if ~isempty(L);L=L(1);end 
        TriNo(i,J)=(j);  
        P1P2FreeFlg(j,:)=P1P2FreeFlg(j,:)+1; 
        [J]=find(TriNo(i,:)==0);if ~isempty(J);J=J(1);end 
    end
    if isequal(SixPntsP2P3(i,:),SixPntsP2P3(j,:))
        TriNo(j,L)=(i);  
        P2P3FreeFlg(i,:)=P2P3FreeFlg(i,:)+1; 
        [L]=find(TriNo(j,:)==0);if ~isempty(L);L=L(1);end 
        TriNo(i,J)=(j);
        P2P3FreeFlg(j,:)=P2P3FreeFlg(j,:)+1; 
        [J]=find(TriNo(i,:)==0);if ~isempty(J);J=J(1);end 
    end
    if isequal(SixPntsP3P1(i,:),SixPntsP3P1(j,:))
        TriNo(j,L)=(i);  
        P3P1FreeFlg(i,:)=P3P1FreeFlg(i,:)+1; 
        [L]=find(TriNo(j,:)==0);if ~isempty(L);L=L(1);end 
        TriNo(i,J)=(j); 
        P3P1FreeFlg(j,:)=P3P1FreeFlg(j,:)+1; 
        [J]=find(TriNo(i,:)==0);if ~isempty(J);J=J(1);end     
    end

	end
end

TriNoNw=zeros(size(TriNo));
for i=1:size(TriNo,1)
    TriNoNw(i,1:numel(unique(TriNo(i,:))))=unique(TriNo(i,:));
end
%Only 3 connections (removing zeros at front)
TriNo=TriNoNw(:,2:4);

%Making flags
P1P2FreeFlg=P1P2FreeFlg==0;
P2P3FreeFlg=P2P3FreeFlg==0;
P3P1FreeFlg=P3P1FreeFlg==0;

end

function [SixPnts]=CreateSortedEdgeVec(Pa,Pb)
SixPnts=zeros(numel(Pa(:,1)),6);
for i=1:numel(Pa(:,1))
    %First sort Pa and Pb so the lowest distance from centre is always first:
    [PaAz,PaEl,PaR] = cart2sph(Pa(i,1),Pa(i,2),Pa(i,3));
    [PbAz,PbEl,PbR] = cart2sph(Pb(i,1),Pb(i,2),Pb(i,3));

    %Only in very rare cases will such a sorting be problamatic 
    [~,I]=sort([(PaAz+PaEl+PaR),(PbAz+PbEl+PbR)],2); %Sort so we have the smallest dist from body centre first
    if I(1)==1
        SixPnts(i,:)=[Pa(i,:),Pb(i,:)];
    else 
        SixPnts(i,:)=[Pb(i,:),Pa(i,:)];
    end
end

% %rounding if needed
% roundV=10;
% SixPnts=round(SixPnts,roundV);

end