function [PaShrd,PbShrd] = FindTheSharedEdgePnt(FreeFlgPaPb,ShrdPntFlg,Pa,Pb)

%Now check where these points sit. Compare if the points XYZ is the same as
%that of the edge triangle points. We want to remove the shared edge point
%and collapse 2 to one tri.
A=ismember(Pa(FreeFlgPaPb,:),Pa(ShrdPntFlg,:),'rows');
B=ismember(Pa(FreeFlgPaPb,:),Pb(ShrdPntFlg,:),'rows');
C=ismember(Pb(FreeFlgPaPb,:),Pa(ShrdPntFlg,:),'rows');
D=ismember(Pb(FreeFlgPaPb,:),Pb(ShrdPntFlg,:),'rows');

Indx=find(FreeFlgPaPb);

PaShrd=[];
PbShrd=[];

%Then grab the actual indx
if sum(A)~=0
    PaShrd=Indx(A);
elseif sum(B)~=0
    PbShrd=Indx(B);
elseif sum(C)~=0
    PaShrd=Indx(C);
elseif sum(D)~=0
    PbShrd=Indx(D);
end

end