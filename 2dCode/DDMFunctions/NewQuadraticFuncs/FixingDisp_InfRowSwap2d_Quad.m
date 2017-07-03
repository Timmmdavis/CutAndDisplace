function [ DsTn,DsTs,DnTn,DnTs,Tn,Ts,NUM ]...
    = FixingDisp_InfRowSwap2d_Quad( DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy,NUM,Tn,Ts,Fdisp)
%FixingDisp_InfRowSwap This swaps rows of the coefficient matrices for
%stress to that of disp. This is be used to fix displacement of chosen mid
%points that have been flagged with 'Fdisp'. 
%Rows in the boundary condition stress vectors are switched for 0
%representing 0 displacement

%   Copyright 2017, Tim Davis, The University of Aberdeen
 %Replacing fixed disp rows with disp inf matrices
    %DsTn,DsTs,Ds_Ux,Ds_Uy,DnTn,DnTs,Dn_Ux,Dn_Uy
    FD=logical(Fdisp);
    DsTn(FD,:)=Ds_Ux(FD,:);
    DsTs(FD,:)=Ds_Uy(FD,:);
    DnTn(FD,:)=Dn_Ux(FD,:);
    DnTs(FD,:)=Dn_Uy(FD,:);
    %Now removing any fixed disp element cols
    DsTn=DsTn(:,~FD);
    DsTs=DsTs(:,~FD);
    DnTn=DnTn(:,~FD);
    DnTs=DnTs(:,~FD);
    %Now fixing the elements to 0 disp, you could put any disp here you
    %want if you wanted a displacement boundary condition
    Tn(FD)=0;
    Ts(FD)=0;
    %Removing any of the fixed elements from the size value
    Sm = sum(Fdisp);
    NUM=NUM-(Sm/3);


end

