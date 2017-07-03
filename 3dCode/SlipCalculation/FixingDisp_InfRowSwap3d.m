function [ DnTn,DnTss,DnTds,DssTn,DssTss,DssTds,DdsTn,DdsTss,DdsTds,NUM,Tnn,Tss,Tds ]...
    = FixingDisp_InfRowSwap3d( DnTn,DnTss,DnTds,DssTn,DssTss,DssTds,DdsTn,...
    DdsTss,DdsTds,Dn_dx,Dn_dy,Dn_dz,Dss_dx,Dss_dy,Dss_dz,Dds_dx,Dds_dy,Dds_dz,NUM,Tnn,Tss,Tds,Fdisp)
%FixingDisp_InfRowSwap This swaps rows of the Coefficient matrices for
%stress to that of disp. This is be used to fix displacement of chosen mid
%points that have been flagged with 'Fdisp'. 
%Rows in the boundary condition stress vectors are switched for 0
%representing 0 displacement

%   Copyright 2017, Tim Davis, The University of Aberdeen
    %Replacing fixed disp rows with disp inf matrices
    FD=logical(Fdisp);
    DnTn(FD,:)=Dn_dx(FD,:);
    DnTss(FD,:)=Dn_dy(FD,:);
    DnTds(FD,:)=Dn_dz(FD,:);
    DssTn(FD,:)=Dss_dx(FD,:);
    DssTss(FD,:)=Dss_dy(FD,:);
    DssTds(FD,:)=Dss_dz(FD,:);
    DdsTn(FD,:)=Dds_dx(FD,:);
    DdsTss(FD,:)=Dds_dy(FD,:);
    DdsTds(FD,:)=Dds_dz(FD,:);
    %Now removing any fixed disp element cols
    DnTn=DnTn(:,~FD);
    DnTss=DnTss(:,~FD);
    DnTds=DnTds(:,~FD);
    DssTn=DssTn(:,~FD);
    DssTss=DssTss(:,~FD);
    DssTds=DssTds(:,~FD);
    DdsTn=DdsTn(:,~FD);
    DdsTss=DdsTss(:,~FD);
    DdsTds=DdsTds(:,~FD);
    %Now fixing the triangles to 0 disp, Replacing tractions with the disp
    %boundary condition you could put any disp here you want
    Tnn(FD)=0;
    Tss(FD)=0;
    Tds(FD)=0;
    %Removing any of the fixed triangles from the size
    Sm = sum(Fdisp);
    NUM=NUM-Sm;


end

