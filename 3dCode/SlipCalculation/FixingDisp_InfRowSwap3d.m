function [ StressInf,NUM,Tn,Tss,Tds ]...
    = FixingDisp_InfRowSwap3d( StressInf,DispInf,NUM,Tn,Tss,Tds,Fdisp)
% FixingDisp_InfRowSwap3d: This swaps rows of the Coefficient matrices for
%                       stress to that of disp. This is be used to fix displacement
%                       of chosen midpoints that have been flagged with 'Fdisp'.
%                       Rows in the boundary condition stress vectors are switched
%                       for 0 representing 0 displacement.
%
% usage #1:
%  [ StressInf,NUM,Tn,Tss,Tds ]...
%    = FixingDisp_InfRowSwap3d( StressInf,DispInf,NUM,Tn,Tss,Tds,Fdisp)
%
% Arguments: (input)
% StressInf          -Structure containing:
%                     DnTn,DnTss,DnTds
%                     DssTn,DssTss,DssTds
%                     DdsTn,DdsTss,DdsTds
%                     Square influence matricies of much a displacement
%                     of one element (first part of name) effects the traction
%                     on another element (last part of name).
%
% DispInf            -Structure containing:
%                     Dn_dx,Dn_dy,Dn_dz
%                     Dss_dx,Dss_dy,Dss_dz
%                     Dds_dx,Dds_dy,Dds_dz 
%                     Square influence matricies of much a displacement
%                     of one element (first part of name) effects the
%                     displacement at the midpoint of another
%                     element (not the elements displacement itself). 
%
%
%        Fdisp        - Flag telling the user if any elements are going to
%                       be fixed.
%             
%  Tn,Tss,Tds         - Traction vector inputs that are the vector B of the
%                       linear equation system D = A\B.
%
%       NUM            - Size of the vectors Tn etc and the edge of the
%                        influence matricies.
%
% Arguments: (output)
% StressInf          -Structure containing:
%                     DnTn,DnTss,DnTds
%                     DssTn,DssTss,DssTds
%                     DdsTn,DdsTss,DdsTds
%                     Square influence matricies of much a displacement
%                     of one element (first part of name) effects the traction
%                     on another element (last part of name).
%
%
%  Tn,Tss,Tds         -  Traction vector inputs that are the vector B of the
%                    linear equation system D = A\B but with rows that are
%                    0 for the locations that we have fixed displacements.
%
%       NUM         - New size of the rectangular matrix (column no)
%
% Example usage:
%
%  [ StressInf,NUM,Tn,Tss,Tds ]...
%     = FixingDisp_InfRowSwap3d( StressInf,DispInf,NUM,Tn,Tss,Tds,Fdisp)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
    

%Replacing fixed disp rows with disp inf matrices
FD=logical(Fdisp);
StressInf.DnTn(FD,:)  =DispInf.DnUx(FD,:);
StressInf.DnTss(FD,:) =DispInf.DnUy(FD,:);
StressInf.DnTds(FD,:) =DispInf.DnUz(FD,:);
StressInf.DssTn(FD,:) =DispInf.DssUx(FD,:);
StressInf.DssTss(FD,:)=DispInf.DssUy(FD,:);
StressInf.DssTds(FD,:)=DispInf.DssUz(FD,:);
StressInf.DdsTn(FD,:) =DispInf.DdsUx(FD,:);
StressInf.DdsTss(FD,:)=DispInf.DdsUy(FD,:);
StressInf.DdsTds(FD,:)=DispInf.DdsUz(FD,:);
%Now removing any fixed disp element cols
StressInf.DnTn  =StressInf.DnTn(:,~FD);
StressInf.DnTss =StressInf.DnTss(:,~FD);
StressInf.DnTds =StressInf.DnTds(:,~FD);
StressInf.DssTn =StressInf.DssTn(:,~FD);
StressInf.DssTss=StressInf.DssTss(:,~FD);
StressInf.DssTds=StressInf.DssTds(:,~FD);
StressInf.DdsTn =StressInf.DdsTn(:,~FD);
StressInf.DdsTss=StressInf.DdsTss(:,~FD);
StressInf.DdsTds=StressInf.DdsTds(:,~FD);
%Now fixing the triangles to 0 disp, Replacing tractions with the disp
%boundary condition you could put any disp here you want
Tn(FD)=0;
Tss(FD)=0;
Tds(FD)=0;
%Removing any of the fixed triangles from the size
Sm = sum(Fdisp);
NUM=NUM-Sm;


end

