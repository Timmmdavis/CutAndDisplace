function [StressInf,Tn,Ts,NUM]...
    = FixingDisp_InfRowSwap2d(StressInf,DispInf,NUM,Tn,Ts,Fdisp)
% FixingDisp_InfRowSwap2d: This swaps rows of the Coefficient matrices for
%               stress to that of disp. This is be used to fix displacement
%               of chosen midpoints that have been flagged with 'Fdisp'.
%               Rows in the boundary condition stress vectors are switched
%               for 0 representing 0 displacement.
%
% usage #1:
% [DsTn,DsTs,DnTn,DnTs,Tn,Ts,NUM]...
% = FixingDisp_InfRowSwap2d(DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy,NUM,Tn,Ts,Fdisp )
%
% Arguments: (input)
% StressInf          -Structure containing:
% 					DsTn,DsTs,DnTn,DnTs
%					Square influence matricies of much a displacement
%                   of one element (first part of name) effects the traction
%                   on another element (last part of name).
%
% DispInf            -Structure containing:
% 					DsUx,DsUy,DnUx,DnUy
%					Square influence matricies of much a displacement
%                   of one element (first part of name) effects the
%                   displacement at the midpoint of another
%                   element (not the element displacement itself). 
%
%        Fdisp        - Flag telling the user if any elements are going to
%                      be fixed
%             
%       Tn,Ts         - Traction vector inputs that are the vector B of the
%                    linear equation system D = A\B.
%
%       NUM            - Size of the vectors Tn and Ts and the edge of the
%                       influence matricies
%
% Arguments: (output)
% DsTn,DsTs,DnTn,DnTs - Rectangular influence matricies of much a displacement
%                   of one element (first part of name) effects the traction
%                   on another element (last part of name) but now for
%                   elements that have fixed displacement the rows
%                   represent how much a displacement
%                   of one element (first part of name) effects the
%                   displacement at a point at the midpoint of another
%                   element (not the elements displacement itself). Cols
%                   for these effected elements have been removed as we
%                   dont want to calculate the dispalcement at these as we
%                   presume these are locked.
%
%       Tn,Ts         - Traction vector inputs that are the vector B of the
%                    linear equation system D = A\B but with rows that are
%                    0 for the locations that we have fixed displacements.
%
%       NUM         - New size of the rectangular matrix (column no)
%
% Example usage:
%
% [DsTn,DsTs,DnTn,DnTs,Tn,Ts,NUM]...
% = FixingDisp_InfRowSwap2d(DsTn,DsTs,DnTn,DnTs,Ds_Ux,Ds_Uy,Dn_Ux,Dn_Uy,NUM,Tn,Ts,Fdisp )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Replacing fixed disp rows with disp inf matrices
FD=logical(Fdisp);
StressInf.DsTn(FD,:)=DispInf.DsUx(FD,:);
StressInf.DsTs(FD,:)=DispInf.DsUy(FD,:);
StressInf.DnTn(FD,:)=DispInf.DnUx(FD,:);
StressInf.DnTs(FD,:)=DispInf.DnUy(FD,:);
%Now removing any fixed disp element cols
StressInf.DsTn=StressInf.DsTn(:,~FD);
StressInf.DsTs=StressInf.DsTs(:,~FD);
StressInf.DnTn=StressInf.DnTn(:,~FD);
StressInf.DnTs=StressInf.DnTs(:,~FD);
%Now fixing the elements to 0 disp, you could put any disp here you
%want if you wanted a displacement boundary condition
Tn(FD)=0;
Ts(FD)=0;
%Removing any of the fixed elements from the size value
Sm = sum(Fdisp);
NUM=NUM-(Sm);


end

