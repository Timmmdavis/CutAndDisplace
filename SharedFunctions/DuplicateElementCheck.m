function DuplicateElementCheck(MidPoint,Fdisp)
% MidPointCreate: Check to see if there are duplicate points on the line
%                   and calls error if so. By duplicate points the
%                   midpoints of the segments are in the same place.
%               
% usage #1:
% DuplicateElementCheck(MidPoint,Fdisp)
%
% Arguments: (input)
% MidPoint          - The midpoints n*3, (XYZ) of each triangle. 
%
% Fdisp             - Typically a flag that says which elements are locked.
%                     For the inhomogenous problem it is different:
%                       0=free boundary E1
%                       1=fixed bits of the free boundary of E1
%                       2=E1-E2 interface, E1 elastic properties, normals point towards E2
%                       3=E2-E1 interface, E2 elastic properties, normals point towards E1
%                       4=If existed would be free boundary E2 
%                       5=fixed bits of the free boundary of E2
%                      Interfaces have duplicate elements. 
%
%
%
% Arguments: (output)
% N/A
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



%Returns midpoint with only unique rows (rounding slightly).
[C] = unique(round(MidPoint,8),'rows'); 
%Number of interface els (duplicate points by definition): 
IFEls=sum(Fdisp == 3| Fdisp == 2)/2;
%Check if the number of unique rows matches that of MidPoint
if size(C,1)+IFEls~=size(MidPoint,1) 
    disp('Number of duplicate elements:')
    disp((size(C,1)+IFEls)-size(MidPoint,1))
    error('Duplicate elements exist, linear equations will fail')
end



