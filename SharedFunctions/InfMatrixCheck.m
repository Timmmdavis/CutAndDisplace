function InfMatrixCheck( InfluenceMatrix )
% InfMatrixCheck: Checks the influence matrix for nans and rank
%           number. Nans mean elements are probably overlapping so a
%           midpoint of another element lies on the plane of one of the
%           elements. 
%           The Rank number gives an idea of the number of
%           additional constraints required. For example in 2D for a single
%           closed body problem this will be '2' if there are no additional
%           constraints already added. This suggests the problem needs to
%           be constrained in some manner, the simple solution suggested by
%           Crouch and Starfield is to add displacement constraints in X
%           and Y in the part of the body you are not interested in, for
%           example the interior of the closed body if you are intersted in
%           evalulating the exterior. 
%
% usage #1:
% [ InfluenceMatrix ] = InfMatrixCheck( InfluenceMatrix )
%
% Arguments: (input)
% InfluenceMatrix   - The dense influence matrix for the boundary element
%                    problem that we check for rank and nans.
%
%
% Arguments: (output)
% N/A
%
% Example usage (1):
%
% %Check influence matrix 'A' is OK.
% InfMatrixCheck( A )
% %No errors reported? Continue with linear eq system:
% D = A\B;
%
% Example usage (1):
%
% %Show this reports nans/inf values
% A=rand(10);
% A(2,4)=inf;
% InfMatrixCheck( A )
% % You should get an error. 
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Check for nans/inf values
bad = isnan(InfluenceMatrix) | isinf(InfluenceMatrix);
if any(bad(:))
    warning('NaN or inf values exist in your influence matrix, are there overlapping segments?')
    return %Leave function without computing rank, doing this with nans throws an error. 
end

%Number of Rows
Rows=size(InfluenceMatrix,2);
%Calculate Rank with MATLAB's func
Rank = rank(InfluenceMatrix);

if Rank<Rows 

    %Number of additional constraints required to make your matrix non
    %singular.
    warning('Number of additional constraints required:')
    disp(Rows-Rank)
    
end    


end

