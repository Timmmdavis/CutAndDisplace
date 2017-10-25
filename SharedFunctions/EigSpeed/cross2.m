function c = cross2(a, b, idA, idB)
%CROSS2  Vector cross product.
%   This function exploits the BAXFUN engine (MATLAB Central, file #23084),
%   which enables multiple products and array expansion (AX).
%
%   Copyright (c) 2009, Paolo de Leva
%   All rights reserved.
%
%   When A and B are 3-element vectors (e.g. 3×1, 1×3, or 1×1×3 arrays):
%
%       C = CROSS2(A, B) returns their cross product. That is, C = A x B.
%       If NDIMS(A) == 2 and NDIMS(B) == 2, CROSS2(A, B) is equivalent to
%       CROSS(A, B).
%
%   More generally, when A and B are arrays of any size containing one or
%   more 3-element vectors:
%
%       C = CROSS2(A, B) is equivalent to C = CROSS2(A, B, IDA, IDB), where
%       IDA and IDB are the first dimensions of A and B whose length is 3.
%
%       C = CROSS2(A, B, DIM) is equivalent to C = CROSS2(A, B, IDA, IDB),
%       where IDA = IDB = DIM. If A and B have the same size, it is also
%       equivalent to C = CROSS(A, B, DIM).
%
%       C = CROSS2(A, B, IDA, IDB) returns the cross products between the
%       vectors contained in A along dimension IDA and those contained in B
%       along dimension IDB. These vectors must have 3 elements. A and B
%       are viewed as "block arrays". IDA and IDB are referred to as their
%       "internal dimensions" (IDs). For instance, a 3×6×2 array may be
%       viewed as an array containing twelve 3-element blocks. In this
%       case, its size is denoted by (3)×6×2, and its ID is 1. Since AX is
%       enabled, A and B may have different size, and IDA may not coincide
%       with IDB (see MULTIPROD).
%
%       Input and output format:
%           Array     Block size     Internal dimension
%           ---------------------------------------------------------------
%           A         3  (1-D)       IDA
%           B         3  (1-D)       IDB
%           C         3  (1-D)       MAX(IDA, IDB)
%           ---------------------------------------------------------------
%           If SIZE(A)==SIZE(B) and IDA==IDB, then SIZE(C)=SIZE(A)=SIZE(B).
%
%   Examples:
%       If A and B are both   (3)×6×2 arrays of vectors,
%       C = CROSS2(A, B) is a (3)×6×2 array  of vectors.
%
%       A single vector B multiplies thirty vectors contained in A: 
%       If  A is ................ a 5×6×(3) array of 30 vectors,
%       and B is ................ a (3)×1   vector,
%       C = CROSS2(A, B, 3, 1) is a 5×6×(3) array of 30 vectors.
%
%   See also DOT2, CROSS, CROSSDIV, OUTER, MAGN, UNIT, PROJECTION,
%            REJECTION, TESTCROSS2.

% $ Version: 2.0 $
% CODE      by:                 Paolo de Leva (IUSM, Rome, IT) 2009 Feb 2
%           optimized by:       Code author                    2009 Feb 12
% COMMENTS  by:                 Code author                    2009 Feb 26
% OUTPUT    tested by:          Code author                    2009 Feb 26
% -------------------------------------------------------------------------

% Allow 2 to 4 input arguments
 narginchk(2, 4) ; 

% Setting IDA and/or IDB
sizeA = size(a);
sizeB = size(b);
switch nargin
    case 2        
        idA = find(sizeA==3, 1, 'first'); % First dim. of length 3
        idB = find(sizeB==3, 1, 'first');
        if isempty(idA) || isempty(idB)
            error('CROSS2:InvalidSize',...
                  'A and B must have at least one dimension of length 3.');
        end        
    case 3
        idB = idA;
end

% Check block size
if (sizeA(idA)~=3) || (sizeB(idB)~=3),
    error('CROSS2:InvalidBlockSize',...
         ['A and B must be of length 3 in the dimensions\n'...
          'in which the cross product is taken.'])
end

% Dimension shift (first step of array expansion)
%     NOTE: The BAXFUN engine is implemented here without calling
%           BAXFUN, to avoid repeated shift of A and B
diff = idB - idA;
if diff >= 0
    id = idB;
    sizeA = [ones(1, diff) sizeA];
    a = reshape(a, sizeA);
else
    id = idA;
    sizeB = [ones(1,-diff) sizeB];
    b = reshape(b, sizeB);
end

% Vectorized indices to rearrange A and B
idxA = ivector(sizeA);
idxB = ivector(sizeB);

% Calculate cross product
%     c = [a(2,:).*b(3,:) - a(3,:).*b(2,:)
%          a(3,:).*b(1,:) - a(1,:).*b(3,:)
%          a(1,:).*b(2,:) - a(2,:).*b(1,:)];
idxA{id} = [2 3 1]; 
idxB{id} = [3 1 2];
c =     bsxfun( @times, a(idxA{:}), b(idxB{:}) );
idxA{id} = [3 1 2]; 
idxB{id} = [2 3 1];
c = c - bsxfun( @times, a(idxA{:}), b(idxB{:}) );


%--------------------------------------------------------------------------
function indices = ivector(sizeA)
%IVECTOR   Vectorizing the indices of an array. 

Ndims = length(sizeA);
indices = cell(1,Ndims); % preallocating
for d = 1 : Ndims
   indices{d} = 1:sizeA(d);
end    