%How to call H matrices approx in code
%When down to the linear EQ -     D = A\B; in the func 'SlipCalculator'

%   Copyright 2017, Tim Davis, The University of Aberdeen

tic;
%if square
[EPS,MaxRank,minn,n1,n2,n3,n4] = InitDefaultHmatVars( A );

%%%Short way
AA = HMatrix(A, [n1,n2],[n3,n4], 'S', [0,0], [0,0],  EPS, MaxRank, minn);
%Running to find unknowns 'DD'
DD = AA\B; %mldivide(AA,B);
%Or like so with LU factorisation (does this anyway)
%[L,U] = lu( AA );
%DD_lu= U\(L\B); %DD_lu= (L*(U*B)); %%could do this if we want inverse of A

%%%Long way:
%%if size(A) == size(A')
%%%Flipping A as we can only do matrix multiplication
%%A=inv(A); %This is actually C now
%%else
%%%or if not square
%%A=pinv(A); %This is actually C now
%%end
%%%Making Hmat coeff matrix
%%A = HMatrix(A, [n3,n4],[n1,n2], 'S', [0,0], [0,0],  EPS, MaxRank, minn);
%%%solving for displacement. 
%%C = mtimes( A, B );

toc

%having tested the above for simple 2D approximation it works ok. Would be good to evaluate each parameter.

%Need to work out a way of populating AH during building of coeff matrix. 
%its not too much use if I do this once the matrix is fully sized.

%If you want to try others - https://github.com/gchavez2/awesome-hierarchical-matrices
