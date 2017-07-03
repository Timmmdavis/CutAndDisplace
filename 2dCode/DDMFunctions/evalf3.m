function F3 = evalf3 (A,X,Y,CON)
%   From Steve Martels BEM scripts from:
%   http://www.soest.hawaii.e.du/martel/Martel.BEM_dir/

%   Copyright 2017, Tim Davis, The University of Aberdeen

% function that evaluates arctangent terms in coeffhs
% checks for observation points that are colinear with
% actual and image element endpoints as described in
% Crouch and Starfield (pgs 50-51)
% i finds all nonzero values of Y
% j finds all points 'within' the element
% k finds all points in the plane of, but outside of each element
% Note for index j the solution is +pi

i = find(Y);
j = find( Y==0 & abs(X) < A );
k = find( Y==0 & abs(X) > A );
F3(i) = atan2(Y(i), (X(i) - A)) - atan2(Y(i), (X(i)+A));
F3(j) = -pi*ones(size(j));
F3(k) = zeros(size(k));
F3 = -CON.*(F3)';
