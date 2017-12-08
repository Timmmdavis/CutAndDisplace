function [External] = InsideEllipsoid(x, y, z, rx, ry, rz, tol_dist)
%   InsideEllipsoid.m 
%
%   Is this cartesian point outside the ellipsoidal inclusion?
%   
%   Arguments:
%       x, y, z = coordinates 
%       a, b, c = semi-axes of ellipsoid
%
%   David Healy
%   October 2008 
%
%   Please let me know about any bugs or errors: 
%       d.healy@curtin.edu.au
%       david.healy@mac.com 
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% Slightly updated by Tim Davis (Aberdeen/Potsdam University). 


%Getting original size of X and Y (so we can return these to vectors
%afterwards).
Xorig=size(x);
%Turning into col vecs
x=x(:);y=y(:);z=z(:);
%Size of X. 
sz=numel(x);
%Our flag, assuming all are inside to start with.
External = zeros(sz,1) ; 

%Doing for every point.
for i=1:sz
    if ( ( ( x(i)^2 / rx^2 ) + ( y(i)^2 / ry^2 ) + ( z(i)^2 / rz^2 ) ) > 1.0+tol_dist )
        External(i) = true ; 
    end  
end

%Reshaping back to the original size. 
External=reshape(External,Xorig);
