function cmap = colormap_cpt(cpt,varargin)
% COLORMAP_CPT   builds a MATLAB colormap from a cpt file 
%
%    colormap_cpt(name,<nSteps>)
%
% name   = the name of a colormap from the cpt-city website, see:
%          <a href="http://soliton.vm.bytemark.co.uk/pub/cpt-city/">http://soliton.vm.bytemark.co.uk/pub/cpt-city/</a>
% nSteps = number of colorsteps (optional)
%
% Note:
%  make sure to have downloaded the cpt-city package from:
%  <a href="http://soliton.vm.bytemark.co.uk/pub/cpt-city/pkg/"
%     >http://soliton.vm.bytemark.co.uk/pub/cpt-city/pkg/</a>
%  Unpack and add folder to the MATLAB search path
%
% Example:
%
% cmap = colormap_cpt('temperature');
% subplot(3,1,[1 2])
%   colormap(cmap)
%   [x,y,z] = cylinder(100:-1:10);surf(x,y,z);
%   camlight(300,20); view([-60 32]); shading flat;
%   material shiny; lighting gouraud; colorbar
% subplot(3,1,3)
%   plot(cmap); axis([1 size(cmap,1) -.01 1.01])
%
% See also: colormaps

% This tools is part of <a href="http://OpenEarth.Deltares.nl">OpenEarthTools</a>.
% OpenEarthTools is an online collaboration to share and manage data and 
% programming tools in an open source, version controlled environment.
% Sign up to recieve regular updates of this function, and to contribute 
% your own tools.

%   --------------------------------------------------------------------
%   Copyright (C) 2009 Deltares for Building with Nature
%       Thijs Damsma
%
%       Thijs.Damsma@deltares.nl	
%
%       Deltares
%       P.O. Box 177
%       2600 MH Delft
%       The Netherlands
%
%   This library is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this library.  If not, see <http://www.gnu.org/licenses/>.
%   --------------------------------------------------------------------

% $Id: colormap_cpt.m,v 1.1 2011/05/26 18:23:16 jjg Exp jjg $
% $Date: 2011/05/26 18:23:16 $
% $Author: jjg $
% $Revision: 1.1 $
% $HeadURL: https://repos.deltares.nl/repos/OpenEarthTools/trunk/MATLAB/general/color_fun/colormaps/colormap_cpt.m $
% $Keywords: $

%% adjust of *.cpt filename so you can copy paste it from the website
cpt = strrep(cpt,' ','_');
if ~strcmpi(cpt(end-3:end),'.cpt')
    cpt = [cpt '.cpt'];
end

%% open cpt file
fid = fopen(cpt);

%% skip the first bit
S = fgetl(fid);
while strcmp(S(1),'#')
    S = fgetl(fid);
end

%% read color data
cdata = sscanf(S,'%f');
S = fgetl(fid);
while S(1) ~= -1 && ~isletter(S(1)) 
    cdata = [cdata sscanf(S,'%f')];
    S = fgetl(fid);
end

%% close cpt file
fclose(fid);

%% set nSteps
if nargin == 2;
    nSteps = varargin{1};
else
    nSteps = size(cdata,2);
end

%% adjust points to allow for interpolation without duplicate points
cdata(5,1:end-1) = cdata(5,1:end-1)-abs(cdata(5,1:end-1)).*eps;
cdata(5,cdata(5,:)==0) = -eps;

%% interpolation values
xi = linspace(cdata(1,1),cdata(5,end),nSteps);

%% reshape cdata
cdata = reshape(cdata,size(cdata).*[.5 2]);

%% interp RGB values
R = interp1(cdata(1,:),cdata(2,:),xi);
G = interp1(cdata(1,:),cdata(3,:),xi);
B = interp1(cdata(1,:),cdata(4,:),xi);

%% collect data
cmap = [R' G' B']/255;