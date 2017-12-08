%% Copyright (C) 2007-2012 David Bateman
%%
%% This file is part of Octave.
%%
%% Octave is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or (at
%% your option) any later version.
%%
%% Octave is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with Octave; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.


%% -*- texinfo -*-
%% @deftypefn  {Function File} {} trisurf (@var{tri}, @var{x}, @var{y}, @var{z})
%% @deftypefnx {Function File} {@var{h} =} trisurf (@dots{})
%% Plot a triangular surface in 3D@.  The variable @var{tri} is the triangular
%% meshing of the points @code{(@var{x}, @var{y})} which is returned
%% from @code{delaunay}.  The variable @var{z} is value at the point
%% @code{(@var{x}, @var{y})}.
%%
%% The optional return value @var{h} is a graphics handle to the created plot.
%% @seealso{triplot, trimesh, delaunay3}
%% @end deftypefn

%   Copyright 2017, Tim Davis, The University of Aberdeen
%   Used here as an override to the trisurf func,
%   Octave fails to plot my triangulated surfaces correctly

function h = trisurf (tri, x, y, z, varargin)

  if (is_octave)
    % catch Octave transparency and ignore 
    aa=cellstr(varargin);
    C=strfind(aa,'FaceAlpha');
    C = ~cellfun(@isempty,C); %Replacing [] for 0's, should become logical
    if any(C) %If any are true it exists
      varargin(:, [C]) = []; %Remove the alpha line
      varargin(:, [C]) = []; %Remove the alpha value
    else
    end  
  else
    % do it MATLAB way
  end



 
ax = gca;% axescheck(varargin{:});
ax = newplot(ax);
start = 1;

  if (nargin < 3)
    print_usage ();
  end

  if (nargin == 3)
    triplot (tri, x, y);
  elseif (ischar (z))
    triplot (tri, x, y, z, varargin{:});
  else
    if (nargin > 4 && isnumeric (varargin{1}))
      c = varargin{1};
      varargin(1) = [];
    else
      c = z;
    end
    if ( any (strcmpi (varargin, 'EdgeColor'))...
        && strcmpi (varargin{nfc+1}, 'interp'))
      varargin(end+(1:2)) = {'EdgeColor', 'none'};
    end
    newplot ();
    handle = patch ('Faces', tri, 'Vertices', [x(:), y(:), z(:)],...
                    'facevertexcdata', c(:),...
                    'facecolor',get(ax,'DefaultSurfaceFaceColor'), ...
                    'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),'parent',ax,...
                    varargin{start:end});
    if (nargout > 0)is_octave
      h = handle;
    end

    if ~ishold(ax), view(ax,3), grid(ax,'on'), end
  end

end

function r = is_octave ()
  persistent x;
  if (isempty (x))
    x = exist ('OCTAVE_VERSION', 'builtin');
  end
  r = x;
end
