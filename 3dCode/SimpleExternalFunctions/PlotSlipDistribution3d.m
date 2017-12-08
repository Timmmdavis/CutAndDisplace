function PlotSlipDistribution3d(Triangles,Points,cmap,varargin )
% PlotSlipDistribution3d: Plots multiple figures of meshes colourmapped by 
%                   the arguments in varargin.
%               
% usage #1:
% PlotSlipDistribution3d(Triangles,Points,cmap,varargin )
%
% Arguments: (input)
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                     locations in "Points" which contains the XYZ location
%                     of each corner of the triangle.
%
% cmap              -  A colourmap that MATLAB can use. See func
%                   "colormap_cpt.m" to produce one. 
%
% varargin          -  Vectors of values the same length as input
%                     "Triangles".
%
% Arguments: (additional inputs)
%
% ,'Flag',Val,      - Specifying 'Flag' and then a value within the varargin
%                    inputs puts 0 at the centre of the colourmap.
%
% Example usage:
% 
% PlotSlipDistribution3d(Triangles,Points,cmap,Dn,Dss,Dds,'Flag',1 )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Finding if additional argument 'Flag' str is in the input args, if it is
%we use the value after as the flag. 
Default=0;
[ varargin,Flag ] = AdditionalArgsInVaragin( 'Flag',varargin,Default );

%Starting drawing
for i=1:numel(varargin)
    %Grabbing the input variables name, this is a string
    VarName=inputname(i+3); 
    %Get the argument for this loop.
    Data=varargin{i};    
    %Making sure its a col.
    [ Data ] = RowVecToCol( Data ); 
    %Draw figure
	figure;
	trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Data);
    %Draw cmap with the imported cmap if its around. 
    if isempty(cmap); colormap('default'); else; colormap(cmap); end 
    %Colourbar and X Y labels. 
    colorbar; 
    xlabel('x'); ylabel('y'); axis('equal');

    %Put 0 at the centre of the Cmap? 1=yes
	if Flag==1
        DivergingCentre( Data )
	end	
	
	%Setting titles
    if isequal(VarName,'Dss')
        title({'\fontsize{14}SsDisp','\fontsize{8}Positive is left lateral movement'})
	elseif isequal(VarName,'Dds')
        title({'\fontsize{14}DsDisp','\fontsize{8}Positive is dipslip faulting if surface normals point up'})
	elseif isequal(VarName,'Dn') 
        title({'\fontsize{14}TsDisp','\fontsize{8}Positive is an opening across the crack walls'})
    else
        %Just use the input variables name. 
        title(VarName);
    end
end

end
