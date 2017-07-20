function PlotSlipDistribution3d(Triangles,Points,cmap,varargin )
% Plots data on mesh surface but can be used with multiple arguments for multiple figs

%Example usage
%PlotSlipDistribution3d(Triangles,Points,cmap2,StrikeSlipDisp,DipSlipDisp,TensileSlipDisp)

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Put 0 at the centre of the Cmap? 1=yes
flag=0;

for i=1:nargin-3
    VarName=inputname(i+3); %grabbing the input variables name, is a string
    Data=varargin{i};    [ Data ] = RowVecToCol( Data ); %making sure its a col
	figure;
	trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),Data);
    colormap(cmap);	
    colorbar; %caxis([min(Data) max(Data)]);
    xlabel('x'); ylabel('y'); axis('equal');

	if flag==1
	DivergingCentre( Data )
	end	
	
	%setting title
	if isequal(VarName,'StrikeSlipDisp')
	title({'\fontsize{14}SSDisp','\fontsize{8}Positive is left lateral movement'})
	elseif isequal(VarName,'DipSlipDisp')
	title({'\fontsize{14}DSDisp','\fontsize{8}Positive is dipslip faulting if surface normals point up'})
	elseif isequal(VarName,'TensileSlipDisp') 
	title({'\fontsize{14}TSDisp','\fontsize{8}Positive is an opening across the crack walls'})
	else
    title(VarName);
	end

end

end
