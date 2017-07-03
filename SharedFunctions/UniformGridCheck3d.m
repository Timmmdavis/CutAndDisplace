function [ Uniform ] = UniformGridCheck3d( X,Y,Z )
%UniformGridCheck2d Checks X and Y, to see if these have a uniform linear
%spacing. Can be used to check before plotting with different functions or
%for interpolation. 

%Can be used to check before plotting with different functions or
%for using in interpolation. Returns 1 if true. Can use in if statements 
%i.e.. if uniform spaced points plot with isocontour else use scatter3. 

%   Copyright 2017, Tim Davis, The University of Aberdeen


if iscolumn(X)
    flagx=1;
    flagy=1;
    flagz=1;
else    
    %checking that grids are uniformly spaced
    flagx=zeros(1,numel(X(:,1))-1);
    for i=1:numel(X(:,1))-1
        FirstSpace=round(X(1,1)-X(1,2),8);
        spacing=round(X(:,i)-X(:,i+1),8);
        LogEq=spacing==FirstSpace;
        if any(~LogEq)
            flagx(i)=1;
        else
            flagx(i)=0;
        end
    end    
    %Grid data check Y
    flagy=zeros(1,numel(Y(:,1))-1);
    for i=1:numel(Y(:,1))-1
        FirstSpace=round(Y(1,1)-Y(1,2),8);
        spacing=round(Y(:,i)-Y(:,i+1),8);
        LogEq=spacing==FirstSpace;
        if any(~LogEq)
            flagy(i)=1;
        else
            flagy(i)=0;
        end
    end    
        %Grid data check Z
    flagz=zeros(1,numel(Z(:,1))-1);
    for i=1:numel(Z(:,1))-1
        FirstSpace=round(Z(1,1)-Z(1,2),8);
        spacing=round(Z(:,i)-Z(:,i+1),8);
        LogEq=spacing==FirstSpace;
        if any(~LogEq)
            flagz(i)=1;
        else
            flagz(i)=0;
        end
    end    
end


Uniform= any(~flagx) || any(~flagy) || any(~flagz); %if the data is a uniform grid this==1


end

