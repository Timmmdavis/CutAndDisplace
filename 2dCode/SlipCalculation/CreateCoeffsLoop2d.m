function [infmatrix]=CreateCoeffsLoop2d(infmatrix,...
    NUM,x,y,xe,ye,a,Beta,Ds,Dn,nu,E,halfspace,StringHS,StringFS)
%Loop that calls the DDM functions of C&S and fills large column
%coeff matrices
if halfspace==1
    % Creating a progress bar that completes during loop
    progressbar(StringHS)
    for i=1:NUM
        %Creating size and space that will be filled in the influence
        %matrix in each loop
        first = (i-1)*NUM+1;
        last = i*NUM;
        %Calling TWODD function
        infmatrix(first:last,:) = coeffhs_func(x,y,xe(i),ye(i),a(i),Beta(i),Ds,Dn,nu,E);
        % Updating progress bar
        progressbar(i/NUM) 
    end
    else 
    progressbar(StringFS)
    for i=1:NUM
        first = (i-1)*NUM+1;
        last = i*NUM;
        infmatrix(first:last,:) = coeff_func(x,y,xe(i),ye(i),a(i),Beta(i),Ds,Dn,nu,E);
        progressbar(i/NUM) 
    end
end

end