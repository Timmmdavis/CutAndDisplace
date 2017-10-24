function [infmatrix,DisplacementXYZ]=CreateCoeffsLoop3d(infmatrix,DisplacementXYZ,...
    NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD,StringHS,StringFS)
%Loop that calls the TDE functions of Mehdi Nikkhoo and fills large column
%coeff matrices
NUM2=size(P1,1);

if halfspace==1
	progressbar(StringHS) % Create figure and set starting time
	for i=1:size(P1,1)
		first = (i-1)*NUM+1;
		last = i*NUM;
		infmatrix(first:last,:) = TDstrain_stressHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,mu,lambda);
        if FD==1
        DisplacementXYZ(first:last,:) = TDdispHS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,nu);
        end
		progressbar(i/NUM2) % Update figure 
	end
else 
	progressbar(StringFS) % Create figure and set starting time
	for i=1:size(P1,1)
		first = (i-1)*NUM+1;
		last = i*NUM;
		infmatrix(first:last,:) = TDstrain_stressFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,mu,lambda);
        if FD==1
		DisplacementXYZ(first:last,:) = TDdispFS(X,Y,Z,P1(i,:),P2(i,:),P3(i,:),Dss,Dds,Dn,nu); 
        end
		progressbar(i/NUM2) % Update figure 
	end
end
end