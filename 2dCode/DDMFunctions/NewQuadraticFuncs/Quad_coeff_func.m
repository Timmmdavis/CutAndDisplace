function [StressDisp] = ...
Quad_coeff_func(x,y,xe,ye,a,Beta,Dx,Dy,pr,G)

%From 'A time domain boundary element method for modeling the quasi-static
%viscoelastic behaviour of asphalt pavements' 
%Ref in the code as "Wang"

%Also can use the ref 'Viscoelastic DDM for Analysis of Pavements' 
%Ref in the code as "Birgisson". This paper has an error in the derivatives
%effecting the displacements. Using the paper above is a safer bet. 

%   Copyright 2017, Tim Davis, Aberdeen University

%Supplies stress and displacement due to a planar discontinuity lying along the x
%axis with length 'a'. This has quadratic shear/normal dispalcement. 

%x - x location of obs points in relation to the frac
%y - y location of obs points in relation to the frac
%xe - centre point of the fracture in cart coords 
%ye - centre point of the fracture in cart coords. 
%a - half length of dislocation, see fig 3 of "Birgisson"
%Beta - angle of element away from X axis. 
%Dy (Nd) - normal disp, 1x3 mat
%Dx (Sd) - Shear disp, 1x3 mat
%Poisson's ratio - Pr
%Shear mod - G



%Things to do with this:
%1. Wrap with coordinate transforms, elements of any orientation *done
%2. Building coeff matrices. Similar to constant, mag of 1 at each loc
%3. Add tip elements (see "Birgisson") 
%4. Add curved displacement between touching angled elements (See "Shou & Crouch 1995" ref in Birgisson)
%5. Drawing of the quadratic displaced crack surface after eq's are solved
%6. Check accuracy. 
%7. Speed up all of the above (this func could be optimised too)
%8. Need to flag points at end of dislocation that get nan values. Ones on
%the dislocation are fine

%Rotation back to flat, Beta is angle away
Beta=-Beta;

%Collation displacement. These points lie in the locations on the
%elements described in fig 3 & eq [8] of "Birgisson"
%Displacement Convention as in : fig 2 of "Birgisson"
%Sd and Nd must be 1x3 vecs: ie.=[0,0,0]

% %Dy at the three collation points along the element (Ds)
% Dy=Nd; %Negative is opening, Twodd conv 
% %Dx at the three collation points (Dn)
% Dx=Sd;  %Negative is right lat, Twodd conv

%local coords so element is at the centre
x=x(:)-xe;
y=y(:)-ye;

%Rotating so element is flat in respect to obs points 
[x,y] = RotateObject2d(x,y,Beta);


%For func that is called
rsq3=1/sqrt(3); %reciprocal of sqrt 3
two3rds=2/3;
fr3rds=2*two3rds;

%Eqs [15] of "Birgisson". Elastic constants pr(Poisson's ratio) and G(shr
%mod)
k1=1/(1-pr);
k2=(1-(2*pr))/(1-pr);
k3=(2*G)/(1-pr);
kDisp=(0.5-(0.5*pr)).*(G); %Scales displacements correctly (Added by Tim Davis)

%Scaling the DD values so the output stress and disps are then correct
%relative to these.
%Found these values by comparing results to P&S analytical solution for a
%line crack. 
Dy=Dy./((pi*(2-(2*pr))))/G;
Dx=Dx./(37.7/3);


%Some constants we introduce to speed up the maths in appendix [B] eqs of "Birgisson".
xaa=x+a; %'a' represents add
xma=x-a; %'m' represents minus
a2=a.^2;  %'2' represents squared
x2=x.^2;  %'2' represents squared
y2=y.^2;  %'2' represents squared

%%Eqs [A5] of "Birgisson" 
%r1 - Distance of point from left hand tip
r1=sqrt((xma.^2)+y2);
%r2 - Distance of point from right hand tip
r2=sqrt((xaa.^2)+y2);
%thet1 - Angle of point away from left hand tip (measured from x) 
thet1=atan2(y,xma); %atan(y/(x-a));
%thet2 - Angle of point away from right hand tip (measured from x) 
thet2=atan2(y,xaa);  %atan(y/(x+a));

%Some additional constants we introduce to speed up the maths in appendix [B] eqs of "Birgisson".
r12=r1.^2;
r22=r2.^2;
r14=r1.^4;
r24=r2.^4;
thet1m2=thet1-thet2;
thet2m1=thet2-thet1;
rr22=1./r22; %'r' at front represents reciprocal 
rr12=1./r12; 
lnr2dr1=log(r2./r1); %'ln' represents ln of values, 'd' for divided
lnr1dr2=log(r1./r2); 


%A constant used in the front of Eqs [B1]
%%Eqs [A1] of "Wang" 
I0  =	lnr2dr1;
%%And IZer's derivatives
I0X =	((xaa./r22)-(xma./r12));
I0Y =	y.*(rr22-rr12);
I0YY=	(((r22-(2.*y2))./r24)-(((r12)-(2.*y2))./r14));

%Second constant introduced for front of [B2]
%%Eqs [A2] of "Wang". IOne and its derivatives
I1   =	((y.*thet1m2)./a)+(x./a).*lnr2dr1-2;
I1X  =	(1./a).*lnr2dr1-(xma./r12)-(xaa./r22);
I1Y  =	(thet1m2./a)-y.*(rr12+rr22);
I1YY =	(x./a).*(rr12-rr22)-(2.*((xma./r12).^2))-(2.*((xaa./r22).^2));

%Third constant introduced for front of [B3]
%%Eqs [A3] of "Wang" 
I2   =	(((2.*x.*y).*thet1m2)./a2)+((y2-x2)./a2).*lnr1dr2-((2.*x)./a);
I2X  =	(((2.*y).*thet1m2)./a2)+((2.*x)./a2).*lnr2dr1+((xaa./r22)-(xma./r12))-(4./a);
I2Y  =	(((2.*x).*thet1m2)./a2)+((2*y)./a2).*lnr1dr2+y.*(rr22-rr12);
I2YY =	(2./a2).*lnr1dr2+(1./a).*((((2.*x)-a)./r12)+(((2.*x)+a)./r22))+2.*((xaa./r22).^2)-2.*((xma./r12).^2);

%%Eqs [A4] of "Wang" 
J0=thet1m2;
J0Y=-I0X;
J0YY=(2.*y).*((xaa./r24)-(xma./r14));

%%Eqs [A5] of "Wang" 
J1=((x.*thet1m2)./a)+(y./a).*lnr1dr2;
J1Y=-I1X;
J1YY=(y./a).*(rr12-rr22)-(2.*y).*((xma./r14)+(xaa./r24));

%%Eqs [A6] of "Wang" 
J2=(((x2-y2).*thet1m2)./a2)+((2.*x.*y)./a2).*lnr1dr2+((2.*y)./a);
J2Y=-I2X;
J2YY=((2.*thet2m1)./a2)+((2.*y)./a).*(rr12+rr22)-(2.*y).*((xma./r14)-(xaa./r24));


%%Equation [12] of "Birgisson", creating the vectors.
I=IJVec(I0,I1,I2,rsq3,two3rds,fr3rds);
Ix=IJVec(I0X,I1X,I2X,rsq3,two3rds,fr3rds);
Iy=IJVec(I0Y,I1Y,I2Y,rsq3,two3rds,fr3rds);
Iyy=IJVec(I0YY,I1YY,I2YY,rsq3,two3rds,fr3rds);
J=IJVec(J0,J1,J2,rsq3,two3rds,fr3rds);
Jy=IJVec(J0Y,J1Y,J2Y,rsq3,two3rds,fr3rds);
Jyy=IJVec(J0YY,J1YY,J2YY,rsq3,two3rds,fr3rds);


Dx=repmat(Dx,numel(x),1);
Dy=repmat(Dy,numel(x),1);

%%For Uy, Eq [13] of "Birgisson" 
%Ux=(sum((-2*J+k1*y*Ix).*Dx+(k2*I+k1*y*Iy).*Dy))*kDisp;  %This value "kDisp" makes it match the an sol 
%Uy=(sum((-k2*I+k1*y*Iy).*Dx+(-2*J+k1*y*Jy).*Dy))*kDisp;
Ux1 = ((bsxfun(@plus,bsxfun(@times,y*k1,Ix),(-2*J))).*Dx);
Ux2 =  bsxfun(@times,(bsxfun(@plus,bsxfun(@times,y*k1,Iy),(k2*I))),Dy);
Ux=sum((Ux1+Ux2),2).*kDisp;

Uy1 = ((bsxfun(@plus,bsxfun(@times,y*k1,Iy),(-k2*I))).*Dx);
Uy2 =  bsxfun(@times,(bsxfun(@plus,bsxfun(@times,y*k1,Jy),(-2*J))),Dy);
Uy=sum((Uy1+Uy2),2).*kDisp;

% %%For Sxx, Eq [14] of "Birgisson" 
% Sxx=sum(-k3*((2*Iy+y*Iyy).*Dx)+(Jy+y*Jyy).*Dy);
% Syy=sum(-k3*-y*(Iyy.*Dx)+(Jy-y*Jyy).*Dy);
% Sxy=sum(-k3*((Jy+y*Jyy).*Dx)-y*(Iyy.*Dy));
Sxx1 = ((bsxfun(@plus,bsxfun(@times,y,Iyy),(2*Iy))).*Dx).*-k3;
Sxx2 = (bsxfun(@plus,(Jy),bsxfun(@times,y,Jyy))).*Dy;
Sxx=sum((Sxx1+Sxx2),2);

Syy1 = bsxfun(@times,(-k3.*-y),(Iyy.*Dx));
Syy2 = (bsxfun(@plus,(Jy),bsxfun(@times,-y,Jyy))).*Dy;
Syy=sum((Syy1+Syy2),2);

Sxy1 =  ((bsxfun(@plus,bsxfun(@times,y,Jyy),(Jy))).*Dx).*-k3;
Sxy2 =  bsxfun(@times,(-y),(Iyy.*Dy));
Sxy=sum((Sxy1+Sxy2),2);


%Converting stresses and displacements back to cart coords.  
%Rotating so element is back in respect to obs points 
Betam=-Beta;
[Ux,Uy] = RotateObject2d(Ux,Uy,Betam);

ax=Betam; 
nx=cos(ax);
ny=cos((pi/2)-ax);

[ Sxx,Syy,Sxy ] = StressTensorTransformation2d(Sxx,Syy,Sxy,nx,ny );

 StressDisp=[Sxx(:),Syy(:),Sxy(:),-Ux(:),-Uy(:)];
 
function [IJVect]=IJVec(p0,p1,p2,rsq3,two3rds,fr3rds)
%%Equation [12] of "Birgisson", getting the vectors.
%Could put this as a function inside in later versions on MATLAB

IJVect=[(-rsq3*p1)+(two3rds*p2),p0-(fr3rds*p2),(rsq3*p1)+(two3rds*p2)];


