function output = Esh_sol_NOPARA(inp,varargin)
%poisson ratio of matrix
vm=inp.vm;
%elastic modulus of matrix
Em=inp.Em;
%hetergeneity poisson ratio
vh=inp.vh;
%hetergeneity elastic modulus
Eh=inp.Eh;
%converts isotropic constants into a matrix
Cm=Ctensord(Em,vm);
Ch=Ctensord(Eh,vh);

%dimensiona of the ellipsoid. Must always be a1>=a2>=a3 and ortation angles
%[alpha beta sigam]
dim = inp.dim;
ang = inp.ang;

% fake a1>=a2>=a3
exh = zeros(3,3);
for i=1:2
    for j=2:3
        if dim(i)<dim(j)
            exh(i,j) = 1;
            tmp = dim(i);
            dim(i) = dim(j);
            dim(j) = tmp;
        end
    end
end
% pre-rotation in order of [z,y,x]
ang_init = pi/2*[exh(2,3) exh(1,3) exh(1,2)];
Rx = [1 0 0;0 cos(ang_init(1)) -sin(ang_init(1));0 sin(ang_init(1)) cos(ang_init(1))];
Ry = [cos(ang_init(2)) 0 sin(ang_init(2));0 1 0;-sin(ang_init(2)) 0 cos(ang_init(2))];
Rz = [cos(ang_init(3)) -sin(ang_init(3)) 0;sin(ang_init(3)) cos(ang_init(3)) 0;0 0 1];
R_init = Rx*Ry*Rz;
angb_init = -ang_init;
Rx = [1 0 0;0 cos(angb_init(1)) -sin(angb_init(1));0 sin(angb_init(1)) cos(angb_init(1))];
Ry = [cos(angb_init(2)) 0 sin(angb_init(2));0 1 0;-sin(angb_init(2)) 0 cos(angb_init(2))];
Rz = [cos(angb_init(3)) -sin(angb_init(3)) 0;sin(angb_init(3)) cos(angb_init(3)) 0;0 0 1];
Rb_init = Rz*Ry*Rx;

delC=Cm-Ch;
%stress ordering is: sigma11, sigma12, sigma13, sigma22,sigma23,sigma33
%this is the applied stress
stressvec=inp.stressvec;
eigp=inp.eigp;
% rotation matrices w.r.t the ellipsoid
Rx = [1 0 0;0 cos(ang(1)) -sin(ang(1));0 sin(ang(1)) cos(ang(1))];
Ry = [cos(ang(2)) 0 sin(ang(2));0 1 0;-sin(ang(2)) 0 cos(ang(2))];
Rz = [cos(ang(3)) -sin(ang(3)) 0;sin(ang(3)) cos(ang(3)) 0;0 0 1];
R = Rz*Ry*Rx;
Rx = [1 0 0;0 cos(-ang(1)) -sin(-ang(1));0 sin(-ang(1)) cos(-ang(1))];
Ry = [cos(-ang(2)) 0 sin(-ang(2));0 1 0;-sin(-ang(2)) 0 cos(-ang(2))];
Rz = [cos(-ang(3)) -sin(-ang(3)) 0;sin(-ang(3)) cos(-ang(3)) 0;0 0 1];
Rb = Rx*Ry*Rz;
% rotate stress against oblique ellipsoid
stressten_init = [stressvec(1:3)';stressvec([2 4 5])'; stressvec([3 5 6])'];
stressten_rot = R_init*Rb*stressten_init*(R_init*Rb)';
stressvec = [stressten_rot(1,:) stressten_rot(2,[2,3]) stressten_rot(3,3)]';
%correspondingly, the applied strain
epsvec=Cm\stressvec;
%call the internal eshelby tensor.
S4=Eshint(vm,dim);
eigen=(delC*S4-Cm)\(-delC*epsvec-Ch*eigp);
% for all the obersvations grids
Ng = length(inp.grid);

for n=1:Ng
    x = inp.grid{n}{1};
    y = inp.grid{n}{2};
    z = inp.grid{n}{3};
    [X Y Z] = meshgrid(x,y,z);
    Nx = size(x,2);
    Ny = size(y,2);
    Nz = size(z,2);
    u = [];
    rD4 = [];
    stress = [];
    strain = [];
    if Nx == max([Nx Ny Nz])
        for k=1:Nz;
            for i=1:Ny
                if length(varargin)>=2 && ismember('disp',lower(varargin))
                    for j=1:Nx
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        [D4 disp] = Esh_D4_disp(vm,dim,pos,eigen);
                        rD4(i,j,k,:,:)=Cmatrix(D4);
                        u(i,j,k,:) = R*Rb_init*disp';
                    end
                elseif length(varargin)>=1 && ~ismember('disp',lower(varargin))
                    for j=1:Nx
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        D4 = Esh_D4(vm,dim,pos);
                        rD4(i,j,k,:,:)=Cmatrix(D4);
                    end
                else
                    for j=1:Nx
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        u(i,j,k,:) = R*Rb_init*Esh_disp(vm,dim,pos,eigen)';
                    end
                end
            end
        end
    elseif Ny == max([Nx Ny Nz])
        for k=1:Nz;
            for j=1:Nx
                if length(varargin)>=2 && ismember('disp',lower(varargin))
                    for i=1:Ny
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        [D4 disp] = Esh_D4_disp(vm,dim,pos,eigen);
                        rD4(i,j,k,:,:)=Cmatrix(D4);
                        u(i,j,k,:) = R*Rb_init*disp';
                    end
                elseif length(varargin)>=1 && ~ismember('disp',lower(varargin))
                    for i=1:Ny
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        D4 = Esh_D4(vm,dim,pos);
                        rD4(i,j,k,:,:)=Cmatrix(D4);
                    end
                else
                    for i=1:Ny
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        u(i,j,k,:) = R*Rb_init*Esh_disp(vm,dim,pos,eigen)';
                    end
                end
            end
        end
    else
        for j=1:Nx
            for i=1:Ny
                if length(varargin)>=2 && ismember('disp',lower(varargin))
                    for k=1:Nz
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        [D4 disp] = Esh_D4_disp(vm,dim,pos,eigen);
                        rD4(i,j,k,:,:)=Cmatrix(D4);
                        u(i,j,k,:) = R*Rb_init*disp';
                    end
                elseif length(varargin)>=1 && ~ismember('disp',lower(varargin))
                    for k=1:Nz
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        D4 = Esh_D4(vm,dim,pos);
                        rD4(i,j,k,:,:)=Cmatrix(D4);
                    end
                else
                    for k=1:Nz
                        pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                        u(i,j,k,:) = R*Rb_init*Esh_disp(vm,dim,pos,eigen)';
                    end
                end
            end
        end
    end
    fstress = ismember('stress',lower(varargin));
    fstrain = ismember('strain',lower(varargin));
    if ~isempty(rD4)
        for j=1:Nx
            for i=1:Ny
                for k=1:Nz
                    pos = (R_init*Rb*[X(i,j,k) Y(i,j,k) Z(i,j,k)]')';
                    if  pos(1)^2/dim(1)^2+pos(2)^2/dim(2)^2+pos(3)^2/dim(3)^2<=1 % for interior points
                        if fstrain
                            strainr = epsvec+S4*eigen;
                            strainten = [strainr(1:3)';strainr([2 4 5])'; strainr([3 5 6])'];
                            strainten = R*Rb_init*strainten*(R*Rb_init)';
                            strain(i,j,k,:) = [strainten(1,:) strainten(2,[2,3]) strainten(3,3)]';
                        end
                        if fstress
                            stressr = stressvec+Cm*(S4*eigen-eigen);
                            stressrten = [stressr(1:3)';stressr([2 4 5])'; stressr([3 5 6])'];
                            stressrten = R*Rb_init*stressrten*(R*Rb_init)';
                            stress(i,j,k,:) = [stressrten(1,:) stressrten(2,[2,3]) stressrten(3,3)]';
                        end
                    else % exterior points
                        if fstress
                            stressr = stressvec+Cm*squeeze(rD4(i,j,k,:,:))*eigen;
                            stressrten = [stressr(1:3)';stressr([2 4 5])'; stressr([3 5 6])'];
                            stressrten = R*Rb_init*stressrten*(R*Rb_init)';
                            stress(i,j,k,:) = [stressrten(1,:) stressrten(2,[2,3]) stressrten(3,3)]';
                        end
                        if fstrain
                            strainr = epsvec+squeeze(rD4(i,j,k,:,:))*eigen;
                            strainten = [strainr(1:3)';strainr([2 4 5])'; strainr([3 5 6])'];
                            strainten = R*Rb_init*strainten*(R*Rb_init)';
                            strain(i,j,k,:) = [strainten(1,:) strainten(2,[2,3]) strainten(3,3)]';
                        end
                    end
                end
            end
        end
    end
    output.grid{n,1} = X;
    output.grid{n,2} = Y;
    output.grid{n,3} = Z;
    output.u{n} = u;
    output.stress{n} = stress;
    output.strain{n} = strain;   
end