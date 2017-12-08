function [ x, err, iter, flag, convergence, msg] = fischer_newton( A, b, x0, max_iter, tol_rel, tol_abs, solver, profile, lambda )
%Note this is an edited version of the code NUM4LCP.m
%Its recommended you use the original found at: 
%(Niebe 2016) = Niebe, S. and Erleben, K., 2015. Numerical methods for linear complementarity problems in physics-based animation. Synthesis Lectures on Computer Graphics and Animation, 7(1), pp.1-159.
% Copyright 2017, Tim Davis, The University of Aberdeen/Potsdam
% Copyright 2011, Kenny Erleben, DIKU
% Copyright 2012, Michael Andersen, DIKU

%	Original copyright below (Tim Davis, 2017): 
%%
%Copyright (c) 2011 Kenny Erleben code.google.com/p/num4lcp/
%
%This software is provided 'as-is', without any express or implied
%warranty. In no event will the authors be held liable for any damages
%arising from the use of this software.
%
%Permission is granted to anyone to use this software for any purpose,
%including commercial applications, and to alter it and redistribute it
%freely, subject to the following restrictions:
%
%    1. The origin of this software must not be misrepresented; you must not
%    claim that you wrote the original software. If you use this software
%    in a product, an acknowledgment in the product documentation would be
%    appreciated but is not required.
%
%    2. Altered source versions must be plainly marked as such, and must not be
%    misrepresented as being the original software.
%
%    3. This notice may not be removed or altered from any source
%    distribution.
%%

%Speedup
if (~issparse(A))
    %If not sparse they might be single, this func is not going to work if so.
    if (isa(A,'single'))
        A=double(A);
    end
    if (isa(b,'single'))
        b=double(b);
    end
  A = sparse(A);	% Make sure M is sparse
end
b = full(b(:)); 	% Make sure q is a column vector

% Just a list of human readable text strings to convert the flag return
% code into something readable by writing msg(flag) onto the screen.
msg = {'preprocessing';  % flag = 1
       'iterating';      % flag = 2
       'relative';       % flag = 3
       'absolute';       % flag = 4
       'stagnation';     % flag = 5
       'local minima';   % flag = 6
       'nondescent';     % flag = 7
       'maxlimit'        % flag = 8
      };

if nargin < 2
    error('Too few arguments');
end

N    = length(b); % Number of variables
flag = 1;

%--- Make sure we got good working default values -------------------------
if nargin<3
    x0 = zeros(N,1);
end
if nargin<4
    max_iter = floor(N/2);
end
if nargin<5
    tol_rel = 0.0001;
end
if nargin<6
    tol_abs = 10*eps; % Order of 10th of numerical precision seems okay
end
if nargin<7
    solver = 'random';
end
if nargin<8
    profile = false;
end
if nargin < 9
    switch lower(solver)
        case 'penalized'
            lambda = 0.95;
        otherwise
            lambda = 1.00; % No support for penalization in other solvers
    end
end

%--- Make sure all values are valid ---------------------------------------
max_iter = max(max_iter,1);
tol_rel  = max(tol_rel,0);
tol_abs  = max(tol_abs,0);
x0       = max(0,x0);

%--- Here comes a bunch of magic constants --------------------------------

%doesn't make much diff
h       = 1e-7;    % Fixed constant used to evaluate the directional detivative
%
alpha   = 0.5;%0.8;%0.5     % Step reduction parameter for projected Armijo backtracking line search
beta    = 0.001;   % Sufficent decrease parameter for projected Armijo backtracking line search
gamma   = 1e-28;   % Perturbation values used to fix near singular points in derivative
rho     = eps;     % Descent direction test parameter used to test if the Newton direction does a good enough job.

%--- Setup values need while iterating ------------------------------------

convergence = []; % Used when profiling to measure the convergence rate

err     = Inf;         % Current error measure
x       = x0;          % Current iterate
iter    = 1;           % Iteration count
display  =1;      
    
    
if display        
        
    %figure;
    figure('units','normalized','outerposition',[0.2 0.2 1/1.5 1/2])  %full screen  
    h1=subplot(1,2,1);
    title('LCP Solver Convergence','HandleVisibility','off');
    ylabel('psi (fisher burmister objective function), log scale','HandleVisibility','off');
    xlabel('time(s)','HandleVisibility','off');
    set(gca,'NextPlot','replacechildren') ; 
    
    h2=subplot(1,2,2); 
    title('Scatter of Fischer Burmeister values after each iteration','HandleVisibility','off');    
    ylabel('y=InputVectorSize','HandleVisibility','off');
    xlabel('x=Z','HandleVisibility','off');
    zlabel('z=FischerBurmeisterValue','HandleVisibility','off');    
    set(gca,'NextPlot','replacechildren') ;     

    %Setting up some more vars before the loop
    timed=0;
    psy=0;
    f_k=0;
    phi_k=zeros(N,1);
end   


solver='zero';   

% Init J
Zers    =  zeros(N,N); 
J       =  sparse(Zers) ;    
 
    
while (iter <= max_iter )
    
    y = A*x + b;
    
    %--- Test all stopping criteria used ------------------------------------
    phi     = phi_lambda(y, x, lambda);         % Calculate fischer function
    old_err = err;
    err     = 0.5*(phi'*phi);       % Natural merit function
    
    if profile
        convergence = [convergence err];
    end
    if (abs(err - old_err) / abs(old_err)) < tol_rel  % Relative stopping criteria
        flag = 3;
        break;
    end
    if err < tol_abs   % Absolute stopping criteria
        flag = 4;
        break;
    end
    
    %--- Solve the Newton system --------------------------------------------
    %--- First we find the Jacobian matrix. We wish to avoid singular points
    %--- in Jacobian of the Fischer function
    dx = zeros(N,1);                     % Allocate space for computing the Newton direction
    S       = abs(phi)<gamma & abs(x)<gamma;  % Bitmask for singular indices
    I       = find(S==0);                     % Bitmask for non-singular indices

    
    switch lower(solver)
        
        case 'random'    % works on full system
            
            q       = rand(N,1) ./ sqrt(2)  - 1;
            p       = rand(N,1) ./ sqrt(2)  - 1;
            J       = sparse( zeros(N,N) );
            q(I)    = (y(I)./((y(I).^2+x(I).^2).^0.5))-1;
            p(I)    = (x(I)./((y(I).^2+x(I).^2).^0.5))-1;
            J(I,I)  = diag(p(I))*eye(length(I)) + diag(q(I))*A(I,I);
            
            [dx, ~]  = gmres( J, (-phi), restart);
            
            
        case 'zero'   
        % works on reduced system would be similar to random-case if all
        % random values were set to 1
            
            %If comp is satisfied yi=xi=0
            %then pi(x)=ai-1
            %then qi(x)=bi-1
            p       = (x(I)./((y(I).^2+x(I).^2).^0.5))-1; %equation 2.165a (Niebe 2016)
            q       = (y(I)./((y(I).^2+x(I).^2).^0.5))-1; %equation 2.165b (Niebe 2016)
            
            %Tim Davis, adapted for the large dense matrices of my BEM
            %model. Faster in this situation than the original
            %implementation

            %%%%%%
            %A) Index then convert to sparse
            J       =  zeros(N,N) ;  %preallocating  
            J(I,I) = diag(kron(p,1)) + bsxfun(@times,(q),A(I,I));  %%equation 2.169 (Niebe 2016)
            J=sparse(J);               %sparse for later calcs  
            
            % %B) Correctly index the sparse each time. 
            % Data= diag(p) + bsxfun(@times,(q),A(I,I));
            % [indxb,indxa]=meshgrid(I,I); %IndxLocs
            % J=sparse(indxa,indxb,Data,N,N);
            %%%%%%
            
            %dx(I)=J(I,I)\-phi(I);                                 
            %equation ? (Niebe 2016)
                      
            
            if min(size(A,1),50)/2<30
                restart=min(size(A,1),50)/2;  
            else
            restart = 30; 
            end

            %Some different methods, '\' always returns a solution.
            %https://de.mathworks.com/help/matlab/math/systems-of-linear-equations.html
            %[dx(I), ~] = bicg      ( J(I,I), (-phi(I)));
            %[dx(I), ~] = bicgstab  ( J(I,I), (-phi(I)));
            %[dx(I), ~] = bicgstabl ( J(I,I), (-phi(I)));
            %[dx(I), ~] = cgs       ( J(I,I), (-phi(I)));
            [dx(I), ~] = gmres     ( J(I,I), -phi(I),restart);
            %[dx(I), ~] = lsqr      ( J(I,I), (-phi(I)));
            %[dx(I), ~] = minres    ( J(I,I), (-phi(I)));
            %[dx(I), ~] = qmr       ( J(I,I), (-phi(I)));
            %[dx(I), ~] = symmlq    ( J(I,I), (-phi(I)));
            %[dx(I), ~] = tfqmr     ( J(I,I), (-phi(I)));  
            %[dx(I), ~] = pcg     ( J(I,I), (-phi(I)));     
            %dx(I)=J(I,I)\-phi(I); 
            %dx(I)=Jfull(I,I)\-phi(I); %full Eq
            
            %doing some allocation for drawing
            pw= zeros(N,1);qw=pw;
            for i=1:numel(p)
                pw(I,1)=p(i);
                qw(I,1)=q(i);
            end
            
%             %Original NUM4LCP code, works but not as fast when working with the dense
%             %matrices coming out of BEM solution, (Comment: Tim).
%             J       = sparse( zeros(N,N) );            
%             J(I,I)  = diag(p)*eye(length(I)) + diag(q)*A(I,I);
%             restart = min(size(A,1),50);%10           
%             The number of iterations done before GMRES should restart
%             [dx(I), ~] = gmres( J(I,I), (-phi(I)), restart);
            

    otherwise
            disp('Unknown solver method for Newton subsystem')
    end
    
   %displaying convergence 
   if display        
       [timed,psy]=drawprep(f_k,psy,timed,h1,iter,max_iter);
       subplot(h2);scatter3(x,1:N,phi_k);
   end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Test if the search direction is smaller than numerical precision. 
    % That is if it is too close to zero.
    if max(abs(dx)) < eps
        flag = 5;
        % Rather than just giving up we may just use the gradient direction
        % instead. However, I am lazy here!
        %  dx = nabla_phi'
        break;
    end
    
    % Test if we have dropped into a local minimia if so we are stuck
    nabla_phi = phi'*J;
    if norm(nabla_phi) < tol_abs
        flag = 6;
        break;
    end
    
    % Test if our search direction is a 'sufficient' descent direction
    if  nabla_phi*dx  > -rho*(dx'*dx)
        flag = 7;
        % Rather than just giving up we may just use the gradient direction
        % instead. However, I am lazy here!
        %  dx = nabla_phi'
        break;
    end
    
    %--- Armijo backtracking combined with a projected line-search ---------
    tau     = 1.0;                  % Current step length
    f_0     = err;
    grad_f  = beta*(nabla_phi*dx);
    x_k     = x;
    
    while true
        
        x_k   = max(0,x + dx*tau);
        y_k   = A*x_k + b;
        phi_k = phi_lambda( y_k, x_k, lambda );
        f_k   = 0.5*(phi_k'*phi_k);
        
        % Perform Armijo codition to see if we got a sufficient decrease
        if ( f_k <= f_0 + tau*grad_f)
            break;
        end
        
        % Test if time-step became too small
        if tau*tau < gamma
            break;
        end
        
        tau = alpha*tau;
    end
    
    % Update iterate with result from Armijo backtracking
    x = x_k;
    
    % Increment the number of iterations
    iter = iter + 1;
    
end %end of while lp. 

if iter >= max_iter
    flag = 8;
    iter = iter - 1;
end


% Just return the message string, not entire cell.
msg = msg{flag};


end

function [ phi ] = fischer(y,x)
% Auxiliary function used by the Fischer-Newton method
% Copyright 2011, Kenny Erleben

phi  = (y.^2 + x.^2).^0.5 - y - x;

end

function phi_l = phi_lambda(a,b,lambda)
%% CCK NCP-function, a convex composition of Fischer-Burmeister and CCK NCP
%
%   Input
%       a -> A column vector size = (n,1)
%       b -> A column vector size = (n,1)
%       l -> A fixed lambda value used to weight the input.
%
%   Output
%       phi_l -> A column vector with the result of the Fischer-Burmeister
%                NCP function, with size = (n,1)

phi_l = lambda*fischer(a,b)+(1-lambda)*(max(0,a).*max(0,b));
%when calculated pretty much just the value from fischer(a,b)

end


function [timed,psy]=drawprep(f_k,psy,timed,h1,iter,max_iter)

%grabbing the times
timer=toc;
timed=[timed,timer];
psy=[psy,f_k];
      
subplot(h1)
draw2(timed,psy,max_iter)
    
if iter==1
    disp('iter = (iteration count), psi= Theta(z), p 449 in the LCP cottle pang stone')
end
%fprintf('iter = %, psi = %3.0e\n',iter,f_k);  
fprintf('iter = %2d, psi = %3.0e\n',iter,f_k); 
    
end


function draw2(times,psy,max_iter)
%drawing a the convergence on a log scale

%Nice curves
tt=linspace(min(times),max(times),max_iter*3);
yy=interp1(times,psy,tt,'PCHIP');
%clearing fig for each loop and redrawing with new vars
cla
plot(tt,yy,'g');
hold on
scatter(times,psy);
set(gca,'yscale','log')
set(gca,'xscale','linear')
drawnow;
hold on

end