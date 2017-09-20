function [ x err iter flag convergence msg] = fischer_newton2d( A, b, x0, max_iter, tol_rel, tol_abs, solver, profile, lambda )
%EDITED BY TIM DAVIS
%IF INSPIRED PLEASE USE THE ORIGINAL NUM4LCP CODE, NOT THIS. 
%Parts is the how many parts I am solving, different bits converge
%differently. I draw these seperatly. 
%(Niebe 2016) = Niebe, S. and Erleben, K., 2015. Numerical methods for linear complementarity problems in physics-based animation. Synthesis Lectures on Computer Graphics and Animation, 7(1), pp.1-159.

%   Copyright 2017, Tim Davis, The University of Aberdeen

% Copyright 2011, Kenny Erleben, DIKU
% Copyright 2012, Michael Andersen, DIKU

%Speedup?
if (~issparse(A))
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
iter    = 1;           % Iteration counter

flag = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TIM ADDED DRAWING CONVERGENCE
    display  =1;      
    if display        
        
    %figure;
    figure('units','normalized','outerposition',[0 0 1 1])  %full screen  
    h1=subplot(5,2,1);
    title({'Newtonian convergence';'looking for the root of this function';'3 parts shown, mean avg of each subsection of the solution vec'}),
    ylabel('DnElDisp'); 
    h2=subplot(5,2,3);
    ylabel('DsElDisp+'); 
    h3=subplot(5,2,5);
    ylabel({'Avg Calculated FB function in comp formulation';'Ts-'}); 
    xlabel('Avg Calculated Z, in comp formulation, see Ritz 2012'); 
    h6=subplot(5,2,[2,4]);
    title('LCP Solver Convergence','HandleVisibility','off');
    ylabel('psi (fisher burmister objective function), log scale','HandleVisibility','off');
    xlabel('time(s)','HandleVisibility','off');
    set(gca,'NextPlot','replacechildren') ; 
    h7=subplot(5,2,[8,10]); 
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
    yD1=[];
    xD1=[];
    yD2=[];
    xD2=[];        
    yD3=[];
    xD3=[];
    end   
    solver='zero';    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
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
            
            
        case 'zero'   % works on reduced system would be similar to random-case if all radom values were set to 1
            
            %If comp is satisfied yi=xi=0
            %then pi(x)=ai-1
            %then qi(x)=bi-1
            p       = (x(I)./((y(I).^2+x(I).^2).^0.5))-1; %equation 2.165a (Niebe 2016)
            q       = (y(I)./((y(I).^2+x(I).^2).^0.5))-1; %equation 2.165b (Niebe 2016)
            
            %Tim Davis, adapted for the large dense matrices of my BEM
            %model. Faster in this situation than the original
            %implementation
            %The calculation of J takes the longest in this function 
            J       =  zeros(N,N) ;  %preallocating  
            J(I,I) = diag(kron(p,1)) + bsxfun(@times,(q),A(I,I)); %equation 2.169 (Niebe 2016)
%            dx(I)=J(I,I)\-phi(I);                                 %equation ? (Niebe 2016)
            J=sparse(J);               %sparse for later calcs
            restart = min(size(A,1),50);%10           % The number of iterations done before GMRES should restart
            %http://www.ms.uky.edu/~carl/ma123/kob98/kob98htm/chap15e.html
            %see dx on first diagram
            [dx(I), ~] = gmres( J(I,I), (-phi(I)), restart);
            %when calculated pretty much just the value from
            %dx(I)=J(I,I)\-phi(I); On big problems the gmres approach is far
            %better at finding the minimum 
            
            %doing some allocation for drawing
            pw= zeros(N,1);qw=pw;
            for i=1:numel(p)
                pw(I,1)=p(i);
                qw(I,1)=q(i);
            end
            
%             %old code, works but not as fast when working with the dense
%             matrices coming out of BEM solution, Tim.
%             J       = sparse( zeros(N,N) );            
%             J(I,I)  = diag(p)*eye(length(I)) + diag(q)*A(I,I);
%             restart = min(size(A,1),50);%10           % The number of iterations done before GMRES should restart
%             [dx(I), ~] = gmres( J(I,I), (-phi(I)), restart);
            

        otherwise
            disp('Unknown solver method for Newton subsystem')
    end
    
   %displaying convergence 
   if display        
   [yD1,xD1,yD2,xD2,yD3,xD3,timed,psy]=drawprep(x,phi_k,f_k,psy,timed,h1,h2,h3,h6,iter,max_iter,yD1,xD1,yD2,xD2,yD3,xD3,N);
   subplot(h7);scatter3(x,1:N,phi_k);
   end    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Test if the search direction is smaller than numerical precision. That is if it is too close to zero.
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
end

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





function [yD1,xD1,yD2,xD2,yD3,xD3,timed,psy]=drawprep(x,y,f_k,psy,timed,h1,h2,h3,h6,iter,max_iter,yD1,xD1,yD2,xD2,yD3,xD3,sz);
  
%setting up variables for the sub functions 'draw' and 'draw2' below


    
    sz=sz/3;
    %mean of each vector
    %1
    ydraw1=mean(y(1:sz));
    xloc1=mean(x(1:sz));
    %2
    ydraw2=mean(y(sz+1:sz*2));
    xloc2=mean(x(sz+1:sz*2));
    %3
    ydraw3=mean(y(sz*2+1:sz*3));
    xloc3=mean(x(sz*2+1:sz*3));

%     %first point from each vector
%     %1
%     ydraw1=mean(y(1));
%     xloc1=mean(x(1));
%     %2
%     ydraw2=mean(y(sz+1));
%     xloc2=mean(x(sz+1));
%     %3
%     ydraw3=mean(y(sz*2+1));
%     xloc3=mean(x(sz*2+1));
%     %4
%     ydraw4=mean(y(sz*3+1));
%     xloc4=mean(x(sz*3+1));    
%     %5
%     ydraw5=mean(y(sz*4+1));
%     xloc5=mean(x(sz*4+1));

    %6 grabbing the times
    timer=toc;
    timed=[timed,timer];
    psy=[psy,f_k];
    
    %conc vects these are all the points from previous loops 
    yD1=[yD1,ydraw1];
    xD1=[xD1,xloc1];
    yD2=[yD2,ydraw2];
    xD2=[xD2,xloc2];
    yD3=[yD3,ydraw3];
    xD3=[xD3,xloc3];    
  
    
    %now drawing each plot
    %if you do not want these as subplots just change the figure handles at
    %the top of this func
    subplot(h1)
    draw(xD1,yD1,iter)
    subplot(h2)
    draw(xD2,yD2,iter)
    subplot(h3)
    draw(xD3,yD3,iter)

    subplot(h6)
    draw2(timed,psy,max_iter)
    
    if iter==1
    disp('iter = (iteration count), psi= Theta(z), p 449 in the LCP cottle pang stone')
    end
    disp(sprintf('iter = %2d, psi = %3.0e',iter,f_k));   
end

function draw(xD,yD,iter)

%drawing the Newton Raphson diagrams explained here- 
%https://www.comsol.com/blogs/solving-nonlinear-static-finite-element-problems/
%the blue dotted line is the curve. we are trying to find where this
%crosses the x axis for every point.
%the greenline is the most recent point the gradient has been calculated
%at.
%the black lines are previous estimations. 

%vectors and the figure number
    %only starting drawing in loop if we have two sets of xy points
    if iter>1.
    %clearing and putting hold on
    cla
    %clf('reset')
    hold on 
    %drawing axis lines through origin
    ax = gca;
    ax.XAxisLocation = 'origin';
    %this loop draws a single vector every loop. After finished all
    %approxmimations should show    
    for i=1:iter
    if i>1 
    %drawing all the lines from xy to new x & y=0, would good if these had
    %arrow heads like in comsol diagram
    line('XData',[xD(i-1);xD(i)],'YData',[yD(i-1);0])
    %on the last bit of this loop updating the most recent line so its 
    %green and drawing a line connecting all points
    if i==iter
    %drawing a single STRAIGHT line connecting the points, we want to know where
    %this passes the x axis (its 'root')
    DataXY= [xD((1):(i-1))',yD((1):(i-1))'];
    DataXsorted = sortrows(DataXY);
    h=line('XData',DataXsorted(:,1),'YData',DataXsorted(:,2),'color','b');
    h.LineStyle = ':';     
    %making the most recent line green so you know what its latest
    %approximation is
    line('XData',[xD(i-1);xD(i)],'YData',[yD(i-1);0],'color','g')
    if iter<=3
    %do nothing with axes    
    elseif iter<5
    biggest=max(DataXY(end-(iter-2):end,1));
    smallest=min(DataXY(end-(iter-2):end,1));
    if biggest~=smallest %stops error
    xlim([smallest biggest])
    end
    %Yaxis
    biggest=max(DataXY(end-(iter-2):end,2));
    smallest=min(DataXY(end-(iter-2):end,2));
    if biggest~=smallest
    ylim([smallest biggest])
    end
    else
    %setting axis, around latest 4 lines
    %Xaxis
    biggest=max(xD(1,end-3:end)); %max(DataXY(end-3:end,1)); 
    smallest=min(xD(1,end-3:end));
    xlength=biggest-smallest;
    if biggest~=smallest
    xlim([smallest-(xlength/2) biggest+(xlength/2)])
    end
    %Yaxis
    biggest=max(yD(1,end-3:end)); %max(DataXY(end-3:end,1)); 
    smallest=min(yD(1,end-3:end));
    ylength=biggest-smallest;
    if biggest~=smallest
    ylim([smallest-(ylength/2) biggest+(ylength/2)])
    end
    end
    end    
    end
    end
    drawnow;        
    end
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
