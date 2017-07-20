function [ x err iter flag convergence msg] = fischer_newton3d( A, b, x0, max_iter, tol_rel, tol_abs, solver, profile, lambda )
%EDITED BY TIM DAVIS
%IF INSPIRED PLEASE USE THE ORIGINAL NUM4LCP CODE, NOT THIS. 
%Parts is the how many parts I am solving, different bits converge
%differently. I draw these seperatly. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

% Copyright 2011, Kenny Erleben, DIKU
% Copyright 2012, Michael Andersen, DIKU

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
    display  =0;      
   if display        
        
    figure;
    h1=subplot(5,2,1);
    title({'Newtonian convergence';'looking for the root of this function';'5 parts shown, mean avg of each subsection of the solution vec'}),
    ylabel('DnElDisp'); 
    h2=subplot(5,2,3);
    ylabel('DssElDisp+'); 
    h3=subplot(5,2,5);
    ylabel({'Avg Calculated f(Z) in comp formulation';'DdsElDisp+'}); 
    h4=subplot(5,2,7);
    ylabel('tssElDisp-');
    h5=subplot(5,2,9);  
    ylabel('tdsElDisp-');
    xlabel('Avg Calculated Z, in comp formulation, see Ritz 2012'); 
    h6=subplot(5,2,[2,4,6,8,10]);
    title('Fric Solver Convergence');
    xlabel('time(s)'); ylabel('psi'); 
    
   end   
   %Setting up some more vars before the loop
    timed=0;
    psy=0;
    f_k=0;
    solver='zero';
    %lambda = 0.95;    
            %random
            %perturbation
            %zero
            %approximation
            %penalized
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
    restart = min(size(A,1),10);              % The number of iterations done before GMRES should restart
    
    switch lower(solver)
        
        case 'random'    % works on full system
            
            q       = rand(N,1) ./ sqrt(2)  - 1;
            p       = rand(N,1) ./ sqrt(2)  - 1;
            J       = sparse( zeros(N,N) );
            q(I)    = (y(I)./((y(I).^2+x(I).^2).^0.5))-1;
            p(I)    = (x(I)./((y(I).^2+x(I).^2).^0.5))-1;
            J(I,I)  = diag(p(I))*eye(length(I)) + diag(q(I))*A(I,I);
            
            [dx, ~]  = gmres( J, (-phi), restart);
            
        case 'perturbation'     % works on full system
            
            px          = x;
            dir         = sign( x );
            dir(dir==0) = 1;
            px(S==1)    = gamma * dir(S==1);
            
            p = (px./((y.^2+px.^2).^0.5))-1;
            q = ( y./((y.^2+px.^2).^0.5))-1;
            J = diag(p)*eye(N) + diag(q)*A;
            
            % Billups solution actually pick the gradient direction as the descent
            % direction if the Newton system can not be solved. Look in his
            % thesis at page 85.
            %
            % %if rcond(J) < eps
            % %  dx = -J'*phi;          % Too bad we pick gradient direction
            % %end
            %
            %Instead one may use the pseudo inverse if badly conditioned
            %
            if rcond(J) < eps*2
                dx = pinv(J) * (-phi); % badly conditioned
            else
                dx = J \ (-phi);       % well conditioned
            end
            %
            % Or one can hope for the best and try GMRES
            %
            %[dx ~]  = gmres( J, (-phi), restart);
            %
            
        case 'zero'   % works on reduced system would be similar to random-case if all radom values were set to 1
            
            J       = spalloc( N,N,numel(I) );            
            q       = (y(I)./((y(I).^2+x(I).^2).^0.5))-1;
            p       = (x(I)./((y(I).^2+x(I).^2).^0.5))-1;
            FM = diag(kron(p,1)) + bsxfun(@times,diag(q),A(I,I));
            
            %putting into sparse mat
            for lp=1:numel(I)
            curr=I(lp)    ;
            J(curr,curr)=FM(lp,lp);
            end
            
            %old code, not good for my dense matrices, Tim.             
            %J(I,I)  = diag(p)*eye(length(I)) + diag(q)*A(I,I);

            [dx(I), ~] = gmres( J(I,I), (-phi(I)), restart);
            

        case 'approximation'  % works on reduced system assuming zero Jacobian for all bad values
            
            J       = sparse( zeros(N,N) );
            q       = (y(I)./((y(I).^2+x(I).^2).^0.5))-1;
            p       = (x(I)./((y(I).^2+x(I).^2).^0.5))-1;
            J(I,I)  = diag(p)*eye(length(I)) + diag(q)*A(I,I);
            
            fun       = @(dx) ( fischer( A(I,I)*(x(I)+ dx*h) + b(I), x(I) + dx*h) - phi(I) ) / h;
            [dx(I), ~] = gmres( fun, (-phi(I)), restart);
            
        case 'penalized'
            % Produces an element of V_{k} \in \partial_C \Phi_{\lambda}
            % according to algorithm 2.4 of A Penalized Fischer-Burmeister
            % Ncp-Function: Theoretical Investigation And Numerical Results
            % (see p. 8).
            SMALL = 1e-10; % Kanzow uses 1e-8
            z = zeros(N,1);
            S1 = abs(x) <= SMALL & abs(y) <= SMALL;
            S2 = x > SMALL & y > SMALL;
            NS = ~(S1 | S2);
            Da = zeros(N,1);
            Db = zeros(N,1);
            %%% Does not really seem to matter whether we use rand or ones!
            z(S1) = ones(length(find(S1)),1);
            % z(S1) = 0.5*ones(length(find(S1)),1);
            % z(S1) = rand(length(find(S1)),1);
            
            denomS1 = sqrt(z(S1).^2+(A(S1,S1)*z(S1)).^2);
            Da(S1) = lambda*(z(S1)./denomS1-1);
            Db(S1) = lambda*((A(S1,S1)*z(S1))./denomS1-1);
            
            denomS2 = sqrt(x(S2).^2+y(S2).^2);
            Da(S2) = lambda*(x(S2)./denomS2-1)+(1-lambda)*y(S2);
            Db(S2) = lambda*(y(S2)./denomS2-1)+(1-lambda)*x(S2);
            
            denomNS = sqrt(x(NS).^2+y(NS).^2);
            Da(NS) = lambda*(x(NS)./denomNS-1);
            Db(NS) = lambda*(y(NS)./denomNS-1);
            
            % J = Vk
            J = diag(Da) + diag(Db)*A;
            
            dx = J \ (-phi);
            


            
        otherwise
            disp('Unknown solver method for Newton subsystem')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TIM ADDED DRAWING CONVERGENCE
   if display        
    %3d, so splitting vector into 5 parts and calculating means
    sz=numel(y)/5; %only works for 3d!!!
    %1
    ydraw1=mean(y(1:sz));
    xloc1=mean(x(1:sz));
    %2
    ydraw2=mean(y(sz+1:sz*2));
    xloc2=mean(x(sz+1:sz*2));
    %3
    ydraw3=mean(y(sz*2+1:sz*3));
    xloc3=mean(x(sz*2+1:sz*3));
    %4
    ydraw4=mean(y(sz*3+1:sz*4));
    xloc4=mean(x(sz*3+1:sz*4));    
    %5
    ydraw5=mean(y(sz*4+1:sz*5));
    xloc5=mean(x(sz*4+1:sz*5));
    %6 grabbing the times
    timer=toc;
    timed=[timed,timer];
    psy=[psy,f_k];
    
    %intitiating vectorss
    YesNo = any(strcmp(who,'yD1'));
    if YesNo==0
        yD1=[];
        xD1=[];
        yD2=[];
        xD2=[];        
        yD3=[];
        xD3=[];
        yD4=[];
        xD4=[];        
        yD5=[];
        xD5=[];        
    end    
    %conc vects these are all the points from previous loops 
    yD1=[yD1,ydraw1];
    xD1=[xD1,xloc1];
    yD2=[yD2,ydraw2];
    xD2=[xD2,xloc2];
    yD3=[yD3,ydraw3];
    xD3=[xD3,xloc3];    
    yD4=[yD4,ydraw4];
    xD4=[xD4,xloc4];
    yD5=[yD5,ydraw5];
    xD5=[xD5,xloc5];    
    
    %now drawing each plot
    %if you do not want these as subplots just change the figure handles at
    %the top of this func
    subplot(h1)
    draw(xD1,yD1,iter)
    subplot(h2)
    draw(xD2,yD2,iter)
    subplot(h3)
    draw(xD3,yD3,iter)
    subplot(h4)
    draw(xD4,yD4,iter)
    subplot(h5)
    draw(xD5,yD5,iter)
    
    subplot(h6)
    draw2(timed,psy,max_iter)
    if iter==1
    disp('iter = (iteration count), psi= Theta(z), p 449 in the LCP cottle pang stone')
    end
    disp(sprintf('iter = %2d, psi = %3.0e',iter,f_k));    
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
% Copyright 2011, Kenny Erleben, DIKU
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

end

function psi_l = psi_lambda(a,b,l)
%% Natural merit function of phi_lambda. psi_lambda : R^n -> R
%
%   Input:
%       a -> A column vector of size = (n,1)
%       b -> A column vector of size = (n,1)
%       l -> Real value describing weight of penalization term
%
%   Output:
%       psi_l -> Merit values at x^(k) based on phi_l(x^(k))

phi_l = phi_lambda(a,b,l);
psi_l = 0.5*(phi_l'*phi_l);

end

function draw(xD,yD,iter)
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
    %Nice curves
    tt=linspace(min(times),max(times),max_iter*3);
    yy=interp1(times,psy,tt,'PCHIP');
    %clearing fig for each loop and redrawing with new vars
    cla
    plot(tt,yy,'g');
    hold on
    scatter(times,psy);
    drawnow;
    hold on
end
