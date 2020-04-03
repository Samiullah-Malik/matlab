
%% Question # 3


% 1. Input Data
a = 0.; b = 1;  % Domain [a,b]
tfinal = 1.;   
c = 0.5;        % Wave speed
cfl = 0.9;


ms = 'pulse';

% 2. Define Manufactured Polynomial and Traveling pulse.
% For manufactured polynomial solution
if( strcmp(ms, 'poly') )
    b0=.5; b1=.7; b2=.9; % time coefficients
    c0=1.; c1=.8; c2=.6; % space coefficients
    
    ue = @(x,t) (b0+b1*t+b2*t^2)*(c0+c1*x+c2*x.^2);
    uet = @(x,t) (b1+2.*b2*t)*(c0+c1*x+c2*x.^2);
    uett = @(x,t) (2.*b2)*(c0+c1*x+c2*x.^2);
    uex = @(x,t) (b0+b1*t+b2*t^2)*(c1+2.*c2*x);
    uexx = @(x,t) (b0+b1*t+b2*t^2)*(2.*c2);
    f = @(x,t) uett(x,t) - (c^2)*uexx(x,t);
end

% For traveling pulse solution
if( strcmp(ms, 'pulse') )
    beta = 20; x0 = 0.25;
    
    ue = @(x,t) exp( -(beta*(x - x0 - c*t) ).^2 );
    uet = @(x,t) (2.*c*beta^2 * (x - x0 - c*t)) .*ue(x,t);
    %uett = @(x,t) 
    uex = @(x,t) -2.*(beta^2)*(x-x0-ct).*ue(x,t);
    uexx = @(x,t) -2.*beta^2 * (-2.*ue(x,t) * beta^2 * (x-x0-c*t)^2 ...
                    + ue(x,t));
    uexx = @(x,t) ( -(2.*beta^2) + (-2.*beta^2*(x-x0-c*t)).^2 ).* ue(x,t);
    f = @(x,t) 0;  %% used symbolab.com for symbolic computation
end    
   
% 3. Define BCs and ICs
ga = @(t) ue(a,t);      % BC at x = a (LHS)
gb = @(t) ue(b,t);      % BC at x = b (RHS)
u0 = @(x) ue(x,0);      % IC for u0(x)
u1 = @(x) uet(x,0);     % IC for u1(x)


% 4(a). Grid Resolutions (1,2,3,4,5)
numResolutions=2;
for m=1:numResolutions
    
    % (b). Setup Grid
    Nx=10*2^m;                  % no. of grid intervals
    dx=(b-a)/Nx;                % size of each interval
    ia=1;                       % index of LHS boundary point at x=a
    ib=ia+Nx;                   % index of RHS boundary point at x=b
    i1=ia+1;                    % first interior index 
    i2=ib-1;                    % last interior index
    I = i1:i2;                  % interior points (PDE is applied here)
    x = linspace(a,b,Nx+1)';    % grid points

    % (c). Allocating space for the solution. (3-levels).
    unm1 = zeros(Nx+1,1);
    un = zeros(Nx+1,1);
    unp1 = zeros(Nx+1,1);
    
    % (d). Adjust time step
    dt = cfl*dx/c;              
    Nt = round(tfinal/dt);      % no. of time-steps
    dt = tfinal/Nt;             % adjust dt to reach tfinal exactly
    
    % (e). Apply Initial Conditions
    unm1 = u0(x);               % at t=0
    un = unm1 + dt*u1(x) + ((dt^2)/2.) * ((c^2) * uexx(x,0) + f(x,0)); % at t=dt
    t=dt;
    
    % (f). Apply time-stepping loop    
    for n = 2:Nt       
        t_current = dt*(n-1);
        t_new = dt*n;        
                
        unp1(I) = 2.*un(I) - unm1(I) + (c*dt/dx)^2*( un(I+1)-2*un(I)+un(I-1) ) + (dt^2)*f(x(I),t_current);
        unp1(ia) = ga(t_new);       % apply BCs on LHS
        unp1(ib) = gb(t_new);       % apply BCs on RHS
        
        unm1 = un;           % set unm1 un and un for next time step
        un = unp1;  

    end
    
    % (g). Compute errors and print
    u_exact_solution = ue(x,tfinal);
    error = un - u_exact_solution;
    errMax(m) = max(abs(error(ia:ib)));
    fprintf('ms=%s: t=%8.2e: m=%d Nx=%3d Nt=%4d dt=%8.2e maxErr=%8.2e',ms,t,m,Nx,Nt,dt,errMax(m));
    if( m==1 ) fprintf('\n'); else fprintf(' ratio=%5.2f rate=%5.2f\n', ... 
        errMax(m-1)/errMax(m),log2(errMax(m-1)/errMax(m))); end

    % (h). Plot for pulse
    if( strcmp(ms,'pulse') )
        subplot(2,1,1)
        plot( x(ia:ib) ,un, 'r-', x, u_exact_solution, 'b-')
        legend('computed', 'exact')
        title('Wave Equation Pulse - Computed vs Exact Solution - m=2')
        
        subplot(2,1,2)        
        plot( x(ia:ib) ,error, 'r-')
        legend('error')
        title('Wave Equation Pulse - Error - m=2')
    end
    
end


