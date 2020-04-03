
%% Question 2

% 1. Input Data
a = 0.; b = 2*pi;       % Domain [a,b]
tfinal = pi/2.;   
c = 1.;                 % Wave speed
kx = 6.;
alpha = 1.; beta = 3.;  % top_hat parameters

ms = 'rk3-upw2';
solution = 'top_hat';
%ms = 'upw1-case2';

% 2. Set cfl for different schemes.
if( strcmp(ms, 'upw1-case1') )
    cfl = 1.;
elseif ( strcmp(ms, 'upw1-case2') )
    cfl = .9;
elseif ( strcmp(ms, 'rk3-upw2') )
    cfl = .5;
end

% 3. Define functions
if ( strcmp(solution, 'sine') )
    uexact = @(x,t) sin(kx*(x-c*t));
elseif ( strcmp(solution, 'top_hat') )
    uexact = @(x,t) periodicTopHat(x-c*t, alpha, beta);
end
u0 = @(x) uexact(x,0);    % Given initial condition function



% 4(a). Grid Resolutions (1,2,3,4,5,6)
numResolutions=6;
for m=1:numResolutions
    
    % (b). Setup Grid
    Nx=10*2^m;           % no. of grid intervals
    dx=(b-a)/Nx;         % size of each interval
    lgp = 2;             % no. of left ghost points
    ia=lgp+1;            % index of LHS boundary point at x=a
    ib=ia+Nx;            % index of RHS boundary point at x=b
    i1=ia;               % first interior index (due to left ghost points)
    i2=ib-1;             % last interior index
    I = i1:i2;           % interior points (PDE is applied here)
    x = linspace(a-(lgp*dx),b,Nx+ia)';    % grid (coordinate values)

    % (c). Allocating space for the solution. (3-levels).
    un = zeros(Nx+ia,1);
    unp1 = zeros(Nx+ia,1);
    if ( strcmp(ms,'rk3-upw2') )
        k1 = zeros(Nx+ia,1);
        k2 = zeros(Nx+ia,1);
        k3 = zeros(Nx+ia,1);
        temp = zeros(Nx+ia,1);
    end
    
    % (d). Adjust time step
    dt = cfl*dx/c;              
    Nt = round(tfinal/dt);      % no. of time-steps
    dt = tfinal/Nt; %           adjust dt to reach tfinal exactly
    
    % (e). Define difference operators
    Dp = @(U,I) ( U(I+1) - U(I) )/ dx;
    Dm = @(U,I) ( U(I) - U(I-1) )/ dx;
    Dm2 = @(U,I) ( 3*U(I) - 4*U(I-1) + U(I-2) )/(2*dx); %upwind scheme
    
    % (f). Apply Initial Conditions
    un = u0(x);               % at t=0
    
    % (g). Apply time-stepping loop    
    for n = 1:Nt       
        t_current = dt*(n-1);
        t_new = dt*n;       
        
        if ( strcmp(ms, 'upw1-case1') || strcmp(ms, 'upw1-case2')  )        
            unp1(I) = un(I) - (c*dt/dx)*( un(I)-un(I-1));
            unp1 = applyBCs(lgp,ia,Nx, unp1);        
        elseif ( strcmp(ms, 'rk3-upw2') )
            k1(I) = -c*Dm2(un,I);               % Stage 1
            temp(I) = un(I) + dt*k1(I);
            temp = applyBCs(lgp,ia,Nx,temp);
            
            k2(I) = -c*Dm2(temp,I);             % Stage 2
            temp(I) = un(I) + (dt/4.)*(k1(I) + k2(I));
            temp = applyBCs(lgp, ia, Nx, temp);
            k3(I) = -c*Dm2(temp,I);
            
            
            unp1(I) = un(I) + (dt/6.)*(k1(I) + k2(I) + 4*k3(I));
            unp1 = applyBCs(lgp, ia, Nx, unp1);         
        end   
        un = unp1; 
    end % time-stepping loop
    
    % (h). Compute errors and print
    u_exact_solution = uexact(x,tfinal);
    error = un - u_exact_solution;
    maxNormError(m) = max(abs(error(ia:ib))); 
    l1NormError(m) = norm(error,1)*dx;
    
    fprintf('ms=%s: solution=%s: t=%6.3f: Nx=%3d Nt=%4d dt=%8.2e errMax=%8.2e errL1=%8.2e',...
            ms, solution, tfinal,Nx,Nt,dt,maxNormError(m),l1NormError(m));
    if( m==1 ) fprintf('\n'); else fprintf(' ratio=[%3.2f,%3.2f] (max,L1)\n',...
               maxNormError(m-1)/maxNormError(m),l1NormError(m-1)/l1NormError(m)); end
    
end % loop grid resolutions







% Apply the Boundary Conditions
function [U] = applyBCs(lgp,ia,Nx,U)  % solution: U, left bdry index: ia
    U(ia+Nx) = U(ia);   % RHS bdry point
    for n=1:lgp         % left ghost points
        U(ia-n) = U(ia-n+Nx);
    end
end

