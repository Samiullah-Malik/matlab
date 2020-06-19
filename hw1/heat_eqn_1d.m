%% Q1. 1D Heat Equation - IBVP 
clear all; clc;

for m = 1:3
k = 2;
kappa = 0.1;
tfinal = 1.0;
a = 0; b = 1; % domain

Nx = 10 * 2^m; 
dx = (b - a) / Nx;
x = linspace(a,b,Nx+1)';      % grid points

uexact = @(x,t) sin(k * pi * x) * exp(-kappa * (k * pi)^2 * t);  
ga = @(t) uexact(a,t);       % ga(t)   - BC
gb = @(t) uexact(b,t);       % gb(t)   - BC
uphi = @(x) uexact(x, 0);    % uphi(x) - IC

% Allocate space
un = zeros(Nx+1, 1);         % U_i^n
unp1 = zeros(Nx+1, 1);       % U_i^(n+1)
cfl = 0.25;                  % cfl parameter
dt = cfl * (dx^2) / kappa; 
Nt = round(tfinal / dt);     % Number of timesteps
dt = tfinal / Nt;            % adjust dt to reach tfinal


%  |-----------|-------------|----------|------------|  
% ia                                                 ib
ia = 1;        % index of bdry point @x = a
ib = Nx + 1;   % index of bdry point @x = b
i1 = ia + 1;   % first interior point
i2 = ib - 1;   % last interior point
I = i1:i2;     % I = i1, i1+1, i1+2, ... , i2

t = 0;
un = uphi(x);    % u_i^n = u

% start time stepping loop
for(n = 1: Nt)
    unp1(I) = un(I) + (kappa * dt / (dx ^ 2)) * (un(I+1) - 2*un(I) + un(I-1));
    
    t = n * dt;
    
    unp1(ia) = ga(t);      % BC @ x = a
    unp1(ib) = gb(t);      % BC @ x = b
    un = unp1;             % set un <- unp1 for next step
      
end

% compute errors
ue = uexact(x,t);  % exact solution
errMax(m) = max(abs(un-ue));   % max norm error

% [Q1.(c)] Print Results 
fprintf('m = %1d t=%10.4e: Nx=%3d Nt=%4d dt=%9.3e maxErr=%8.2e', m, t,Nx,Nt,dt,errMax(m));

if(m == 1) fprintf('\n');
else fprintf(' ratio=%8.2e, rate=%4.2f\n',errMax(m-1)/errMax(m),...
log2(errMax(m-1)/errMax(m)));
end


% [Q1.(d)]
figure(1)
if (m == 1)   
    % Plot exact vs computed solution
    plot(x, un)  % Plot computed solution
    hold on
    plot(x, ue)  % Plot exact solution
    legend('Computed Solution', 'Exact Solution')
    xlabel('x')
    ylabel('u(x,t=1)')
    title('Exact vs Computed Solution')
    hold off   
end

% All errors on one graph
figure(2)
plot(x, un-ue)
%legend('m=1','','')
hold on
%

end

figure(2)
legend('m = 1','m = 2','m = 3')
title('Errors')
xlabel('x')
ylabel('error')

