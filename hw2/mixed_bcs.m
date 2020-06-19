%
% Solve the heat equation with mixed BC's and manufactured solutions
%
% u_t = kappa * u_xx + f(x,t), a < x < b , 0 < t <= tFinal
% alphaA*u(a,t) + betaA*ux(a,t) = ga(t),
% alphaB*u(b,t) + betaB*ux(b,t) = gb(t)
% u(x,0) = u_0(x)
%
clear; clf; fontSize=14; lineWidth=2;

kappa=.1; % coefficient of diffusion
a=0.; b=1.; % space interval interval
tFinal=1.; % final time
cfl=1.;

alphaA=.5; betaA=-1.; % coefficients of the mixed BC
alphaB=.25; betaB=+2.;

% Time-stepping
ts = 'FE'; % forward Euler
ts = 'AB2'; % Adams-Bashforth order 2

% Choose the manufactured solution:
ms = 'poly';
ms = 'trig'; % un-comment for trig MS

% Choose BC option
bcOption = 'centered'; bcOrder=2;
bcOption = 'oneSided'; bcOrder=2;
bcOption = 'oneSided'; bcOrder=1;

if( strcmp(ms,'poly') )
% Polynomial manufactured solution:
 % b0=.5; b1=.7; b2=.0; % time coefficients
 b0=.5; b1=.7; b2=.9; % time coefficients

 c0=1.; c1=.8; c2=.6; % space coefficients
 % c0=1.; c1=.8; c2=.0; % space coefficients
 ue = @(x,t) (b0+b1*t+b2*t^2)*(c0+c1*x+c2*x.^2);
 uet = @(x,t) (b1+2.*b2*t )*(c0+c1*x+c2*x.^2);
 uex = @(x,t) (b0+b1*t+b2*t^2)*(c1+2.*c2*x);
 uexx = @(x,t) (b0+b1*t+b2*t^2)*(2.*c2);
 else
 % Trigonometric manufactured solution:
 kt = pi; kx=7*pi;
 ue = @(x,t) cos(kt*t)*cos(kx*x);
 uet = @(x,t) -kt*sin(kt*t)*cos(kx*x);
 uex = @(x,t) -kx*cos(kt*t)*sin(kx*x);
 uexx = @(x,t) (-kx^2)*ue(x,t);
 end

 ga = @(t) alphaA*ue(a,t) + betaA*uex(a,t); % BC RHS at x=a
 gb = @(t) alphaB*ue(b,t) + betaB*uex(b,t); % BC RHS at x=b
 u0 = @(x) ue(x,0); % initial condition function
 f = @(x,t) uet(x,t) - kappa*uexx(x,t); % forcing


 numResolutions=5;
 for m=1:numResolutions
 Nx=10*2^m; % number of space intervals
 dx=(b-a)/Nx; % grid spacing
 numGhost=1;
 Nd = Nx + 1 + 2*numGhost; % total number of points
 x = linspace(a-numGhost*dx,b+numGhost*dx,Nd)'; % spatial grid
 if( abs(x(2)-x(1) -dx) > 1.e-12*dx ) fprintf('ERROR in dx!\n'); pause; pause; end;

 % allocate space for the solution at two levels
 um = zeros(Nd,1); % holds U_i^{n-1}
 un = zeros(Nd,1); % holds U_i^n
 unp1 = zeros(Nd,1); % holds U_i^{n+1}

 dt = cfl*.25*dx^2/kappa; % time step (adjusted below)
 Nt = round(tFinal/dt); % number of time-steps
 % fprintf('dt=%9.3e, adjusted=%9.3e, diff=%8.2e\n',dt,tFinal/Nt,abs(dt-tFinal/Nt));
 dt = tFinal/Nt; % adjust dt to reach tFinal exactly

 ia=numGhost+1; % index of boundary point at x=a
 ib=ia + Nx; % index of boundary point at x=b
 i1=ia+1; % first interior pt
 i2=ib-1; % last interior pt
 I = i1:i2; % interior points
 Ib = ia:ib; % interior + boundary points

 % Define D+D- operator
 DpDm = @(u,I) (u(I+1)-2.*u(I)+u(I-1))/(dx^2);

 t=0.;
 un = u0(x); % initial conditions
 um = ue(x,t-dt); % set old time to exact
 % --- Start time-stepping loop ---
 for( n=1:Nt )

 tm = (n-1)*dt; % old time
 t = n*dt; % new time

 if( strcmp(ts,'FE') )
 % Forward Euler in time:
 unp1(Ib) = un(Ib) + (kappa*dt)*DpDm(un,Ib) + dt*f(x(Ib),tm);
 else
 % AB2
 unp1(Ib) = un(Ib) + dt*( 1.5*( kappa*DpDm(un,Ib) + f(x(Ib),t-dt ) ) ...
 -.5*( kappa*DpDm(um,Ib) + f(x(Ib),t-2*dt) ) );
 end

 if( strcmp(bcOption,'centered') )
 % Centered BC: assign ghost:
 unp1(ia-1)=unp1(ia+1) + (2.*dx/betaA)*( alphaA*unp1(ia) - ga(t) ); % set left ghost from BC at x=a
 unp1(ib+1)=unp1(ib-1) - (2.*dx/betaB)*( alphaB*unp1(ib) - gb(t) ); % set right ghost from BC at x=b
 else
% one-sided BC
 if( bcOrder==1 )
 % first-order one-sided
 unp1(ia) = ( ga(t) - betaA*( unp1(ia+1) )/dx )/( alphaA - betaA/dx );
 unp1(ib) = ( gb(t) - betaB*( -unp1(ib-1) )/dx )/( alphaB + betaB/dx );
 else
 % second-order one-sided
 unp1(ia) = ( ga(t) - betaA*( 4*unp1(ia+1) - unp1(ia+2) )/(2.*dx) )/( alphaA - 3.*betaA/(2.* dx) );
 unp1(ib) = ( gb(t) - betaB*( -4*unp1(ib-1) + unp1(ib-2) )/(2.*dx) )/( alphaB + 3.*betaB/(2.* dx) );
 end
 % extrapolate ghost
 unp1(ia-1)=3*unp1(ia)-3*unp1(ia+1)+unp1(ia+2);
 unp1(ib+1)=3*unp1(ib)-3*unp1(ib-1)+unp1(ib-2);
 end;


 um=un; % set um <- un for next step
 un=unp1; % Set un <- unp1 for next step

 end;
 % --- End time-stepping loop ---

 uexact = ue(x,t); % eval exact solution:
 err = un-uexact; % error
 errMax(m) = max(abs(err)); % max-norm error
 fprintf('Heat: %s, %s%d, %s: t=%4.2f: Nx=%3d Nt=%5d dt=%9.3e maxErr=%8.2e',ts,bcOption,bcOrder,ms,t,Nx,Nt,dt,errMax(m));
 if( m==1 ) fprintf('\n'); else fprintf(' ratio=%8.2e\n',errMax(m-1)/errMax(m)); end;

 % plot results
 figure(1);
 plot( x,un,'b-o', x,uexact,'k-','Linewidth',lineWidth);
 legend('computed','true');
 title(sprintf('Heat-Eqn MS=%s t=%9.2e (Nx=%d) dt=%8.2e, bc=%s',ms,t,Nx,dt,bcOption));
 xlabel('x'); ylabel('u'); set(gca,'FontSize',fontSize); grid on;
 if( m==1 )
 print('-depsc2',sprintf('heatMixedMS%s%s%d.eps',ms,bcOption,bcOrder)); % save as an eps file
 end;

 % plot error
 figure(2);

 plot(x,err,'b-x','Linewidth',lineWidth);
 legend('Error');
 title(sprintf('Heat Equation Error t=%9.3e bc=%s',t,bcOption));
 xlabel('x'); ylabel('u'); set(gca,'FontSize',fontSize); grid on;

 % save grid and errors in a data structure
 data{m}.x=x; data{m}.err=err;

 pause(1);

 end; % end for m


 % plot errors
 figure(2);
 % plot(x1,err1,'r-x', x2,err2,'g-o',x3,err3,'b-+','Linewidth',lineWidth);
 plot(data{1}.x,data{1}.err,'r-x', data{2}.x,data{2}.err,'g-o',data{3}.x,data{3}.err,'b-+','Linewidth',lineWidth);
 legend('error Nx=20','error Nx=40','error Nx=80','Location','best');
 title(sprintf('Heat Equation Errors MS=%s t=%9.3e',ms,t));
 xlabel('x'); ylabel('u'); set(gca,'FontSize',fontSize); grid on;
 fileName=sprintf('heatMixedMS%s%s%dErrors.eps',ms,bcOption,bcOrder);
 print('-depsc2',fileName); % save as an eps file
 fprintf('Wrote file=[%s] with errors\n',fileName);