%
% Solve the heat equation with the theta-method
%
% u_t = kappa * u_xx -alpha*u a < x < b , 0 < t <= tFinal
% u(a,t) = ga(t), u(b,t)=gb(t)
% u(x,0) = u_0(x)
%
clear; clf; fontSize=14; lineWidth=2;

theta=.5; % BE: theta=1, FE: theta=0, CN theta=.5
theta=1;
theta=.25;
% method Name : (replace '.' by 'p' so we can use in file names)
methodName = regexprep(sprintf('Theta%3.2f',theta),'\.','p'); % theta0p5 or theta1p0 etc.

kappa=.05; % coefficient of diffusion
alpha=.1;
kx=4.*pi; % wave number in the IC and exact solution
a=0.; b=1.; % space interval interval
tFinal=.5; % final time

ga = @(t) 0 ; % BC RHS at x=a
gb = @(t) 0 ; % BC RHS at x=b
u0 = @(x) sin(kx*x); % initial condition function
% exact solution function:
uexact = @(x,t) sin(kx*x)*exp((-kappa*kx^2-alpha)*t);

fprintf('Solve the modified equation heat: theta=%4.2f\n',theta);
numResolutions=4;
for m=1:numResolutions
 Nx=10*2^m; % number of space intervals
 dx=(b-a)/Nx; % grid spacing
 x = linspace(a,b,Nx+1)'; % spatial grid

 % allocate space for the solution at two levels
 un = zeros(Nx+1,1); % holds U_i^n
 unp1 = zeros(Nx+1,1); % holds U_i^{n+1}

 dt=dx; % time step (adjusted below)
 Nt = round(tFinal/dt); % number of time-steps
 dt = tFinal/Nt; % adjust dt to reach tFinal exactly

 ia=1; % index of boundary point at x=a
 ib=Nx+1; % index of boundary point at x=b
 i1=ia+1; % first interior pt
 i2=ib-1; % last interior pt
 I = i1:i2;

 % Form the implicit matrix
 rhs=zeros(Nx,1);
 A=sparse(Nx+1,Nx+1);
 for i=i1:i2
  A(i,i-1)= -theta*dt*( kappa*( 1./dx^2) );
  A(i,i) = 1 -theta*dt*( kappa*(-2./dx^2) -alpha );
  A(i,i+1)= -theta*dt*( kappa*( 1./dx^2) );
 end;
 A(ia,ia)=1.; % BC at x=a
 A(ib,ib)=1.; % BC at x=b

 % Define D+D- operator
 DpDm = @(u,I) (u(I+1)-2.*u(I)+u(I-1))/(dx^2);

 t=0.;
 un = u0(x); % initial conditions
 % --- Start time-stepping loop ---
 for( n=1:Nt )

 t = n*dt; % new time
 % Theta-method
 rhs(I) = un(I) + ((1.-theta)*dt)*( kappa*DpDm(un,I) - alpha*un(i1:i2) );
 rhs(ia)=ga(t); % BC at x=a
 rhs(ib)=gb(t); % BC at x=b

 unp1 = A\rhs;

 un=unp1; % Set un <- unp1 for next step
 end;
 % --- End time-stepping loop ---

 ue = uexact(x,t); % eval exact solution:
 err = un-ue;
 errMax(m) = max(abs(err)); % max-norm error
 fprintf('theta=%3.1f: t=%10.4e: Nx=%3d Nt=%4d dt=%9.3e maxErr=%8.2e',theta,t,Nx,Nt,dt,errMax(m));
 if( m==1 ) fprintf('\n'); else fprintf(' ratio=%8.2e\n',errMax(m-1)/errMax(m)); end;

 % plot results
 figure(1);
 plot( x,un,'b-o', x,ue,'k-','Linewidth',lineWidth);
 legend('computed','true');
 title(sprintf('Mod. Heat: \\theta=%4.2f t=%8.2e Nx=%d dt=%7.1e',theta,t,Nx,dt));
 xlabel('x'); ylabel('u'); set(gca,'FontSize',fontSize); grid on;

 if( m==1 )
 print('-depsc2',sprintf('heatEquation%s.eps',methodName)); % save as an eps file
 end;

 % plot error
 figure(2);
 plot(x,err,'b-x','Linewidth',lineWidth);
 legend('Error');
 title(sprintf('Mod. Heat: Error: \\theta=%3.1f t=%8.2e Nx=%d dt=%7.1e',theta,t,Nx,dt));
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
 title(sprintf('Mod. Heat. Errors \\theta=%3.1f t=%9.3e',theta,t));
 xlabel('x'); ylabel('u'); set(gca,'FontSize',fontSize); grid on;

 fileName = sprintf('heatEquation%sErrors.eps',methodName);
 print('-depsc2',fileName); % save as an eps file
 fprintf('Wrote file=[%s]\n',fileName);