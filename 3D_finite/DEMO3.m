clear, close all

N = 10;                 % number of resonators

%%% Material parameters
high = 5000;

rho0 = high;             % density of background
kappa0 = high;           % bulk modulus of background
v = sqrt(kappa0/rho0);  % speed of sound in background

rho_b = 1;            % density of resonators  
kappa_b = 1;          % bulk modulus of resonators
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air

% High contrast parameters \delta
delta=rho_b/rho0;

% Define size of resonators
R = ones(1,N);

% Define positions of resonators
cx = 5*linspace(1,N,N);

%%% Plot the geometry
cy = zeros(1,N); cz = zeros(1,N);
figure, hold on
t = linspace(0,2*pi);
for n = 1:length(R)
    plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
end
daspect([1 1 1])
hold off

%% Compute the resonances using the multipole expansion method and Muller's method
% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
disp('Full multipole method')
tic
N_multi = 2;

range = [0.01,0.04];
multires = multipoleres(cx,R,N_multi,rho0,rho_b,kappa0,kappa_b,delta,range);
toc

%% Compute the resonances using the capacitance matrix, which is approximated using the multipole expansion method
disp('Capacitance matrix approximation')
tic
C = real(capacitance(cx,R,rho0,rho_b,kappa0,kappa_b,delta));
Cvol = 3/4/pi*diag(R.^-3)*C;

[V,D] = eig(Cvol);
Vol = 4/3*pi*diag(R.^3);

capres = zeros(1,N);
for n = 1:N
    denom = V(:,n)'*Vol*V(:,n);
    tau = v_b^2/v/8/pi*V(:,n)'*Vol*ones(N)*Vol*V(:,n)/denom;
    capres(n) = sqrt(delta*v_b^2*D(n,n)) - 1i*delta*tau*D(n,n)^2;
end
toc

%% Plot the resonances
% Real parts only
figure
scatter(real(multires),1.1*ones(1,N),'xb')
hold on
scatter(real(capres),ones(1,N),'sb')
ylim([0.95, 1.2]); set(gca,'ytick',[]); h = gca; h.YAxis.Visible = 'off';
legend('full multipole method','capacitance matrix')

% Complex plane
figure
sz = 100;
scatter(real(multires),imag(multires),sz,'xb')
hold on
scatter(real(capres),imag(capres),sz,'sb')
legend('full multipole method','capacitance matrix','Location','SouthEast','interpreter','latex')
set(gca,'FontSize',17,'TickLabelInterpreter','latex','XAxisLocation','top')
set(gcf, 'Position', [0, 0, 750, 300])
xlabel('Real part','interpreter','latex')
ylabel('Imaginary part','interpreter','latex')