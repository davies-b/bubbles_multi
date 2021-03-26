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
cx = 3*linspace(1,N,N);

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
N_multi = 2;

range = [0.01,0.04];
multires = multipoleres(cx,R,N_multi,rho0,rho_b,kappa0,kappa_b,delta,range);

%% Compute the resonances using the capacitance matrix, which is approximated using the multipole expansion method
C = capacitance(cx,R,rho0,rho_b,kappa0,kappa_b,delta);
Cvol = 3/4/pi*diag(R.^-3)*C;
capres = sqrt(delta*v_b^2*eig(real(Cvol)));

%% Compute the resonances using the dilute capacitance matrix
Cd = capacitancedilute(cx,cy,cz,R);
Cdvol = 3/4/pi*diag(R.^-3)*Cd;
diluteres = sqrt(delta*v_b^2*eig(real(Cdvol)));

if norm(imag(diluteres)) > 0
    disp('WARNING: Dilute matrix giving a poor approximation, the resonators are too large')
end

%% Plot the resonances
figure
scatter(real(multires),1.1*ones(1,N),'xb')
hold on
scatter(capres,ones(1,N),'ob')
scatter(diluteres,0.9*ones(1,N),'sb')
ylim([0.85, 1.2]); set(gca,'ytick',[]);
legend('full multipole method','capacitance matrix','dilute capacitance matrix')
