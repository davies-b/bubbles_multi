clear, close all

L = 35e-3;              % length of cochlea
N = 10;                 % number of resonators
L_end = L - 0.005;

s = 1.05;
Rad = L_end*(1-s)/(1-s^N)/3;
for i = 1:N
    R(i) = Rad*s^(i-1);
end

%%% Material parameters
rho0 = 1e3;             % density of water
kappa0 = 2e9;           % bulk modulus of water
v = sqrt(kappa0/rho0);  % speed of sound in water

rho_b = 1.2;            % density of resonators  
kappa_b = 1e5;          % bulk modulus of resonators
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air


% High contrast parameters \delta
delta=rho_b/rho0;

% Define positions of resonators
cx = R(1)*ones(1,N);
if N > 1
    for i = 2:N
        cx(i) = cx(i-1) + 2*R(i-1) + R(i);
    end
end

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
N_multi = 3;

range = [1,30000];
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
