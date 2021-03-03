% RUN_resonances.m
%
% Computes the resonant frequencies for an array
% of N resonators in two dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B
%
% 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

%% Define parameters

L = 35e-3;              % length of cochlea
N = 22;                 % number of resonators
L_end = L;

s = 1.05;
a = 0.0001;

for i = 1:N
    R(i) = a*s^(i-1);
end

%%% Material parameters
rho0 = 1e3;                 % density of water
kappa0 = 2e9;               % bulk modulus of water
v = sqrt(kappa0/rho0);      % speed of sound in water

rho_b = 1.2;                % density of resonators/air
kappa_b = 1e5;              % bulk modulus of resonators/air
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air

% High contrast parameter \delta
delta=rho_b/rho0;

% Define the position of the centres of the resonators
cx = linspace(0,L,N+2);
cx = cx(2:end-1);
cy = zeros(1,N);

% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
% If we use higher order, then accuracy improves. Usually 3 is sufficient.
N_multi = 3;

%%% Plot the geometry
figure, hold on
t = linspace(0,2*pi);
for n = 1:length(R)
    plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
end
daspect([1 1 1])
hold off


%% Compute initial guesses for the resonances
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) min(eig(MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy)));

x = linspace(1, 2*pi*22000, 200);
init = [];
y = zeros(1, length(x));
for i = 1:length(x)
    y(i) = abs(f(x(i)));
end
for i = 2:length(x)-1
    if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e0)
        init = [init x(i)];
    end
end

if length(init) < length(R)
    disp('WARNING: fewer than N initial guesses created')
end

init = sort(init);

%% Use Muller's method to find the resonances

distTol = 5e-5; fTol = 1e-5; iterMax = 10;
resonances = [];
n = 1;
for initGuess = init
        
    z0 = initGuess;
    z1 = initGuess - 1i;
    z2 = initGuess - 2i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e0
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi (%.0f Hz) \n'], real(res), imag(res), real(res)/2/pi)
       resonances = [resonances res];
       n = n + 1;
    end
end

figure
scatter(real(resonances), imag(resonances), 'x')
xlabel('$\Re\omega_n$ (Hz)','interpreter','latex')
ylabel('$\Im\omega_n$','interpreter','latex')