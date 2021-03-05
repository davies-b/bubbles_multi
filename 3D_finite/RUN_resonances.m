% RUN_resonances.m
%
% Computes the resonant frequencies for an array
% of N spherical resonators in three dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B
%
% 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
cy = zeros(1,N);
cz = zeros(1,N);


% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
N_multi = 0;


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

x = linspace(1,30000,300);
init = [];
y = zeros(1, length(x));
for i = 1:length(x)
    y(i) = abs(f(x(i)));
end
for i = 2:length(x)-1
    if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e-8)
        init = [init, x(i)];
    end
end

% if length(init) < length(R)
%     disp('WARNING: fewer than N initial guesses created')
% %     pause
% end

init = sort(init);

%% Use Muller's method to compute the resonances

distTol = 5e-5; fTol = 1e-5; iterMax = 10;
resonances = [];
n = 1;
for initGuess = init
    
    z0 = initGuess;
    z1 = initGuess - 0.00001i;
    z2 = initGuess - 0.00002i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e-7
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi \n'], real(res), imag(res))
       resonances = [resonances res];
       n = n + 1;
    end
end

figure
scatter(real(resonances),imag(resonances),'xb')

