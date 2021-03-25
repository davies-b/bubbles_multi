function resonances = multipoleres(cx,R,N_multi,rho0,rho_b,kappa0,kappa_b,delta)

N = length(cx);
cy = zeros(1,N);
cz = zeros(1,N);

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

if length(init) < length(R)
    disp('WARNING: fewer than N initial guesses created by multipole method')
end

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