function C = capacitance(cx,R,rho0,rho_b,kappa0,kappa_b,delta)

N = length(cx);
cy = zeros(1,N);
cz = zeros(1,N);

N_multi = 0;

z = 0.000001;
A = MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);

B = zeros(N);
for i = 1:N
    B(i,i) = A(2*i-1,2*i-1);
end
for i = 1:N
    for j = 1:N
        if i ~= j
            B(i,j) = -A(2*i-1,2*j);
        end
    end
end

B = -B/4/pi/R(1)^2;

C = B^-1;