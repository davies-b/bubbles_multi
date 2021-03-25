function Cd = capacitancedilute(cx,cy,cz,R)

N = length(R);
Cd = zeros(N);

for i = 1:N
    for j = 1:N
        if i == j
            Cd(i,j) = 4*pi*R(i);
        else
            d = sqrt((cx(i)-cx(j))^2+(cy(i)-cy(j))^2+(cz(i)-cz(j))^2);
            Cd(i,j) = -4*pi*R(i)*R(j)/d;
        end
    end
end