%**********************
% 2022 rewrite
% Density of PS: 1050 kg/m3
% 1 Dalton = 1.66e-27 kg
close all

k = 2.3e-4;     % Young's Modulus of PEG (N/m)
q_total = -6*1.6e-19;  % charge of molecule (e)
Hdiv = 10;      % Differential number
E = 7e8;        % electric field intensity at surface (V/m)
d = 10e-9;      % size of molecule (m)
w = 2*pi*40;    % augular frequency (rad/s)
% m = 5.498e-19;   % molecule mass (kg) of 100 nm PSNP
m = 2.49e-22;   % molecule mass (kg) of 150 kDa
c = 8.9e-4;     % viscosity of water (N s/ m2)
De = 7.89e-9;   % Debye length (m)

H = linspace(0, d, Hdiv+1);
H = H(1:end-1);
qi = linspace(q_total/Hdiv, q_total/Hdiv, Hdiv);

%%%%% ODE solver
qESum = @(x) 0;
for i = 1 : Hdiv
    qESum = @(x) qESum(x) + qi(i) * E*exp(-(x+H(i))/De); 
end

dz2 = @(t, z)[z(2); (qESum(z(1))*(cos(w*t))*1 - k*z(1) - 6*pi*c*d*z(2))/m];
z10 = 0e-9;
z20 = 0;
step = 0.001;
flag = 0;
lastwarn('', '');
figure(1)
for i = 0:step:1-step
    if isempty(lastwarn)
        if i == 0 
            [t_sol, z_sol] = ode15s(dz2, [i, i+step], [z10, z20]);
        else
            if flag == 0
                [t_sol, z_sol] = ode15s(dz2, [i, i+step], [z_sol(end, 1), z_sol(end, 2)]);
            else
                flag = 0;
                [t_sol, z_sol] = ode15s(dz2, [i, i+step], [0, 0]);
            end
        end  
    else
        lastwarn('', '');
        for j = 1:length(z_sol)
            if z_sol(j, 1) < 0
                z_sol(j, 1) = 0;
            end
        end
        flag = 1;
        continue
    end
    plot(t_sol, z_sol(:, 1), 'b')
    hold on
end
ylim([0 6e-8])