% Insurance Model! SIHR

clear all
close all


contact_matrix = [25 10 5 10 5 5 5 5 5;
    10 40 10 5 15 10 5 5 5;
    5 10 20 15 10 10 5 5 5;
    10 5 15 20 10 10 5 5 5;
    5 15 10 10 10 10 5 5 5;
    5 10 10 10 10 5 5 5 5;
    5 5 5 5 5 5 5 10 5;
    5 5 5 5 5 5 10 10 10;
    5 5 5 5 5 5 5 10 10;]; % the average number of daily contacts between age groups % let's have 9 age groups as in https://arxiv.org/pdf/2004.04144.pdf

issymmetric(contact_matrix) % make sure it's symmetric...just checking

contact_matrix = contact_matrix./sum(contact_matrix,1) % we want to make it so that sum of each row and column is 1...how to do this??...close enough for now

S_u_0_total = 26e6;  %26 million...this is the initial number of uninsured susceptible people (we will use US data here, for now)
S_i_0_total = 300e6;  % 300 million

uninsured_fraction_by_age = [.05 .09 .16 .21 .18 .16 .13 .01 .01];  %sum should be 1

insured_fraction_by_age = 1-uninsured_fraction_by_age

S_u_0 = S_u_0_total * uninsured_fraction_by_age

sum(S_u_0)

S_i_0 = S_i_0_total * insured_fraction_by_age
sum(S_i_0)

M_S_u_0 = diag(S_u_0);

M_S_i_0 = diag(S_i_0);



N = sum(S_u_0) + sum(S_i_0) % total population remains constant


I_u_0 = 20; I_i_0 = 20;
H_u_0 = 0; H_i_0 = 0;
R_u_0 = 0; R_i_0 = 0; 
D_u_0 = 0; D_i_0 = 0; % other initial conditions


t0 = 0;
tf = 200;
time_steps=100;
tee=linspace(t0,tf,time_steps);

p_i = S_i_0/(S_i_0+S_u_0); % percentage of total population that is insured
p_u = S_u_0/(S_i_0+S_u_0);

h = 0.05; % this the percentage of symptomatic infected people who need ICU hospitilization...should be taken from literature or data analysis!
% h = p_u*c_u + p_i*c_i

c_i = 0.045; % we want to ensure that c_u > c_i so we start with a c_i that is slightly smaller than h

c_u = (h - p_i*c_i)/p_u % probability that uninsured symptomatic infected will need ICU hospitalization


d_u = 0.6; % probability that uninsured ICU patient will die % we will need to find a systematic way to determine this
d_i = 0.55; % make sure d_u > d_i

alpha_u = 1/14; % rate at which ICU Hosptializations rocever
alpha_i = 1/14; 
delta_u = 1/14; % rate at which infected recover (without need for ICU hospitalization)
delta_i = 1/14;
gamma_u = 1/5;  % rate which infected go to ICU
gamma_i = 1/5;
ksi_u = 1/3;    % death rate from ICU
ksi_i = 1/3;





[t,y] = ode45(@(t,y) sihr(t, y, N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i), tee, [S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0]); 

figure; plot(t,y)
legend %labels lines

%myplot(t, y, t0, tf)

function myplot(t,y,t0,tf)

% below we plot the results
color = get(gca,'colororder'); % different colors for plotting

figure(1)
subplot(5,2,1)
plot(t,y(:,1),'-o','Color',color(1,:))
hold on
title('Susceptible (Uninsured)')
xlim([t0 tf])

subplot(5,2,2)
plot(t,y(:,2),'-o','Color',color(2,:))
hold on
title('Susceptible (Insured)')
xlim([t0 tf])

subplot(5,2,3)
plot(t,y(:,3),'-o','Color',color(3,:))
hold on
title('Infected (Uninsured)')
xlim([t0 tf])

subplot(5,2,4)
plot(t,y(:,4),'-o','Color',color(4,:))
hold on
title('Infected (Insured)')
xlim([t0 tf])

subplot(5,2,5)
plot(t,y(:,5),'-o','Color',color(5,:))
hold on
title('ICU Hospitalized (Uninsured)')
xlim([t0 tf])

subplot(5,2,6)
plot(t,y(:,6),'-o','Color',color(6,:))
hold on
title('ICU Hospitalized (Insured)')
xlim([t0 tf])

subplot(5,2,7)
plot(t,y(:,7),'-o','Color',color(5,:))
hold on
title('Recovered (Uninsured)')
xlim([t0 tf])

subplot(5,2,8)
plot(t,y(:,8),'-o','Color',color(6,:))
hold on
title('Recovered (Insured)')
xlim([t0 tf])

 subplot(5,2,9)
plot(t,y(:,9),'-o','Color',color(5,:))
hold on
title('Dead (Uninsured)')
xlim([t0 tf])

subplot(5,2,10)
plot(t,y(:,10),'-o','Color',color(6,:))
hold on
title('Dead (Insured)')
xlim([t0 tf])
end

function aprime = sihr(t, y, N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i)

S_u = y(1); 
S_i = y(2); 
I_u = y(3); 
I_i = y(4); 
H_u = y(5); 
H_i = y(6); 
R_u = y(7); 
R_i = y(8); 
D_u = y(9); 
D_i = y(10);

I =  I_u+I_i;

beta = Beta(t);  % contact rate. currently an arbitrarily chosen value.  when we include age structuring, we will abandon this in favor of a contact matrix



aprime = [-beta * S_u * I / N ; % dS_u/dt
    -beta * S_i * I / N ; % dS_i/dt
    beta * S_u * I / N - (gamma_u * c_u * I_u) - delta_u * (1-c_u) * I_u; % dI_u/dt
    beta * S_i * I / N - (gamma_i * c_i * I_i) - delta_i * (1-c_i) * I_i; % dI_i/dt
    gamma_u * c_u * I_u - (ksi_u * d_u * H_u) - alpha_u * (1 - d_u) * H_u; % dH_u/dt
    gamma_i * c_i * I_i - (ksi_i * d_i * H_i) - alpha_i * (1 - d_i) * H_i; % dH_i/dt
    delta_u * (1-c_u) * I_u + alpha_u * (1 - d_u) * H_u; % dR_u/dt
    delta_i * (1-c_i) * I_i + alpha_i * (1 - d_i) * H_i; % dR_i/dt
    ksi_u * d_u * H_u; % dD_u/dt
    ksi_i * d_i * H_i;]; % dD_i/dt 
end



function beta = Beta(t)
% this function returns the time-varying beta


beta = 0.25; % default value
end