% Insurance Model! SIHR

%change
% make a change

% different branch


S_u_0 = 26e6  %26 million...this is the initial number of uninsured susceptible people (we will use US data here, for now)
S_i_0 = 300e6  % 300 million

I_u_0 = 1; I_i_0 = 0;
H_u_0 = 0; H_i_0 = 0;
R_u_0 = 0; R_i_0 = 0; 
D_u_0 = 0; D_i_0 = 0; % other initial conditions


t0 = 0;
tf = 200;
time_steps=100;
tee=linspace(t0,tf,time_steps);
[t,y] = ode45(@sihr, tee, [S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0]); 

myplot(t,y, t0,tf)

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

function aprime = sihr(t,y)

S_u = y(1); S_i = y(2); I_u = y(3); I_i = y(4); H_u = y(5); H_i = y(6); R_u = y(7); R_i = y(8); D_u = y(9); D_i = y(10);

I =  I_u+I_i;

S_u_0 = 26e6;  %26 million...this is the initial number of uninsured susceptible people (we will use US data here, for now)
S_i_0 = 300e6;  % 300 million


N = S_u_0 + S_i_0; % total population remains constant

beta = Beta(t); % contact rate. currently an arbitrarily chosen value.  when we include age structuring, we will abandon this in favor of a contact matrix

p_i = S_i_0/(S_i_0+S_u_0); % percentage of total population that is insured
p_u = S_u_0/(S_i_0+S_u_0);

h = 0.05; % this is the rate at which symptomatic infected people need ICU hospitilization...should be taken from literature or data analysis!
% h = p_u*c_u + p_i*c_i

c_i = 0.045; % we want to ensure that c_u > c_i so we start with a c_i that is slightly smaller than h

c_u = (h - p_i*c_i)/p_u % probability that uninsured symptomatic infected will need ICU hospitalization



d_u = 0.6; % probability that uninsured ICU patient will die % we will need to find a systematic way to determine this
d_i = 0.55; % make sure d_u > d_i

l = L(d_u*H_u, d_i*H_i, t);
g = G(d_u*H_u, d_i*H_i, t);

aprime = [-beta * S_u * I / N + l * S_i - g * S_u; % dS_u/dt
    -beta * S_i * I / N - l * S_i + g * S_u; % dS_i/dt
    beta * S_u * I / N - c_u * I_u - (1/14)*(1-c_u) * I_u; % dI_u/dt
    beta * S_i * I / N - c_i * I_i - (1/14)*(1-c_i) * I_i; % dI_i/dt
    c_u * I_u - d_u * H_u - (1 - d_u) * H_u; % dH_u/dt
    c_i * I_i - d_i * H_i - (1 - d_i) * H_i; % dH_i/dt
    (1 - c_u) * I_u + (1 - d_u) * H_u; % dR_u/dt
    (1 - c_i) * I_i + (1 - d_i) * H_i; % dR_i/dt
    d_u * H_u; % dD_u/dt
    d_i * H_i;]; % dD_i/dt 
end

function l = L(r1,r2,t)
% this function returns the rate at which insured susceptible -> uninsured
% susceptible

l = 0; % dummy value
end

function g = G(r1,r2,t)
% this function returns the rate at which uninsured susceptible -> insured
% susceptible


g = 0; % dummy value
end

function beta = Beta(t)
% this function returns the time-varying beta


beta = 0.25; % default value
end