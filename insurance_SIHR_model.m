% Insurance Model! SIHR

clear all
close all
clc

% -----------------------initial population------------------------------

S_u_0 = 27.5e6;  %26 million...this is the initial number of uninsured susceptible people (we will use US data here, for now)
S_i_0 = 296e6;  % 300 million

I_u_0 = 1; I_i_0 = 1;
H_u_0 = 0; H_i_0 = 0;
R_u_0 = 0; R_i_0 = 0; 
D_u_0 = 0; D_i_0 = 0; % other initial conditions

N = S_u_0 + S_i_0; % total population remains constant


%---------setting and computing parameters ----------------------------%
beta = 0.5; % contact rate % taken from literature...mask or not to mask paper here

global beta_values;
beta_values = [];

p_i = S_i_0/(S_i_0+S_u_0); % percentage of total population that is insured
p_u = S_u_0/(S_i_0+S_u_0);

h = 0.05; % this the percentage of symptomatic infected people who need ICU hospitilization...should be taken from literature or data analysis!
% h = p_u*c_u + p_i*c_i

c_i = 0.045; % we want to ensure that c_u > c_i so we start with a c_i that is slightly smaller than h

c_u = (h - p_i*c_i)/p_u % probability that uninsured symptomatic infected will need ICU hospitalization

k = 0.23 % probability that ICU patient will die.  taken from literature.  k = p_u*d_u + p_i * d_i

d_i = 0.225; % set d_i slightly lower than k
d_u = (k - p_i*d_i)/p_u;


alpha_u = 1/14; % rate at which ICU Hosptializations rocever
alpha_i = 1/14; 
delta_u = 1/14; % rate at which infected recover (without need for ICU hospitalization)
delta_i = 1/14;
gamma_u = 1/5;  % rate which infected go to ICU
gamma_i = 1/5;
ksi_u = 1/3;    % death rate from ICU
ksi_i = 1/3;
eta = 1/30;  % rate at which insured susceptible <-> uninsured susceptible  (30 days <- we need to figure out if this is viable number!)


global loss_values;
loss_values = [];


global gain_values;
gain_values = [];

%----------gain of coverage--------
t_start_coverage  = 150; % day on which universal coverage is passed  % TINKER WITH ME

fraction_each_time_step_that_gains_coverage = 1/20; % TINKER WITH ME! % for step function

delta_period = 10;  % relevant if using periodic delta gain function


coverage_implementation_type = 1; % 1 is step function at time t, 2 is periodic delta





%----------turn off or on features-------------
universal_coverage_feature = 1; % 0 is off

unemployment_feature = 1; % on or off, 0 is off ...this helps us simulate baseline resutlts

time_varying_beta = 1; % 0 is off



%----------Let's solve this thing!-----------------
t0 = 1;
tf = 600; % unit = days
time_steps = 600;
tee=linspace(t0,tf,time_steps);

[t,y] = ode45(@(t,y) sihr(t, y, N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i, eta, unemployment_feature, time_varying_beta, beta, universal_coverage_feature, t_start_coverage, coverage_implementation_type, fraction_each_time_step_that_gains_coverage, delta_period), tee, [S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0]); 


%-------------Let's plot the results-----------------
S_u = y(:,1); 
S_i = y(:,2); 
I_u = y(:,3); 
I_i = y(:,4); 
H_u = y(:,5); 
H_i = y(:,6); 
R_u = y(:,7); 
R_i = y(:,8); 
D_u = y(:,9); 
D_i = y(:,10);

sprintf("peak infections %d", getPeakInfections(I_u,I_i))
sprintf("peak ICU hosp %d", getPeakICUHospitalizations(H_u,H_i))
sprintf("total deaths %d", getTotalDeaths(D_u, D_i))


%compareCoverageStartToSpeed(N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i, eta, unemployment_feature, time_varying_beta, beta, universal_coverage_feature, coverage_implementation_type, tee, S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0)




plotCompartmentsSeparately(t, y, t0, tf)

plotLossOfCoverage();

plotGainOfCoverage();

plotBeta();

plotAsFractions(t,y, t0, tf);


%----------function declarations/definitions below this line------------%




function dydt = sihr(t, y, N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i, eta, unemployment_feature, time_varying_beta_on, default_beta, universal_coverage_feature, t_start_coverage, coverage_implementation_type, fraction_each_time_step_that_gains_coverage, delta_period)

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


I = I_u + I_i;

beta = Beta(default_beta, time_varying_beta_on, H_i+H_u);  % contact rate. currently an arbitrarily chosen value.  when we include age structuring, we will abandon this in favor of a contact matrix

l = 0; 
g = 0;

u = unemployment_val(H_i+H_u, unemployment_feature);

if (t > t_start_coverage) & (universal_coverage_feature ~= 0)
    g = eta*gain_func(round(t), coverage_implementation_type, fraction_each_time_step_that_gains_coverage, delta_period);
elseif unemployment_feature ~= 0  % we don't allow unemployment data once universal coverage is in play
    l = eta*u;
end


dydt = [-beta * S_u * I / N + l * S_i - g * S_u; % dS_u/dt
    -beta * S_i * I / N - l * S_i + g * S_u; % dS_i/dt
    beta * S_u * I / N - (gamma_u * c_u * I_u) - delta_u * (1-c_u) * I_u; % dI_u/dt
    beta * S_i * I / N - (gamma_i * c_i * I_i) - delta_i * (1-c_i) * I_i; % dI_i/dt
    gamma_u * c_u * I_u - (ksi_u * d_u * H_u) - alpha_u * (1 - d_u) * H_u; % dH_u/dt
    gamma_i * c_i * I_i - (ksi_i * d_i * H_i) - alpha_i * (1 - d_i) * H_i; % dH_i/dt
    delta_u * (1-c_u) * I_u + alpha_u * (1 - d_u) * H_u; % dR_u/dt
    delta_i * (1-c_i) * I_i + alpha_i * (1 - d_i) * H_i; % dR_i/dt
    ksi_u * d_u * H_u; % dD_u/dt
    ksi_i * d_i * H_i;]; % dD_i/dt 
end


function loss = unemployment_val(H, on)
alpha=1;
k = 1e5; % threshold
extreme_val = 0.005; % this amount to approx the maximum percent of people that we have seen unemployed in a single day, on average (taken from monthly unemployment data in 2020)

if on == 1
    loss = extreme_val*(alpha*H/(H+k));
else
    loss = 0;
end

global loss_values;
loss_values = [loss_values loss];

end


function y = gain_func(x, implementation_type, fraction_each_time_step_that_gains_coverage, delta_period)

if implementation_type == 1
    y = fraction_each_time_step_that_gains_coverage; % we don't need to use heaviside func since we already check t > t_start in sihr func * heaviside(x-t_start_coverage);
elseif implementation_type == 2
    y = dirac(mod(x, delta_period));  % NEED TO MAKE SURE THIS IS WORKING CORRECTLY!
end

global gain_values;
gain_values = [gain_values y];

end

function beta = Beta(default_val,on, H)
% this function returns the time-varying beta
alpha = 1;% confidence
k = 1e5; % threshold

if on == 1
    beta = default_val*(1-alpha*H/(H+k));
else
    beta = default_val;
end

global beta_values;
beta_values = [beta_values beta];

end



%-----------below this line are plotting and measure finding functions-----------

function val = getPeakInfections(uninsured_vec,insured_vec)

combined_vec = uninsured_vec + insured_vec;
val = max(combined_vec);
 
end

function val = getPeakICUHospitalizations(uninsured_vec,insured_vec)
combined_vec = uninsured_vec + insured_vec;
val = max(combined_vec);
    
end

function val = getTotalDeaths(uninsured_vec,insured_vec)
combined_vec = uninsured_vec + insured_vec;
val = max(combined_vec);
    
end

function compareCoverageStartToSpeed(N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i, eta, unemployment_feature, time_varying_beta, beta, universal_coverage_feature, coverage_implementation_type, tee, S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0)
    total_start_days=150;
    k_end = 100;
    peak_hosp = zeros(total_start_days, k_end);
    peak_deaths = zeros(total_start_days, k_end);
    if coverage_implementation_type == 1 %step func
        for start_day = 1:total_start_days
            for k = 1:k_end % we will take 1/k to be the fraction gaining coverage on each time step
                global beta_values;
                beta_values = [];
                global loss_values;
                loss_values = [];
                [t,y] = ode45(@(t,y) sihr(t, y, N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i, eta, unemployment_feature, time_varying_beta, beta, universal_coverage_feature, start_day, coverage_implementation_type,1/k,1), tee, [S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0]); 
                peak_deaths(start_day,k) = getTotalDeaths(y(:,9),y(:,10));
                peak_hosp(start_day,k) = getPeakICUHospitalizations(y(:,5),y(:,6));
            end
        end
        figure()
        [C,h] = contourf(peak_deaths)
        clabel(C,h)
        title('Total Death')
        xlabel('k (note: 1/k is the fraction of uninsured that gains coverage each day)') 
        ylabel('Start of Coverage (# days after start of pandemic)') 
        figure()
        [C,h] = contourf(peak_hosp)
        clabel(C,h)
        title('Peak Hosp')
        xlabel('k (note: 1/k is the fraction of uninsured that gains coverage each day)') 
        ylabel('Start of Coverage (# days after start of pandemic)')
    elseif coverage_implementation_type == 2 %periodic delta
        for start_day = 1:total_start_days
            for k = 1:30 % we will take 1/k to be the fraction gaining coverage on each time step
                [t,y] = ode45(@(t,y) sihr(t, y, N, d_u, d_i, c_u, c_i, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i, unemployment_vector, eta, unemployment_feature, time_varying_beta, beta, universal_coverage_feature, start_day, coverage_implementation_type,1 , k), tee, [S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0]); 
                mat(start_day,k) = getPeakDeaths(y(:,9),y(:,10));
            end
        end
        [C,h] = contourf(mat)
        clabel(C,h)
        title('Total Death')
        xlabel('k (period of periodic delta)') 
        ylabel('Start of Coverage (# days after start of pandemic)') 
    end
        

end

function plotBeta()
    global beta_values;
    figure()
    plot(beta_values);
    title('Beta(t)')
    xlabel('time (ODE solution time steps)') 
    ylabel('Beta(t)')

end

function plotLossOfCoverage()
    global loss_values;
    figure()
    plot(loss_values);
    title('L(t)')
    xlabel('time (ODE solution time steps)') 
    ylabel('L(t)')

end

function plotGainOfCoverage()
    global gain_values;
    figure()
    plot(gain_values);
    title('G(t)')
    xlabel('time (ODE solution time steps)') 
    ylabel('G(t)')

end

function plotCompartmentsSeparately(t,y,t0,tf)

% below we plot the results
color = get(gca,'colororder'); % different colors for plotting
figure()
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



function plotAsFractions(t, y, t0, tf)

fntsz=16;
labels={'S_u','S_i','I_u','I_i','H_u','H_i','R_u','R_i','D_u','D_i'};
pltcolors={'blue','blue','red','red','green','green','magenta','magenta','cyan','cyan'};
pltstyles={'--','-'};
figure('Position',[100,100,800,400]);
hold on;
set(gca,'FontSize',fntsz,'defaultTextInterpreter','latex')
Nu=sum(y(:,1:2:end),2); %compute total uninsured at each time
Ni=sum(y(:,2:2:end),2); %compoute total insured at each time
Ncompartments=size(y,2)/2;
for jj=1:Ncompartments
    %plot fraction uninsured in each compartment
    plot(t,y(:,2*jj-1)./Nu,'LineWidth',2,'Color',pltcolors{2*jj-1},'LineStyle','-')
    %plot fraction insured in each compartment
    plot(t,y(:,2*jj)./Ni,'LineWidth',2,'Color',pltcolors{2*jj},'LineStyle','--')
end
legend(labels); %legend
xlim([t0 tf])
xlabel('$t$ (days)')
ylabel('fraction of insured or uninsured')

end