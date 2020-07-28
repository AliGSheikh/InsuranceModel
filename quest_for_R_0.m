%--------------------------------------------------- Finding R_0

%-----------------------------------------Punit's method
%----------------- Inputing the parameter values
beta = 0.5;
N = 326000040;
S_0 = 326e6;
delta_u = 1/14;
delta_i = 1/14;
gamma_u = 1/5;
gamma_i = 1/5;

t0 = 1;
tf = 150; % unit = days
time_steps=150;
tee=linspace(t0,tf,time_steps);


% ----------------- Calculating c_i
h = 0.05;
S_u_0 = 26e6;
S_i_0 = 300e6;
p_i = S_i_0/(S_i_0+S_u_0); % percentage of total population that is insured
p_u = S_u_0/(S_i_0+S_u_0);
beta = 0.045;
c_u = (h - p_i*beta)/p_u;

%------------------ Making the matrix
M = [(beta*S_0)/N-delta_u*(1-c_u)-gamma_u*c_u (beta*S_0)/N; (beta*S_0)/N (beta*S_0)/N-delta_i*(1-beta)-gamma_i*beta];


%--------------------------------------Next Generation method
F = [beta*S_u_0/N beta*S_u_0/N; beta*S_i_0/N beta*S_i_0/N];
V = [delta_u*(1-c_u)+gamma_u*c_u 0; 0 delta_i*(1-beta)+gamma_i*beta];
inv(V);
det(F*inv(V))
trace(F*inv(V))
x = eig(F*inv(V));
R_0_ng = max(x)

%-------------------- Another Next Generation method
T = [beta*S_u_0/N beta*S_u_0/N; beta*S_i_0/N beta*S_i_0/N];
S = [-delta_u*(1-c_u)-gamma_u*c_u 0; 0 -delta_i*(1-beta)-gamma_i*beta];
det(-T*inv(S))
trace(-T*inv(S))
y = eig(-T*inv(S));
R_0_ng2 = max(y)


%----------------------- Plotting to find R_0 using Punit's method

%{
plot(-beta*S_u_0*d_i*(1-c_i)/N - beta*S_u_0*gamma_i*c_i/N - beta*S_i_0*delta_u*(1-c_u)/N + 
        delta_u*(1-c_u)*delta_i*(1-c_i) + delta_u*(1-c_u)*gamma_i*c_i - beta*S_i_0*gamma_u*c_u/N +
        delta_i*(1-c_i)*gamma_u*c_u + gamma_u*c_u*gamma_i*c_i, 
%}
%-beta*S_u_0*delta_i*(1-c_i)/N - beta*S_u_0*gamma_i*c_i/N - beta*S_i_0*delta_u*(1-c_u)/N + delta_u*(1-c_u)*delta_i*(1-c_i) + delta_u*(1-c_u)*gamma_i*c_i - beta*S_i_0*gamma_u*c_u/N + delta_i*(1-c_i)*gamma_u*c_u + gamma_u*c_u*gamma_i*c_i


%%
Nscan=11;%number of parameter values to scan over
beta_vals=linspace(0.25,1.5,Nscan); %set up a range of values for c_I to scan over
 
R_0_equation=zeros(Nscan,1); %initialize vector to store peak uninsured hospitalizations
 
for jj=1:Nscan %add a loop to run the simulation for each value
    
    beta=beta_vals(jj); %update the value of c_u
  
    %c_u = (h - p_i*beta)/p_u % varying c_i but avg of hosp rate is the same 
    
    %c_u_vals(jj)= c_u; 
    %run simulation with new value of c_u
    [t,y] = ode45(@(t,y) sihr(t, y, N, d_u, d_i, c_u, beta, alpha_u, alpha_i, delta_u, delta_i, gamma_u, gamma_i, ksi_u, ksi_i), tee, [S_u_0, S_i_0, I_u_0, I_i_0, H_u_0, H_i_0, R_u_0, R_i_0, D_u_0,D_i_0]); 
    
    %store peak number of uninsured hopitatlizations
    R_0_equation(jj)=beta\N*(S_u_0*(delta_i*(1-c_i)-gamma_i*c_i) + S_i_0*(delta_u*(1-c_u) - gamma_u*c_u))/((delta_i*(1-c_i) - gamma_i*c_i)*(delta_u*(1-c_u) - gamma_u*c_u));  %get max of 5th column of y, which stores H_u
   %myplot(t, y, t0, tf)
end
 
%plot the peak hospitalization of uninsured for the different values of c_u
figure;
plot(beta_vals,R_0_equation,'o-')
xlabel('beta')
ylabel('R_0')
 







