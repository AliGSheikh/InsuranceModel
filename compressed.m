%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AIM CTMC SEIR Group Project %%
%&      June 2020              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%close all
set(0,'DefaultAxesFontSize', 18)
set(gca,'fontsize',18);
%initial values
initi_i=1; %I(0)=initi_i;
initi_u=1;
%inite=0; %E(0)=inite

%Parameter values
alpha=1/4; %recovery rate
gamma=1/4; %infectious -> hospitalized rate
xsi=.01; %disease-related death rate

beta=2*(alpha+xsi); %%%CHANGE to 1.5, 2,  4
R0=beta/(alpha+xsi);
N=100; %poulation size
sim=3; %number of sample paths
time=100;

  for k=1:sim
    clear s i h r d
    cas=initi_i;
    tot=N;
    i(1)=initi_i + initi_u;
    s(1)=N-initi_i-initi_u;
    r(1)=0;
    h(1)-0;
    j=1;
    sumdeath=0;
    while e(j)+i(j)>0 %when E+I=0, then epidemic stops
		u1=rand;
	    u2=rand;
        totev=(beta/tot)*i(j)*s(j)+gamma*e(j)+xsi*i(j)+alpha*i(j);%interevent time
		t(j+1)=-log(u1)/totev+t(j); %Interevent time
        ev1=(beta/tot)*i(j)*s(j)/totev;
        ev2=ev1+gamma*e(j)/totev;
        ev3=ev2+xsi*i(j)/totev; 
        ev4=ev3+alpha*i(j)/totev; % Equals one
		if (u2<= ev1) % transmission
		   e(j+1)=e(j)+1;
		   s(j+1)=s(j)-1;
           i(j+1)=i(j);
           r(j+1)=r(j); 
        elseif (u2> ev1 &&  u2<=ev2) % Exposed->Infectious
			e(j+1)=e(j)-1;
			i(j+1)=i(j)+1;
            s(j+1)=s(j);
            r(j+1)=r(j);
            cas=cas+1; %count number of infectious cases 
        elseif  (u2>ev2  &&  u2<=ev3) %Disease-related death
           i(j+1)=i(j)-1;
           s(j+1)=s(j);
           e(j+1)=e(j);
           r(j+1)=r(j);
           sumdeath=sumdeath+1;
        elseif (u2>ev3  && u2<= ev4) %Recovery
            i(j+1)=i(j)-1;
            r(j+1)=r(j)+1;
            s(j+1)=s(j);
            e(j+1)=e(j);
        end
        j=j+1;
        tot=s(j)+e(j)+i(j)+r(j);
    end
    figure(3)
    if k==1
      subplot(1,2,2)
      stairs(t,i,'r-','Linewidth',2);
      xlabel('Day')
      ylabel('Infectious Individuals')
      hold on
      subplot(1,2,1)
      stairs(t,e,'r-','Linewidth',2)
      xlabel('Day')
      ylabel('Exposed Individuals')
      hold on
    elseif k==2
      subplot(1,2,2)
      stairs(t,i,'g-','Linewidth',2);
      hold on
      subplot(1,2,1)
      stairs(t,e,'g-','Linewidth',2);
      hold on
    elseif k==3
      subplot(1,2,2)
      stairs(t,i,'b-','Linewidth',2);
      hold on
      subplot(1,2,1)
      stairs(t,e,'b-','Linewidth',2)
      hold on
    end
    casetot(k)=cas; %total number of cases
    Tend(k)=t(j); %time epidemic ends
    Deaths(k)=sumdeath; %total number of deaths
  end
 x(1)=N-inite-initi;
 u(1)=inite;
 y(1)=initi;
 z(1)=0;
 tot=N;
 dt=.005;
 sumcases=initi;
 sumdeaths=0;
for k=1:time/dt  %% Euler's method to solve ODE
    x(k+1)=x(k)+dt*(-beta*x(k)*y(k)/tot);
    u(k+1)=u(k)+dt*(beta*x(k)*y(k)/tot-gamma*u(k));
    y(k+1)=y(k)+dt*(gamma*u(k)-alpha*y(k)-xsi*y(k));
    z(k+1)=z(k)+dt*(alpha*y(k));
    tot=x(k+1)+u(k+1)+y(k+1)+z(k+1);
    sumcases=sumcases+dt*(gamma*u(k));
    sumdeaths=sumdeaths+dt*(xsi*y(k));
end
subplot(1,2,2)
plot([0:dt:time],y,'k-','LineWidth',2);
hold off
subplot(1,2,1)
plot([0:dt:time],u,'k-','LineWidth',2);
hold off
 str = sprintf( 'CTMC SEIR, R_0 = %.2f',R0);
 sgtitle(str);

 disp((sprintf('---------------------------------------------------------------------------')))
 disp((sprintf('N= %5.0f, gamma= %5.2f, alpha=%5.2f, delta=%5.2f, beta= %5.2f,', N, alpha, xsi, gamma, beta)))
 disp((sprintf('I(0)= %2.0f,  E(0)=%2.0f, R0= %2.2f, ODECase=%3.2f, ODEDeath=%3.2f', initi, inite,  R0, sumcases,sumdeaths)));
 disp((sprintf('Case1=%3.0f,  Case2=%3.0f, Case3=%3.0f ', casetot(1), casetot(2),casetot(3))))
 disp((sprintf('Tend1=%5.2f, Tend2=%5.2f, Tend3=%5.2f ', Tend(1), Tend(2), Tend(3))))
 disp((sprintf('Death1=%2.0f, Death2=%2.0f, Death3=%2.0f', Deaths(1), Deaths(2), Deaths(3))))
 disp((sprintf('----------------------------------------------------------------------------')))
