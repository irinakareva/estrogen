function Kareva_Brown_estrogen_Liebig
% This code reproduces key figures from the manuscript 
%"Estrogen as an essential resources and the coexistence of ER+ and ER- 
% cancer cells" by Irina Kareva and Joel Brown, (c) 2021

% System of Equations is given in System (3) of the manuscript, and is derived in the Appendix

clc
options=odeset('RelTol',1.0e-10);

tmax = 2000; %simulation runs until time tmax

%% initial conditions

p.N0     = 1;       %initial population size
p.R1init = 1.1; %initial condition R1 (estrogen)
p.R2init = 3.9; %initial condition R2 (glucose)

p.q0     = 0;  % initial condition for auxiliary variable q(t)
p.pp0    = 0; % initial condition for auxiliary variable p(t)
p.tt0    = 0; 

%% parameter values

p.a11  = 0.9; %growth rate ER+ (strategy 1) from R1 (estrogen)
p.a21  = 3; %growth rate ER+ (strategy 1) from R2 (glucose)
p.a22  = 0.25; %growth rate ER- (strategy 2) from R2 (glucose)
%b_{ji} is conversion rate of resource j into proliferation when in state i
%b_{ji} is the probability of taking up resource j when in state i

p.b11  = 0.9; %handling time of R1 for growth of ER+
p.b21  = 1; %handling time of R2 for growth of ER+
p.b22  = 1; %handling time of R2 for growth of ER

p.d    = 0.01; %death rate all cells

p.k1   = 0.01; % clearance of R1 (estrogen)
p.R01  = p.R1init*p.k1; %inflow R1 (estrogen)

p.k2   = 0.01; % clearance of R2 (glucose)
p.R02  = p.R2init*p.k2; %inflow R2 (glucose)

p.mu = -2;
p.a=0; p.b=1; %boundaries of truncated initial distribution

%% solving system of equations

time=[0;tmax]; %time interval over which we integrate
in_cond=[p.q0 p.pp0 p.R1init p.R2init p.tt0]; %initial conditions q(t)=0, t=0;

[T, Y]=ode23s(@estrogen,time,in_cond,options,p);
%ode23 is a solver that will numerically integrate function "estrogen"
%over time period "time" with initial conditions "in_cond", with precision
%defined in "options" and parameter values defined in "p"
% The results will be given in vector Y

q    = Y(:,1); %auxiliary "keystone" variable from the solution
pp   = Y(:,2); %auxiliary "keystone" variable from the solution
R1_E = Y(:,3);
R2_G = Y(:,4);
tt   = Y(:,5); 

%% statistical characteristics
% in this section we calculate moment generating function, the expected
% value of parameter alpha and its variance as they change over time
% here, we formulas are written for truncated initial distribution, and
% should be modified for other distributions
for i=1:length(Y)
    m01(i)      = (p.mu*(exp(p.a*(q(i)-p.mu))-exp(p.b*(q(i)-p.mu))))./((p.mu-q(i))*(exp(-p.mu*p.a)-exp(-p.mu*p.b))); %moment generating function for truncated exponential distribution on interval [0 1]
    exp_val(i)  = p.b+1/(p.mu-q(i))+(p.a-p.b)*exp(p.b*p.mu+p.a*q(i))./(exp(p.b*p.mu+p.a*q(i))-exp(p.a*p.mu+p.b*q(i))); %expected value of alpha for truncated exponential distribution on interval [0 1]
    variance(i) = 1./(p.mu-q(i))^2-(((p.a-p.b)^2)*exp((p.a+p.b)*(p.mu+q(i))))./((exp(p.b*p.mu+p.a*q(i))-exp(p.a*p.mu+p.b*q(i))).^2); %variance of alpha for truncated exponential distribution on interval [0 1]
    Nt(i)       = p.N0*exp(pp(i))*m01(i);% total population size 
end
    
% here, we calclate the change in distribution over clones within the
% population over time
mu1=linspace(0,1,10); % defining the space for paramter mu on interval [0 1] 
for i=1:length(Y)
    for j=1:length(mu1)
        p0al(i)    = p.mu*exp(-p.mu*mu1(j))/(1-exp(-p.mu));
        pdist(i,j) = p0al(i)*(exp(mu1(j)*(q(i))))/m01(i);
    end
end


%% Plotting the results
figure
subplot(2,3,1)
hold on
plot(T,Nt,'LineWidth',2)
ylabel('total population size N(t)')
xlabel('time')

subplot(2,3,2)
hold on
plot(T,R1_E,'LineWidth',2)
ylabel('estrogen resource')
xlabel('time')

subplot(2,3,3)
hold on
plot(T,R2_G,'LineWidth',2)
ylabel('glucose resource')
xlabel('time')
ylim([0.02 0.1])

subplot(2,3,4)
hold on
plot(T,exp_val,'LineWidth',2)
ylabel('expected value E^t[\alpha]')
xlabel('time')

subplot(2,3,5)
hold on
plot(T,variance,'LineWidth',2)
ylabel('variance Var^t[\alpha]')
xlabel('time')

subplot(2,3,6)
surf(mu1,T,pdist,'EdgeAlpha',0)  % The 3D Graph
xlabel('\alpha')
ylabel('Time t')
zlabel('distribution P^t[\alpha]')
ylim([0 tmax])
grid on, shading faceted %or 'shading interp'
colormap(sqrt(hsv))  %We can modify the color matrix by exponentiation ;-]
axis xy



%defining the equations that are used in ode23 solver
function F1=estrogen(t,y,p)

q  = y(1); %auxiliary "keystone" variable from the solution
pp = y(2);
R1 = y(3);
R2 = y(4);
tt = y(5); %time from the solution

m01     = (p.mu*(exp(p.a*(q-p.mu))-exp(p.b*(q-p.mu))))/((p.mu-q)*(exp(-p.mu*p.a)-exp(-p.mu*p.b))); %calculating the moment generation function
nt      = p.N0*exp(pp)*m01; %calculating total population size
exp_val = p.b+1/(p.mu-q)+(p.a-p.b)*exp(p.b*p.mu+p.a*q)./(exp(p.b*p.mu+p.a*q)-exp(p.a*p.mu+p.b*q));

min4q   = 1*min(p.a11*p.b11.*R1,p.a21*p.b21.*R2);
min4R1  = 1*min(p.a11*R1,p.a21*p.b21.*R2/p.b11);
min4R2  = 1*min(p.a11*p.b11.*R1/p.b21,p.a21.*R2);


% these are the equations that are actually being solved
% we then calculate all the necessary characteristics in the main program
dq  = min4q-p.a22*p.b22*R2;
dp  = p.a22*p.b22*R2-p.d;
dR1 = p.R01-p.k1*R1-exp_val.*nt.*min4R1;
dR2 = p.R02-p.k2*R2-exp_val.*nt*min4R2-(1-exp_val).*nt*p.a22*R2;
dtt = 1;

F1=[dq;dp;dR1;dR2;dtt];
