% Lambda = [0.2,-0.2,10,-1.2];
Lambda = Lambda(1,:);

%%

%%%%%%%%% Initialize Variables %%%%%%%%%

%%%%%%%%% Unpack Lambda %%%%%%%%%
theta1_dot = Lambda(1);
theta2_dot = Lambda(2);
theta3_dot = Lambda(3);
theta_dots = [theta1_dot;theta2_dot;theta3_dot];
dvd = [0;Lambda(4);0];  %Relative extrusion velocity (m/s)

%%% Time Variables %%%
dt = 0.001; %Time step size (s)
T = 3;      %Total simulation time (s)

%%% Droplet and Relevent Properties %%%
g = 9.81;                %Gravitational acceleration (m/s^2)
R = 0.001;               %Droplet Radius (m)
V = (4/3)*pi*R^3;        %Droplet volume (m^3)
nu_2 = 0.25;             %Phase 2 volume fraction
rho_1 = 2000;            %Phase 1 density (kg/m^3)
rho_2 = 7000;            %Phase 2 density (kg/m^3)

rho_star = (1-nu_2)*rho_1 + nu_2*rho_2; %Effective density (kg/m^3)
mi = V*rho_star %Droplet mass (kg)

L1 = 0.3; L2 = 0.2; L3 = 0.08; %Rod lengths (m)

%%%%%%%%% Initial Positions %%%%%%%%%
theta10 = pi/2;
theta20 = 0;
theta30 = 0;
theta0s = [theta10;theta20;theta30];
t = linspace(0,T,T/dt);
thetas = theta0s + t.*theta_dots;
r0 = [0;0.5;0];
xd = L1*cos(thetas(1,:))+L2*cos(thetas(2,:))+L3*sin(thetas(3,:));
yd = L1*sin(thetas(1,:))+L2*sin(thetas(2,:));
zd = L3*cos(thetas(3,:));
rd = r0+[xd;yd;zd];
xd_dot = -L1*theta1_dot*sin(thetas(1,:))-L2*theta2_dot*sin(thetas(2,:))+L3*theta3_dot*cos(thetas(3,:));
yd_dot = L1*theta1_dot*cos(thetas(1,:))+L2*theta2_dot*cos(thetas(2,:));
zd_dot = -L3*theta3_dot*sin(thetas(3,:));
vd = [xd_dot;yd_dot;zd_dot];
r0_drops = rd;
v0_drops = vd+dvd;
a0_drops = zeros(size(r0_drops))-[0;9.8;0];
% flight_times = sqrt(2*r0_drops(2,:)/g);
a = 0.5*a0_drops(2,:);
b = v0_drops(2,:);
c = r0_drops(2,:);
flight_times = (-b-sqrt(b.^2-4.*a.*c))./(2.*a);
r_drops = r0_drops+v0_drops.*flight_times+0.5.*a0_drops.*flight_times.^2;

plot(r_drops(3,:),-r_drops(1,:));
xlabel('x');
ylabel('z');
%%
num = 0;
for i = 1:length(r_des)
    num = num + norm(r_des(:,i)-r_gen(:,i));
end
denom = 0;
for i = 1:length(r_des)-1
    denom = denom + norm(r_des(i)-r_des(i+1));
end
num/denom
%% Genetic Algorithm without reevaluating parents
data = load('robotprint_data.mat');
r_des = data.ri;
cf = @(r_gen) sum(vecnorm(r_des-r_gen,1),2)./sum(vecnorm(diff(r_des)),2);
lambda_lower = [15,15,6,-3.5];
lambda_upper = [16,16,7,-3];
parents = 10;
TOL_GA = 2*10^-2;
G = 1000;
S = 100;
dv = 4;

% Generate empty output arrays.
PI = zeros(G,S);
Orig = zeros(G,S);
% Randomly generate first generation.
Lambda = lambda_lower + rand(S,dv).*(lambda_upper-lambda_lower);

t0 = tic;
for g = 1:G
    myProgressBar(toc(t0), g, G);
    % Test fitness of members of first generation.
    if (g == 1)
        for i = 1:S
            lambda = Lambda(i,:);
            r_gen = simplePrinterSim(lambda);
            PI(g,i) = cf(r_gen);
        end
    else
        PI(g,1:parents) = PI(g-1,1:parents);
        for i = parents+1:S
            lambda = Lambda(i,:);
            r_gen = simplePrinterSim(lambda); 
            PI(g,i) = cf(r_gen);
        end
    end
    
    % Rank genetic strings.
    [PI(g,:),Orig(g,:)] = sort(PI(g,:));
    Lambda = Lambda(Orig(g,:)',:);
    % Break if tolerance is met
    if (PI(g,1) <= TOL_GA)
        break;
    end
    % Mate top pairs.
    children = zeros(parents,dv);
    for p = linspace(1, parents-1, parents/2)
        parent1 = Lambda(p,:);
        parent2 = Lambda(p+1,:);
        phi1 = rand;
        phi2 = rand;
        children(p,:)   = parent1.*phi1 + parent2.*(1-phi1);
        children(p+1,:) = parent1.*phi2 + parent2.*(1-phi2);
%         Lambda(parents+p,:) = parent1.*phi1 + parent2.*(1-phi1);
%         Lambda(parents+p+1,:) = parent1.*phi1 + parent2.*(1-phi2);
    end
%     Lambda(2*parents+1:end,:) = lambda_lower + rand(S-2*parents,dv).*(lambda_upper-lambda_lower);
    Lambda(2*parents+1:end,:) = lambda_lower + rand(S-2*parents,dv).*(lambda_upper-lambda_lower);
end
%%
r_gen = simplePrinterSim(Lambda(1,:));
figure
plot(r_gen(3,:),r_gen(1,:))
hold on
plot(r_des(3,:),r_gen(1,:))
hold off

figure
plot(linspace(0,G,G),PI(:,1));
hold on
plot(linspace(0,G,G),mean(PI(:,1:parents),2));
plot(linspace(0,G,G),mean(PI,2));
legend('Best designs','Mean of parents','Mean of population');
title('Cost vs Generation')
xlabel('Generation');
ylabel('Cost');
hold off

figure
familyTree(Orig, parents, 100)

figure
comparepattern(Lambda(1,1),Lambda(1,2),Lambda(1,3),Lambda(1,4))
%%
i = 4;
a = Lambda(i,1);
b = Lambda(i,2);
c = Lambda(i,3);
d = Lambda(i,4);
comparepattern(a,b,c,d)