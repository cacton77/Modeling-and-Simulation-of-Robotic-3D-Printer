function [r_drops, flight_times] = simplePrinterSim(Lambda)

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
mi = V*rho_star; %Droplet mass (kg)

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
a = 0.5*a0_drops(2,:);
b = v0_drops(2,:);
c = r0_drops(2,:);
flight_times = (-b+sqrt(b.^2-4.*a.*c))./(2.*a);
r_drops = r0_drops+v0_drops.*flight_times+0.5.*a0_drops.*flight_times.^2;
end
