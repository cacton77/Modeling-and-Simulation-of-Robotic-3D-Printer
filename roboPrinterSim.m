function [r_drops, flight_times] = roboPrinterSim(Lambda,animate)

%Animated GIF Initialization
if (animate ~= 0)
    close all;
    init = false;
    first = true;
    last = false;
    f = figure;
    f.WindowState = 'maximized';
    axis manual tight;
    filename = '1ElFieldOn.gif';
end

%%%%%%%%% Initialize Variables %%%%%%%%%

%%% Time Variables %%%
dt = 0.001; %Time step size (s)
T = 3;      %Total simulation time (s)

%%% Droplet and Relevent Properties %%%
g = 9.81;                %Gravitational acceleration (m/s^2)
elec_perm = 8.854*10^12; %Electric permittivity (F/m)
R = 0.001;               %Droplet Radius (m)
AiD = pi*R^2;            %Drag reference area (m^2)
V = (4/3)*pi*R^3;        %Droplet volume (m^3)
nu_2 = 0.25;             %Phase 2 volume fraction
rho_1 = 2000;            %Phase 1 density (kg/m^3)
rho_2 = 7000;            %Phase 2 density (kg/m^3)
q_1 = 0;                 %Phase 1 charge capacity (C/m3)
q_2 = 1*10^-3;           %Phase 2 charge capacity (C/m^3)
qp = -8*10^-5;           %Grid pixel charge (C)
rho_a = 1.225;           %Surrounding medium density (kg/m^3)
mu_f = 1.8*10^-5;        %Surrounding medium viscosity (Pa/s)

rho_star = (1-nu_2)*rho_1 + nu_2*rho_2; %Effective density (kg/m^3)
q_star = (1-nu_2)*q_1 + nu_2*q_2;        %Effective charge capacity (C/m^3)

mi = V*rho_star; %Droplet mass (kg)
qi = V*q_star;   %Droplet charge (C)

vf  = [0.5;0;0.5]; %Surrounding medium velocity (m/s)

L1 = 0.3; L2 = 0.2; L3 = 0.08; %Rod lengths (m)

%%%%%%%%% Creating the charged grid %%%%%%%%%
L_bed = 0.8; %Charged grid side length (m)
range = linspace(-L_bed/2,L_bed/2,10);
bed = zeros(10,3,10);
[bed(:,1,:),bed(:,2,:),bed(:,3,:)] = meshgrid(-range,0,range);

%%%%%%%%% Initial Positions %%%%%%%%%
theta10 = pi/2;
theta20 = 0;
theta30 = 0;
theta0s = [theta10;theta20;theta30];
thetas = theta0s;
r0 = [0;0.5;0];
xd = L1*cos(theta0s(1))+L2*cos(theta0s(2))+L3*sin(theta0s(3));
yd = L1*sin(theta0s(1))+L2*sin(theta0s(2));
zd = L3*cos(theta0s(3));
rd = r0+[xd;yd;zd];

r_drops = Inf(3,T/dt);
v_drops = Inf(3,T/dt);
stopped = zeros(1,T/dt);
start_times = zeros(T/dt,1);
flight_times = zeros(T/dt,1);
Nd = 0;

%%%%%%%%% Unpack Lambda %%%%%%%%%
theta1_dot = Lambda(1);
theta2_dot = Lambda(2);
theta3_dot = Lambda(3);
theta_dots = [theta1_dot;theta2_dot;theta3_dot];
dvd = [0;Lambda(4);0];  %Relative extrusion velocity (m/s)

%%%%%%%%% Start the Time Loop %%%%%%%%%
t = dt;
while sum(stopped) < length(r_drops)
    %%%%%%%%% Animate %%%%%%%%%
    if (t == T)
        last = true;
    end
    if (mod(t/dt,10) < 1 && animate ~= 0)
        time = strcat('time = ', num2str(fix(t)));
        %Plot current Positions, and capture image
        r1 = [0,L1*cos(thetas(1));0.5,0.5+L1*sin(thetas(1));0,0];
        r2 = r1(:,2) + [0,L2*cos(thetas(2));0,L2*sin(thetas(2));0,0];
        r3 = [r2(:,2),rd];
        if (animate == 1)
            subplot(1,2,1)
            plot3(r1(3,:),r1(1,:),r1(2,:));
            pbaspect([1,1,1]);
            hold on
            plot3(r2(3,:),r2(1,:),r2(2,:));
            plot3(r3(3,:),r3(1,:),r3(2,:));
            plot3([0,r0(3)],[0,r0(1)],[0,r0(2)]);
            c = flight_times;
            scatter3(r_drops(3,:),r_drops(1,:),r_drops(2,:),15,c);
            scatter3(reshape(bed(:,3,:),[100,1]),...
                     reshape(bed(:,1,:),[100,1]),...
                     reshape(bed(:,2,:),[100,1]),'*','g');
            axis([-L_bed/2 L_bed/2 -L_bed/2 L_bed/2 -2 2]);
            xlabel('z')
            ylabel('x')
            zlabel('y')
            cb = colorbar('eastoutside');
            cb.Label.String = 'Flight time (s)';
            title('Electric Field On')
            grid on
            hold off
            subplot(2,2,2)
            plot(r1(3,:),r1(1,:));
            pbaspect([1,1,1]);
            hold on
            plot(r2(3,:),r2(1,:));
            plot(r3(3,:),r3(1,:));
            scatter(r_drops(3,:),r_drops(1,:),15,c);
            scatter(reshape(bed(:,3,:),[100,1]),...
                     reshape(bed(:,1,:),[100,1]),'*','g');
            axis([-L_bed/2 L_bed/2 -L_bed/2 L_bed/2]);
            xlabel('z')
            ylabel('x')
            hold off
            subplot(2,2,4)
            plot(r1(1,:),r1(2,:));
            pbaspect([1,1,1]);
            hold on
            plot(r2(1,:),r2(2,:));
            plot(r3(1,:),r3(2,:));
            plot([0,r0(3)],[0,r0(2)]);
            scatter(r_drops(1,:),r_drops(2,:),15,c);
            scatter(reshape(bed(:,1,:),[100,1]),...
                     reshape(bed(:,2,:),[100,1]),'*','g');
            axis([-L_bed/2 L_bed/2 -0.1 2]);
            xlabel('x')
            ylabel('y')
            title(time)
            hold off
        elseif (animate == 2)
            if (first == true)
                first = false;
                scatter3(rd(3),rd(1),rd(2),'p');
                text(rd(3),rd(1),rd(2)+.1,'Start');
                hold on
            elseif (last == true)
                scatter3(rd(3),rd(1),rd(2),'p');
                text(rd(3),rd(1),rd(2)+.1,'End');
            else
                scatter3(rd(3),rd(1),rd(2),'.','r');
            end
            pbaspect([1,1,1]);
            axis([-L_bed/2 L_bed/2 -L_bed/2 L_bed/2 -2 2]);
            scatter3(reshape(bed(:,3,:),[100,1]),...
                     reshape(bed(:,1,:),[100,1]),...
                     reshape(bed(:,2,:),[100,1]),'*','g');
            xlabel('z')
            ylabel('x')
            zlabel('y')
            cb = colorbar('eastoutside');
            cb.Label.String = 'Flight time (s)';
        end
        
        drawnow
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if init == false 
            init = true;
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.01); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01); 
        end
    end    
    if (t <= T)
        % Update thetas
        thetas = theta0s+theta_dots.*t;
        % Update arm position
        xd = L1*cos(thetas(1))+L2*cos(thetas(2))+L3*sin(thetas(3));
        yd = L1*sin(thetas(1))+L2*sin(thetas(2));
        zd = L3*cos(thetas(3));
        rd = r0+[xd;yd;zd];
        % Update arm velocity
        xd_dot = -L1*theta1_dot*sin(thetas(1))-L2*theta2_dot*sin(thetas(2))+L3*theta3_dot*cos(thetas(3));
        yd_dot = L1*theta1_dot*cos(thetas(1))+L2*theta2_dot*cos(thetas(2));
        zd_dot = -L3*theta3_dot*sin(thetas(3));
        vd = [xd_dot;yd_dot;zd_dot];

        %When the arm deposits a droplet we initialize the droplet's position and
        %velocity as this:
        ri0 = rd;
        vi0 = vd+dvd;
        Nd = Nd + 1;
        r_drops(:,Nd) = ri0;
        v_drops(:,Nd) = vi0;
        start_times(Nd) = t;
    end
    %%%%%%% Droplet dynamics %%%%%%
    for i = 1:Nd
        ri = r_drops(:,i);
        vi = v_drops(:,i);
        %%%Check to see if droplet has reached substrate%%%
        if (stopped(i) == true)
            continue;
        end
        %flight_times(i) = t-flight_times(i);
        %%%Force of Gravity%%%
        Fi_grav = [0;-mi*g;0];
        %%%Electric Field Force%%%
        bed_ri = -bed+permute(ri, [2 1 3]);
        norm_bed_ri = vecnorm(bed_ri,2,2);
        bed_1_2_3 = bed_ri./norm_bed_ri.^2;
        bed_sum_2_3 = sum(bed_1_2_3,1);
        bed_sum_2 = sum(bed_sum_2_3,3);
        bed_sum_1 = permute(bed_sum_2,[2 3 1]);
        Fi_elec = (qp.*qi)./(4.*pi.*elec_perm).*bed_sum_1;
        %Fi_elec = [0;0;0];
        %%%Force of Drag%%%
        %First we find the Reynolds Number.
        Re = 2*R*rho_a*norm(vf-vi)/mu_f;
        %Then we use that to calculate the drag coefficient.
        if (Re > 0 && Re <= 1)
            CDi = 24/Re;
        elseif (Re > 1 && Re <= 400)
            CDi = 24/Re^0.0646;
        elseif (Re > 400 && Re <= 3*10^5)
            CDi = 0.5;
        elseif (Re > 3*10^5 && Re < 2*10^6)
            CDi = 0.000366*Re^0.4275;
        elseif (Re > 2*10^6)
            CDi = 0.18;
        else
            ME = MException('ReNum:InvalidNum','Invalid Reynolds Number');
            throw(ME);
        end
        Fi_drag = (1/2)*rho_a*CDi*norm(vf-vi).*(vf-vi)*AiD;
        %Sum all forces on i.
        Psii_tot = Fi_grav + Fi_elec + Fi_drag;
        %Calculate new 
        ai = Psii_tot/mi; %Acceleration of i caused by forces at time t.
        vi_new = vi+dt*ai; %Velocity of i at t + dt due to accel at t.
        ri_new = ri+dt*vi; %Position of i at t + dt due to vel at t.
        flight_times(i) = t-start_times(i);
        if (ri_new(2) <= 0)
            ri_new(2) = 0;
            vi_new = [0;0;0];
            stopped(i) = true;
        end
        v_drops(:,i) = vi_new; %Update velocity of i to next time step.
        r_drops(:,i) = ri_new; %Update position of i to next time step.
    end
    t = t + dt;
end
end

