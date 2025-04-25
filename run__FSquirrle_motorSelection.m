%%
clearvars

v = 1; % desired speed
a = 2; % desired acceleration
F = 20; % magnitude of the desired force

T = 0; % external rotating moment on robot;
Omega = 0; % rotational velcoity of the robot

AnchorCircleRadius = 0.6;
A1 = AnchorCircleRadius*[cosd(90); sind(90); 0];
A2 = AnchorCircleRadius*[cosd(210); sind(210);0];
A3 = AnchorCircleRadius*[cosd(-30); sind(-30);0];

RobotCircleRadius = 0.1;
C1_u = [RobotCircleRadius*[cosd(90);  sind(90)] ; 0.1]; % the anchor points relative to robot's handle
C2_u = [RobotCircleRadius*[cosd(210); sind(210)]; 0.1];
C3_u = [RobotCircleRadius*[cosd(-30); sind(-30)]; 0.1];

R = [0;0;0.1]; % robot's handle position. 

r_pulley = 0.005; % radius of of pulley, in meters
p_screw = 0.010; % lead screw pitch, in meters
eff_screw = 0.5; % efficiency of lead screw

hf999 = figure(999); 
clf(hf999)
plot_FSquirrleBases(R,C1_u, C2_u, C3_u, A1, A2, A3, figure(999))


%%
m = 5;

gearRatio = 1;

phiRange = 0:0.05:2*pi;

xRange = -0.3:0.05:0.3;
yRange = -0.3:0.05:0.1;
zRange = 0.1:0.1:0.5;

numberOfVectors = 1000;
unitVectors = fibonacci_sphere(numberOfVectors);


[tau_F_only_line,tau_F_only_screw,omega_pulley_line,omega_pulley_screw] = ...
    deal(nan(length(xRange),length(yRange),length(zRange)));

for i = 1:length(xRange)
    for j = 1:length(yRange)
        for k = 1:length(zRange)
            R = [xRange(i); yRange(j); zRange(k)]

            J_inv = FSquirrle_Jacobian_inverse(R(1),R(2),R(3),A1(1),A1(2),A2(1),A2(2),A3(1),A3(2),C1_u(1),C1_u(2),C1_u(3),C2_u(1),C2_u(2),C2_u(3),C3_u(1),C3_u(2),C3_u(3),r_pulley,p_screw);
            J = FSquirrle_Jacobian(R(1),R(2),R(3),A1(1),A1(2),A2(1),A2(2),A3(1),A3(2),C1_u(1),C1_u(2),C1_u(3),C2_u(1),C2_u(2),C2_u(3),C3_u(1),C3_u(2),C3_u(3),r_pulley,p_screw);
            [B,~,B_l] = FSquirrle_bases(R(1),R(2),R(3),A1(1),A1(2),A2(1),A2(2),A3(1),A3(2),C1_u(1),C1_u(2),C1_u(3),C2_u(1),C2_u(2),C2_u(3),C3_u(1),C3_u(2),C3_u(3));
            
            if min([dot(B_l(:,1),B_l(:,2)) , dot(B_l(:,1),B_l(:,3)) , dot(B_l(:,2),B_l(:,3))])<cosd(165)
                continue
            end

            omega_pulley = [];
            tau_F_only = [];
            % clf(hf999)
            % plot_FSquirrleBases(R,C1_u, C2_u, C3_u, A1, A2, A3, hf999)
            for u = unitVectors'
                handForce = -u*F;
                
                lineForces = lsqnonneg(B,handForce);
                tau_F_only = [tau_F_only, lineForces.*[r_pulley/gearRatio; r_pulley/gearRatio; r_pulley/gearRatio; p_screw/2/pi/eff_screw]];

                % ax = a*cos(phi+phaseshift);
                % ay = a*sin(phi+phaseshift);
                % robotAcc = [ax;ay;0];
                % tau_a_only = [tau_a_only, J' * (M*robotAcc)];
                %

                Vel = u*v;
                omega_pulley = [omega_pulley, (J_inv*Vel).*[gearRatio; gearRatio; gearRatio; 1]];
            end

            omega_pulley_line(i,j,k) = max(omega_pulley(1:3,:),[],'all');
            omega_pulley_screw(i,j,k) = max(omega_pulley(4,:),[],'all');
            tau_F_only_line(i,j,k) = max(tau_F_only(1:3,:),[],'all');
            tau_F_only_screw(i,j,k) = max(tau_F_only(4,:),[],'all');
        end
    end
end

%%
% tau_combined = tau_a_only + tau_F_only;

figure(1); clf
[xx,yy] = meshgrid(xRange,yRange);

for i = 1:size(tau_F_only_line,3)
    subplot(1,2,1)
    surface(xx, yy, tau_F_only_line(:,:,i)')
    hold all
    title('max torque - line motors')

    subplot(1,2,2)
    surface(xx, yy, tau_F_only_screw(:,:,i)'/gearRatio)
    hold all
    title('max torque - screw motor')

end


figure(2); clf
for i = 1:size(omega_pulley_line,3)
    subplot(1,2,1)
    surface(xx, yy, omega_pulley_line(:,:,i)'*gearRatio * 60/2/pi)
    hold all
    title('max speed - line motors')
    % 
    subplot(1,2,2)
    surface(xx, yy, omega_pulley_screw(:,:,i)' * 60/2/pi)
    hold all
    title('max speed - screw motor')

end

