%%
clearvars

v = 1; % desired speed
a = 2; % desired acceleration
F = 20; % magnitude of the desired force

T = 0; % external rotating moment on robot;
Omega = 0; % rotational velcoity of the robot

phaseshift = pi; % the angle between robot's acceleration and hand force. maximum torques at pi

AnchorCircleRadius = 0.6;
A1 = AnchorCircleRadius*[cosd(90); sind(90); 0];
A2 = AnchorCircleRadius*[cosd(210); sind(210);0];
A3 = AnchorCircleRadius*[cosd(-30); sind(-30);0];

RobotCircleRadius = 0.1;
C1_u = [RobotCircleRadius*[cosd(90);  sind(90)] ; 0.1]; % the anchor points relative to robot's handle
C2_u = [RobotCircleRadius*[cosd(210); sind(210)]; 0.1];
C3_u = [RobotCircleRadius*[cosd(-30); sind(-30)]; 0.1];

R = [0;0;0.1]; % robot's handle position. 

r_pulley = 0.005; % radius of of pulley
screwAngle = 0.01; % some sort of pitch, like tangent of screw angle.



hf999 = figure(999); 
clf(hf999)
plot_FSquirrleBases(R,C1_u, C2_u, C3_u, A1, A2, A3, figure(999))


%%
m = 5;

gearRatio = 1;

phiRange = 0:0.05:2*pi;

xRange = -0.3:0.02:0.3;
yRange = -0.3:0.02:0.1;
zRange = 0.1:0.1:0.5;


[tau_F_only_max,tau_F_only_min,omega_pulley_max,omega_pulley_min] = ...
    deal(nan(length(xRange),length(yRange),length(zRange)));

for i = 1:length(xRange)
    for j = 1:length(yRange)
        for k = 1:length(zRange)
            R = [xRange(i); yRange(j); zRange(k)];

            J_inv = FSquirrle_Jacobian_inverse(R(1),R(2),R(3),A1(1),A1(2),A2(1),A2(2),A3(1),A3(2),C1_u(1),C1_u(2),C1_u(3),C2_u(1),C2_u(2),C2_u(3),C3_u(1),C3_u(2),C3_u(3),r_pulley);
            J = FSquirrle_Jacobian(R(1),R(2),R(3),A1(1),A1(2),A2(1),A2(2),A3(1),A3(2),C1_u(1),C1_u(2),C1_u(3),C2_u(1),C2_u(2),C2_u(3),C3_u(1),C3_u(2),C3_u(3),r_pulley);
            [B,~,B_l] = FSquirrle_bases(R(1),R(2),R(3),A1(1),A1(2),A2(1),A2(2),A3(1),A3(2),C1_u(1),C1_u(2),C1_u(3),C2_u(1),C2_u(2),C2_u(3),C3_u(1),C3_u(2),C3_u(3));
            
            if min([dot(B_l(:,1),B_l(:,2)) , dot(B_l(:,1),B_l(:,3)) , dot(B_l(:,2),B_l(:,3))])<cosd(165)
                continue
            end

            omega_pulley = [];
            tau_F_only = [];
            % clf(hf999)
            % plot_FSquirrleBases(R,C1_u, C2_u, C3_u, A1, A2, A3, hf999)
            for phi = phiRange
                Fx = F*cos(phi);
                Fy = F*sin(phi);
                Fz = 0;

                handForce = -[Fx;Fy;Fz];
                
                lineForces = lsqnonneg(B,handForce);
                tau_F_only = [tau_F_only, lineForces.*[r_pulley; r_pulley; r_pulley; screwAngle]];

                % ax = a*cos(phi+phaseshift);
                % ay = a*sin(phi+phaseshift);
                % robotAcc = [ax;ay;0];
                % tau_a_only = [tau_a_only, J' * (M*robotAcc)];
                %
                vx = v*cos(phi);
                vy = v*sin(phi);
                vz = 0;

                V = [vx;vy;vz];
                omega_pulley = [omega_pulley, J_inv*V];
            end

            omega_pulley_max(i,j,k) = max(omega_pulley(:));
            omega_pulley_min(i,j,k) = min(omega_pulley(:));
            tau_F_only_max(i,j,k) = max(tau_F_only(1:3,:),[],'all');
            tau_F_only_min(i,j,k) = min(tau_F_only(1:3,:),[],'all');
        end
    end
end

%%
% tau_combined = tau_a_only + tau_F_only;

figure(1); clf
[xx,yy] = meshgrid(xRange,yRange);

for i = 1:size(tau_F_only_max,3)
    subplot(1,2,1)
    surface(xx, yy, tau_F_only_max(:,:,i)'/gearRatio)
    hold all
    zlim([0,1])
    title('max torque')

    subplot(1,2,2)
    surface(xx, yy, tau_F_only_min(:,:,i)'/gearRatio)
    hold all
    zlim([-1,0])
    title('min torque')

end


%%
figure(2); clf
for i = 1:size(omega_pulley_max,3)
    subplot(1,2,1)
    surface(xx, yy, omega_pulley_max(:,:,i)'/gearRatio)
    hold all
    title('max motor speed')

    subplot(1,2,2)
    surface(xx, yy, omega_pulley_min(:,:,i)'/gearRatio)
    hold all
    title('min motor speed')

end

