%%
clearvars

v = 1; % desired speed
a = 2; % desired acceleration
F = 20; % magnitude of the desired force

T = 0; % external rotating moment on robot;
Omega = 0; % rotational velcoity of the robot

phaseshift = pi; % the angle between robot's acceleration and hand force. maximum torques at pi

D = 0.15;
r = 0.04;
theta = 0;

J = HAMSTER_Jacobian(D,r,theta);

m = 5;
I = 0.4;
M = diag([m,m,I]);

gearRatio = 1;

tau_F_only = [];
tau_a_only = [];
omega_wheel = [];

phiRange = 0:0.1:2*pi;

for phi = phiRange
    Fx = F*cos(phi);
    Fy = F*sin(phi);
    
    handForce = -[Fx;Fy;T];
    tau_F_only = [tau_F_only, J' * handForce];
    
    ax = a*cos(phi+phaseshift);
    ay = a*sin(phi+phaseshift);
    robotAcc = [ax;ay;0];
    tau_a_only = [tau_a_only, J' * (M*robotAcc)];

    vx = v*cos(phi);
    vy = v*sin(phi);

    V = [vx;vy;Omega];
    omega_wheel = [omega_wheel, J\V];
    

end

tau_combined = tau_a_only + tau_F_only;

figure(1); clf
subplot(3,1,1)
plot(rad2deg(phiRange), tau_F_only/gearRatio)
title('motor torque for hand force only')

subplot(3,1,2)
plot(rad2deg(phiRange), tau_a_only/gearRatio)
title('motor torque for linear acceleration only')

subplot(3,1,3)
plot(rad2deg(phiRange), tau_combined/gearRatio)
title('combined effect of acceleration and force')

figure(2); clf
plot(rad2deg(phiRange), omega_wheel*gearRatio)
title('motor speed')