clearvars

%% 
M = 5; % actual mass of the robot
A = [0 1 ; 0 0];
B = [0; 1/M];

M_ad = 0.5; % admitance mass of the virtual system
Kd_ad = 0; % damping in virtual system
A_ad = [0 1 ; 0 -Kd_ad];
B_ad = [0; 1/M_ad];

F_des = 0; % desired force we want to track

dt = 0.001;

Kp_hand = 50; % stiffness of virtual hand
Kd_hand = 50; % damping of virtual hand

hand_frequency = .5*pi; % frequency of hand movements
hand_amp = 0.3; % Amplitude of the hand motion

Q = 1e3*diag([1, 1]);
R = 1; 

K_lqr = lqr(A,B,Q,R);

gain_p = 10;
gain_d = 10;

tEnd = 10;


robotNoise = 0.001;
forceNoise = 0.05;

delay = round(0.01/dt);


N = round(tEnd/dt);
IC = [0;hand_frequency*hand_amp];
X = zeros(length(IC),N+1);
X(:,1) = IC;
X_ad = X;

F_measured_log = zeros(N,1);

hand_log = zeros(N,2);
robot_log = zeros(N,2);
model_log = zeros(N,2);



for i = 1:N
    t = (i-1)*dt;

    hand_p = hand_amp*sin(hand_frequency*t); % random-looking motion of the "user"
    hand_v = hand_amp*hand_frequency*cos(hand_frequency*t);

    robot_p = X(1,i) + robotNoise*randn; % robot's position
    robot_v = X(2,i) + robotNoise*randn; %robot's velocity



    model_p = X_ad(1,i); % virtual system's position
    model_v = X_ad(2,i); % virtual system's velocity
    % model_p = hand_p;
    % model_v = hand_v;
    
    
    F_int = ((hand_p-robot_p)*Kp_hand + (hand_v-robot_v)*Kd_hand); % force from environment on robot

    F_measured = F_int + forceNoise*randn; % assume we measure force without knowing where it came from

    F_measured_log(i) = F_measured;
    hand_log(i,:) = [hand_p, hand_v];
    robot_log(i,:) = [robot_p, robot_v];
    model_log(i,:) = [model_p, model_v];


    if i<=delay
        F_feedback = 0;
        robot_feedback = [0;0];
    else
        F_feedback = F_measured_log(i-delay);
        robot_feedback = robot_log(i-delay,:)';
    end
        
    Xdot_ad = A_ad*X_ad(:,i) + B_ad*(F_des+F_feedback); % vritual admitance model's dynamics. calculated using feedback force
    X_ad(:,i+1) = X_ad(:,i) + dt*Xdot_ad;
    


    u_pid = gain_p*(model_p-robot_feedback(1)) + gain_d*(model_v-robot_feedback(2));

    u_lqr = -K_lqr*(robot_feedback-X_ad(:,i)); % tracking the motion of the virtual system. 

    u_k = -12 * (F_des - F_feedback);

    Xdot = A*X(:,i) + B*(u_pid+F_measured); % physical system's dynamics. calculated using real (current) force value


    X(:,i+1) = X(:,i) + dt*Xdot;


end

figure(1); clf

subplot(1,2,1)
plot(linspace(0,tEnd+dt,N),[robot_log(:,1),model_log(:,1),hand_log(:,1)], LineWidth=2)
title('states')
legend('robot x','model x','hand x')


subplot(1,2,2)
plot(linspace(0,tEnd,N), F_measured_log)
ylabel('f measured')