clearvars

%% 
M = 10; % actual mass of the robot
A = [0 1 ; 0 0];
B = [0; 1/M];

M_ad = 1; % admitance mass of the virtual system
Kd_ad = 0.5; % damping in virtual system
A_ad = [0 1 ; 0 -Kd_ad];
B_ad = [0; 1/M_ad];

F_des = 10; % desired force we want to track

dt = 0.001;

Kp_hand = 1; % stiffness of virtual hand
Kd_hand = 0; % damping of virtual hand

hand_frequency = 2*pi; % frequency of hand movements

Q = 1e5*diag([1, 1]);
R = 1; 

K_lqr = lqr(A,B,Q,R);

tEnd = 100;

N = round(tEnd/dt);
IC = [0;2*pi];
X = zeros(length(IC),N+1);
X(:,1) = IC;
X_ad = X;

F_measured_log = zeros(N,1);

for i = 1:N
    t = (i-1)*dt;
    
    robot_p = X(1,i);
    robot_v = X(2,i);

    hand_p = sin(hand_frequency*t); % random-looking motion of the "user"
    hand_v = hand_frequency*cos(hand_frequency*t);
    
    F_int = ((hand_p-robot_p)*Kp_hand + (hand_v-robot_v)*Kd_hand); % force of hand on robot

    F_measured = F_int; % assume we measure force without knowing where it came from
    F_measured_log(i) = F_measured;

    Xdot_ad = A_ad*X_ad(:,i) + B_ad*(F_des+F_measured); % vritual admitance model's dynamics
    X_ad(:,i+1) = X_ad(:,i) + dt*Xdot_ad;

    u_lqr = -K_lqr*(X(:,i)-X_ad(:,i)); % tracking the motion of the virtual system. 
    
    Xdot = A*X(:,i) + B*u_lqr + B*F_measured; % actual system's dynamics

    X(:,i+1) = X(:,i) + dt*Xdot;


end

subplot(1,2,1)
plot(linspace(0,tEnd+dt,N+1),X')

subplot(1,2,2)
plot(linspace(0,tEnd,N), F_measured_log)