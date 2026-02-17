clearvars


syms R [3,1] % robot's absolute position
syms dR [3,1] % robot's differential position change
syms A1 A2 A3 [3,1] % three anchor points
syms C1_u C2_u C3_u [3,1] % three attachment points to the robot's upper section. These are relative to R

syms t
syms r_p p_s % pully radius and screw pitch
%% This section calculates the Jacobian (actually its inverse)


%neglecting dimensions of the pullies at anchor points.

[A1(3), A2(3), A3(3)] = deal(0); % assuming anchors' heigh are global zero height 

C1_u_global = C1_u + R;
C2_u_global = C2_u + R;
C3_u_global = C3_u + R

C1_l_global  = C1_u_global .* [1;1;0]; % assuming the lower section is level with the height of anchor poitns
C2_l_global  = C2_u_global .* [1;1;0];
C3_l_global  = C3_u_global .* [1;1;0]

l1_u = A1 - C1_u_global;
l2_u = A2 - C2_u_global;
l3_u = A3 - C3_u_global

l1_l = A1 - C1_l_global;
l2_l = A2 - C2_l_global;
l3_l = A3 - C3_l_global


l1_u_dir = l1_u./sum(l1_u.^2)^0.5;
l2_u_dir = l2_u./sum(l2_u.^2)^0.5;
l3_u_dir = l3_u./sum(l3_u.^2)^0.5

l1_l_dir = l1_l./sum(l1_l.^2)^0.5;
l2_l_dir = l2_l./sum(l2_l.^2)^0.5;
l3_l_dir = l3_l./sum(l3_l.^2)^0.5

dl1_u = sum(l1_u_dir.*dR);
dl2_u = sum(l2_u_dir.*dR);
dl3_u = sum(l3_u_dir.*dR)

dl1_l = sum(l1_l_dir.*dR);
dl2_l = sum(l2_l_dir.*dR);
dl3_l = sum(l3_l_dir.*dR)

dl1 = dl1_l + dl1_u;
dl2 = dl2_l + dl2_u;
dl3 = dl3_l + dl3_u

dl = [dl1;dl2;dl3];

dq = dl/r_p;
dq(end+1) = dR(3)/(p_s/2/pi);

J_inv = simplify(equationsToMatrix(dq, dR))

% symvar(J_inv);

matlabFunction(J_inv, Vars=[R1,R2,R3,A11,A12,A21,A22,A31,A32,C1_u1,C1_u2,C1_u3,C2_u1,C2_u2,C2_u3,C3_u1,C3_u2,C3_u3,r_p,p_s],...
    Optimize=true,File='FSquirrle_Jacobian_inverse')

%%
%% This section calculates the basis vectors

B_u = [l1_u_dir, l2_u_dir, l3_u_dir, zeros(3,1)]
B_l = [l1_l_dir, l2_l_dir, l3_l_dir, [0;0;1]]

B = B_l+B_u;

matlabFunction(B,B_u,B_l, Vars=[R1,R2,R3,A11,A12,A21,A22,A31,A32,C1_u1,C1_u2,C1_u3,C2_u1,C2_u2,C2_u3,C3_u1,C3_u2,C3_u3],...
    Optimize=true,File='FSquirrle_bases')