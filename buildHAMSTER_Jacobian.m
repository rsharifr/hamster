%% 
% 

%%
clearvars

syms('phi', [3,1]);
syms theta % angle of each wheel relative to global X
syms('omega',[3,1]);
syms v_x v_y Omega
syms r D % r: radius of wheel. D: distance from robot center to wheel


phi1 = pi+theta;
phi2 = 5*pi/3 + theta;
phi3 = pi/3 + theta;

omegaEqs = 1/r * (v_x*cos(phi) + v_y*sin(phi) + D*Omega) ;
omegaEqs = subs(omegaEqs)

J_inv = equationsToMatrix(omegaEqs, [v_x; v_y; Omega]);

J_inv = simplify(J_inv)

J = simplify(inv(J_inv))

matlabFunction(J,Optimize=true,File='HAMSTER_Jacobian')