function J = HAMSTER_Jacobian(D,r,theta)
%HAMSTER_Jacobian
%    J = HAMSTER_Jacobian(D,R,THETA)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    01-Apr-2025 14:06:45

t2 = 1.0./D;
t3 = pi./3.0;
t4 = t3+theta;
t5 = (r.*t2)./3.0;
J = reshape([r.*cos(theta).*(-2.0./3.0),r.*sin(theta).*(-2.0./3.0),t5,r.*cos(-t3+theta).*(2.0./3.0),r.*cos(theta+pi./6.0).*(-2.0./3.0),t5,r.*cos(t4).*(2.0./3.0),r.*sin(t4).*(2.0./3.0),t5],[3,3]);
end
