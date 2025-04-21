function J = FSquirrle_Jacobian(R1,R2,R3,A11,A12,A21,A22,A31,A32,C1_u1,C1_u2,C1_u3,C2_u1,C2_u2,C2_u3,C3_u1,C3_u2,C3_u3,r_p,p_s)

J_inv = FSquirrle_Jacobian_inverse(R1,R2,R3,A11,A12,A21,A22,A31,A32,C1_u1,C1_u2,C1_u3,C2_u1,C2_u2,C2_u3,C3_u1,C3_u2,C3_u3,r_p,p_s);
J = pinv(J_inv);