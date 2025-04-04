function plot_FSquirrleBases(R,C1_u, C2_u, C3_u, A1, A2, A3, hf)
%%
figure(hf)
hold all
axis equal

plot(A1(1),A1(2),'*')
plot(A2(1),A2(2),'*')
plot(A3(1),A3(2),'*')

plot3(R(1),R(2),R(3),'s')

plot3(R(1)+C1_u(1),R(2)+C1_u(2),R(3)+C1_u(3),'o')
plot3(R(1)+C2_u(1),R(2)+C2_u(2),R(3)+C2_u(3),'o')
plot3(R(1)+C3_u(1),R(2)+C3_u(2),R(3)+C3_u(3),'o')

plot3(R(1)+C1_u(1),R(2)+C1_u(2),0,'o')
plot3(R(1)+C2_u(1),R(2)+C2_u(2),0,'o')
plot3(R(1)+C3_u(1),R(2)+C3_u(2),0,'o')

%%
[B,B_u,B_l]= FSquirrle_bases(R(1),R(2),R(3),A1(1),A1(2),A2(1),A2(2),A3(1),A3(2),C1_u(1),C1_u(2),C1_u(3),C2_u(1),C2_u(2),C2_u(3),C3_u(1),C3_u(2),C3_u(3));

for i = 1:size(B_u,2)-1
    c = eval("C"+i+"_u");
    plot3(R(1)+c(1)+[0,B_u(1,i)],R(2)+c(2)+[0,B_u(2,i)],R(3)+c(3)+[0,B_u(3,i)])
    plot3(R(1)+c(1)+[0,B_l(1,i)],R(2)+c(2)+[0,B_l(2,i)],[0,B_l(3,i)])
end
plot3(R(1)+[0,B_l(1,end)],R(2)+[0,B_l(2,end)],R(3)+[0,B_l(3,end)])
%%

