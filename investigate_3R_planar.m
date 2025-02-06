zv = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

kin.H = [ez ez ez];
kin.P = [zv ex ex 0.25*ex];
kin.joint_type = [0 0 0];

%%
syms x y z

[~, P] = fwdkin(kin, [x; y; z])
% P = simplify(P)

%%
E = P - 1.2*ex;
E = E(1:2);

solve(E)
%%
fimplicit(subs(E, z, 0),[-pi, pi])
%%
% q = [0;0;0]
q = [0.2; 0.2; 0.2]

[R,P,P_inter] = fwdkin_inter(kin, q, [1 2 3]);
P_inter = [P_inter P]
plot(P_inter(1,:), P_inter(2,:), '-o'); hold on
plot(1.2*ex(1), 1.2*ex(2), 'x');
hold off
axis equal

%% Consider 2R + distance constraint

zv = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

kin.H = [ez ez];
kin.P = [zv ex ex];
kin.joint_type = [0 0];

syms x y real

[~, P] = fwdkin(kin, [x; y])
d_3T = P-1.2*ex
E = d_3T'*d_3T == 0.25^2

fimplicit(E, [-pi, pi])


%% Variable EE

zv = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

kin.H = [ez ez];
kin.P = [zv 3*ex ex];
kin.joint_type = [0 0];
L_3T = 0.25

syms q1 q2 X real

[~, P] = fwdkin(kin, [q1; q2])
d_3T = P-X*ex
E = d_3T'*d_3T == L_3T^2

syms x y z
fimplicit3(subs(E,[q1 q2 X], [z x y]), [-pi pi  0 5 -pi pi], MeshDensity=85, LineWidth=0.1)
zlabel("q_1")
xlabel("q_2")
ylabel("\rho")
view(107, 13)