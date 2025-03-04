% Test cases from RobotStudio when e_r = e_y
% Controller > Configuration > Motion > Robot > Arm-Angle Reference Direction

% J is joint angle from RobotStudio where joint 3 <-> Joint 7 


q_1 = zeros([7 1]);
p_0T_1 = [341.5; 0; 598];
RPY_1 = deg2rad([0; 90; 0]);
psi_1_z = deg2rad(0);
psi_1_y = deg2rad(90);

q_2 = deg2rad(20)*ones([7 1]);
p_0T_2 = [292.46; 234.73; 348.75];
RPY_2 = deg2rad([148.67; 16.13; 175.45]);
psi_2_z = deg2rad(6.86);
psi_2_y = deg2rad(105.54);

q_3 = deg2rad(30)*ones([7 1]);
p_0T_3 = [165.49; 290.86; 255.37];
RPY_3 = deg2rad([157.34; -19.80; -175.99]);
psi_3_z = deg2rad(14.50);
psi_3_y = deg2rad(99.26);

q_4 = deg2rad(-40)*ones([7 1]);
p_0T_4 = [-219.17; 12.62; 796.12];
RPY_4 = deg2rad([13.67; -3.39; -152.29]);
psi_4_z = deg2rad(106.67);
psi_4_y = deg2rad(15.68);

q_5 = deg2rad(-80)*ones([7 1]);
p_0T_5 = [-126.60; 499.23; 441.88]; 
RPY_5 = deg2rad([5.30; 17.08; 34.27]);
psi_5_z = deg2rad(99.01);
psi_5_y = deg2rad(-26);

%% Compare p_0T and RPY
kin = define_yumi;
kin.P(:,end) = [36;0;0];

[R_t_1, p_t_1] = fwdkin(kin, q_1);
[R_t_2, p_t_2] = fwdkin(kin, q_2);
[R_t_3, p_t_3] = fwdkin(kin, q_3);
[R_t_4, p_t_4] = fwdkin(kin, q_4);
[R_t_5, p_t_5] = fwdkin(kin, q_5);

P_t = [p_t_1 p_t_2 p_t_3 p_t_4 p_t_5];
P_robotstudio = [p_0T_1 p_0T_2 p_0T_3 p_0T_4 p_0T_5];
round(P_t - P_robotstudio, 2)

%%
R_7T = round(rot([0;1;0], pi/2));
[RPY_1A, RPY_1B] = R_to_RPY(R_t_1*R_7T);
[RPY_2A, RPY_2B] = R_to_RPY(R_t_2*R_7T);
[RPY_3A, RPY_3B] = R_to_RPY(R_t_3*R_7T);
[RPY_4A, RPY_4B] = R_to_RPY(R_t_4*R_7T);
[RPY_5A, RPY_5B] = R_to_RPY(R_t_5*R_7T);

RPY_t = [RPY_1A RPY_2A RPY_3A RPY_4A RPY_5A]
RPY_robotstudio = [RPY_1 RPY_2 RPY_3 RPY_4 RPY_5]
round(rad2deg(RPY_t-RPY_robotstudio), 2)

%% Test RPY function
% R = rot(e_3, Y)*rot(e_2, P)*rot(e_1, Y)
e_1 = [1;0;0];
e_2 = [0;1;0];
e_3 = [0;0;1];

R = rot(rand_normal_vec, rand_angle);
[RPY_A, RPY_B] = R_to_RPY(R);
rot(e_3, RPY_A(3))*rot(e_2, RPY_A(2))*rot(e_1, RPY_A(1))-R
rot(e_3, RPY_B(3))*rot(e_2, RPY_B(2))*rot(e_1, RPY_B(1))-R
%% Compare SEW angles using e_r = e_z

kin = define_yumi;
SEW = yumi.sew_abb([0;0;1]);

psi_vec = [SEW.fwd_kin(kin, q_1)
SEW.fwd_kin(kin, q_2)
SEW.fwd_kin(kin, q_3)
SEW.fwd_kin(kin, q_4)
SEW.fwd_kin(kin, q_5)];

psi_vec_robotstudio = [psi_1_z
psi_2_z
psi_3_z
psi_4_z
psi_5_z];

rad2deg([psi_vec psi_vec_robotstudio])
round(rad2deg(psi_vec - psi_vec_robotstudio), 2)

%% Compare SEW angles using e_r = e_y

kin = define_yumi;
SEW = yumi.sew_abb([0;1;0]);

psi_vec = [SEW.fwd_kin(kin, q_1)
SEW.fwd_kin(kin, q_2)
SEW.fwd_kin(kin, q_3)
SEW.fwd_kin(kin, q_4)
SEW.fwd_kin(kin, q_5)];

psi_vec_robotstudio = [psi_1_y
psi_2_y
psi_3_y
psi_4_y
psi_5_y];

rad2deg([psi_vec psi_vec_robotstudio])
round(rad2deg(psi_vec - psi_vec_robotstudio), 2)

function [RPY_1, RPY_2] = R_to_RPY(R)
% R = rot(e_3, Y)*rot(e_2, P)*rot(e_1, R)
e_1 = [1;0;0];
e_2 = [0;1;0];
e_3 = [0;0;1];

% R*e_1 = rot(e_3, Y)*rot(e_2, P)*e_1

[Y_t, P_t] = subproblem.sp_2(R*e_1, e_1, -e_3, e_2);

Y_1 = Y_t(1); P_1 = P_t(1);
Y_2 = Y_t(end); P_2 = P_t(end);

% R*e_2 = rot(Y, e_3)*rot(e_2, P)*rot(e_1, R)*e_2
%  rot(e_1, R)*e_2 = (rot(e_3, Y)*rot(e_2, P))' * R*e_2
R_1 = subproblem.sp_1(e_2, (rot(e_3, Y_1)*rot(e_2, P_1))' * R*e_2, e_1);
R_2 = subproblem.sp_1(e_2, (rot(e_3, Y_2)*rot(e_2, P_2))' * R*e_2, e_1);

RPY_1 = [R_1; P_1; Y_1];
RPY_2 = [R_2; P_2; Y_2];
end