% From "Singularities of ABB’s YuMi 7-DOF robot arm"
clc
h_1 = [0;0;1];
h_4_0 = rand_normal_vec;
p_0W = rand_vec;


e_0 = h_1; % Unit vector along axis 1
w = p_0W; % Position vector of WCP: Intersection of axis 7 and common normal of axes 6 and 7
e_3 = h_4_0; % Unit vector along axis 4

n = cross(e_0, w) / norm(cross(e_0, w));

psi_A = atan2( ...
    sign(dot(e_3, e_0))*norm(cross(n, e_3)),...
    dot(n, e_3))
% norm([sign(dot(e_3, e_0))*norm(cross(n, e_3)); dot(n, e_3)]) % always 1
sign(dot(e_3, e_0))*acos(dot(e_3, n))
% SEW = sew_conv(h_1);
% psi_B = SEW.fwd_kin([0;0;0], h_4_0, p_0W)


% psi_B = subproblem.sp_1(e_3, w, e_0)
% rad2deg(psi_A - psi_B)
% rad2deg(psi_A + psi_B)


%% Try out with RobotStudio
% TODO: Arm angle is supposed to be affected by q_6
clc
% Single arm YuMi
kin = define_yumi;
kin.P(:,1) = [0;0;140];
kin.P(:,end) = [36;0;0];

% q = zeros([7 1]);


q = deg2rad(-40)*ones([7 1]);
% q(6) =deg2rad(138);
% N = 10000;
% psi_A_list = NaN(N, 1);
% psi_B_list = NaN(N, 1);
% 
% for i = 1:N

% q = rand_angle([7 1]);
% [R, T] = fwdkin(kin, q)
% R_6T = round(rot([0;1;0], pi/2));

[R, T, P_SEW] = fwdkin_inter(kin, q, [1 4 7]);


h_1 = [0;0;1];
h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
p_0W = P_SEW(:,3);

e_0 = h_1; % Unit vector along axis 1
w = p_0W; % Position vector of WCP: Intersection of axis 7 and common normal of axes 6 and 7
e_3 = h_4_0; % Unit vector along axis 4

n = cross(e_0, w) / norm(cross(e_0, w));

psi_A = atan2( ...
    sign(dot(e_3, e_0))*norm(cross(n, e_3)),...
    dot(n, e_3))
% psi_A_deg = rad2deg(psi_A)


SEW = yumi.sew_abb();
SEW.fwd_kin(kin, q)
%%


SEW = sew_conv([0;0;1]);
% psi_B = SEW.fwd_kin([0;0;0], h_4_0, P_SEW(:,3));
psi_B = SEW.fwd_kin([0;0;0], P_SEW(:,2), P_SEW(:,3));
psi_B_deg = rad2deg(psi_B)

% psi_A_list(i) = psi_A;
% psi_B_list(i) = psi_B;
% 
% end
% 
% plot(psi_B_list, psi_A_list,'.')
% axis image

%%
diagrams.setup(); hold on
diagrams.robot_plot(kin, q, auto_scale=true)
diagrams.arrow([0;0;0], P_SEW(:,3));
diagrams.redraw(); hold off


%% Double check Jacobian

e_r = rand_normal_vec;
p_0W = rand_vec;

delta = 1e-10;
direction = rand_normal_vec;
delta_vec = delta*direction;

(e_x_fwd(e_r, p_0W+delta_vec) - e_x_fwd(e_r, p_0W)) / delta

e_x = e_x_fwd(e_r, p_0W);
Jacobian_x_W(e_r, p_0W, e_x)*direction


function e_x = e_x_fwd(e_r, p_0W)
    e_x = cross(e_r, p_0W) / norm(cross(e_r, p_0W));
end

function J_x_W = Jacobian_x_W(e_r, p_0W, e_x)
    num = cross(e_r, e_x)* e_x';
    den = norm(cross(e_r, e_x)'* p_0W);
    J_x_W = num/den;
end

%%
q = rand_angle([7 1]);
q(3) = 0;
direction = rand([7 1]);
direction = direction / norm(direction);
delta = 1e-15;
delta_vec = delta*direction;

SEW = yumi.sew_abb();
(SEW.fwd_kin(kin, q+delta_vec) - SEW.fwd_kin(kin, q)) / delta
SEW.jacobian(kin, q) * direction

%%
clc
q = rand_angle([7 1]);
% q(3) = 0;
q(2) = pi
% q(3) = 1e-10;
SEW.fwd_kin(kin, q)
[J_psi, sign_term] = SEW.jacobian(kin, q)

%% TODO put into RS correctly (separate axis for 7)
% q = deg2rad([0
% -43.52
% -6.71
% 7.04
% 285
% 138
% 116
% ])

q = rand_angle([7 1])

q(2) = -1e-15

rad2deg(SEW.fwd_kin(kin, q))
[R, T] = fwdkin(kin, q)

%% Investigate q_A
clc
q = rand_angle([7 1]);
q(6)=0;

J = robotjacobian(kin,q);
[J_psi, sign_term_i] = SEW.jacobian(kin, q)
J_A = [J; J_psi];

det_i = det(J_A)
psi = SEW.fwd_kin(kin, q)
sign_term = sign_term_i