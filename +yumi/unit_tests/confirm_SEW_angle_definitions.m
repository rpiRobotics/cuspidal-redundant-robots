%% Make sure definitions given for SEW angle are equivalent

kin = define_yumi;
e_r = rand_normal_vec;
SEW = yumi.sew_conv_h4(e_r);
q = rand_angle([7 1]);

% Conventional SEW angle (using h_4 as shoulder-elbow direction)
psi_1 = SEW.fwd_kin_q(q, kin);

[~, ~, P_SEW] = fwdkin_inter(kin,q, [1 4 7]);
h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
p_0S = P_SEW(:,1);
p_0W = P_SEW(:,3);
p_SW = p_0W - p_0S;
e_SW = p_SW  / norm(p_SW);

% ATAN2 definition for psi_conv
psi_2 = atan2(cross(e_SW,e_r)' * h_4_0, -cross(e_SW,cross(e_SW,e_r))'*h_4_0);

% ATAN2 definition for psi_ABB
psi_3 = atan2(-cross(e_SW,cross(e_SW,e_r))'*h_4_0, -cross(e_SW,e_r)' * h_4_0);
psi_4 = atan2(-e_r'*hat(e_SW)^2*h_4_0, e_r'*hat(e_SW) * h_4_0);

% These should all be the same
wrapToPi([psi_1+pi/2; psi_2+pi/2; psi_3; psi_4])