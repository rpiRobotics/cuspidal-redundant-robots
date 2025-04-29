kin  = define_yumi;
SEW = yumi.sew_conv_h4([0;0;1]);
%% Baseline

q = rand_angle([7 1]);
J_A = SEW.J_aug(q, kin);
det(J_A)
%% q_2 singularity
% Joints 1 and 3 become collinear

q = rand_angle([7 1]);
q(2) = 0;
J_A = SEW.J_aug(q, kin);
det(J_A)

%% q_6 singularity
% Joints 5 and 7 become collinear

q = rand_angle([7 1]);
q(6) = 0;
J_A = SEW.J_aug(q, kin);
det(J_A)

%% Investigate what self-motion looks like using 6x7 matrix
J = robotjacobian(kin, q)
null(J)