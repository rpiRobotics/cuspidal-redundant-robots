kin = define_yumi;
SEW = sew_conv([0;0;1]);

%%

% Random pose

q = rand_angle([7 1]);
[R, p, P_SEW] = fwdkin_inter(kin, q, [1 4 7]); % TODO which joints?
psi = SEW.fwd_kin(P_SEW(:,1),P_SEW(:,2),P_SEW(:,3));

%%
% All IK solns
Q = SEW_IK.IK_gen_7_dof_mex(R, p, SEW, psi, kin)

%%

% sgn(det(J)) for each soln
signs = NaN([1 width(Q)]);
for i = 1:numel(signs)
    J = robotjacobian(kin, Q(:,i));
    J_psi = full_J_psi(kin, SEW, Q(:,i));
    J_aug = [J; J_psi];
    signs(i) = sign(det(J_aug));
end

idx_pos = find(signs>0);
idx_neg = find(signs<0);

% Path between two random poses with same sign(det(J))
q_A = Q(:,idx_pos(randperm(numel(idx_pos), 1)));
q_B = Q(:,idx_pos(randperm(numel(idx_pos), 1)));


N = 100;
lambda = linspace(0, 1,  N);
q_path = lambda.*q_B + (1-lambda).*q_A;
det_path = NaN(1,N);

for i = 1:N
    J = robotjacobian(kin, q_path(:,i));
    J_psi = full_J_psi(kin, SEW, q_path(:,i));
    J_aug = [J; J_psi];
    det_path(i) = det(J_aug);
end

plot(lambda, det_path)
yline(0);




function J_psi = full_J_psi(kin, SEW, q)
kin_E = kin;
kin_E.joint_type = [0 0 0];
kin_W = kin;
kin_W.joint_type = [0 0 0 0 0 0];

    [~, ~, P_SEW] = fwdkin_inter(kin,q, [1 4 7]);
    [J_psi_E, J_psi_W] = SEW.jacobian(P_SEW(:,1), P_SEW(:,2), P_SEW(:,3));
    
    J_E = robotjacobian(kin_E, q);
    J_E = [J_E(4:6, :) zeros(3, 4)];
    
    J_W = robotjacobian(kin_W, q);
    J_W = [J_W(4:6, :) zeros(3,1)];
    
    J_psi = J_psi_E*J_E + J_psi_W*J_W;
end