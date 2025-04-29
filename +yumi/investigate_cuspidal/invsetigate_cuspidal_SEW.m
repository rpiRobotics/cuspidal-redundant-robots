[kin, q_min, q_max] = define_yumi;
SEW = yumi.sew_conv_h4([0;0;1]);

%%

% Random pose

q = q_min + rand([7 1]).*(q_max - q_min);
[R, p] = fwdkin(kin, q);
psi = SEW.fwd_kin_q(q, kin);


% All IK solns
Q_nonfilter = yumi.IK_SEW_mex(R, p, SEW, psi, kin, false, 500)
Q = unique_q_tol(yumi.filter_Q_joint_limits([q Q_nonfilter], q_min, q_max, mode='remove'), deg2rad(0.1), "infinity")


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

% Paths for all positive and all negative solutions
N = 100;
lambda = linspace(0, 1,  N);

q_A_list = [];
q_B_list = [];

if numel(idx_pos) == 4
q_A_list = [Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_neg(1)) Q(:,idx_neg(1)) Q(:,idx_neg(1))];
q_B_list = [Q(:,idx_pos(2)) Q(:,idx_pos(3)) Q(:,idx_pos(4)) Q(:,idx_neg(2)) Q(:,idx_neg(3)) Q(:,idx_neg(4))];
elseif numel(idx_pos) == 3
q_A_list = [Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_neg(1)) Q(:,idx_neg(1))];
q_B_list = [Q(:,idx_pos(2)) Q(:,idx_pos(3)) Q(:,idx_neg(2)) Q(:,idx_neg(3))];
elseif numel(idx_pos) == 2
q_A_list = [Q(:,idx_pos(1))  Q(:,idx_neg(1))];
q_B_list = [Q(:,idx_pos(2))  Q(:,idx_neg(2))];
end


det_path_mat = NaN(width(q_A_list),N);
for i_pair = 1:width(q_A_list)
    q_A = q_A_list(:,i_pair);
    q_B = q_B_list(:,i_pair);
    q_path = lambda.*q_B + (1-lambda).*q_A;
    for i = 1:N
        J = robotjacobian(kin, q_path(:,i));
        J_psi = full_J_psi(kin, SEW, Q(:,i));
        J_aug = [J; J_psi];
        det_path_mat(i_pair, i) = det(J_aug);
    end
end

% disp(attempt)
if any(all(det_path_mat'>1e-2)) || any(all(det_path_mat'<-1e-2))
    disp("found it!")

    plot(lambda, det_path_mat')
    xlabel("\lambda")
    ylabel("det(J)")
    yline(0);
    % break
end



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