[kin, q_min, q_max] = define_yumi;
kin.P = kin.P/100;
SEW = yumi.sew_conv_h4([0;0;1]);

%%
for attempt = 1:1000

% random e_R
e_r = rand_normal_vec;
SEW = yumi.sew_conv_h4(e_r);

% Random pose
q = q_min + rand([7 1]).*(q_max - q_min);
% q = rand_angle([7 1]);
[R, p] = fwdkin(kin, q);
psi = SEW.fwd_kin_q(q, kin);


% All IK solns
Q_nonfilter = yumi.IK_SEW_mex(R, p, SEW, psi, kin, false, 700);
Q = unique_q_tol(yumi.filter_Q_joint_limits([q Q_nonfilter], q_min, q_max, mode='remove'), deg2rad(0.1), "infinity");
% Q = unique_q_tol(Q_nonfilter, deg2rad(0.1), "infinity");

% sgn(det(J)) for each soln
signs = NaN([1 width(Q)]);
for i = 1:numel(signs)
    J_aug = SEW.J_aug(Q(:,i), kin);
    signs(i) = sign(det(J_aug));
end

idx_pos = find(signs>0);
idx_neg = find(signs<0);

% Paths for all positive and all negative solutions
N = 100;
lambda = linspace(0, 1,  N);

q_A_list = [];
q_B_list = [];

if numel(idx_pos) > 4 || numel(idx_neg) > 4
    disp("more than 4!")
    break
end

if numel(idx_pos) == 4
q_A_list = [Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_pos(2)) Q(:,idx_pos(2)) Q(:,idx_pos(3))];
q_B_list = [Q(:,idx_pos(2)) Q(:,idx_pos(3)) Q(:,idx_pos(4)) Q(:,idx_pos(3)) Q(:,idx_pos(4)) Q(:,idx_pos(4))];
elseif numel(idx_pos) == 3
q_A_list = [Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_pos(2)) ];
q_B_list = [Q(:,idx_pos(2)) Q(:,idx_pos(3)) Q(:,idx_pos(3))];
elseif numel(idx_pos) == 2
q_A_list = [Q(:,idx_pos(1))];
q_B_list = [Q(:,idx_pos(2))];
end

if numel(idx_neg) == 4
q_A_list = [q_A_list Q(:,idx_neg(1)) Q(:,idx_neg(1)) Q(:,idx_neg(1)) Q(:,idx_neg(2)) Q(:,idx_neg(2)) Q(:,idx_neg(3))];
q_B_list = [q_B_list Q(:,idx_neg(2)) Q(:,idx_neg(3)) Q(:,idx_neg(4)) Q(:,idx_neg(3)) Q(:,idx_neg(4)) Q(:,idx_neg(4))];
elseif numel(idx_neg) == 3
q_A_list = [q_A_list Q(:,idx_neg(1)) Q(:,idx_neg(1)) Q(:,idx_neg(2)) ];
q_B_list = [q_B_list Q(:,idx_neg(2)) Q(:,idx_neg(3)) Q(:,idx_neg(3))];
elseif numel(idx_neg) == 2
q_A_list = [q_A_list Q(:,idx_neg(1))];
q_B_list = [q_B_list Q(:,idx_neg(2))];
end



det_path_mat = NaN(width(q_A_list),N);
for i_pair = 1:width(q_A_list)
    q_A = q_A_list(:,i_pair);
    q_B = q_B_list(:,i_pair);
    q_path = lambda.*q_B + (1-lambda).*q_A;
    for i = 1:N
        J_aug = SEW.J_aug(q_path(:,i), kin);
        det_path_mat(i_pair, i) = det(J_aug);
    end
end

disp(attempt)
if any(all(det_path_mat'>1e-2)) || any(all(det_path_mat'<-1e-2))
    disp("found it!")

    plot(lambda, det_path_mat')
    xlabel("\lambda")
    ylabel("det(J)")
    yline(0);
    break
end

end

%%

all_pos = all(det_path_mat > 1e-2, 2);
all_neg = all(det_path_mat < -1e-2, 2);
success_idx = find(all_pos | all_neg)
plot(det_path_mat(success_idx,:));
yline(0);
q_A = q_A_list(:,success_idx)
q_B = q_B_list(:,success_idx)

%% Confirm same FK
[R_A, p_A] = fwdkin(kin, q_A)
psi_A = SEW.fwd_kin_q(q_A, kin)
[R_B, p_B] = fwdkin(kin, q_B)
psi_B = SEW.fwd_kin_q(q_B, kin)

[R_A p_A] - [R_B p_B]
psi_A - psi_B
