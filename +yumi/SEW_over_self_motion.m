[kin, q_min, q_max] = define_yumi;

% Pick a fixed EE pose
% p = [300;0;500];
% R = eye(3);
q = rand_angle([7 1]);
[R, p] = fwdkin(kin, q);

% IK over all choices of q_1
N = 100;
q1_list = linspace(-pi, pi, N);
Q_path = NaN([6 16 N]);

for i = 1:N

    kin_i = fwdkin_partial(kin, q1_list(i), 1);

    Q = IK.IK_4_6_intersecting_mex(R, p, kin_i); %% TODO make sure to use R_7T
    if ~isempty(Q)
        Q_path(:, 1:width(Q), i) = Q;
    end

end
%
Q_path_filter = yumi.filter_Q_joint_limits(Q_path, q_min(2:7), q_max(2:7));
q2_list = squeeze(Q_path_filter(1,:,:)); % since q is q2 through q7
plot(q1_list, q2_list', '.')
xline(q_min(1));
xline(q_max(1));
yline(q_min(2));
yline(q_max(2));
%
% q2_list = squeeze(Q_path(1,:,:)); % since q is q2 through q7
% plot(q1_list, q2_list', '.')

% Extract one self-motion manifold

% [~, Q_path_idx] =  min(wrapTo2Pi(q2_list));
[~, Q_path_idx] =  min((q2_list));
Q_path_extract = NaN([6 N]);
for i = 1:N
    Q_path_extract(:,i) = Q_path_filter(:,Q_path_idx(i), i); %%%%%
end

q2_list_extract = Q_path_extract(1,:);
plot(q1_list, q2_list', 'k.'); hold on
plot(q1_list, q2_list_extract', 'xr'); hold off

% Calculate sign term and ABB SEW angle over path
Q7_path = NaN([7 N]);
Q7_path(1,:) = q1_list;
Q7_path(2:end, :) = Q_path_extract;

sign_term_path = NaN([1 N]);
psi_path = NaN([1 N]);
psi_conv_path = NaN([1 N]);

e_r = [0;0;1];
% e_r = [0;1;0];
% e_r = rand_normal_vec;

SEW = yumi.sew_abb(e_r);

SEW_conv = sew_conv(e_r);

for i = 1:N
    psi_path(i) = SEW.fwd_kin(kin, Q7_path(:,i));
    [~, sign_term_path(i)] = SEW.jacobian(kin, Q7_path(:,i));

    [~, ~, P_SEW] = fwdkin_inter(kin, Q7_path(:,i), [1 4 7]);
    psi_conv_path(i) = SEW_conv.fwd_kin(P_SEW(:,1), P_SEW(:,2), P_SEW(:,3));
end


%% Plot parameterizations over path

plot(q1_list, q1_list); hold on
plot(q1_list, q2_list_extract);
plot(q1_list, psi_path, 'k.', LineWidth=2);
plot(q1_list, sign_term_path, 'k:');
plot(q1_list, psi_conv_path, '.');
yline(0);
hold off

xline(q_min(1));
xline(q_max(1));

legend('q_1', 'q_2', '\psi^{ABB}', 'sign term', '\psi^{conv}','','', Location='southeast')
xlabel("\lambda")
ylim([-pi pi])
xlim([-pi, pi])
xline(q1_list(85));

%% Print out one of the joint angles
i_disp=75;
q_disp = Q7_path(:,i_disp);

q_ABB_deg = yumi.q2ABB(rad2deg(q_disp));
q_ABB_deg(1:6)
q_ABB_deg(7)
%%

[R_t, T_t] = fwdkin(kin, Q7_path(:,85))

% rad2deg(SEW.fwd_kin(kin, q_disp))
% rad2deg(psi_conv_path(:,i_disp))