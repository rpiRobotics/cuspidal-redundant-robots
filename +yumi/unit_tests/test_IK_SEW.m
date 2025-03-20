for i = 1:1e4

[kin, q_min, q_max] = define_yumi;
% e_r = rand_normal_vec;
e_r = [0;0;1];
SEW = yumi.sew_conv_h4(e_r);

q = rand_angle([7 1]);
% q = (1:7)'/10
% q = q_min + rand([7 1]).*(q_max - q_min);

[R_07, p_0T] = fwdkin(kin, q);


psi = SEW.fwd_kin_q(q, kin);

[Q, is_LS_vec] = yumi.IK_SEW_mex(R_07, p_0T, SEW, psi, kin);
% Q_filter = unique_q_tol(yumi.filter_Q_joint_limits(Q, q_min, q_max, mode='remove'), 1e-4);
Q_filter = unique_q_tol(Q, 1e-4);
if width(Q_filter) > 15
    beep;
    disp("found it")
    Q_filter
    break
end

disp(i)
end
%%
[Q, is_LS_vec] = yumi.IK_SEW(R_07, p_0T, SEW, psi, kin)
%%
codegen +yumi/IK_SEW.m -args {R_07, p_0T, SEW, psi, kin}
%%
[Q, is_LS_vec] = yumi.IK_SEW_mex(R_07, p_0T, SEW, psi, kin)
%%
chi = [R_07(:); p_0T; psi];
chi_Q = NaN(13, width(Q_filter));

for i = 1:width(Q_filter)
    q_i = Q_filter(:,i);
    [R_07_i, p_0T_i] = fwdkin(kin, q_i);
    psi_i = SEW.fwd_kin_q(q_i, kin);
    chi_i = [R_07_i(:); p_0T_i; psi_i];
    chi_Q(:,i) = chi_i;
end
disp(chi_Q-chi)

%%
unique_q_tol(Q)

%%
Q_disp = Q_filter;

diagrams.setup(); hold on

opts = {'auto_scale', true, 'show_arrows', false, 'show_joint_labels', false, 'show_arrow_labels', false};

for i = 1:width(Q_disp)
    diagrams.robot_plot(kin, Q_disp(:,i), opts{:});
end

diagrams.redraw(); hold off

