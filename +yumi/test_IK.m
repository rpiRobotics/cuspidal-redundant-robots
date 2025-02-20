kin = define_yumi;
kin.P(:,1) = [0;0;140];
kin.P(:,end) = [36;0;0];
SEW = yumi.sew_abb();

% q = rand_angle([7 1]);
q = (1:7)'/10

[R_07, p_0T] = fwdkin(kin, q)
psi = SEW.fwd_kin(kin, q)
%% 
[Q, is_LS_vec] = yumi.IK_SEW_ABB(R_07, p_0T, SEW, psi, kin)
%%
codegen +yumi/IK_SEW_ABB.m -args {R_07, p_0T, SEW, psi, kin}
%%
[Q, is_LS_vec] = yumi.IK_SEW_ABB_mex(R_07, p_0T, SEW, psi, kin)
%%
chi = [R_07(:); p_0T; psi];
chi_Q = NaN(13, width(Q));

for i = 1:width(Q)
    q_i = Q(:,i);
    [R_07_i, p_0T_i] = fwdkin(kin, q_i);
    psi_i = SEW.fwd_kin(kin, q_i);
    chi_i = [R_07_i(:); p_0T_i; psi_i];
    chi_Q(:,i) = chi_i;
end
disp(chi_Q-chi)

%%
unique_q_tol(Q)

%%
diagrams.setup(); hold on

opts = {'auto_scale', true, 'show_arrows', false, 'show_joint_labels', false, 'show_arrow_labels', false};

for i = 1:width(Q)
    diagrams.robot_plot(kin, Q(:,i), opts{:});
end

diagrams.redraw(); hold off

