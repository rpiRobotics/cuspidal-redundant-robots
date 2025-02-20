
kin = define_yumi;
kin.P = kin.P / 100; % fix scaling for det(J)
SEW = yumi.sew_abb();

%% Random pose
for attempt = 1:100
q = rand_angle([7 1]);
[R, p] = fwdkin(kin, q);
psi = SEW.fwd_kin(kin, q);

% All IK solns
Q = yumi.IK_SEW_ABB_mex(R, p, SEW, psi, kin);
Q = unique_q_tol(Q, 1e-4);

% sgn(det(J)) for each soln
signs = NaN([1 width(Q)]);
for i = 1:numel(signs)
    J = robotjacobian(kin, Q(:,i));
    J_psi = SEW.jacobian(kin, Q(:,i));
    J_A = [J; J_psi];
    signs(i) = sign(det(J_A));
end

idx_pos = find(signs>0);
idx_neg = find(signs<0);

% Paths for all positive and all negative solutions
N = 100;
lambda = linspace(0, 1,  N);

[q_A_list_pos, q_B_list_pos] = generate_pairs(Q(:,idx_pos));
[q_A_list_neg, q_B_list_neg] = generate_pairs(Q(:,idx_neg));
q_A_list = [q_A_list_pos q_A_list_neg];
q_B_list = [q_B_list_pos q_B_list_neg];


det_path_mat = NaN(width(q_A_list),N);
for i_pair = 1:width(q_A_list)
    q_A = q_A_list(:,i_pair);
    q_B = q_B_list(:,i_pair);
    q_path = lambda.*q_B + (1-lambda).*q_A;
    for i = 1:N
        J = robotjacobian(kin, q_path(:,i));
        J_psi = SEW.jacobian(kin, q_path(:,i));
        J_A = [J; J_psi];
        det_path_mat(i_pair, i) = det(J_A);
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

%% Which was succesful?
all(det_path_mat'>1e-2)
all(det_path_mat'<-1e-2)

%% Path between 2 solns
q_A = q_A_list(:,1); % [ 2.6542    2.8648    0.3867    2.3326    0.8902    0.9457    2.3849]
q_B = q_B_list(:,1); % [-1.1841    2.6043    2.8249    1.8852    1.1459    0.8405    1.9855]


% 2.0137    3.1221   -3.0802   -0.2020    1.5407    2.9876    0.4203
% -1.8861    2.6673   -0.5224    0.2350    1.2883    2.8092    0.4334

N = 100;
lambda = linspace(0, 1,  N);
q_path = lambda.*q_B + (1-lambda).*q_A;
det_path = NaN(1,N);
for i = 1:N
    J = robotjacobian(kin, q_path(:,i));
    J_psi = SEW.jacobian(kin, q_path(:,i));
    J_A = [J; J_psi];
    det_path(i) = det(J_A);
end

plot(lambda, det_path)
xlabel("\lambda")
ylabel("det(J)")
yline(0)


%% Plot SEW angle over path
q_A = q_A_list(:,1);
q_B = q_B_list(:,1);

N = 100;
lambda = linspace(0, 1,  N);
q_path = lambda.*q_B + (1-lambda).*q_A;
psi_path = NaN(1,N);
for i = 1:N
    psi_path(i) = SEW.fwd_kin(kin, q_path(:,i));
end

plot(lambda, psi_path)
xlabel("\lambda")
ylabel("\psi")

%% Confirm q_A and q_B lead to same EE pose


[R_07_A, p_0T_A] = fwdkin(kin, q_A);
psi_A = SEW.fwd_kin(kin, q_A);
chi_A = [R_07_A(:); p_0T_A; psi_A];

[R_07_B, p_0T_B] = fwdkin(kin, q_B);
psi_B = SEW.fwd_kin(kin, q_B);
chi_B = [R_07_B(:); p_0T_B; psi_B];
chi_A-chi_B

%%
diagrams.setup(); hold on

opts = {'auto_scale', true, 'show_arrows', false, 'show_joint_labels', false, 'show_arrow_labels', false};


diagrams.robot_plot(kin, q_A, opts{:}, link_color=diagrams.colors.red);
diagrams.robot_plot(kin, q_B, opts{:});

diagrams.redraw(); hold off





function [q_A_list, q_B_list] = generate_pairs(Q)
    n = width(Q);
    
    if n < 2
        q_A_list = [];
        q_B_list = [];
        return;
    end
    
    q_A_list = [];
    q_B_list = [];
    
    for i = 1:n
        for j = i+1:n
            q_A_list = [q_A_list, Q(:,i)];
            q_B_list = [q_B_list, Q(:,j)];
        end
    end
end
