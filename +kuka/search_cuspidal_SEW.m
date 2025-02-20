kin = robot_kin.kuka;
% SEW = sew_conv([0;0;1]);
SEW = sew_conv([1;0;0]);
kin.P(:,end) = 0;
%%

% Random pose
q = rand_angle([7 1]); 
[R, p, P_SEW] = fwdkin_inter(kin, q, [1 3 5]);
psi = SEW.fwd_kin(P_SEW(:,1), P_SEW(:,2), P_SEW(:,3));
% All IK solns
% Q = SEW_IK.IK_2R_2R_3R(R, p, SEW, psi, kin)
Q_A = SEW_IK.IK_2R_2R_3R(R, p, SEW, psi, kin);
% Q_B = SEW_IK.IK_2R_2R_3R(R, p, SEW, psi+rand_angle(), kin);
Q_B = SEW_IK.IK_2R_2R_3R(R, p, SEW, psi, kin); % Same SEW angle

%

% sgn(det(J)) for each soln
signs_A = NaN([1 width(Q_A)]);
for i = 1:numel(signs_A)
    J = robotjacobian(kin, Q_A(:,i));
    J_psi = full_J_psi(kin, SEW, Q_A(:,i));
    signs_A(i) = sign(det([J; J_psi]));
end
signs_B = NaN([1 width(Q_B)]);
for i = 1:numel(signs_B)
    J = robotjacobian(kin, Q_B(:,i));
    J_psi = full_J_psi(kin, SEW, Q_B(:,i));
    signs_B(i) = sign(det([J; J_psi]));
end


idx_pos_A = find(signs_A>0);
idx_neg_A = find(signs_A<0);
idx_pos_B = find(signs_B>0);
idx_neg_B = find(signs_B<0);
signs_A
signs_B

% Path between 2 solns
q_A = Q_A(:,idx_pos_A(randperm(numel(idx_pos_A), 1)));
q_B = Q_B(:,idx_pos_B(randperm(numel(idx_pos_B), 1)));

% Connect *ANY* two poses. Sign of det(J_i) shouldn't matter
% q_A = Q_A(:, randi(width(Q_A)));
% q_B = Q_B(:, randi(width(Q_B)));


N = 100;
lambda = linspace(0, 1,  N);
q_path = lambda.*q_B + (1-lambda).*q_A;
det_path_1 = NaN(1,N);
det_path_2 = NaN(1,N);
det_path_3 = NaN(1,N);
det_path_4 = NaN(1,N);
det_path_5 = NaN(1,N);
det_path_6 = NaN(1,N);
det_path_7 = NaN(1,N);
det_path_aug = NaN(1,N);



for i = 1:N
    J = robotjacobian(kin, q_path(:,i));
    J_psi = full_J_psi(kin, SEW, q_path(:,i));


    det_path_1(i) = det(J(:, [  2 3 4 5 6 7]));
    det_path_2(i) = det(J(:, [1   3 4 5 6 7]));
    det_path_3(i) = det(J(:, [1 2   4 5 6 7]));
    det_path_4(i) = det(J(:, [1 2 3   5 6 7]));
    det_path_5(i) = det(J(:, [1 2 3 4   6 7]));
    det_path_6(i) = det(J(:, [1 2 3 4 5   7]));
    det_path_7(i) = det(J(:, [1 2 3 4 5 6  ]));
    det_path_aug(i) = det([J; J_psi]);
    

    % J = robotjacobian(kin_7, [q_path(1:2,i); pi/2 ;q_path(3:end,i)]);
    % det_path(i) = min((svd(J)));
end

plot(lambda, det_path_1); hold on
plot(lambda, det_path_2); 
plot(lambda, det_path_3); 
plot(lambda, det_path_4); 
plot(lambda, det_path_5); 
plot(lambda, det_path_6); 
plot(lambda, det_path_7); 
plot(lambda, det_path_aug, 'k', LineWidth=3)
hold off
xlabel("\lambda")
ylabel("det(J_i)")
yline(0)

legend([string(1:7) "aug"])
%% Double check FK
[R, p, P_SEW] = fwdkin_inter(kin, q_A, [1 3 5])
[R, p, P_SEW] = fwdkin_inter(kin, q_B, [1 3 5])

%% Paths for all positive and all negative solutions
N = 100;
lambda = linspace(0, 1,  N);

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
        det_path_mat(i_pair, i) = det(J(:, [1 2 4 5 6 7]));
    end
end


plot(lambda, det_path_mat')
xlabel("\lambda")
ylabel("det(J)")
yline(0);

%% Iterate through psi = q_3

N_1 =50;
N_2 = 170;
q3_path_1 = linspace(pi/2, pi, 50+1);
q3_path_1 = q3_path_1(1:end-1);
q3_path_2 = linspace(pi/2, -pi,N_2+1);
q3_path_2 = q3_path_2(1:end-1);
Q_path_1 = NaN([6, N_1]);
Q_path_1(:,1) = q_B;
Q_path_2 = NaN([6, N_2]);
Q_path_2(:,1) = q_B;
det_path_1 = NaN(1,N_1);
det_path_2 = NaN(1,N_2);

for i = 2:N_1
    disp(i);
    kin_i = fwdkin_partial(kin_7, q3_path_1(i), 3);
    Q_i =  IK.IK_4_6_intersecting(R, p, kin_i);
    Q_path_1(:,i) = closest_q(Q_i, Q_path_1(:,i-1));

    % J = robotjacobian(kin_i, Q_path(:,i));
    % det_path(i) = det(J);
    J = robotjacobian(kin_7, [Q_path_1(1:2,i); q3_path_1(i) ;Q_path_1(3:end,i)]);
    det_path_1(i) = min((svd(J)));
end

for i = 2:N_2
    disp(i);
    kin_i = fwdkin_partial(kin_7, q3_path_2(i), 3);
    Q_i =  IK.IK_4_6_intersecting(R, p, kin_i);
    Q_path_2(:,i) = closest_q(Q_i, Q_path_2(:,i-1));

    % J = robotjacobian(kin_i, Q_path(:,i));
    % det_path(i) = det(J);
    J = robotjacobian(kin_7, [Q_path_2(1:2,i); q3_path_2(i) ;Q_path_2(3:end,i)]);
    det_path_2(i) = min((svd(J)));
end
%%
q3_path = [q3_path_1 q3_path_2];
Q_path = [Q_path_1, Q_path_2];
plot(q3_path, Q_path', 'x')

%%
plot(q3_path, [det_path_1, det_path_2])
xlabel("\lambda")
ylabel("det(J)")

%% what is sign(det(J)) when we fix each joint?

J = robotjacobian(kin_7, [q_A(1:2); pi/2 ; q_A(3:end)]);
sign_vec = NaN([7 1]);
for i = 1:7
    % sign_vec(i) = sign(det(J(:, ~(1:7 ==i))));
    sign_vec(i) = (det(J(:, ~(1:7 ==i))));
end
sign_vec

%%
kin_i = fwdkin_partial(kin_7, -pi+0.01, 3);
    Q_i =  IK.IK_4_6_intersecting(R, p, kin_i);

%% Fully plot of all IK solutions showing all self-motion manifolds

N = 1000;
psi_path = linspace(-pi, pi, N);
Q_SM = [];

for i = 1:N
    disp(i);
    Q_i =  SEW_IK.IK_2R_2R_3R(R, p, SEW, psi_path(i), kin);
    Q_SM = [Q_SM Q_i];
end
%%
scatter( (Q_SM(1, :)), (Q_SM(5, :)), [], Q_SM(3,:), '.');
% hold on
% plot(q_path(1,:), q_path(5,:), 'k', LineWidth=3);
% hold off
colormap hsv
xlabel("q_1")
ylabel("q_5")
axis tight




function J_psi = full_J_psi(kin, SEW, q)
kin_E = kin;
kin_E.joint_type = [0 0 0];
kin_W = kin;
kin_W.joint_type = [0 0 0 0 0];

    [~, ~, P_SEW] = fwdkin_inter(kin,q, [1 3 5]);
    [J_psi_E, J_psi_W] = SEW.jacobian(P_SEW(:,1), P_SEW(:,2), P_SEW(:,3));
    
    J_E = robotjacobian(kin_E, q);
    J_E = [J_E(4:6, :) zeros(3, 4)];
    
    J_W = robotjacobian(kin_W, q);
    J_W = [J_W(4:6, :) zeros(3,2)];
    
    J_psi = J_psi_E*J_E + J_psi_W*J_W;
end