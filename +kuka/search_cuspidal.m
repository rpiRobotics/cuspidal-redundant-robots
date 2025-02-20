%% Lock joint 3
kin_7 = robot_kin.kuka;
q_3 = 0;
kin = fwdkin_partial(kin_7, q_3, 3);
% 
[R, p] = fwdkin(kin, rand_angle([6,1]))

Q = IK.IK_spherical_2_intersecting(R, p, kin)

[R_t, p_t] = fwdkin(kin, Q(:,1))
%% Lock joint 1
% Move joint origins so p_12 =0
kin_7 = robot_kin.kuka;
kin_7.P(:,4) = kin_7.P(:,3);
kin_7.P(:,3) = 0;

q_1 = pi/3;
kin = fwdkin_partial(kin_7, q_1, 1);
% 
[R, p] = fwdkin(kin, rand_angle([6,1]))

Q = IK.IK_spherical_2_intersecting(R, p, kin)

[R_t, p_t] = fwdkin(kin, Q(:,1))
%%
% Random pose

q = rand_angle([6 1]); 
[R, p] = fwdkin(kin, q);
% All IK solns
Q = IK.IK_spherical_2_intersecting(R, p, kin)
%

% sgn(det(J)) for each soln
signs = NaN([1 width(Q)]);
for i = 1:numel(signs)
    J = robotjacobian(kin_7, [Q(1:2,i); q_3; Q(3:6,i)]);
    J_psi = [0 0 1 0 0 0 0];
    % J = robotjacobian(kin_7, [q_1; Q(:,i)]);
    % J_psi = [1 0 0 0 0 0 0];
    J_aug = [J; J_psi]; 
    signs(i) = sign(det(J_aug));
end

idx_pos = find(signs>0);
idx_neg = find(signs<0);
signs

% %% Path between 2 solns
% q_A = Q(:,idx_pos(1));
% q_B = Q(:,idx_pos(2));
% 
% 
% N = 100;
% lambda = linspace(0, 1,  N);
% q_path = lambda.*q_B + (1-lambda).*q_A;
% det_path = NaN(1,N);
% for i = 1:N
%     J = robotjacobian(kin, q_path(:,i));
%     det_path(i) = det(J);
% 
%     % J = robotjacobian(kin_7, [q_path(1:2,i); pi/2 ;q_path(3:end,i)]);
%     % det_path(i) = min((svd(J)));
% end
% 
% plot(lambda, det_path)
% xlabel("\lambda")
% ylabel("det(J)")
% yline(0)

% Paths for all positive and all negative solutions
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
% for i_pair = 1
    q_A = q_A_list(:,i_pair);
    q_B = q_B_list(:,i_pair);
    
    % q_path = lambda.*q_B + (1-lambda).*q_A;
    
    q_offset = q_A;
    q_offset(2) = q_B(2);
    q_offset_2 = q_offset;
    q_offset_2(3) = q_B(3);
    
    q_path_part_1 = lambda.*q_offset + (1-lambda).*q_A;
    q_path_part_2 = lambda.*q_offset_2 + (1-lambda).*q_offset;
    q_path_part_3 = lambda.*q_B + (1-lambda).*q_offset_2;
    q_path = [q_path_part_1 q_path_part_2 q_path_part_3];
    


    for i = 1:width(q_path)
        J = robotjacobian(kin_7, [q_path(1:2,i); q_3; q_path(3:6,i)]);
        J_psi = [0 0 1 0 0 0 0];
        % J = robotjacobian(kin_7, [q_1; q_path(:,i)]);
        % J_psi = [1 0 0 0 0 0 0];
        J_aug = [J; J_psi];
        det_path_mat(i_pair, i) = det(J_aug);
    end
end


% plot(lambda, det_path_mat')
plot(det_path_mat')
xlabel("\lambda")
ylabel("det(J)")
yline(0);
%%
plot(det_path_mat(end,:))
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
% Iterate through psi = q_3

N = 1000;
q3_path = linspace(-pi, pi, N);
Q_SM = [];

for i = 1:N
    disp(i);
    kin_i = fwdkin_partial(kin_7, q3_path(i), 3);
    Q_i =  IK.IK_4_6_intersecting(R, p, kin_i);
    if ~isempty(Q_i)
        Q_i_7 = [Q_i(1:2,:); q3_path(i)*ones([1 width(Q_i)]);Q_i(3:6,:)];
        Q_SM = [Q_SM Q_i_7];
    end
end
%%
scatter( (Q_SM(1, :)), (Q_SM(2, :)), [], Q_SM(3,:), '.'); hold on
scatter(Q(1,idx_pos), Q(2,idx_pos), 200, 'rx'); 
scatter(Q(1,idx_neg), Q(2,idx_neg), 200, 'kx');
plot(q_path(1,:), q_path(2,:), 'k');
hold off
colormap hsv
xlabel("q_1")
ylabel("q_2")
axis tight

%% Plot robot in 3D (along with self-motion direction)
i=50
J = robotjacobian(kin, q_path(:,i));
det(J);
[U,S,V] = svd(J);
delta_q = V(:,end);

diagrams.setup(); hold on
diagrams.robot_plot(kin, q_path(:,i), auto_scale=true);
diagrams.robot_plot(kin, q_path(:,i)+0.2*delta_q, auto_scale=true);
diagrams.robot_plot(kin, q_path(:,i)-0.2*delta_q, auto_scale=true);


diagrams.redraw(); hold off