
kin_7 = define_yumi;

kin = fwdkin_partial(kin_7, 0, 1);
kin.P = kin.P / 100; % fix scaling for det(J)

[R, p] = fwdkin(kin, zeros([6 1]));

Q = IK.IK_4_6_intersecting(R, p, kin)

[R_t, p_t] = fwdkin(kin, Q(:,1));
[R p] - [R_t p_t]


%% Random pose
for attempt = 1:1000
q = rand_angle([6 1]);
% q = []'
[R, p] = fwdkin(kin, q);

% All IK solns
Q = IK.IK_4_6_intersecting_mex(R, p, kin);
%

% sgn(det(J)) for each soln
signs = NaN([1 width(Q)]);
for i = 1:numel(signs)
    J = robotjacobian(kin, Q(:,i));
    signs(i) = sign(det(J));
end

idx_pos = find(signs>0);
idx_neg = find(signs<0);

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
    q_A = q_A_list(:,i_pair);
    q_B = q_B_list(:,i_pair);
    q_path = lambda.*q_B + (1-lambda).*q_A;
    for i = 1:N
        J = robotjacobian(kin, q_path(:,i));
        det_path_mat(i_pair, i) = det(J);
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



%% Path between 2 solns
% q_A = Q(:,idx_pos(1)); % [0.4407   -2.6237    1.4665    0.8085   -2.5726   -0.8051]'
% q_B = Q(:,idx_pos(2)); % [-2.6952   -1.7465    1.6249   -2.0826   -0.4787    2.3917]'

% pos: 
% [-1.7731   -1.5194    1.6451   -1.7093   -1.3842    2.0765]
% [3.0542   -1.8400    1.6153    2.4274   -0.1745   -2.1431]

N = 100;
lambda = linspace(0, 1,  N);
q_path = lambda.*q_B + (1-lambda).*q_A;
det_path = NaN(1,N);
for i = 1:N
    J = robotjacobian(kin, q_path(:,i));
    det_path(i) = det(J);

    J = robotjacobian(kin_7, [q_path(1:2,i); pi/2 ;q_path(3:end,i)]);
    det_path(i) = min((svd(J)));
end

plot(lambda, det_path)
xlabel("\lambda")
ylabel("det(J)")
yline(0)

%% Confirm q_A and q_B lead to same EE pose

[R_A, p_A] = fwdkin(kin, q_A)
[R_B, p_B] = fwdkin(kin, q_B)

[R_A p_A] - [R_B p_B]