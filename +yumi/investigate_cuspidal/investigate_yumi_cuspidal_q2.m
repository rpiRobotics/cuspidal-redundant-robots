
kin_7 = define_yumi;

kin = fwdkin_partial(kin_7, pi/3, 2);
kin.P = kin.P / 100; % fix scaling for det(J)

[R, p] = fwdkin(kin, zeros([6 1]));
% [R, p] = fwdkin(kin, rand_angle([6 1]));

Q = IK.IK_4_6_intersecting(R, p, kin)

[R_t, p_t] = fwdkin(kin, Q(:,1));
[R p] - [R_t p_t]


%% Random pose
for attempt = 1:1000
q = rand_angle([6 1]);
% q = []'
[R, p] = fwdkin(kin, q);

% All IK solns
Q = IK.IK_4_6_intersecting(R, p, kin);
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

if numel(idx_pos) >= 4 && numel(idx_neg) >= 4
q_A_list = [Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_neg(1)) Q(:,idx_neg(1)) Q(:,idx_neg(1))];
q_B_list = [Q(:,idx_pos(2)) Q(:,idx_pos(3)) Q(:,idx_pos(4)) Q(:,idx_neg(2)) Q(:,idx_neg(3)) Q(:,idx_neg(4))];
elseif numel(idx_pos) >= 3 && numel(idx_neg) >= 3
q_A_list = [Q(:,idx_pos(1)) Q(:,idx_pos(1)) Q(:,idx_neg(1)) Q(:,idx_neg(1))];
q_B_list = [Q(:,idx_pos(2)) Q(:,idx_pos(3)) Q(:,idx_neg(2)) Q(:,idx_neg(3))];
elseif numel(idx_pos) >= 2 && numel(idx_neg) >= 2
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
% q_A = Q(:,idx_pos(1)); % [-1.2353595643770409751027727907058, 0.83995490411487372384158334170934, 1.6667289900385404699534319661325, -0.17112461375420415232717630260595, -2.3093090236793156755368272570195, 1.1043509798935793320140419382369]
 
% q_B = Q(:,idx_pos(2)); % [1.7921170938638899539085969081498, -0.97088919297651199435961189010413, 1.8132481598686596147018690317054, 0.81539702404301772631356470810715, -1.0769669541497495224291469639866, 3.1255250151859299556633686734131]



N = 100;
lambda = linspace(0, 1,  N);
q_path = lambda.*q_B + (1-lambda).*q_A;
det_path = NaN(1,N);
for i = 1:N
    J = robotjacobian(kin, q_path(:,i));
    det_path(i) = det(J);
end

plot(lambda, det_path)
xlabel("\lambda")
ylabel("det(J)")
yline(0)

%% Confirm q_A and q_B lead to same EE pose

[R_A, p_A] = fwdkin(kin, q_A)
[R_B, p_B] = fwdkin(kin, q_B)

[R_A p_A] - [R_B p_B]