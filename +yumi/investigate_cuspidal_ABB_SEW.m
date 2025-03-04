
[kin, q_min, q_max] = define_yumi;
kin.P = kin.P / 100; % fix scaling for det(J)
SEW = yumi.sew_abb();

%% Random pose
for attempt = 1:1000
q = rand_angle([7 1]);
[R, p] = fwdkin(kin, q);
psi = SEW.fwd_kin(kin, q);

% All IK solns
Q = yumi.IK_SEW_ABB_mex(R, p, SEW, psi, kin);
Q = unique_q_tol(Q, 1e-3);
Q = expand_q_with_2pi(Q);
Q = q_within_range(Q, q_min, q_max);

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
    % skip if the same by 2*pi apart
    if norm(wrapToPi(q_A - q_B)) < 1e-6
        continue
    end

    % Skip if q_6=0 is passed, as that's an augmentation singularity
    if q_A(6)*q_B(6) < 0
        continue
    end

    q_path = lambda.*q_B + (1-lambda).*q_A;
    for i = 1:N
        J = robotjacobian(kin, q_path(:,i));
        J_psi = SEW.jacobian(kin, q_path(:,i));
        J_A = [J; J_psi];
        det_path_mat(i_pair, i) = det(J_A);
    end
end

disp(attempt)
if any(all(det_path_mat'>0)) || any(all(det_path_mat'<0))
% if any(all(det_path_mat'>-1e-1)) || any(all(det_path_mat'<1e-1))
% if ~isempty(det_path_mat(~isnan(det_path_mat)))
    disp("found it!")
    beep; pause(.5); beep

    plot(lambda, det_path_mat')
    xlabel("\lambda")
    ylabel("det(J)")
    yline(0);
    break
end


end

%% Which was succesful?
find(all(det_path_mat'<-1e-1) | all(det_path_mat'>1e-1))

%% Path between 2 solns %  2     5    23    27    31    41
q_A = q_A_list(:,6)
q_B = q_B_list(:,6)

%%
q_A = [2.0137    3.1221   -3.0802   -0.2020    1.5407    2.9876  0.4203];
q_B = [-1.8861    2.6673   -0.5224    0.2350    1.2883    2.8092    0.4334];

%%
% -1.3752    2.7595   -2.6584   -2.1838    2.1299    1.7472   -1.3041
 % 0.3107    3.1244   -0.0105   -2.0362    1.1804    1.9455   -1.1798

%% Great example - very close q
q_A = [2.4187    2.3280   -0.5604    0.7826   -2.0013   -3.0373   -0.9235]'
q_B = [2.3610    2.2931   -0.9101    0.8240    2.2905   -3.0097   -2.6589]'

%% Check if within RobotStudio limits
[q_min < q_A , q_A < q_max]
[q_min < q_B , q_B < q_max]
%%
N = 1000;
lambda = linspace(0, 1,  N);
q_path = lambda.*q_B + (1-lambda).*q_A;
det_path = NaN(1,N);
sign_term_path = NaN(1,N);
for i = 1:N
    J = robotjacobian(kin, q_path(:,i));
    [J_psi, sign_term_i] = SEW.jacobian(kin, q_path(:,i));
    J_A = [J; J_psi];
    det_path(i) = det(J_A);
    sign_term_path(i) = sign_term_i;
end

plot(lambda, det_path)
xlabel("\lambda")
ylabel("det(J)")
yline(0)
%%
plot(lambda, sign_term_path)
xlabel("\lambda")
ylabel("`e_r^T R_{03} h_4")
yline(0)

%% Plot SEW angle over path


q_path = lambda.*q_B + (1-lambda).*q_A;
psi_path = NaN(1,N);
for i = 1:N
    psi_path(i) = SEW.fwd_kin(kin, q_path(:,i));
end

plot(lambda, psi_path)
xlabel("\lambda")
ylabel("\psi")
yline(0);

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


%% Optimize a midpoint - adjust straight line path
% q_A = [ -0.3951    0.2511   -2.3602   -1.1889    0.1224    0.1595    0.3801]';
% q_B = [1.1792   -0.2441    2.3177   -1.8662    0.0365    0.4881    0.4695]';
% q_M_i = [2.2307   -0.2262    1.4681   -1.5702    0.3516   -0.0109    0.7228]';
q_M = q_A*0.5 + q_B*0.5;

delta = rand_angle([7 1]) /10;

q_M_i = q_M + delta;

N = 100;
lambda = linspace(0, 1,  N);
lambda_half = linspace(0, 1,  N/2);

q_path_AM = lambda_half.*q_M_i + (1-lambda_half).*q_A;
q_path_MB = lambda_half.*q_B + (1-lambda_half).*q_M_i;
q_path = [q_path_AM q_path_MB];
det_path = NaN(1,N);
psi_path = NaN(1,N);
sign_term_path = NaN(1,N);
for i = 1:N
    J = robotjacobian(kin, q_path(:,i));
    [J_psi, sign_term_i] = SEW.jacobian(kin, q_path(:,i));
    J_A = [J; J_psi];

    det_path(i) = det(J_A);
    psi_path(i) = SEW.fwd_kin(kin, q_path(:,i));
    sign_term_path(i) = sign_term_i;
end

tiledlayout(3,1)
nexttile
plot(lambda, det_path);
yline(0)
xlabel("\lambda")
ylabel("det(J)")



nexttile
plot(lambda, psi_path);
yline(0)
xlabel("\lambda")
ylabel("\psi")

nexttile
plot(lambda, sign_term_path);
yline(0)
xlabel("\lambda")
ylabel("e_r^T R_{03} h_4")


%%

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

function Q = q_within_range(Q, q_min, q_max)
    % Filters Q to keep only values within (q_min, q_max)
    Q = Q(:, all(Q > q_min & Q < q_max, 1));
end

function Q_expanded = expand_q_with_2pi(Q)
    % Generates additional q_i values by increasing or decreasing each element by 2*pi
    [dim, num] = size(Q);
    Q_expanded = Q;
    
    for i = 1:num
        for k = 1:dim
            q_plus = Q(:, i);
            q_minus = Q(:, i);
            q_plus(k) = q_plus(k) + 2*pi;
            q_minus(k) = q_minus(k) - 2*pi;
            Q_expanded = [Q_expanded, q_plus, q_minus];
        end
    end
end
