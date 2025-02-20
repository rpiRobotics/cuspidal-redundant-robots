%% Proof that YuMi parameterized by any joint angle is cuspidal
%  If we parameterize by q_i, we can show the subrobot formed by locking
%  joint i (with some joint angle) is cuspidal. We just need to find some
%  q_A and q_B where the MoveJ in between is nonsingular

kin_7 = define_yumi;
kin_7.P = kin_7.P / 100; % Better scaling for det(J)

%% Joint 1

kin = fwdkin_partial(kin_7, 0, 1);
q_A = [0.4407   -2.6237    1.4665    0.8085   -2.5726   -0.8051]';
q_B = [-2.6952   -1.7465    1.6249   -2.0826   -0.4787    2.3917]';

fwd_kin_error(kin, q_A, q_B)
plot_det_J_path(kin, q_A, q_B);

%% Joint 2
kin = fwdkin_partial(kin_7, pi/3, 2);
q_A = [-1.2353595643770409751027727907058, 0.83995490411487372384158334170934, 1.6667289900385404699534319661325, -0.17112461375420415232717630260595, -2.3093090236793156755368272570195, 1.1043509798935793320140419382369]';
q_B = [1.7921170938638899539085969081498, -0.97088919297651199435961189010413, 1.8132481598686596147018690317054, 0.81539702404301772631356470810715, -1.0769669541497495224291469639866, 3.1255250151859299556633686734131]';
fwd_kin_error(kin, q_A, q_B)
plot_det_J_path(kin, q_A, q_B);
 
%% Joint 3
kin = fwdkin_partial(kin_7, pi/2, 3);
q_A = [-1.6884    1.2167   -2.9637   -0.1269   -0.9153    0.9605]';
q_B = [1.7871   -1.3096    0.6505    2.9181   -0.4242    1.2237]';

fwd_kin_error(kin, q_A, q_B)
plot_det_J_path(kin, q_A, q_B);

%% Joint 4
kin = fwdkin_partial(kin_7, pi/3, 4);

q_A = [0.9580   -1.1292   -2.4903    0.2235   -2.1057    2.4094]';
q_B = [-0.1976   -3.0522   -2.1780   -0.9074   -1.3382    0.2224]';

fwd_kin_error(kin, q_A, q_B)
plot_det_J_path(kin, q_A, q_B);

%% Joint 5
kin = fwdkin_partial(kin_7, pi/3, 5);
q_A = [3.0726   -2.6636    2.4202    0.2566   -0.0434   -1.6936]';
q_B = [2.8386   -3.1323   -1.2856   -2.8371   -0.3598    1.8212]';

fwd_kin_error(kin, q_A, q_B)
plot_det_J_path(kin, q_A, q_B);

%% Joint 6
kin = fwdkin_partial(kin_7, pi/3, 6);

q_A = [1.7299   -1.5880    2.6435    1.6695    2.1005   -2.6614]';
q_B = [2.0352   -0.2238   -1.7780    1.9442    2.7318   -0.8155]';

fwd_kin_error(kin, q_A, q_B)
plot_det_J_path(kin, q_A, q_B);

%% Joint 7
kin = fwdkin_partial(kin_7, pi/3, 7);

q_A = [1.0639   -0.7255   -3.0816    1.9216   -1.9533   -1.1198]';
q_B = [-0.2112   -2.3409    0.7570    1.9369    1.7178    0.6800]';
fwd_kin_error(kin, q_A, q_B)
plot_det_J_path(kin, q_A, q_B);

%%

function err_mat = fwd_kin_error(kin, q_A, q_B)
    [R_A, p_A] = fwdkin(kin, q_A);
    [R_B, p_B] = fwdkin(kin, q_B);
    err_mat = [R_A p_A] - [R_B p_B];
end

function plot_det_J_path(kin, q_A, q_B)
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
end