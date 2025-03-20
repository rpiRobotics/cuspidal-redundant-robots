function [Q_in_limits, Q_out_limits] = filter_Q_joint_limits(Q, q_min, q_max, NameValueArgs)
% For a joint vector trajectory Q = [q^(1) q^(2) ... q^(N)]
% And lower limit q_min and upper limit q_max
% Replaces q^(i) with NaN if any element of q^(i) falls outside the range
% or removes q^(i) entirely based on the 'Mode' option.
% (Opposite for Q_out_limits)
%
% Also works on 3D matrices where 
% Q(:,:,j) = [q^(1, j) q^(2, j) ... q^(N, j)]

arguments
    Q
    q_min
    q_max
    NameValueArgs.Mode (1,1) string {mustBeMember(NameValueArgs.Mode, ["replace", "remove"])} = "replace"
end

dims = ndims(Q);

if dims == 2
    % 2D matrix: Q is m x n
    in_limit_mask = all(Q >= q_min & Q <= q_max, 1);
    
    if NameValueArgs.Mode == "replace"
        Q_in_limits = Q;
        Q_in_limits(:, ~in_limit_mask) = NaN;
        Q_out_limits = Q;
        Q_out_limits(:, in_limit_mask) = NaN;
    else  % mode == 'remove'
        Q_in_limits = Q(:, in_limit_mask);
        Q_out_limits = Q(:, ~in_limit_mask);
    end
    
elseif dims == 3
    % 3D matrix: Q is m x n x p. Process each slice independently.
    if NameValueArgs.Mode == "replace"
        Q_in_limits = Q;
        Q_out_limits = Q;
        for j = 1:size(Q,3)
            in_limit_mask = all(Q(:,:,j) >= q_min & Q(:,:,j) <= q_max, 1);
            Q_in_limits(:, ~in_limit_mask, j) = NaN;
            Q_out_limits(:, in_limit_mask, j) = NaN;
        end
    else  % mode == 'remove'
        Q_in_limits = {};
        Q_out_limits = {};
        for j = 1:size(Q,3)
            in_limit_mask = all(Q(:,:,j) >= q_min & Q(:,:,j) <= q_max, 1);
            Q_in_limits{j} = Q(:, in_limit_mask, j);
            Q_out_limits{j} = Q(:, ~in_limit_mask, j);
        end
        Q_in_limits = cat(3, Q_in_limits{:});
        Q_out_limits = cat(3, Q_out_limits{:});
    end
else
    error('Input Q must be a 2D or 3D matrix.');
end

end
