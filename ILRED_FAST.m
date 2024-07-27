function [Q_tilt,S,Q_hat,Q_hatPlus] = ILRED_FAST(Q_tilt,S,Q_hat,Q_hatPlus,a,d,tol_alpha)
% FUNCTION: ILRED_FAST algorithm in "Efficient Eigen-Decomposition for 
% Low-Rank Symmetric Matrices in Graph Signal Processing: An Incremental 
% Approach"

QaLpha = Q_hat'*(Q_tilt'*a);
e = a - Q_tilt*(Q_hat*QaLpha); 
p = norm(e);
k = size(S,1);

if p < tol_alpha % if error is small, use (III.15) & (III.18)
    Y = [S  QaLpha
         QaLpha' d];
    [QY,SY] = eigs(Y,k+1);
    if abs(SY(k+1,k+1)) < tol_alpha
        W = QY(1:k,1:k);
        w = QY(k+1,1:k);
        Q_hat = Q_hat*W;
        if 1-w*w' > realmin
            Wplus = W' + (w'*(w*W'))/(1-w*w');
        else
            Wplus = W';
        end
        Q_hatPlus = Wplus*Q_hatPlus;
        Q_tilt_add = w*Q_hatPlus;
        Q_tilt = [Q_tilt;Q_tilt_add];
        S = SY(1:k,1:k);
    else
        Q_hat(end+1,end+1) = 1;
        Q_hat = Q_hat*QY;
        Q_hatPlus(end+1,end+1) = 1;
        Q_hatPlus = QY'*Q_hatPlus;
        Q_tilt(end+1,end+1) = 1;
        S = SY;
    end
else % if error is large, use (III.21)
    e = e/p;
    Y = [    S         zeros(k,1)  QaLpha
         zeros(1,k)       0        p
             QaLpha'           p        d];
    [QY,SY] = eigs(Y,k+2);
    if abs(SY(k+1,k+1)) < tol_alpha  
        % if the rank is unchanged, multiply directly, faster
        Q_tilt = [Q_tilt*(Q_hat*QY(1:k,1:k))+e*QY(1+k,1:k);QY(k+2,1:k)];
        Q_hat = eye(k);
        Q_hatPlus = eye(k);
        S = SY(1:k,1:k);
    else  % if the rank is increased
        W = QY(1:k+1,1:k+1);
        w = QY(k+2,1:k+1);
        Q_hat(end+1,end+1) = 1;
        Q_hat = Q_hat*W;
        if 1-w*w' > realmin
            Wplus = W' + (w'*(w*W'))/(1-w*w');
        else
            Wplus = W';
        end
        Q_hatPlus(end+1,end+1) = 1;
        Q_hatPlus = Wplus*Q_hatPlus;
        Q_tilt = [[Q_tilt e];w*Q_hatPlus];
        S = SY(1:k+1,1:k+1);
    end
end
end