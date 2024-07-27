function [Q,S] = ILRED_BASIC(Q,S,a,d,tol_alpha)
% FUNCTION: ILRED_BASIC algorithm in "Efficient Eigen-Decomposition for 
% Low-Rank Symmetric Matrices in Graph Signal Processing: An Incremental 
% Approach"
QaLpha = Q'*a;
e = a - Q*QaLpha; 
p = norm(e);
k = size(S,1);

if p < tol_alpha % if error is small, use (III.6) or (III.7)
    Y1 = [S  QaLpha
         QaLpha' d];
    [Q1,S1] = eigs(Y1,k+1);
    if abs(S1(k+1,k+1)) < tol_alpha % if the eigenvalue is small
        aa = Q*Q1(1:k,1:k); % use (III.6)
        Q = [aa;Q1(k+1,1:k)];
        S = S1(1:k,1:k);
    else
        aa = Q*Q1(1:k,:);
        Q = [aa;Q1(k+1,:)];
        S = S1;
    end
else % if error is large, use (III.9) or (III.10)
    e = e/p;
    Y2 = [    S         zeros(k,1)  QaLpha
         zeros(1,k)       0        p
             QaLpha'      p        d];
    [Q2,S2] = eigs(Y2,k+2);
    if abs(S2(k+1,k+1)) < tol_alpha
        aa = [Q e]*Q2(1:k+1,1:k);
        Q = [aa;Q2(k+2,1:k)];
        S = S2(1:k,1:k);
    else  
        aa = [Q e]*Q2(1:k+1,1:k+1);
        Q = [aa;Q2(k+2,1:k+1)];
        S = S2(1:k+1,1:k+1);
    end
end
