function [Q,S,tailThres] = ILRED_FAST_MultiDimension(initMtx,finalMtx,k)
% FUNCTION: use ILRED_FAST alg. to solve the multi-frame adding problem, 
% tol_alpha and tailthres are decided based on the initial eigenvalues

% parameters
dim_init = size(initMtx,1);
dim_final = size(finalMtx,1);
tol_alpha = 1e-10;

% calculate initial Q and S using eigs
rank_init = min(dim_init,k+10);
[Q_init,S] = eigs(initMtx,rank_init,'largestabs','Tolerance',tol_alpha);
S_diag = diag(S);
tailThres = S_diag(k);
tol_alpha = 0.1*tailThres;
rank_init = k;
S = S(1:rank_init,1:rank_init);
Q_tilt = Q_init(:,1:rank_init);
Q_hat = eye(rank_init);
Q_hatPlus = eye(rank_init);

% main loop
for indxDim = dim_init+1:dim_final
    a = finalMtx(1:indxDim-1,indxDim);
    d = finalMtx(indxDim,indxDim);
    [Q_tilt,S,Q_hat,Q_hatPlus] = ILRED_FAST(Q_tilt,S,Q_hat,Q_hatPlus,a,d,tol_alpha);
end

% convert EVD result to Q
Q = Q_tilt*Q_hat;

% cut the tail of S (eigenvalues) using tailThres
S_diag = diag(S);
S_diag_keep = S_diag(S_diag > tailThres);
num_s_keep = length(S_diag_keep);
S = diag(S_diag_keep);
Q = Q(:,1:num_s_keep);
end