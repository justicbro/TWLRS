function [X, G, Core, Out] = TW_LRG_sparsityGc_TC(F, Omega, opts)
if isfield(opts, 'tol');         tol   = opts.tol;              end
if isfield(opts, 'maxit');       maxit = opts.maxit;            end
if isfield(opts, 'R');           max_R = opts.R;                end
if isfield(opts, 'lambda');      lambda= opts.lambda;           end
if isfield(opts, 'rho');         rho   = opts.rho;              end

%
R = max(max_R, 2);
Ndim = ndims(F); 
Nway = size(F);
X = F;
XT = opts.Xtrue;
p = 1;
q = 1;

%
Factors_dims = factor_dims(Nway, R);
rng('default')
%
G = cell(Ndim,1);
Gunfold = cell(Ndim,1);
for i=1:Ndim
    G{i}=rand(Factors_dims(i,:));
    tempG{i} = my_Unfold(G{i},size(G{i}),2);
    [n1, n2] = size(tempG{i});
    W{i} = tempG{i};
    H{i} = zeros(n1,n2);
end
Core = rand(R(2,:));
Out.RSE = [];
Out.RE = [];
Out.PSNR = [];

%
for k = 1:maxit
    X_old = X;
    % Update G_k, k=1,2,...,N.
    for num = 1:Ndim
        GCrest = tenremat(circ_tnprod_rest(G,Core,num), Ndim); 
        tempA = (rho*W{num} - H{num} + tenmat_sb(X,num)*GCrest');
        tempB = (GCrest*GCrest') + rho*eye(size(GCrest,1));
        Gunfold{num} = tempA*pinv(tempB);
        G{num} = my_Fold(Gunfold{num},size(G{num}),2);
        
%       Update W_k
        W{num} = myshrink2(Gunfold{num} + H{num}/rho, lambda(1)/rho);
    end

    % Update X 
    X = my_Core_contract_G_approximateX(Core, G);   
    X(Omega) = F(Omega);

    % Update Core
    if k==1  || ( k>p && mod(k, q)==0 )
        Core = my_X_contract_G_approximateGc2(X, G);
        Core = soft(Core, lambda(2));
    end

    %Update Hk
    for num = 1:Ndim
        H{num} = H{num} + rho *(Gunfold{num} - W{num});
    end
    
    %% check the convergence
    re = norm(X(:)-X_old(:)) / norm(X_old(:));
    rse = norm(XT(:)-X(:)) / norm(XT(:));
    PSNR     = psnr_ybz(XT*255,X*255);
    Out.RSE = [Out.RSE, rse];
    Out.PSNR = [Out.PSNR, PSNR];
    Out.RE = [Out.RE, re];
    
    if k == 1 || mod(k, 20) == 0
        fprintf('TWLRS_TC: iter = %d   PSNR=%f   RE=%.10f   \n', k, PSNR, re);
    end

    if re < tol 
        break;
    end
    
end
end
