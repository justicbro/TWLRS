function [X] = my_Core_contract_G_approximateX(Core,G)
    Ndim = ndims(Core);
    for num=1:Ndim
        Nway(num) = size(G{num},2);
    end

    GCrest = tenremat(circ_tnprod_rest(G,Core,1), Ndim);
    G1_m = my_Unfold(G{1},size(G{1}),2);
    X_m = G1_m*GCrest;

    X = my_Fold(X_m,Nway,1);
end