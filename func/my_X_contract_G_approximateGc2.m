function [Core] = my_X_contract_G_approximateGc2(X,G)
    Ndim = ndims(X);
    Nway = size(X);
    
    for num = 1:Ndim
        R(num) = size(G{num},1);
        L(num) = size(G{num},3);
    end

    if Ndim >2
        X_m1 = tenmat_sb(X,1);
        G_m1 = my_Unfold(G{1},size(G{1}),2);
        Y1_m = pinv(G_m1'*G_m1)*G_m1'*X_m1;
        tempNwayX = size(X);
        tempNwayX(1) = [];
        que = 3:Ndim;
        Y1 = reshape(Y1_m,[R(1),L(1),R(2),tempNwayX]);
        Y1 = permute(Y1,[3,4,2,que+2,1]);
        Y1_m = ypl_Unfold(Y1,[1,2]);
    
        Y2_m = Y1_m;
        tempL = [];
        for num = 2:Ndim-1
            G_m2 = reshape(G{num},size(G{num},1)*size(G{num},2),[]);
            pinv_G_m2 = pinv(G_m2'*G_m2)*G_m2';
            Y2 = pinv_G_m2*Y2_m;
            tempNway = Nway;
            tempNway(1:num) = [];
            tempL = [L(num-1), tempL];
            Y2 = reshape(Y2,[L(num),R(num+1),tempL,tempNway,R(1)]);
            que = 1:ndims(Y2);
            que(1:3)=[];que(num-1)=[];
            Y2 = permute(Y2,[2,num+2,1,3,que]);
            Y2_m = ypl_Unfold(Y2,[1,2]);
        end
        GN_m = ypl_Unfold(G{Ndim}, [1,2,4]);
        pinv_GN_m = pinv(GN_m'*GN_m)*GN_m';
        Y2dim = ndims(Y2);
        temp_Y2que = 1:Y2dim;
        Y2que = temp_Y2que;
        Y2que(3) = temp_Y2que(end);
        Y2que(4:end) = temp_Y2que(3:Y2dim-1);
        Y3 = permute(Y2, Y2que);
        Y3_m = ypl_Unfold(Y3, [1,2,3]);

        Core_m = pinv_GN_m*Y3_m;
        tempL = flip(L);
        Core = reshape(Core_m, tempL);
        Core = (permute(Core,flip(1:Ndim)));
    else
        G1_m = my_Unfold(G{1},size(G{1}),2);
        pG1_m = pinv(G1_m'*G1_m)*G1_m';
        Y1_m = pG1_m*X;
        Y1 = reshape(Y1_m, [R(1), L(1), R(2), Nway(2)]);
        G2_m = ypl_Unfold(G{2},[1,2,4]);
        pG2_m = pinv(G2_m'*G2_m)*G2_m';
        Y2 = permute(Y1,[3,4,1,2]);
        Y2_m = ypl_Unfold(Y2,[1,2,3]);
        Core = pG2_m*Y2_m;
        Core = permute(Core,[2,1]);
    
    end


end