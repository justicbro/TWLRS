function [Out] = myshrink2(X,theta)
    Nway = size(X);
    [U,S,V] = svd(X);
    Out = zeros(Nway);
    if ndims(X)>3
        printf("shrink matrix demetion must samller than 2");
        return;
    end
    N = min(Nway);
    for i = 1:N
        if S(i,i) < -theta
            S(i,i) = S(i,i) + theta;
        elseif S(i,i) > theta
            S(i,i) = S(i,i) - theta;
        else
            S(i,i) = 0;
        end
    end
    Out = U*S*V';

end