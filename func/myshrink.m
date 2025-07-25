function [Out] = myshrink(X,theta)
    Nway = size(X);
    Out = zeros(Nway);
    if ndims(X)>3
        printf("shrink matrix demetion must samller than 2");
        return;
    end
    for i = 1:Nway(1)
        for j = 1:Nway(2)
            if X(i,j) < -theta
                Out(i,j) = X(i,j) + theta;
            elseif X(i,j) > theta
                Out(i,j) = X(i,j) - theta;
            else
                Out(i,j) = 0;
            end
        end
    end
end