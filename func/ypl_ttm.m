function [Y] = ypl_ttm(X,D,num,flagT)

Ndim = ndims(X);
if flagT
    for i = 1:Ndim
        D{i} = D{i}';
    end
end

tempY = X;
for i = 1:Ndim
    if ismember(i,num)
        Nway = size(tempY);
        new_Nway = Nway;
        new_Nway(i) = size(D{i},1);
        tempY_m = my_Unfold(tempY,Nway,i);
        newtempY_m = D{i}*tempY_m;
        tempY = my_Fold(newtempY_m,new_Nway,i);
    end
end

Y = tempY;