%The Input W is the tensor need to unfold and i is the unfold of the i-TH dimension. i can be a vector or a scalar.
%The output is the Unfold matrix
function W = ypl_Unfold(W, i)
dim = size(W);
ilength = length(i);
n=1:length(dim); 
n([i])=[]; 
m=zeros(length(dim),1);
m(1:ilength)=i;
m(ilength+1:end)=n;
W = permute(W,m);
Wshape = 1;
for i = 1:ilength
    Wshape = Wshape*size(W,i);
end
W = reshape(W, Wshape, []);
end