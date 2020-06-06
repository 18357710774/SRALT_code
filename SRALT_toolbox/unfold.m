function Xunfold = unfold(X, mode)
% unfold the mode-"mode" of orderX-way array X (tensor)
% N = [n1, n2,... nN];

switch mode 
    case 2
        X = permute(X,[2 3 1]);
    case 3
        X = permute(X,[3 1 2]);
end
N = size(X);
dim2 = N(2)*N(3);
Xunfold = reshape(X,N(1),dim2);
end
