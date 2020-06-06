function X = refold(Xunfold, mode, N)
% refold the mode-"mode" unfolded orderX-way array X (tensor)
% N = [n1, n2,... nN];

if mode == 1
    X = reshape(Xunfold, N);
elseif mode == 2
    N = circshift(N,[1 -1]);
    X = reshape(Xunfold, N);
    X = permute(X,[3 1 2]);
else
    N = circshift(N,[1 -2]);
    X = reshape(Xunfold, N);
    X = permute(X,[2 3 1]);
end   
end