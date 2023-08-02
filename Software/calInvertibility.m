function invert = calInvertibility(H2, X)
Xt = X(:);
nx = size(X,1);
ny = size(X,2);
Y = H2'*Xt;
Xinv = pinv(H2')*Y;
rec = reshape(Xinv,[nx ny]);
invert = ssim(double(rec), X);

end