csratio = 0.05;
z = x;
p = 80;
step = ceil(p/(p*sqrt(csratio)));
zz = zeros(size(z));
zz = imresize(z(1:step:end,1:step:end),[p p],'box');
figure; imshow([zz x],[])