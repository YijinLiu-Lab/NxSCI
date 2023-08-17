function resolution = getFSC(im1, im2)

scale_F = 1;

[corrCoeffs, spatialFrequency, meanIntensity] = FourierShellCorrelate(imresize(im1(1:2:end,1:2:end),scale_F,"bicubic"), imresize(im2(2:2:end,2:2:end),scale_F,"bicubic"),150,1);

corrCoeffs = smooth(corrCoeffs);
% figure, plot(spatialFrequency,corrCoeffs,'k','LineWidth',1)
A=[spatialFrequency; corrCoeffs'];
B=[spatialFrequency; 0.143*ones(numel(spatialFrequency),1)'];
P = InterX(A,B);
if numel(P)>0
    resolution = 1./P(1,1);
else
    resolution = 0;
end
resolution = resolution./scale_F;

end