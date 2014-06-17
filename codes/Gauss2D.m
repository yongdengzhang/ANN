function [out,M]=Gauss2D(x,sigma)
% Gauss2D	apply a 2 dimensional gauss filter
%
%    out = Gauss2D(x,sigma);
%
%    INPUT: x      image
%           sigma  of gauss filter
%
%    OUTPUT: out   filtered image
%            M     gaussian mask
%

% bug fix: AP - 10.07.02

R = ceil(3*sigma);   % cutoff radius of the gaussian kernel
M = zeros(2*R+1); % KJ
for i = -R:R,
   for j = -R:R,
      M(i+R+1,j+R+1) = exp(-(i*i+j*j)/2/sigma/sigma);
   end
end
M = M/sum(M(:));   % normalize the gaussian mask so that the sum is
%                    % equal to 1
                   
% more correct version - and probably a bit faster
% M = GaussMask2D(sigma,2*R+1,[],1);

% Convolute matrices
out=filter2(M,x);

