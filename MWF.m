% the MWF beamformer
% SDWMWF:   h = (PhiX + mu * PhiN)^-1 * PhiSN(:,bin)
% input:    PhiN, (Nch, Nch, Nbin) the noise covariance matrix
%           PhiX, (Nch, Nch, Nbin) the speech covariance matrix
%           PhiSN,(Nch, Nbin) the noise&speech across-covariance matrix
%           mu, the speech distortion/noise reduction trade-off parameter
% output:   h, (Nch, Nbin)  the beamformer coefficients
% author : Xu Changlai,6/1,2019

function h = MWF(PhiX,PhiN,PhiSN,mu)
if nargin < 3
    mu = 1;                 % typical value {0, 1}
end

[Nch, ~, Nbin] = size(PhiX);
h = zeros(Nch, Nbin);

for bin = 1:Nbin
    if rcond(PhiN(:,:,bin)) < eps
     %   disp(['bin ' num2str(bin) ': Noise covariance ill-conditioned.']);
        PhiN(:,:,bin) = PhiN(:,:,bin) + 1e-10 * eye(Nch);
    end
    h(:,bin)  = (PhiX(:,:,bin) + mu * PhiN(:,:,bin)) \ PhiSN(:,bin);
end


