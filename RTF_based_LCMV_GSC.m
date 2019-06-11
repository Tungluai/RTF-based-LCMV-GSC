% LCMV-GSC for speech enhancement
% author : Xu Changlai,6/2,2019

clear all
close all

[speech , fs ] = audioread('male_female_pure_mixture.wav');
speech = speech';
[Nch,Nz] = size(speech);
Nfft =floor( fs*64/1000); % 64 ms per frame
Nbin = floor(Nfft/2+1);
Nfrm = floor(Nz/Nbin)-1;
win = sqrt(hanning(Nfft))';

yout = zeros(1,Nz);
Ybin_nonclosed = zeros(1,Nbin);
 
q = zeros(Nch+1, Nbin);
pest = zeros(Nbin,1);
mu = 0.05;
alphaP = 0.9;
Yfbf = zeros(Nch,Nfft);
phi_x = zeros(Nbin);
phi_n = zeros(Nbin);
PhiN = zeros(Nch+1, Nch+1, Nbin);
PhiS = zeros(Nch+1, Nch+1, Nbin);
PhiSN = zeros(Nch+1,Nbin);

% processing
[RTF,SPP,Mark] = segmentation (speech,fs,75,64);
% resampling C 
C1 = shiftdim(RTF,2);
for nsr = 1 : size(C1,3)
    C2(:,:,nsr) = resample(C1(:,:,nsr),64,64);
end
C = shiftdim(C2,1);

nsrce = size(C,2);
g = [1;zeros(nsrce-1,1)];
g = flipud(g); % change the enhanced person in the case of 2 speakers
enhansp = find(1 == g); 
 for frm = 1 : Nfrm  
    %STFT
    for ch = 1 : Nch
         Y(ch ,:) = fft(win .* speech(ch ,(frm-1)*Nbin+1:(frm-1)*Nbin+Nfft),Nfft);
    end

    for bin=1:Nbin  
        w0(:,bin) = C(:,:,bin)/(C(:,:,bin)'* C(:,:,bin)) * g;   
        B(:,:,bin) = eye(Nch,Nch) - C(:,:,bin) /(C(:,:,bin)'*C(:,:,bin))*C(:,:,bin)';
    % processing      
        % FBF filtering
        Yfbf(bin) = w0(:,bin)' * Y(:,bin)/norm(w0(:,bin));    
        % BM filtering
        u(:,bin) = B(:,:,bin) * Y(:,bin);
        
        Yout(bin) = Yfbf(bin) - q(:,bin)'* [Yfbf(bin);u(:,bin)];

        % SDW-MWF
        S(:,bin) = [Yout(bin);1e-10 * ones(Nch,1)];
        N(:,bin) = [Yfbf(bin)-Yout(bin); u(:,bin)];
        PhiS(:,:,bin) = 0.98 * PhiS(:,:,bin) + 0.02 * S(:,bin) * S(:,bin)';
        PhiN(:,:,bin) = 0.98 * PhiN(:,:,bin) + 0.02 * N(:,bin) * N(:,bin)';
        PhiSN(:,bin) = 0.98 * PhiSN(:,bin) + 0.02 * N(:,bin) * (Yfbf(bin)-Yout(bin))';
        q(:,bin) = MWF(PhiS(:,:,bin), PhiN(:,:,bin),PhiSN(:,bin),1.11); 
    end
    % load all stuff
    yout((frm-1)*Nfft/2+1:(frm-1)*Nfft/2+Nfft) = yout((frm-1)*Nfft/2+1:(frm-1)*Nfft/2+Nfft) + win .* real(ifft([Yout,conj(Yout(end-1:-1:2))]));   
end
audiowrite('RTF.wav', yout,fs);

% plot
figure(2);
subplot(3,1,1);
plot(audioread('male.wav'));
hold on 
plot(Mark(1,:));
subplot(3,1,2);
plot(audioread('female.wav'));
hold on 
plot(Mark(2,:));
subplot(3,1,3);
plot(yout);

% [ scoresbefore ] = pesq( 'male.wav', 'male_female_pure_mixture.wav' );
% [ scoresafter ] = pesq( 'male.wav', 'RTF.wav' );
% [ scoresideal ] = pesq( 'male.wav', 'male.wav' );
% fprintf('scorebefore: %f\n',scoresbefore);
% fprintf('scoreafter: %f\n',scoresafter);
% fprintf('scoreideal: %f\n',scoresideal);
% fprintf(['improved PESQ socre : %f\n'],scoresafter-scoresbefore);