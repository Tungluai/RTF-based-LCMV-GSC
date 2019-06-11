function [RTF,SPP,Mark] = segmentation (speech,fs,ov,t_p_frm)

% reference : Data-Driven Source Separation Based on Simplex Analysis,2018,Bracha
%******input
% speech : source in time domain
% fs :sample rate
% t_p_frm : time delay in per frame (ms)
% ov :  overlap ov%

%******output
% RTF: Nch x nsrce x Nbin   relative transmission fuction
% nsrce : source number
% Mark: nsrce x (length of speech)    to mark segments

%author : Xu Changlai,6/2,2019

[Nch,Nz] = size(speech);
Nfft =floor( fs*t_p_frm/1000); % 64 ms per frame
Nbin = floor(Nfft/2+1);
Nov = 100 / (100-ov) ;
Lbin = floor( Nfft/Nov ); 
Nfrm = floor(Nz/Lbin)-(Nov-1);
win = sqrt(hanning(Nfft))';

for frm = 1 : Nfrm          
    %STFT
    for ch = 1 : Nch
         Y(ch ,frm,:) = fft(win .* speech(ch ,(frm-1)*Lbin+1:(frm-1)*Lbin+Nfft),Nfft);
    end
end

upbin = floor(4.8*1000/fs*Nfft);
lowbin = floor(0.3*1000/fs*Nfft+1);
for frm = 1: Nfrm
  for bin = lowbin : upbin   %  0.2 ~ 4.8 KHz
    if frm == 1
          Am(:,frm,bin-lowbin+1) = (Y(2:end,frm,bin) * Y(1,frm,bin)' + Y(2:end,frm+1,bin) * Y(1,frm+1,bin)') ./...
              (Y(1,frm,bin) * Y(1,frm,bin)' + Y(1,frm+1,bin) * Y(1,frm+1,bin)'); 
    
    elseif frm == Nfrm    
          Am(:,frm,bin-lowbin+1) = (Y(2:end,frm-1,bin) * Y(1,frm-1,bin)' + Y(2:end,frm,bin) * Y(1,frm,bin)') ./...
              (Y(1,frm-1,bin) * Y(1,frm-1,bin)' + Y(1,frm,bin) * Y(1,frm,bin)'); 
    else
          Am(:,frm,bin-lowbin+1) = (Y(2:end,frm-1,bin) * Y(1,frm-1,bin)' + Y(2:end,frm,bin) * Y(1,frm,bin)'+ Y(2:end,frm+1,bin) * Y(1,frm+1,bin)') ./...
              (Y(1,frm-1,bin) * Y(1,frm-1,bin)' + Y(1,frm,bin) * Y(1,frm,bin)'+ Y(1,frm+1,bin) * Y(1,frm+1,bin)'); 
    end
  end
  ac(:,frm) = reshape((reshape(Am(:,frm,:),Nch-1,upbin-lowbin+1))',(Nch-1)*(upbin-lowbin+1),1);
  a(:,frm) = [real(ac(:,frm));imag(ac(:,frm))];
  for n = 1 : frm
     W(frm,n) = a(:,frm)'* a(:,n)/(2*(Nch-1)*(upbin-lowbin+1));
     W(n,frm) = W(frm,n);
  end
end

% EVD on W
[U,D] = eig(W);
norm_eigv = diag(D)/D(end,end);
V = [];
for cnt = length(norm_eigv):-1:1
  if norm_eigv(cnt) < .119  %  0.11 ~ 0.128
      nsrce = length(norm_eigv) - cnt;
      break;
  end
  V = [V,U(:,cnt)];
end
% find probability vector
[~,I1] = max(sum(V.^2,2));
e(1,:) = V(I1,:); 
V1 = V - repmat(e(1,:),size(V,1),1);
[~,I2] = max(sum(V1.^2,2));
e(2,:) = V(I2,:);
if nsrce > 2
  Er = [];
  for r = 3:nsrce      
      er = e(r-1,:)-e(1,:); 
      Er = [Er er'];
      temp = pinv(Er' * Er);%pseudoantique
      P = eye(nsrce)- Er * temp * Er';
      [~,I] = max(sum((P * V1').^2,1));
      e(r,:) = V(I,:);
  end  
 end
Q = e';
SPP = (Q \ V')'; % source present probality
% clustering and estimate the RTF 
for n = 1 : nsrce
    Ydom = zeros(Nch,1,Nbin);
    Yref = zeros(1,1,Nbin);
    mark = zeros(1,Nz);
    L = find (SPP(:,n) > .96);%classic .95
    disp(L);
    for i = 1 : length(L)
      Ydom  = Ydom + Y(: ,L(i),1:Nbin) .* repmat(conj(Y(1 ,L(i),1:Nbin)),Nch,1);
      Yref  = Yref + Y(1 ,L(i),1:Nbin) .* conj(Y(1 ,L(i),1:Nbin));
      mark((L(i)-1)*Lbin+1:(L(i)-1)*Lbin+Nfft) = ones(1,Nfft);
    end 
    RTF(:,n,:) = Ydom ./ repmat(Yref,Nch,1,1);
    % making marks
    Mark(n,:) = mark/10;
end

figure(1);
plot(speech(1,:));
hold on
plot(Mark(1,:));
hold on
plot(Mark(2,:));
end