function c = mfcc(s, fs)
% MFCC Calculate the mel frequencey cepstrum coefficients (MFCC) of a signal
%
% Inputs:
%       s       : speech signal
%       fs      : sample rate in Hz
%
% Outputs:
%       c       : MFCC output, each column contains the MFCC's for one speech frame

N = 256;
M = 100;

blk = zeros(N,1);

i = 1;
while i*M + N < length(s) 
    blk(:,i) = s(i*M:i*M+N-1);
    i = i+1;
end

%
imagesc(blk);

% applying window to frames
for i = 1:(size(blk,2))
    blk(:,i) = blk(:,i).*hamming(N);
end

%
imagesc(blk);

% f is now power spectrum
f = abs(fft(blk)).^2;

% Obtaining the mel-spectrum
p = 30; % number of mel filtersbanks
n2 = 1 + floor(N/2);
m = melfb(p, N, fs);

ms = zeros(p,size(f,2));
for i = 1:(size(f,2))
   ftemp = f(:,i);
   ms(:,i) = m * ftemp(1:n2); 
end

%
imagesc(ms);

% Last step, compute mel-frequency cepstrum coefficients
c = idct(log10(ms));
c(1,:) = [];    % exclude 0'th order cepstral coefficient

%
imagesc(c);

end

