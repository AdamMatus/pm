function c = mfcc(s, fs)

% MFCC Calculate the mel frequencey cepstrum coefficients (MFCC) of a signal
%
% Inputs:
%   s	: speech signal
% 	 fs 	: sample rate in Hz
%
% Outputs:
%	c 	: MFCC output, each column contains the MFCC's for one speech frame

% All previous steps...

% Obtain the mel-spectrum in the variable: ms

% Last step, compute mel-frequency cepstrum coefficients

c = dct(log(ms));

c(1,:) = [];    % exclude 0'th order cepstral coefficient
end

