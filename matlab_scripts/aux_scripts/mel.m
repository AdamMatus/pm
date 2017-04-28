function [ m ] = mel( f )
% freqence to mel scale
m = 2595*log10(1+f/700);

end

