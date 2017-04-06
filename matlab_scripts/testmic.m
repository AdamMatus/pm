function [numDetected, diffDetected] = testmic(code)

fs = 12500;

recObj = audiorecorder(fs,16,1);
disp('Start speaking.')
recordblocking(recObj, 2);
disp('End of Recording.');
y = getaudiodata(recObj);

 v = mfcc(y, fs);            % Compute MFCC's
   
    diffDetected = inf;
    numDetected = 0;
   
    for l = 1:length(code)      % each trained codebook, compute distortion
        d = disteu(v, code{l}); 
        dist = sum(min(d,[],2)) / size(d,1);
      
        if dist < diffDetected
            diffDetected = dist;
            numDetected = l;
        end      
    end
   
    msg = sprintf('Speaker captured from microphone matches with speaker %d', numDetected);
    disp(msg);

end

