function newtrainigdata(number)
% number of speaker

fs = 12500;

recObj = audiorecorder(fs,16,1);
disp('Start speaking.');
recordblocking(recObj, 2);
disp('End of Recording.');
y = getaudiodata(recObj);

file = sprintf('train\\s%d.wav',number);
wavwrite(y,fs,file);

end

