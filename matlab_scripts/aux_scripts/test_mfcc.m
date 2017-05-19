test = load('speak.mat');
test = test.test;

speak_frame = sprintf('static const std::array<double, %d>> speak_test_frame{\n',length(test));
arr = sprintf(',%.9f',test(2:end));
speak_frame = sprintf('%s%.9f%s};',speak_frame,test(1),arr);

fd=fopen('../../cpp/test_inc/speak_frame.h','wt');
fprintf(fd, speak_frame);
speak_frame