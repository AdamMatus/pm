v = load('v.mat');
v = v.v;
mfcc_vq = sprintf('static const std::vector<std::array<double, 29>>mfcc_test_frames{\n');
for i = 1:size(v,2)
    arr = sprintf(',%.9f',v(2:end,i));
    arr = sprintf('{%.9f%s},\n', v(1,i), arr);
    mfcc_vq = sprintf('%s%s',mfcc_vq, arr);
end
mfcc_vq = sprintf('%s};',mfcc_vq);

fd=fopen('mfcc_test.h','wt');
fprintf(fd, mfcc_vq);
mfcc_vq