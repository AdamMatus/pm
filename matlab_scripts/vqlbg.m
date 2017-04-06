function c = vqlbg(d, k)
% VQLBG Vector quantization using the Linde-Buzo-Gray algorithm
%
% Inputs:
%       d contains training data vectors (one per column)
%       k is number of centroids required
%
% Outputs:
%       c contains the result VQ codebook (k columns, one for each centroids)

% 1-vector codebook, computing first centroid
c = zeros(length(d(:,1)),1);
c(:, 1) = mean(d, 2);

plot(d(5,:),d(6,:),'.');
hold on
plot(c(5,:),c(6,:),'*r');
hold off

% 1-vector codebook, computing first centroid
eps = 0.01;

while size(c,2) < k
    for i = 1:size(c,2)
        c = [c c(:,i)*(1+eps)];
        c(:,i) = c(:,i)*(1-eps);
    end

    plot(d(5,:),d(6,:),'.');
    hold on
    plot(c(5,:),c(6,:),'*r');
    hold off

    flag = true;
    Distortion = zeros(1,size(c,2));
    oldDistortion = inf*ones(1,size(c,2));

    while flag == true

        %nearest-neighbor search
        z = disteu(d, c);
        [m, ind] = min(z, [], 2);
        %centroids update
        for i = 1:size(c,2)
            c(:, i) = mean(d(:, find(ind == i)), 2);
            Distortion(1,i) = sum(m(find(ind == i)));
        end
        %computin relative distortion Dold-Dnew/Dnew < eps
        flag = false;
        for i = 1:size(c,2)
            if ((oldDistortion(1,i) - Distortion(1,i)) / Distortion(1,i)) > eps
                flag = true;
            end
            oldDistortion(1,i) = Distortion(1,i);
        end

        plot(d(5,:),d(6,:),'.');
        hold on
        plot(c(5,:),c(6,:),'*r');
        hold off
    end
end

end

