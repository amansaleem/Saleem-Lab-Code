function [F, c1, c2, con, H] = smoothhist2D_corrected(Z, smth_win, M, bins)
if isempty(Z)
    F = nan*ones(M);
    c1 = bins;
    c2 = bins;
    con = F;
    H = F;
    return;
end

X = Z(:,1);
    Y = Z(:,2);
[F, c1, c2, H] = smoothhist2D_AS(Z, smth_win, M, bins);

KX = smoothhist2D([X ones(size(Y))], smth_win, [M(1) 2]);
KY = smoothhist2D([ones(size(Y)) Y], smth_win, [2 M(2)]);
KX = KX(2,:);
KY = KY(:,2);

% con = KY*KX;
con = repmat(KX,size(F,1),1);
con = con ./ sum(con(:)) * sum(F(:));

F = (F ./ con);
H = (H ./ con);
% imagesc(c1, c2, F)
end