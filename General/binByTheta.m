function es_out = binByTheta(es, chn)

if chn==2
    diffAdj = angle(es.theta.B.hill) + circshift(angle(es.theta.B.hill),1);
    gradAng = gradient(angle(es.theta.B.hill));
else
    diffAdj = angle(es.theta.A.hill) + circshift(angle(es.theta.A.hill),1);
    gradAng = gradient(angle(es.theta.A.hill));
end

zeroCross = find(abs(diffAdj)<1 & gradAng<0);
zeroCross(find(diff(zeroCross)==1)+1) = [];

outSize = length(zeroCross)-1;
numCells = size(es.spikeTrain,2);
es_out = initialise_es(outSize, numCells);
es_out.orig = es;

es_out.spikeIDs = es.spikeIDs;
es_out.isolDist = es.isolDist;

for ibin = 1:outSize
    % sum over
    es_out.mua(ibin) = sum(es.mua(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.spikeTrain(ibin,:) = sum(es.spikeTrain(zeroCross(ibin):zeroCross(ibin+1),:),1);
    es_out.lick(ibin) = sum(es.lick(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.reward(ibin) = sum(es.reward(zeroCross(ibin):zeroCross(ibin+1)));
    % take mean
    es_out.sampleTimes(ibin)  = nanmean(es.sampleTimes(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.smthBallSpd(ibin)  = nanmean(es.smthBallSpd(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.traj(ibin)  = nanmean(es.traj(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.trajPercent(ibin)  = nanmean(es.trajPercent(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.trajspeed(ibin)  = nanmean(es.trajspeed(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.ballspeed(ibin)  = nanmean(es.ballspeed(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.smthTrajSpd(ibin)  = nanmean(es.smthTrajSpd(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.contrast(ibin)  = nanmean(es.contrast(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.start(ibin)  = nanmean(es.start(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.gain(ibin)  = nanmean(es.gain(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.blanks(ibin)  = nanmean(es.blanks(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.active(ibin)  = nanmean(es.active(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.rewardPos(ibin)  = nanmean(es.rewardPos(zeroCross(ibin):zeroCross(ibin+1)));
    es_out.roomLength(ibin)  = nanmean(es.roomLength(zeroCross(ibin):zeroCross(ibin+1)));
    % take mean and round
    es_out.iexp(ibin)  = round(nanmean(es.iexp(zeroCross(ibin):zeroCross(ibin+1))));
    es_out.trialID(ibin)  = round(nanmean(es.trialID(zeroCross(ibin):zeroCross(ibin+1))));
    es_out.outcome(ibin)  = round(nanmean(es.outcome(zeroCross(ibin):zeroCross(ibin+1))));
end

    function es_out = initialise_es(outSize, numCells)
        tmp = zeros(outSize, 1);
        es_out.mua         = tmp;
        es_out.smthBallSpd = tmp;
        es_out.sampleTimes = tmp;
        es_out.lick = tmp;
        es_out.reward = tmp;
        es_out.rewardPos = tmp;
        es_out.trialID = tmp;
        es_out.traj  = tmp;
        es_out.trajspeed = tmp;
        es_out.ballspeed = tmp;
        es_out.smthBallSpd = tmp;
        es_out.smthTrajSpd = tmp;
        es_out.contrast = tmp;
        es_out.start = tmp;
        es_out.gain = tmp;
        es_out.roomLength = tmp;
        es_out.blanks = tmp;
        es_out.active = tmp;
        es_out.outcome = tmp;
        es_out.roomLength = tmp;
        es_out.iexp = tmp;
        es_out.mua  = tmp;
        es_out.trajPercent = tmp;
        es_out.spikeTrain = zeros(outSize, numCells);
    end

end