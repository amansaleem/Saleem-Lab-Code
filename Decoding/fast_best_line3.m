function [ bestRes, bestSpd, bestYInterc ] = fast_best_line( pMat )
%FAST_BEST_LINE Find best straight trajectory through postprob matrix
% Decoder returns a posterior probability matrix [pos bins [1cm] by time
% bins [1ms], typically [100,300] i.e. 100cm by 300ms]. Each bin is the
% prob the rat was at that pos at that time. Goal is to find the best
% constant speed trajactory (straight line) path through this matrix i.e.
% the line that is highest probability.
%
% Method is simple - for a given line sum the prob in the bins within a
% certain distance of the line (can be graded - see below), then
% iterate for different lines and find the
% one that includes the highest sum (actually mean) probability.
%
% Implementation is done with conv2 for speed - allows all offsets for a
% given line gradient to be tested at the same time. Note only lines of a
% given range of gradient are tested and only certain offsets are valid. A
% final criteria is that lines must include sufficient bins to be
% considered trustworthy (i.e. don't accept lines that just 'clip' the
% postProb matrix).
%
% DEFINING ROI AROUND LINE
% See inline function (around line 220) for code that defines the line with
% ROI - this can either be a hard cut off or graded. If the latter then the
% yRange needs to be large (100). Comment lines to change behaviour.


% ARGS              [note most args are now hard coded]
% pMat              Posterior prob matrix with dim1 being spatial bins and
%                   dim2 being time. This is output of decoder indicating
%                   for each time bin the prob animal is at given point.
%                   Columns should sum to 1. [nPosBin x nTimeBin]
%

%
% spdBnd            max speeds to consider in cm/s e.g. [200]. This defines
%                   the range of gradients that are tested. Specifically
%                   gradients will run from minus this number to positive
%                   this number e.g. for 200 we test 'speeds' of -200cm/s
%                   to 200cm/s.
%

%
% RETURNS
% bstGrad           gradient of best fit line (see fitLine_CB2 for details
%                   units) but basically is in timeBin/posBin
%
% bstYInt           yIntercept of best fit line (see fitLine_CB2) but
%                   basically is posBin of intercept
%
% bstRes            residual or probability of best fit line
%
% e.g.
%  [ bstGrd, bstInt, bstRes ] = fast_best_line( pMat, 0.001, 1, 100 )


% --- DEFINE HARD CODED VARS --------------------------------------------
% Temporal bin size in s (temporal axis is the y-axis).
tBinSz          =0.001; %[0.001] i.e. 1ms

% Spatial bin isze in cm, likely to be 1cm
sBinSz          =1;

% Min and max speed to consider in cm/s e.g. [200] would be 200cm/s.
% Effecitvly defines the range of gradients that are tested. Specifically
% gradients will run from minus this number to positive this number e.g.
% for 200 we test 'speeds' of -200cm/s to 200cm/s.
spdBnd          =1500;

% Also define the increment step to iterate through speed e.g. should we
% test speeds every 1cm/s (i.e. -200, -199, -198) or every 2cm/s or 5cm/s.
% Larger numbers will give fewer increments (faster) but lower fidelity.
spdInc          =20;

% Define size of 'band' around line that is used to calculated it's
% goodness of fit. This value should be an integer and indicates the number
% of bins above and below (in the y-axis) that are counted. Must be >0
yRange          =100; %[100 of graded] [15 for hard]

% Also have a requirement for final line to be some minimum length - to
% avoid situations where a line just cuts through a small section of the
% pMat and scores a (spurious) high score. Lines of less than this number
% are discarded.
minLgth         =150;



% --- HOUSE KEEPING -----------------------------------------------------
% Use parallel loops (parfor) to speed up calculation - first check if
% matlabpool is started and if not start one
if ~(matlabpool('size') > 0)
    matlabpool local;
end


% --- MAIN CODE ---------------------------------------------------------
% Size of pMat - used several times
% szPM1           =size(pMat,1);
szPM2           =size(pMat,2);

% Deterine the range of gradient that will need to be tested and
% express in terms of sBins per tBin
spd2Test        =-spdBnd:spdInc:spdBnd;
spd2Test        =(spd2Test ./sBinSz) .*tBinSz;
nSpd2Test       =length(spd2Test);
clear           tBinSz sBinSz spdBnd


% Now loop over each gradient to test and in each loop test all possible
% offsets for a line of that gradient. First do some preallocation for
% speed. Each loop returns the mean probability of bins in the pMat covered
% by the line and the band around it as well as the raw number of bins
% covered. These are stored in 2d mats of size [nOffsets tested x nSpds
% tested]
[meanP, bestOffSet]           ...
    =deal (zeros(nSpd2Test,1)); %nOffset x nSpds tested.

pMatOnes        =ones(size(pMat)); %Mat of ones the same size as pMat
parfor nn       =1:nSpd2Test
    
    %Main step is to calculate the sum of pMat bins that overlap with a
    %line of specific gradient for all allowed offsets. Start by defining
    %line and the band around it as a mat.
    
    
    %Use an inline function to define the line - partly to keep code tiday
    %and because it's called in multiple places. The two mats returned are
    % a zeros and ones mat defining just the line (tstLn) and a similar mat
    %defining the line plus band arround it defined by yRange (tstLnBnd).
    %Each has dim [variables x szPM2] i.e. second dim is same as pMat
    [tstLn, tstLnBnd]   =il_define_line(szPM2, spd2Test(nn), yRange);
    

    %Now do the main step - use filter2 to convolve the tstLn with the pMat
    %so effectivly testing the overlap between the line and pMat at all
    %offsets - to make this step faster just use the 'valid' opperator but
    %first zero pad dim1 of pMat
    tmpMeanP        =filter2(tstLnBnd, padarray(pMat,[size(tstLn,1),0]), 'valid');
    
    %Filter2 just gives a sum - so normalise by dividing by the number of
    %bins included in the sum
    tmpNBinOverLap  =filter2(tstLnBnd, padarray(pMatOnes,[size(tstLn,1),0]), 'valid');
    
    %Use similar method to get the actual lenght of the line
    tmpLineLength   =filter2(tstLn, padarray(pMatOnes,[size(tstLn,1),0]), 'valid');
    
    %Now discard data points corresponding to line that is too short and
    %store the details of the line that gives the highest fit.
    tmpMeanP        =tmpMeanP./tmpNBinOverLap; %Normalise by n bins
    tmpMeanP(tmpLineLength<minLgth)  =0; %Set ones that are too short to 0
    [meanP(nn), bestOffSet(nn)]  =max(tmpMeanP); %Store max and ind.
end
clear tmp* tst* nn

%Finally from the result of the loop extract the details of the line that
%gave the best fit i.e. the value of the fit, the gradient and intercept.

%NL. bestRes is the highest mean prob per bin of the overlap between pMat
%and line found. tmpInd is the index into meanP (i.e. loop index) that gave
%this result.
[bestRes, tmpInd]   =max(meanP);
%bestSpd is the spd (still in sBins per tBins) that gave the highest result
bestSpd             =spd2Test(tmpInd); %Speed in sBins per tBin

%Offset is harder to understand - this number is the offset in the y-dim
%between the pMat and tstLnBnd during filter2 that have the highest result.
%Need to conver this to an intercept on the y-axis. To do this need to know
%the size of tstLnBnd - so recreate this
[tstLn, ~]          =il_define_line(szPM2, spd2Test(tmpInd), yRange);
bestYInterc         =bestOffSet(tmpInd) - size(tstLn,1) + find(tstLn(:,1));


% If no line is found - e.g. if the pMat is too short (in time) to support
% any valid lines (all lines are too short) or if there is no data. Then
% bestRes will be 0. Catch these and return nans
if bestRes == 0 %No fit (probably) - so return nan
    [bestRes, bestSpd, bestYInterc]         =deal(nan);
end




% --- DEBUG STUFF --------------------------------------------------------
%Can comment all this out once code is functioning correctly

% subplot(1,3,1) %First show pMat
% imagesc(pMat);
% hold on
% plot([1, szPM2], [bestSpd + bestYInterc, bestSpd * szPM2 + bestYInterc]);
% hold off
% 
% subplot(1,3,2) %Second best fit line
% imagesc(tstLn);
% 
% subplot(1,3,3); %Third the best fit line with the band around it
% imagesc(tstLnBnd);

end



% --- INLINE FUNCTIONS ---------------------------------------------------
function [tstLn, tstLnBnd]   =il_define_line(szPM2, spd2Test, yRange)
%Inline function that returns the logical matrix to define the line
%(trajecotry through pMat) of interest as well as the line plus band.

%Define line and band - y extent of line is equivalent to gradient *
%size(pMat,2) + the  band width. Define that line as a set of 1s in
%a mat of zeros. Then expand to include the band by filter2 with a
%kernel of the band width. Note always build the line as a +vs gradient
%and flip if gradient is -ve

if ceil(abs(spd2Test)) == 0 %If gradient is 0 hard code tstLn
    tstLn           =ones(1, szPM2);
    
else %Otherwise do the normal thing
    tstLn           =zeros(ceil(abs(spd2Test).* szPM2), szPM2);
    %NL for each x bin determin the y bin of the line based on gradient
    tstLnY          =ceil([1:szPM2] .* abs(spd2Test));
    tstLn(sub2ind(size(tstLn), tstLnY', (1:szPM2)'))        =1;
end

% Create the kernel that is used to define the ROI around the line - this
% might be a hard cut off or a graded area. Comment following lines to
% control this.
% line - this is convolved with the line to get the graded matrix.
% COMMENT ONE OF THE NEXT TWO LINES
kern            =[1:yRange, fliplr(1:yRange-1)]'; %Graded
% kern            =ones(2*yRange -1, 1);             %Hard cut off


tstLn           =padarray(tstLn, [yRange,0]);
tstLnBnd        =filter2(kern, tstLn); %Now line with boundary
if sign(spd2Test) == -1;
    tstLn=flipud(tstLn);
    tstLnBnd=flipud(tstLnBnd);
end
end

% imagesc(pMat);
% pause;
% close all;