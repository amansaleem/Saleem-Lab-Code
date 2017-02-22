function es = loadCA1V1data(iseries, istheta, area, smth_win, box)

if nargin<5
    box = 1;
end
if nargin<4
    smth_win = 0;
end

switch area
    case 'CA1'
        loadCA1 = true;
    case 'V1'
        loadCA1 = false;
    otherwise
        disp('Error, no such area!!!')
end

animalList = [{'M130920_BALL'}, {'M130918_BALL'}, ...
    {'M140501_BALL'}, {'M140501_BALL'}, {'M140501_BALL'}, ...
    {'M140501_BALL'}, {'M140502_BALL'}, {'M140502_BALL'}];
seriesList = [1025, 1030, ...
    530, 531, 601,...
    602, 603, 604];

seriesIdx = find(iseries==seriesList);

fileName = ['es_' animalList{seriesIdx} '_' num2str(iseries)];
[~,k] = system('hostname');
if ~istheta
    sampleRate = 60;
    if strcmp('Aman-PC',k(1:end-1))
        dirName = ['C:\Users\Aman\Dropbox\Work\Data\CA1V1' filesep 'norm'];
    else
        dirName = ['E:\Dropbox\Work\Data\CA1V1' filesep 'norm'];
    end
    %     dirName = ['C:\Users\aman\Dropbox\Work\Data\CA1V1' filesep 'norm'];
    if loadCA1
        fileName = [fileName '_CA1'];
    else
        fileName = [fileName '_V1'];
    end
    load([dirName filesep fileName]);
elseif istheta==1
    sampleRate = 7.5;
    if strcmp('Aman-PC',k(1:end-1))
        dirName = ['C:\Users\Aman\Dropbox\Work\Code\Decoding\Data'];
    else
        dirName = ['E:\Dropbox\Work\Code\Decoding\Data'];
    end
    %     dirName = ['C:\Users\aman\Dropbox\Work\Code\Decoding\Data'];
    %     dirName = ['C:\Users\aman\Dropbox\Work\Data\CA1V1' filesep 'theta'];
    if loadCA1
        if iseries==531
            fileName = [fileName '_thetaBinsB'];
        else
            fileName = [fileName '_thetaBinsA'];
        end
        %         fileName = [fileName '_thetaA_CA1'];
    else
        if strcmp('Aman-PC',k(1:end-1))
            dirName = ['C:\Users\Aman\Dropbox\Work\Data\CA1V1' filesep 'theta'];
        else
            dirName = ['E:\Dropbox\Work\Data\CA1V1' filesep 'theta'];
        end
        %         dirName = ['C:\Users\aman\Dropbox\Work\Data\CA1V1' filesep 'theta'];
        fileName = [fileName '_thetaA_V1'];
    end
    load([dirName filesep fileName]);
elseif istheta==2
    sampleRate = 7.5;
    if strcmp('Aman-PC',k(1:end-1))
        dirName = ['C:\Users\Aman\Dropbox\Work\Code\Decoding\Data'];
    else
        dirName = ['E:\Dropbox\Work\Code\Decoding\Data'];
    end
    %     dirName = ['C:\Users\aman\Dropbox\Work\Code\Decoding\Data'];
    if loadCA1
        if iseries==531
            fileName = [fileName '_thetaBinsA'];
        else
            fileName = [fileName '_thetaBinsB'];
        end
        %         fileName = [fileName '_thetaA_CA1'];
    else
        if strcmp('Aman-PC',k(1:end-1))
            dirName = ['C:\Users\Aman\Dropbox\Work\Data\CA1V1' filesep 'theta'];
        else
            dirName = ['E:\Dropbox\Work\Data\CA1V1' filesep 'theta'];
        end
        %         dirName = ['C:\Users\aman\Dropbox\Work\Data\CA1V1' filesep 'theta'];
        fileName = [fileName '_thetaA_V1'];
    end
    load([dirName filesep fileName]);
end

if smth_win>0
    spkRate = zeros(size(es.spikeTrain));
    for icell = 1:size(es.spikeTrain,2)
        if box
            spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), ...
                sampleRate, smth_win, 'box');
        else
            spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), ...
                sampleRate, smth_win);
        end
    end
    es.spikeTrain = spkRate;
    es.sampleRate = sampleRate;
    es.mua = sum(es.spikeTrain,2);
end