function [output]           = special_smooth_3d(input, win, flag_sp, bins1, bins2, bins3, nGrid)
        if nargin<3
            flag_sp = [0 0];
        elseif (nargin <4 | isempty(bins1) | isempty(bins2) | isempty(bins3))
            flag_sp = [0 0 0];
            display('Do not have speed bins to do special smoothing');
        end
        output = zeros(size(input));
        
        nGrid1 = nGrid(1);
        nGrid2 = nGrid(2);
        nGrid3 = nGrid(3);
        
        % do the smoothing
        r1 = (-nGrid1:nGrid1)/nGrid1;
        Smoother1 = exp(-r1.^2/win(1)^2/2);
        
        r2 = (-nGrid2:nGrid2)/nGrid2;
        Smoother2 = exp(-r2.^2/win(2)^2/2);
        
        r3 = (-nGrid3:nGrid3)/nGrid3;
        Smoother3 = exp(-r3.^2/win(3)^2/2);
        Smooth3 = repmat(Smoother3, [size(Smoother1,2)*size(Smoother2,2) 1]);
        Smooth3 = reshape(Smooth3, [size(Smoother1,2) size(Smoother2,2) size(Smooth3,2)]);
        
        smoother = Smoother1'*Smoother2;
        smoother = repmat(smoother,[1 1 size(Smooth3,3)]);
        smoother = smoother.*Smooth3;
        
        if win(1)==1
            trash = nanmean(input,1);
            input = repmat(trash,[size(input,1) 1 1]);
        end
        if win(2)==1
            trash = nanmean(input,2);
            input = repmat(trash,[1 size(input,2) 1]);
        end
        if win(3)==1
            trash = nanmean(input,3);
            input = repmat(trash,[1 1 size(input,3)]);
        end
        
        if ~flag_sp
            output = convn(input, smoother, 'same');
        else
            output = convn(input, smoother, 'same');
        end
    end