%% Discretizing a sample of n survival times into sub-intervals with no. of 'binCount' subjects in each sub-interval
function [ID,binwv] = discrBinNA(T,binCount)
len = length(T); % sample size
binid = 1:binCount:len;
nbins = length( binid );
binedg = T(binid);

if binCount == 1
    binwv = binedg(1:nbins) - [0,binedg(1:(nbins-1))]; %binwv: bin widths
    binedg = [0;binedg];
    binID = 1:n; %binID: bin IDs
else
    t = 0;
    while  t < ( nbins - 1 )
        t = t + 1;
        ntied = sum( binedg(t) == binedg );
        if ( t + ntied ) <= nbins
            binedg = [binedg(1:t); binedg((t + ntied):nbins)];
        else
            binedg = binedg(1:t);
            nbins = length( binedg );
        end

        if binedg(nbins) == max(T)
            nbins = nbins - 1;
            binedg = binedg(1:nbins);
        end
        binedg(1) = min(T);
        binedg(nbins + 1) = max(T) + (1e-2)*( binedg(nbins) - binedg(nbins-1) );
        binedg(nbins) = binedg(nbins)*(1- 1e-5 );
        binwv = binedg(2:(nbins+1)) - binedg(1:nbins);
        binID = zeros(1,len);
        for  i = 1:1:len
            for  j = 1:nbins
                if T(i) < binedg(j+1) && T(i) >= binedg(j)
                    binID(i) = j;
                end
            end
        end
    end
end

ID = binID;
% binedg = binedg;

end



