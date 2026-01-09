function [dur,pres] = percentdurRMS(seg,frac)
for n = 1:length(seg(1,:))
    [ste(:,n),rd(:,n)] = percentdur(seg(:,n),frac);
end
[~,l] = max(rd);

pres = rms(seg(ste(1,l):ste(2,l),:));
dur = [ste(1,l) ste(2,l)];
