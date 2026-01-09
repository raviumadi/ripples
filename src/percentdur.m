function [duration,pres, fmax] = percentdur(seg,frac)
sm = 2*floor((length(seg)/10)/2)+1;
csig = cumsum(seg.^2);
csig(1) = 0;
csig = csig/max(csig);
csig = smooth(smooth(csig,sm),sm);
int = find(csig>frac(1) & csig< frac(2));
if isempty(int)
    duration = [1 2];
    pres = 0;
else
    duration = [int(1) int(end)];
    pres = rms(seg(int));
    mag = 20*log10(abs(fft(seg(int(1):int(end)), 192))); % each mag value is 1kHz
    [~, fmax] = max(mag(1:length(mag)/2));

end
