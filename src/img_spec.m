function img_spec(sig,n,o,F,fs,range)
% img_spec(sig,n,o,F,fs,range)
% sig is the wavfile
% n is the window size
% o is the overlap
% f is the vector
% fs is the samplerate
% range is range
% exmaple Img_spec(sig,128,127,4096,fs,100)

[B,f,t]=spectrogram(sig,n,o,F,fs);
bmin=max(max(abs(B)))/range;
imagesc(t,f,20*log10(max(abs(B),bmin)/bmin));
set(gca,'YDir','normal');
colormap(flipud(hot))