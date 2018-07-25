x1 = [1:0.5:10000];
x2 = [1:0.5:100];
sin1 = sin(x1);
sin2 = sin(x2);

spectrogram(sin1');
spectrogram(sin2');
