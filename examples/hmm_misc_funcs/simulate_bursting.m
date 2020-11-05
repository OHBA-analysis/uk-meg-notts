function [data_sim, Gamma_true]=simulate_bursting(S)

a = -poly(.82);

seconds = 40;
sample_rate = S.sample_rate;
time_vect = linspace(0,seconds,seconds*sample_rate)';

noise = S.noise_amp*randn(size(time_vect));
%noise = noise .* (.3.*sin(2*pi*.1*time_vect)+1);
x = 1*filter(1,a,noise);

starts = [100 500 900 1500 1800 2400 2800 3500 4000 4300 4800];
starts2 = [100 500 900 1500 1800 2400 2800 3500 4000 4300 4800] + 5120;
duration = [150 100 50 100 50 200 50 100 100 200 50]*1.5;
starts3 = starts + 200;

starts=round(sample_rate*starts/512);
starts2=round(sample_rate*starts2/512);
starts3 = round(sample_rate*starts3/512);
duration=round(sample_rate*duration/512);

x2 = zeros(size(x));
x3 = zeros(size(x));
x4 = zeros(size(x));

Gamma_true = zeros(size(x,1),4);
Gamma_true(:,1)=1;

amp_fast=S.high_freq_amp;
amp_slow=S.low_freq_amp;

for ii = 1:length(starts)
    % Add slow burst
    tmp = amp_slow*sin( 2*pi*20*time_vect(starts(ii):starts(ii)+duration(ii)));
    tmp = tmp.*tukeywin(length(tmp),.25);
    x2(starts(ii):starts(ii)+duration(ii)) = x2(starts(ii):starts(ii)+duration(ii)) + tmp;
    Gamma_true(starts(ii):starts(ii)+duration(ii),2)=1;
    Gamma_true(starts(ii):starts(ii)+duration(ii),1)=0;
    
    % Add fast burst
    tmp = amp_fast*sin( 2*pi*40*time_vect(starts2(ii):starts2(ii)+duration(ii)));
    %tmp = randn( size(tmp) ); % broadband noise
    tmp = tmp.*tukeywin(length(tmp),.25);
    x3(starts2(ii):starts2(ii)+duration(ii)) = x3(starts2(ii):starts2(ii)+duration(ii)) + tmp;
    Gamma_true(starts2(ii):starts2(ii)+duration(ii),3)=1;
    Gamma_true(starts2(ii):starts2(ii)+duration(ii),1)=0;
    
    % Add second slow burst
    tmp = amp_slow*sin( 2*pi*30*time_vect(starts3(ii):starts3(ii)+duration(ii)));
    %tmp = randn( size(tmp) ); % broadband noise
    tmp = tmp.*tukeywin(length(tmp),.25);
    x4(starts3(ii):starts3(ii)+duration(ii)) = x4(starts3(ii):starts3(ii)+duration(ii)) + tmp;
    Gamma_true(starts3(ii):starts3(ii)+duration(ii),4)=1;
    Gamma_true(starts3(ii):starts3(ii)+duration(ii),1)=0;
    
end

data_sim = x'+x2'+x3'+x4';

end

