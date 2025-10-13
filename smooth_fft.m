function spec_smooth = smooth_fft(field,spec,N)

%% Correction of linear baseline

    p = round(length(spec)/100);
    x_min = mean(field(1:p-1));
    y_min = mean(spec(1:p-1));
    x_max = mean(field(end-p:end));
    y_max = mean(spec(end-p:end));
    factor = (y_max - y_min)/(x_max - x_min);
    x_lin = linspace(field(1),field(end),length(field))';
    y_lin = factor*(x_lin - x_min) +y_min;
    spec_dif = spec - y_lin;
    
%% FFT
    L = length(field);

    spec_fft_2s = fft(spec_dif);
    spec_fft_1s = spec_fft_2s(1:L/2+1);
    spec_fft_1s(2:end-1) = 2*spec_fft_1s(2:end-1);
    
    
    f_ref = (1/(field(2) - field(1)))/2;
    f = f_ref*(0:(L/2))/L;
    freq = ((0:length(f)-1)*f_ref)';
    cut_freq = f_ref/(2*N);
    
    
    for i=1:length(f)
        spec_fft_parabole(i,1) = spec_fft_1s(i,1)*(1 - ( ((f(1,i)).^2) / ((cut_freq)^2) ));
        if f(1,i) > cut_freq
            spec_fft_parabole(i,1) = 0;
        end
    end

    spec_smooth = real(ifft(spec_fft_parabole,L));