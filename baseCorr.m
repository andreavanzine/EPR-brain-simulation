function [spec,SecInt,slopeSecInt] = baseCorr(field,spec_raw,method,n_end,n_end_int,tol,step)
% e.g.: [ specFt_Corrected, SecondInt, SlopeSecInt ] = LinCorr_2ndInt( data.B,data.subtr_rawbuf_sample,20,120,0.2,0.0001 )
% author: J.A. Labra-Munoz
% This function is based on the script BaselineCorrection_SecondIntegral,
% but instead of manually looking for the optimal linear slope so that the
% last 120 points of the 2nd integral reach a zero-slope plateau, here it
% automatically search for that point.
% It assumes that the input spectrum is the result of Ft-Buffer spectra.
% n_end: Last # of points to consider for linear correction of raw spectrum
% n_end_int: Last # of points to consider for the 2nd int plateu.
%
% Modified: Fabio S. Otsuka
%   - included an option to either select the number of end points, or to
%     set a specified field value and extract all the points according to
%     that field value
%   - this change was made into the variable 'method' which should be
%     called by 'point' or 'range'
%

% defining line
x_init = 0;  %for starting line reference
y_init = 0;   % -0.01  -4
flag = 1;
iter = 1;
slope = 0;

while flag
    
    if strcmp(method,'point')
        x_end = mean(field(end - n_end:end));
        y_end = mean(spec_raw(end - n_end:end));
        x_lin = linspace(field(1),field(end),length(field))';
        factor = (y_end - y_init)/(x_end - x_init);
        y_lin = x_lin*factor -x_init*factor +y_init;
        spec_dif = spec_raw - y_lin;
    elseif strcmp(method,'range')
        x_end = mean(field(field>n_end));
        y_end = mean(spec_raw(field>n_end));
        x_lin = linspace(field(1),field(end),length(field))';
        factor = (y_end - y_init)/(x_end - x_init);
        y_lin = x_lin*factor -x_init*factor +y_init;
        spec_dif = spec_raw - y_lin;
    end

    z1 = cumtrapz(field,spec_dif); %First integral
    z2 = cumtrapz(field,z1); %Second integral

    %Checking slope of last part of 2nd integral
    if strcmp(method,'point')
        xLin2Int = field(end - n_end_int:end); 
        yLin2Int = z2(end - n_end_int:end);
        pLin2Int = polyfit(xLin2Int,yLin2Int,1);
        linfit2Int = pLin2Int(1)*field + pLin2Int(2);
    elseif strcmp(method,'range')
        xLin2Int = field(field>n_end_int); 
        yLin2Int = z2(field>n_end_int);
        pLin2Int = polyfit(xLin2Int,yLin2Int,1);
        linfit2Int = pLin2Int(1)*field + pLin2Int(2);
    end
%     linfit2Int
    if abs(pLin2Int(1))<=tol  
        flag = 0;
        spec = spec_dif;
        SecInt = z2;
        slopeSecInt = pLin2Int(1);
    end
    
    if slope*pLin2Int(1) < 0
        step = -step*1e-1;
    end
    
    if abs(slope) < abs(pLin2Int(1))
        step = -step;
    end
    
    y_init = y_init - step;
    slope = pLin2Int(1);
end

end

