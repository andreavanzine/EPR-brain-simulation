%% Human Brain EPR spectrum processing pipeline
%
% A routine to process the EPR spectrum of human brain samples
% It consists of the following steps:
%   1) pre-processing
%       - spectrum load
%       - EPR cavity signal subtraction and baseline correction
%
%   2) fitting
%       - 3 peaks model (each peak is fitted individually)
%       - estimation of the weight of each peak
%       - perform the simulation with optimized parameters
%
%   3) analysis
%       - calculation of peak-to-peak linewidth and 2nd integral
%
% 1st version: Fabio Seiji Otsuka (2022)
% 2nd version: André Avanzine (2023)

%% 1) pre-processing

% Data loading
data_dir = '';
data_dir_cav = '';

% If using EPR file
[field spec_raw] = eprload(data_dir);
[fieldT1 specT1] = eprload(data_dir_cav);
                                                                                                                                                                                                                                                  
spec_final = spec_raw - specT1;

% make sure that we get the real part of the data
field = real(field)';
spec_final = real(spec_final);

% Obs. if 'field' is in gauss, you'll need to convert to mT by dividing by
% 10:
% field = field/10;

% Baseline correction (Jacqui's method)
[spec, SecondInt, SlopeSecInt] = baseCorr(field,spec_final,'point',20000,20000,0.25,1e-2); 

figure; plot(field, spec,'black'); %hold on; plot(field, spec_final);

% 2) fitting

n = 1; %select which spectra to fit
spec_n = spec(:,n);

% Experimental Parameters
Exp.Range = [min(field) max(field)]; % in mT
Exp.Temperature = 193.15; % in K
Exp.mwFreq = 9.15; % in GHz
Exp.nPoints = length(field); % number of points (not obrigatory)

% 2.1) g = 4.3
%define field range for g = 4.3
field1 = field(field<225);
spec1 = spec_n(field<225);

Exp1.Range = [min(field1) max(field1)];
Exp1.Temperature = Exp.Temperature;
Exp1.mwFreq = Exp.mwFreq;
Exp1.nPoints = length(field1);

Sys1.S = 5/2;
Sys1.g = [1.818 1.955 2.129];
Sys1.gStrain = [0.222 0.042 0.388];
Sys1.D = [2.5326e+04 8.7066e+03]; %in MHz

Vary1.g = [0.05 0.05 0.05]; %0.01 0.01 0.01
Vary1.gStrain = [0.1 0.1 0.1]; %0.05 0.01 0.05
Vary1.D = [5 0.5]*1e3;

FitOpt.Method = 'simplex fcn';
FitOpt.Scaling = 'lsq2';
FitOpt.maxTime = 3; %minutes

%%
esfit('pepper',spec1,Sys1,Vary1,Exp1,[],FitOpt);

%% After fitting, export the result and update Sys1:

Sys1.g = fit1.bestvalues(1:3);
Sys1.gStrain = fit1.bestvalues(4:6);
Sys1.D = fit1.bestvalues(7:8);

%% Get the weighting factor "w1" for the simulated peak:
w_ref = (fit1.expSpec' * fit1.expSpec) \ (fit1.expSpec' * spec1);

[x y] = pepper(Sys1,Exp1);
y = transpose(y);

y_diff = gradient(y)./gradient(field1);
y_diff2 = gradient(y_diff)./gradient(field1);
y_diff_ref = gradient(transpose(fit1.fitSpec))./gradient(field1);
y_diff_ref2 = gradient(y_diff_ref)./gradient(field1);
w_diff = (y_diff2' * y_diff2) \ (y_diff2' * y_diff_ref2);

w1(n) = w_diff*w_ref;

clear x y

% Simulate over the full range:
[x y] = pepper(Sys1,Exp);
Fehs(:,n) = transpose(y);
clear x y

% 2.2) copper

% Subtract g = 4.3 peak from the experimental spectrum
spec_sub1 = spec_n - w1(n)*Fehs(:,n);

% You can also perform the smoothing at OriginLab using the following method:
% Analysis -> Signal Processing -> Smooth -> FFT Filter
% (check "Auto Preview" in order to select the best parameter for
% "Points in Window")
spec_smooth = smooth_fft(field, spec_sub1,5000);

% Check if agrees
%plot(field,spec_sub1); hold on; plot(field,spec_smooth); 
spec_test = spec_sub1 - spec_smooth; % removing the smoothed spectra

% Crop the field range to show only the copper peaks
field2 = field(220<field);
spec2 = spec_test(220<field);
field2 = field2(field2<370);
spec2 = spec2(field2<370);

Exp2.Range = [min(field2) max(field2)];
Exp2.Temperature = Exp.Temperature;
Exp2.mwFreq = Exp.mwFreq;
Exp2.nPoints = length(field2);

Sys2.S = 1/2;
Sys2.Nucs = 'Cu';
Sys2.g = [2.0389 2.2194];
Sys2.A = [38.8917 550.1459];
Sys2.gStrain = [0.0310 0.0527];

Vary2.g = [0.05 0.05];
Vary2.A = [5 20];
Vary2.gStrain = [0.02 0.02];

FitOpt.Scaling = 'lsq0';
FitOpt.maxTime = 3; %minutes

%%
esfit('pepper',spec2,Sys2,Vary2,Exp2,[],FitOpt);

%% After fitting, export the result and update Sys2:
Sys2.g = fit1.bestvalues(1:2);
Sys2.gStrain = fit1.bestvalues(5:6);
Sys2.A = fit1.bestvalues(3:4);

%% Get the weighting factor "w2" for the simulated peak:

w_ref = (fit1.expSpec' * fit1.expSpec) \ (fit1.expSpec' * spec2);

[x y] = pepper(Sys2,Exp2);
y = transpose(y);

y_diff = gradient(y)./gradient(field2);
y_diff2 = gradient(y_diff)./gradient(field2);
y_diff_ref = gradient(transpose(fit1.fitSpec))./gradient(field2);
y_diff_ref2 = gradient(y_diff_ref)./gradient(field2);
w_diff = (y_diff2' * y_diff2) \ (y_diff2' * y_diff_ref2);

w2(n) = w_diff*w_ref;

clear x y

% Simulate over the full range:
[x y] = pepper(Sys2,Exp);
Cu(:,n) = transpose(y);
clear x y

% 2.3) g = 2.01 (broad peak)

% Subtract g = 4.3 and copper peaks from the experimental spectrum
spec3 = spec_sub1 - w2(n)*Cu(:,n);
%spec3 = spec_n - w1(n)*Fehs(:,n);

% Here we use the two spin system's model (Sys3a and Sys3b) to simulate the
% broad peak:

Sys3a.S = 10; 
Sys3a.g = 2.01; 
Sys3a.D = 155.0974; % in MHz 
Sys3a.HStrain = 1509.3; % in MHz 
Sys3a.weight = 0.7225;  

Vary3a.g = 0.05; 
Vary3a.D = 50;
Vary3a.HStrain = 1000;
Vary3a.weight = 0.15;

Sys3b.S = 10;  
Sys3b.g = 2.01; 
Sys3b.D = 244.6983; % in MHz 
Sys3b.HStrain = 3.9319e+03; % in MHz 
Sys3b.weight = 1.86;

Vary3b.g = 0.05;
Vary3b.D = 50;
Vary3b.HStrain = 1000;
Vary3b.weight = 0.15;

FitOpt.Scaling = 'lsq2';
FitOpt.maxTime = 10; %minutes

%%
esfit('pepper',spec3,{Sys3a,Sys3b},{Vary3a,Vary3b},Exp,[],FitOpt);

%% After fitting, export the result and update Sys3a and Sys3b:
Sys3a.g = fit1.bestvalues(1);
Sys3a.D = fit1.bestvalues(2);
Sys3a.HStrain = fit1.bestvalues(3);
Sys3a.weight = fit1.bestvalues(4);

Sys3b.g = fit1.bestvalues(5);
Sys3b.D = fit1.bestvalues(6);
Sys3b.HStrain = fit1.bestvalues(7);
Sys3b.weight = fit1.bestvalues(8);

%% Get the weighting factor "w3" for the simulated peak:
w_ref = (fit1.expSpec' * fit1.expSpec) \ (fit1.expSpec' * spec3);

[x y] = pepper({Sys3a,Sys3b},Exp);
y = transpose(y);

y_diff = gradient(y)./gradient(field);
y_diff2 = gradient(y_diff)./gradient(field);
y_diff_ref = gradient(transpose(fit1.fitSpec))./gradient(field);
y_diff_ref2 = gradient(y_diff_ref)./gradient(field);
w_diff = (y_diff2' * y_diff2) \ (y_diff2' * y_diff_ref2);

w3(n) = w_diff*w_ref;

Fh(:,n) = y;

% Sample dry mass used in the EPR 
mass = 1; % In miligrams

spec_FeHs = w1*Fehs/mass;
spec_Cu = w2*Cu/mass;
spec_Ft = w3*Fh/mass;

% Absorption spectra of each component
first_int_Fehs = cumtrapz(field, spec_FeHs);
first_int_Cu = cumtrapz(field, spec_Cu);
first_int_Ft = cumtrapz(field, spec_Ft);

sec_int_Fehs = cumtrapz(field, first_int_Fehs);
sec_int_Cu = cumtrapz(field, first_int_Cu);
sec_int_Fh = cumtrapz(field, first_int_Ft);

area_Fehs = trapz(field, sec_int_Fehs);
area_Cu = trapz(field, sec_int_Cu);
area_Fh = trapz(field, sec_int_Ft);

% Absorption value of each component
susc_Fehs = trapz(field, cumtrapz(field, spec_FeHs));
susc_Cu = trapz(field, cumtrapz(field, spec_Cu));
susc_Ft = trapz(field, cumtrapz(field, spec_Ft));

table1(n,1:2) = Sys1.D;
table1(n,3) = Sys1.S;
table1(n,4:6) = Sys1.gStrain;
table1(n,7:9) = Sys1.g;
table2(n,1:2) = Sys2.A;
table2(n,3) = Sys2.S;
table2(n,4:5) = Sys2.gStrain;
table2(n,6:7) = Sys2.g;
if exist('Sys3a') && exist('Sys3b')
    table3(n,1) = Sys3a.D;
    table3(n,2) = Sys3a.HStrain;
    table3(n,3) = Sys3a.g;
    table3(n,4) = Sys3a.weight;
    table3(n,5) = Sys3b.D;
    table3(n,6) = Sys3b.HStrain;
    table3(n,7) = Sys3b.g;
    table3(n,8) = Sys3b.weight;
else
    table3(n,1:8) = 0;
end
table4 = [susc_Fehs, susc_Cu, susc_Ft];
    
figure, plot(field,spec_n); hold on; plot(field,w1(n)*Fehs(:,n) + w2(n)*Cu(:,n) + w3(n)*Fh(:,n));

%% Save all the simulation results
save
clc