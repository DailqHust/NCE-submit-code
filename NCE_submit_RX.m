% Copyright (c) 2023, Longquan Dai, Li Di, Lei Deng in Wuhan National Laboratory 
% for Optoelectronics and School of Optical and Electronic Information, 
% Huazhong University of Science and Technology, Wuhan 430074, China
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation 
%    and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its contributors may
%    be used to endorse or promote products derived from this software without 
%    specific prior written permission.
%
% 4. In case results are published that rely on this source code, please cite
% our paper published in <Nature Communications Engineering> entitled 
% "Frequency-dependent Impairment Calibration and Estimation for a 96 GBaud Coherent Optical Transceiver" [1]. 

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
% OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
% [1] https://www.nature.com/articles/...
%% 
clear
close all
% load calibration data
prename = ['TRx_calibration_data'];
load (['D:\simulation_out\',prename,'.mat']);
y_rx_XI = y_Rx(1,:);  
y_rx_XQ = y_Rx(2,:);
y_rx_YI = y_Rx(3,:);
y_rx_YQ = y_Rx(4,:);
% multi-tone parameter
Fdso=120e9;
N_Sample = 512*500;  % signal length
Ntone = 64;  % The number of tones in each tributary
f_baud = 234375*2 ;    % base frequency
sub_f(1) = f_baud * 4*1;   % Original value of sub-frequency interval 
f00 = f_baud*1600 ;  % Original value of frequency interval
df = sub_f:sub_f:sub_f*Ntone;  
for i = 1:Ntone-1
   sub_f(i+1) =  sub_f(i)+ df(i);   % sub-frequency interval 
end
for i = 1:Ntone
   f0(i) = f00*i+sub_f(i); %  Original value of frequency interval
end
delta_f0 = 0*(f00+df)/4;
delta_f1 = 1*(f00+df)/4;
delta_f2 = 2*(f00+df)/4;
delta_f3 = 3*(f00+df)/4;
% frequency of XI/XQ/YI/YQ
for i = 1:n
    fXI(i)=f0(i); 
    fXQ(i)=f0(i)+delta_f2(i);
    fYI(i)=f0(i)+delta_f1(i); 
    fYQ(i)=f0(i)+delta_f3(i); 
end
% signal capture
[S_II_freq,S_II_amp,S_II_pha] = FFt_capture(y_rx_XI,Fdso,fXI); % amplitude and phase of fXI in y_rx_XI          
[S_QI_freq,S_QI_amp,S_QI_pha] = FFt_capture(y_rx_XI,Fdso,fXQ); % amplitude and phase of fXQ in y_rx_XI
[S_IQ_freq,S_IQ_amp,S_IQ_pha] = FFt_capture(y_rx_XQ,Fdso,fXI); % amplitude and phase of fXI in y_rx_XQ
[S_QQ_freq,S_QQ_amp,S_QQ_pha] = FFt_capture(y_rx_XQ,Fdso,fXQ); % amplitude and phase of fXQ in y_rx_XQ

% phase unwrap
rng(1)  
Phase_ini_XI = rand(1,Ntone)*2*pi-pi;
rng(2)
Phase_ini_XQ = rand(1,Ntone)*2*pi-pi; 

result_S_II_pha = S_II_pha-Phase_ini_XI;
result_S_QI_pha = S_QI_pha-Phase_ini_XQ;
result_S_IQ_pha = S_IQ_pha-Phase_ini_XI;
result_S_QQ_pha = S_QQ_pha-Phase_ini_XQ;

result_S_II_pha_unwrap = unwrap_phase(result_S_II_pha);
result_S_QI_pha_unwrap = unwrap_phase(result_S_QI_pha);
result_S_IQ_pha_unwrap = unwrap_phase(result_S_IQ_pha);
result_S_QQ_pha_unwrap = unwrap_phase(result_S_QQ_pha);
% phase after unwrap to calculate skew
result_S_II_pha_urp = result_S_II_pha_unwrap - result_S_II_pha_unwrap(1);
result_S_QI_pha_urp = result_S_QI_pha_unwrap - result_S_QI_pha_unwrap(1);
result_S_IQ_pha_urp = result_S_IQ_pha_unwrap - result_S_IQ_pha_unwrap(1);
result_S_QQ_pha_urp = result_S_QQ_pha_unwrap - result_S_QQ_pha_unwrap(1);
% skew calculate
S_II_QI=result_S_II_pha_urp - result_S_QI_pha_urp;
S_IQ_QQ=result_S_IQ_pha_urp - result_S_QQ_pha_urp;

[curve_fit,~]= polyfit_s(S_II_freq,S_II_QI, 1); 
delay_TX1 = curve_fit(1)*1e12./(2*pi);    
[curve_fit,~]= polyfit_s(S_IQ_freq,S_IQ_QQ, 1); 
delay_TX2 = curve_fit(1)*1e12./(2*pi);    
delay_TX=(delay_TX1+delay_TX2)/2;

S_II_IQ=result_S_II_pha_urp-result_S_IQ_pha_urp;
S_QI_QQ=result_S_QI_pha_urp-result_S_QQ_pha_urp;


[curve_fit,~]= polyfit_s(S_II_freq,S_II_IQ, 1); 
delay_RX1 = curve_fit(1)*1e12./(2*pi)*1;     
[curve_fit,~]= polyfit_s(S_IQ_freq,S_QI_QQ, 1); 
delay_RX2 = curve_fit(1)*1e12./(2*pi)*1;     
delay_RX=(delay_RX1+delay_RX2)/2;

disp(['Tx skew = ',num2str(delay_TX),'Rx skew = ',num2str(delay_RX)])

return