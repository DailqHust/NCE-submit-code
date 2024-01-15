% Copyright (c) 2023, Longquan Dai, Di Li, Lei Deng in Wuhan National Laboratory 
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
% [1] https://www.nature.com/articles/s44172-023-00147-3
%%
clear
close all
% multi-tone parameter
Fs = 120e9; %DAC sample rate
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
% random phase value
rng(1)  
Phase_ini_XI1 = 1*rand(1,Ntone)*2*pi-pi;
rng(2)
Phase_ini_XQ1 = 1*rand(1,Ntone)*2*pi-pi;
rng(3)
Phase_ini_YI1 = 1*rand(1,Ntone)*2*pi-pi;
rng(4)
Phase_ini_YQ1 = 1*rand(1,Ntone)*2*pi-pi;
% multi-tone signal generation
Tx_I = 0;
Tx_Q = 0;
Ty_I = 0;
Ty_Q = 0;
% frequency of XI/XQ/YI/YQ
delta_f0 = 0*(f00+df)/4;
delta_f1 = 1*(f00+df)/4;
delta_f2 = 2*(f00+df)/4;
delta_f3 = 3*(f00+df)/4;
t = (0:N_Sample - 1)/Fs;
for i = 1:Ntone
    Tx_I = Tx_I + cos(2*pi*(f0(i))*t + Phase_ini_XI1(i));   % multi-tone signal of XI 
    Tx_Q = Tx_Q + cos(2*pi*(f0(i)+delta_f2(i))*t + Phase_ini_XQ1(i));   % multi-tone signal of XQ
    Ty_I = Ty_I + cos(2*pi*(f0(i)+delta_f1(i))*t + Phase_ini_YI1(i));   % multi-tone signal of YI
    Ty_Q = Ty_Q + cos(2*pi*(f0(i)+delta_f3(i))*t + Phase_ini_YQ1(i));   % multi-tone signal of YQ
end
%


