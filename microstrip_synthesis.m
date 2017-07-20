clc;
clear all;
close all;


%%
% constant
varepsilon_0  = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c = 1/sqrt(varepsilon_0*mu0);

%%
%electric parameters
Zc = 50 ;% characteristic impedance.
El = 90;%
fc = 35e9; %center frequency :Hz
%%
% substate parameters
er = 2.2 ; % relative permittivity constant
h = 0.508e-3; % substrate height:m
t = 0;%metal thickness :m
mur = 1; % ralative permeability constant
rough = 0 ; %
cond = 5.88e7; % conductivity
%substrate = struct('er',2.2,'h',0.508,'t',0.1,'tand',0.0009,'cond',5.8e7,'rought',0)
mu = mur*4*pi*1.0e-7;%   permeability
%t = float(lineEdit_ms_T.text())/1.0e3

lambda0 = c/fc;% % wavelengt in free space : m
lx = 1000*25.4e-6;
wmin = 0.01*25.4e-6;
wmax = 499.0*25.4e-6;
%impedance convergence tolerance (ohms)
abstol = 1e-6;
reltol = 0.1e-6;
maxiters = 50;
A = ((er - 1)/(er + 1)) * (0.226 + 0.121/er) + (pi/377.0)*sqrt(2*(er+1))*Zc;
w_h = 4/(0.5*exp(A) - exp(-A));
if w_h > 2.0
    B = pi*377.0/(2*Zc*sqrt(er));
    w_h = (2/pi)*(B - 1 - log(2*B - 1) + ((er-1)/(2*er))*(log(B-1) + 0.293 - 0.517/er));
end
wx = h * w_h;

if wx >= wmax
    wx = 0.95*wmax;
end

if wx <= wmin
    wx = wmin;
end
wold = 1.01*wx;
Zold =  microstrip_z_calc(wold,h,lx,t,fc,er);


if Zold < Zc
    wmax = wold;
else
    wmin = wold;
end

iters = 0;
done = 0;

while done == 0
    iters = iters + 1;
    Z0 = microstrip_z_calc(wx,h,lx,t,fc,er);
    if Z0 < Zc
        wmax = wx;
    else
        wmin = wx;
    end
    if abs(Z0-Zc) < abstol
        done = 1;
    elseif abs(wx-wold) < reltol
        done = 1;
    elseif iters >= maxiters
        done = 1;
    else
        dzdw = (Z0 -Zold)/(wx-wold);
        wold = wx;
        Zold = Z0;
        wx = wx -(Z0-Zc)/dzdw;
        if (wx > wmax) || (wx < wmin)
            wx = (wmin+ wmax)/2.0;
        end
    end
end

er_eff = er_eff_calc(wx,h,t,fc,er);
v =  c/sqrt(er_eff);
l = El/360.0*v/fc;
num2str(wx*1.0e3)
num2str(er_eff)
num2str(l*1.0e3)