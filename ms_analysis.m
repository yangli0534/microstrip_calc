function ms = ms_analysis(substrate, w, l,fc)
% MS_ANALYSIS calc electric specification of microstrip line 
% Author : LEON, yangli0534@yahoo.com
% blog : www.cnblogs.com/hiramlee0534
% parameters : 
% substrate is a struct ,('er',varepsilon_r,'h',h,'mur',mu_r, 't',t, 'cond', cond, 'rough',0,'tand',tanD);
% w , width of ms
% l , length of ms
% fc, frequency
% ms was a struct, ms = struct('Zc',ms_Zc, 'er_eff',ms_ereff, 'El',ms_el,'beta',ms_beta, 'loss',loss, 'delay',delay, 'delta',skin_depth);
% Zc: charateristic impdance
% er_eff: effective relative permitivity
% El : electric length
% 
%%
% constant
varepsilon_0  = 8.854187817e-12;
mu0 = 4*pi*1e-7;
c = 1/sqrt(varepsilon_0*mu0);

%%
%electric parameters
%Zc = 50 ;% characteristic impedance.
%El = 90;%
%fc = 35e9; %center frequency :Hz
%%
% substate parameters
ms_er = substrate.er;%2.2 ; % relative permittivity constant
ms_h = substrate.h/1e3;%0.508e-3; % substrate height:m
ms_t = substrate.t/1e3;%;0;%metal thickness :m
mur = substrate.mur;%1; % ralative permeability constant
rough = substrate.rough;%T0 ; %
cond = substrate.cond;%5.88e7; % conductivity
ms_tand = substrate.tand;
%substrate = struct('er',2.2,'h',0.508,'t',0.1,'tand',0.0009,'cond',5.8e7,'rought',0)
mu = mur*4*pi*1.0e-7;%   permeability
%t = float(lineEdit_ms_T.text())/1.0e3
%%
%freq
lambda0 = c/fc;% % wavelengt in free space : m
ms_fc = fc;
ms_w = w/1e3;
ms_l = l/1e3;
%%
%calc
U = ms_w/ms_h; % ratio of trace width to substrate thickness
if ms_t > 0
    T = ms_t/ms_h; %ratio of conductor thickness to substrate thickness
    %(T/PI)*log(1.0 + 4.0*exp(1.0)/(T*power(coth(sqrt(6.517*u)),2.0)))
    U1 = U +(T*log(1.0+4.0*exp(1)/T/power(1.0/tanh(sqrt(6.517*U)),2.0)))/pi; % from Hammerstad and Jensen
    %   0.5*(1.0 + 1.0/cosh(sqrt(er-1.0)))*deltau1
    Ur = U +(U1-U)*(1.0+1.0/(cosh(sqrt(ms_er-1))))/2.0; % from Hammerstad and Jensen
else
    U1 = U;
    Ur = U;
end
Y = ee_HandJ(Ur,ms_er);
%        %Z0 =z0_HandJ(Ur)/sqrt(Y)
%ereff0 = Y*power(Z01_U1/Z01_Ur,2)
ereff0 = Y*power(z0_HandJ(U1)/z0_HandJ(Ur),2.0);

fn = ms_fc*ms_h/1e7;
P1 = 0.27488 + (0.6315 + (0.525 / (power((1 + 0.157*fn),20))) )*U - 0.065683*exp(-8.7513*U);
P2 = 0.33622*(1 - exp(-0.03442*ms_er));
P3 = 0.0363*exp(-4.6*U)*(1 - exp(-power((fn / 3.87),4.97)));
P4 = 1 + 2.751*( 1 -  exp(-power((ms_er/15.916),8)));
P = P1*P2*power(((0.1844 + P3*P4)*10*fn),1.5763);
ms_ereff = (ms_er*P+ereff0)/(1+P); % equavlent ralative dielectric constant
ms_Zc = microstrip_z_calc(ms_w,ms_h,ms_l,ms_t,ms_fc,ms_er);
z1 = 0.434907*((power(ms_ereff,0.81) + 0.26)/(power(ms_ereff,0.81) - 0.189))*(power(U,0.8544) + 0.236)/(power(U,0.8544) + 0.87);
z2 = 1 + (power(U,0.371))/(2.358*ms_er + 1);
z3 = 1 + (0.5274*atan(0.084*(power(U,(1.9413 / z2)))))/(power(ms_ereff,0.9236));
z4 = 1 + 0.0377*atan(0.067*(power(U,1.456)))*(6 - 5*exp(0.036*(1-ms_er)));
z5 = 1 - 0.218*exp(-7.5*U);
deltal = ms_h * z1 * z3 * z5 / z4;
v = c/sqrt(ms_ereff);
length = ms_l/(v/ms_fc);
% delay
delay = 1e9*ms_l/v;
% dielectric losses
ld = pi*ms_fc/v*ms_er/ms_ereff*(ms_ereff-1)/(ms_er-1)*ms_tand ;% unit in nepers/meter
ld = 20.0/log(10)*ld ;% unit in dB/meter
ld = ld * ms_l; % unit in dB
% conduction losses
% skin depth in meters
skin_depth = sqrt(1/(pi*ms_fc*mu*cond));
if skin_depth <= ms_t
    Z2 = microstrip_z_calc(ms_w,ms_h,ms_l,ms_t,ms_fc,1);
    Z1 = microstrip_z_calc(ms_w-skin_depth,ms_h+skin_depth,ms_l,ms_t-skin_depth,ms_fc,1);
    % conduction losser
    lc = pi*ms_fc/c*(Z1-Z2)/ms_Zc;
elseif ms_t > 0
    %resistance per meter = 1/(Area*conductivity)
    res = 1/(ms_w*ms_t*cond);
    %conduction losses, nepers per meter
    lc = res/(2.0*ms_Zc);
    skin_depth = ms_t;
else
    lc = 0;
end
lc = 20.0/log(10)*lc;
lc = lc * ms_l;
lc = lc * (1.0+2.0/pi*atan(1.4*power(rough/skin_depth,2)));
loss = ld + lc;
ms_k0 = 2*pi/c*ms_fc;
ms_beta = ms_k0*sqrt(ms_ereff);
ms_el = sqrt(ms_ereff)*ms_k0*ms_l/1000/pi*180.0;
ms_el = 360*length;
%         lineEdit_ms_Zc.setText(str(ms_Zc))
%         lineEdit_ms_ereff.setText(str(ms_ereff))
%         lineEdit_ms_el.setText(str(ms_el))
%         lineEdit_ms_beta.setText(str(ms_beta))
%         lineEdit_ms_loss.setText(str(loss))
%         label_delay.setText('delay = '+str(delay)+' ns')
%         label_skin_depth.setText('skin depth = '+str(skin_depth)+' m')
ms = struct('Zc',ms_Zc, 'er_eff',ms_ereff, 'El',ms_el,'beta',ms_beta, 'loss',loss, 'delay',delay, 'delta',skin_depth);
end