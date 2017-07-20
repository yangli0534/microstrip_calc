function Zc = microstrip_z_calc(w,h,l,t,fc,er)

U = w/h; % ratio of trace width to substrate thickness
if t > 0
    T = t/h; %ratio of conductor thickness to substrate thickness
    %(T/PI)*log(1.0 + 4.0*exp(1.0)/(T*power(coth(sqrt(6.517*u)),2.0)))
    U1 = U +(T*log(1.0+4.0*exp(1)/T/power(1.0/tanh(sqrt(6.517*U)),2.0)))/pi; % from Hammerstad and Jensen
    %   0.5*(1.0 + 1.0/cosh(sqrt(er-1.0)))*deltau1
    Ur = U +(U1-U)*(1.0+1.0/(cosh(sqrt(er-1))))/2.0; % from Hammerstad and Jensen
else
    U1 = U;
    Ur = U;
end
Y = ee_HandJ(Ur,er);
Z0 =z0_HandJ(Ur)/sqrt(Y);
%ereff0 = Y*power(Z01_U1/Z01_Ur,2)
ereff0 = Y*power(z0_HandJ(U1)/z0_HandJ(Ur),2.0);
fn = fc*h/1e7;
P1 = 0.27488 + (0.6315 + (0.525 / (power((1 + 0.157*fn),20))) )*U - 0.065683*exp(-8.7513*U);
P2 = 0.33622*(1 - exp(-0.03442*er));
P3 = 0.0363*exp(-4.6*U)*(1 - exp(-power((fn / 3.87),4.97)));
P4 = 1 + 2.751*( 1 -  exp(-power((er/15.916),8)));
P = P1*P2*power(((0.1844 + P3*P4)*10*fn),1.5763);
ereff = (er*P+ereff0)/(1+P); % equavlent ralative dielectric constant
%ms_Zc = Z0*power(ereff0/ms_ereff,0.5)*(ms_ereff-1)/(ereff0-1)
fn = fc*h/1e6;
R1 = 0.03891*(power(er,1.4));
R2 = 0.267*(power(U,7.0));
R3 = 4.766*exp(-3.228*(power(U,0.641)));
R4 = 0.016 + power((0.0514*er),4.524);
R5 = power((fn/28.843),12.0);
R6 = 22.20*(power(U,1.92));
R7 = 1.206 - 0.3144*exp(-R1)*(1 - exp(-R2));
R8 = 1.0 + 1.275*(1.0 -  exp(-0.004625*R3*power(er,1.674)*power(fn/18.365,2.745)));
R9 = (5.086*R4*R5/(0.3838 + 0.386*R4))*(exp(-R6)/(1 + 1.2992*R5));
R9 = R9 * (power((er-1),6))/(1 + 10*power((er-1),6));
R10 = 0.00044*(power(er,2.136)) + 0.0184;
R11 = (power((fn/19.47),6))/(1 + 0.0962*(power((fn/19.47),6)));
R12 = 1 / (1 + 0.00245*U*U);
R13 = 0.9408*(power( ereff,R8)) - 0.9603;
R14 = (0.9408 - R9)*(power(ereff0,R8))-0.9603;
R15 = 0.707*R10*(power((fn/12.3),1.097));
R16 = 1 + 0.0503*er*er*R11*(1 - exp(-(power((U/15),6))));
R17 = R7*(1 - 1.1241*(R12/R16)*exp(-0.026*(power(fn,1.15656))-R15));
Zc = Z0*(power((R13/R14),R17));% characteristic impedance
end