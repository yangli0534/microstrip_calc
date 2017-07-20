function er_eff =er_eff_calc(w,h,t,fc,er)
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
er_eff = (er*P+ereff0)/(1+P)% equavlent ralative dielectric constant
end
