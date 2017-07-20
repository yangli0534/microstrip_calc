
function z01 = z0_HandJ(u)
        LIGHTSPEED = 2.99792458e8;
        FREESPACEZ0 = 4.0*pi*1.0e-7*LIGHTSPEED;
        % from Hammerstad and Jensen.  'u' is the normalized width
        F = 6.0 + (2.0*pi - 6.0)*exp(-1*power((30.666/u),0.7528));
        % from Hammerstad and Jensen
        z01 = (FREESPACEZ0/(2*pi))*log(F/u + sqrt(1.0 + power((2/u),2.0)));
end