function Y = ee_HandJ(u,er)
        %Au = 1.0 + (1.0/49.0)*log(power(Ur,4.0) + power(Ur/52.0,2.0)/(power(Ur,4.0)+0.432))+ (1.0/18.7)*log(1.0+ power(Ur/18.1,3.0))
        A = 1.0 + (1.0/49.0)*log((power(u,4.0) + power((u/52.0),2.0))/(power(u,4.0) + 0.432))...
            + (1.0/18.7)*log(1.0 + power((u/18.1),3.0));
        %Ber = 0.564*power((er-0.9)/(er+3),0.053)
        B = 0.564*power(((er-0.9)/(er+3.0)),0.053);
        %Y = (er+1.0)/2.0+(er-1.0)/2.0*powerer(1+10.0/Ur,-(A*B))
        Y= (er+1.0)/2.0 + ((er-1.0)/2.0)*power((1.0 + 10.0/u),(-A*B));
end