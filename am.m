function a=am(v)       % Alpha for activation gating variable m in Na current
     E0m = v + 41;
     deltam = 1e-5; %mV
     if (abs(E0m) < deltam)
        a = 2000;
     else
        a = (200*E0m)/(1-exp(-0.1*E0m));
     end
end


