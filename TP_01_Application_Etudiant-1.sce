clear;
exec("TP_01_Outils_Quadrature.sce", -1);
function y = modele(x)
    y = exp(x);
endfunction
Iref = %e-1;
ns = [4 8 16 32];
ET = zeros(1, length(ns)); ES = ET;
errT = ET; errS = ET;
for k = 1:length(ns)
    n = ns(k);
    [t, y] = echantillonner(modele, 0, 1, n);
    ET(k) = trapezes(t, y);
    ES(k) = simpson(t, y);
    errT(k) = abs(ET(k) - Iref);
    errS(k) = abs(ES(k) - Iref);
    mprintf("n=%d T=%.9f S=%.9f ET=%.3e ES=%.3e\n", ..
            n, ET(k), ES(k), errT(k), errS(k));
end
for k = 1:length(ns)-1
    pT = ordre_observe(errT(k), errT(k+1));
    pS = ordre_observe(errS(k), errS(k+1));
    mprintf("n=%d -> %d : pT=%.3f pS=%.3f\n", ..
            ns(k), ns(k+1), pT, pS);
end
tmin = [0 10 20 30 40];
Pkw = [2.0 2.6 3.0 2.5 1.9];
th = zeros(1, length(tmin));
for i = 1:length(tmin)
    th(i) = tmin(i)/60;
end
Etrap = trapezes(th, Pkw);
Esimp = simpson(th, Pkw);
mprintf("Pompe : trapezes=%.6f kWh ; Simpson=%.6f kWh\n", Etrap, Esimp);
mprintf("Ecart=%.6f kWh ; borne des lectures=%.6f kWh\n", abs(Etrap-Esimp), 0.1*40/60);
tc = [0 1 2 3 4]; pc = 10 + 4*tc - tc.^2;
mprintf("Controle CM1 : T=%.6f S=%.6f 152/3=%.6f\n", trapezes(tc,pc), simpson(tc,pc), 152/3);
mprintf("Iref=%.9f\n", Iref);
