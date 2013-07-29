%# cutoff frequency:
%# should be located between the carrier plasmon peak and the onset of 
%# interband transition where epsinv is approximately zero
wcut = 1.5;
%# head of epsilon for smallest qpoint:
epsinv = load("EpsInvDyn");
%#*******************************************************************

Nfreq = size(epsinv)(1);

%# find carrier plasmon peak position:
[a,b] = min(epsinv(:,3));
w_eV = epsinv(b,1);

%# find weight in carrier plasmon peak by integrating up to wcut:
[a,b] = min(abs(wcut - epsinv(:,1)));
int_eV = trapz( epsinv(1:b,1), epsinv(1:b,3));

%# get ImEpsInv without carrier plasmon peak;
mat = zeros(Nfreq,2);
mat(:,1) = epsinv(:,1);

[a,b] = min(abs(wcut - epsinv(:,1)));
mat(b+1:Nfreq,2) = epsinv(b+1:Nfreq,3);

plot(epsinv(:,1),epsinv(:,3),'sr-',mat(:,1),mat(:,2),'sb-');axis([0 4 -2 0]);
printf("epsinv(q,w) = aq * delta(w-wq) \n");
printf("wq = %f eV \n",w_eV);
printf("aq = %f eV \n",int_eV);

save dIepsI_dat mat w_eV int_eV;