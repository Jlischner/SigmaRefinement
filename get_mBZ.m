df = load("EpsInvDyn2_fine");
dc = load("EpsInvDyn2");

Nfine = size(df)(1);
Ncoarse = size(dc)(1);

%# first stitch together d and d2:
%# fine sampling from 0 to 1 eV

[a,b] = min(abs(1-dc(:,1)));

Nstitch = Nfine + Ncoarse-b;
d = zeros(Nstitch,3);
d(1:Nfine,:) = df;
d(Nfine+1:Nstitch,:) = dc(b+1:Ncoarse,:);

%# integrate IepsI up to 1.5 eV:
[a,b] = min(d(:,3));
w_eV = d(b,1);

[a,b] = min(abs(1.5-d(:,1)));
int_eV = trapz(d(1:b,1),d(1:b,3));

%# write out dIepsI;
[a,b] = min(abs(1.5-dc(:,1)));
mat = zeros(Ncoarse,2);
mat(:,1) = dc(:,1);
mat(b+1:Ncoarse,2) = dc(b+1:Ncoarse,3);

save dIepsI_dat mat w_eV int_eV;
