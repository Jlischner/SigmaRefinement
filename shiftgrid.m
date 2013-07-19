dfull = load("Derek72grid/kshifted_full");
dirr = load("Derek72grid/kshifted_irr");

ks_full = dfull(:,2:4);
ks_irr  = dirr(:,2:4);
%# map from full zone to irr zone
%# note: irr zone kpoints are labeled by their position in full zone
full2irr = dfull(:,6); 

%# map from irr zone to full
%# gives location of irr zone kpts in full zone
irr2full = dirr(:,6);

%# tolerance
tlc = 1/72/4;

Nirr = size(ks_irr)(1);
Nfull = size(ks_full)(1);
partner = zeros(Nfull,1);

for ii = 1:Nfull;
  
  if( full2irr(ii) < tlc);
    indx_full = ii;
  else;
    indx_full = full2irr(ii);
  endif;

  [dis,indx] = min(abs( irr2full-indx_full ));
  partner(ii) = indx;

endfor;

fid = fopen("partner_shifted","w");
fprintf(fid,"%d \n",partner');
fclose(fid);
