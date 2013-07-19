function [BZks,neq,indx_list,partner,partner_ks] = BZk(syms,kfull,k,tol);

  %# first find little group of k
  Nk    = size(kfull)(1);
  Nsyms = size(syms)(1);
  partner = zeros(Nk,1);
  partner_ks = zeros(Nk,1);

  ntranq = 0;
  indq   = [];
  for jj = 1:Nsyms;
    
    S  = reshape(syms(jj,:),3,3);
    dk = transpose( S*k' ) - k;
    dk(dk < -tol)  += 1;
    dk(dk > 1-tol) -= 1;
  
    if( abs(dk) < tol );
      ntranq += 1;
      indq = [indq; jj];
    endif;

  endfor;
  
  %# reduce full BZ using little group of k
  %# assume first kpoint is [0 0 0]
  nrq = 1;
  BZks = [0 0 0];
  neq = 1;
  indx_list = [1];
  partner(1) = 1;
  partner_ks(1) = 1;
  for ii = 2:Nk;
    
    q = kfull(ii,:);
    has_partner = 0;
    for it = 2:ntranq;
      
      indS = indq(it);
      S  = reshape(syms(indS,:),3,3);
      Sq = transpose( S*q' );
      Sq(Sq < -tol)  += 1;
      Sq(Sq > 1-tol) -= 1;
      
      [diff,indx] = min(sqrt(sum( (BZks - ones(nrq,1)*Sq).^2,2)));
      if( diff < tol);
	has_partner = 1;
	neq(indx)  += 1;
	partner(ii) = indx;
	partner_ks(ii) = indx_list(indx); 
	break;
      endif;
      
    endfor;
    
    if( has_partner == 0);
      BZks = [BZks; q];
      indx_list = [indx_list; ii];
      nrq += 1;
      neq_save = neq;
      neq = zeros(nrq,1);
      neq(1:nrq-1) = neq_save;
      neq(nrq) = 1;
      partner_ks(ii) = ii;
      partner(ii) = nrq;
    endif;

  endfor;
    
endfunction;