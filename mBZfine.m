function avg_t = mBZfine(indx1,ws,band,wse,eta,dIepsI,w_eV,int_eV,EF);

  global G;
  global Gmat;
  global d;
  global tlc;
  global lenS;
  global t;
  global u;
  global zc;
  global fact;
  global Nfreq;
  global Nk;
  global partner_irr;
  global energies;
  global k;
  global ks;
  global ks_shifted;
  global partner_shifted;
  global Nkx;
  global BZkf;
  global NBZf;
  global neqf;
  global tlcf;

  rq72 = sqrt(G/Nkx); 
  alpha = w_eV/27.21/rq72; %# a.u.
  C = int_eV/27.21/rq72; 

  avg_cs = zeros(size(ws));
  avg_c  = zeros(size(ws));

  q1 = ks(indx1,:);

  %# find four corners of coarse grid square:
  q2 = bringBZ( q1 + [d 0 0] );
  q3 = bringBZ( q1 + [d d 0] );
  q4 = bringBZ( q1 + [0 d 0] );
  
  %# find kq = k-q
  kq1 = bringBZ( k - q1 );
  kq2 = bringBZ( k - q2 );
  kq3 = bringBZ( k - q3 );
  kq4 = bringBZ( k - q4 );

  [dis1,ind1] = min( sqrt( sum( (ks_shifted-ones(Nk,1)*kq1).^2,2) ) ); 
  if( dis1(1) > tlc);
    printf("could not find kq1=%f %f in list of kpoints! \n",kq1(1),kq1(2));
  else;
    indx1q = ind1(1);
    e1q = energies(partner_shifted(indx1q),band);
  endif;
  
  [dis2,ind2] = min( sqrt( sum( (ks_shifted-ones(Nk,1)*kq2).^2,2) ) ); 
  if( dis2(1) > tlc);
    printf("could not find kq2 in list of kpoints! \n");
  else;
    indx2q = ind2(1);
    e2q = energies(partner_shifted(indx2q),band);
  endif;
  
  [dis3,ind3] = min( sqrt( sum( (ks_shifted-ones(Nk,1)*kq3).^2,2) ) );
  if( dis3(1) > tlc);
    printf("could not find kq3 in list of kpoints! \n");
  else;
    indx3q = ind3(1);
    e3q = energies(partner_shifted(indx3q),band);
  endif;
  
  [dis4,ind4] = min( sqrt( sum( (ks_shifted-ones(Nk,1)*kq4).^2,2) ) );
  if( dis4(1) > tlc);
    printf("could not find kq4 in list of kpoints! \n");
  else;
    indx4q = ind4(1);
    e4q = energies(partner_shifted(indx4q),band);
  endif;


  #% interpolate energies
  Ekqs= (1-t).*(1-u)*e1q + t.*(1-u)*e2q + t.*u*e3q + (1-t).*u*e4q;
  #% Coulomb potential on fine grid
  qclose = q1;
  qclose(qclose>0.5) -= 1;
  qs_check = ones(lenS,1)*q1 + d*[t u zeros(lenS,1)];
  qs = ones(lenS,1)*qclose + d*[t u zeros(lenS,1)];
  aqs = sqrt( sum( (qs*Gmat).^2,2));
  vcoul = 8*pi./aqs.^2 .* (1-exp(-aqs*zc))*fact; 
  ws  /= 27.21;
  wse /= 27.21;
  eta /= 27.21;
  EF  /= 27.21;

  for jj = 1:lenS;

    qjj = qs_check(jj,:);
    [djj,ijj] = min( sqrt( sum( (BZkf-ones(NBZf,1)*qjj).^2,2) ) );
    if( djj < tlcf );
      
      %# avoid singularity at q=0
      if(jj == 1);
	jj = 2;
      endif;

      tc = t(jj);
      uc = u(jj);
      vc = vcoul(jj);
      ekq= Ekqs(jj)/27.21;
      csign = sign(ekq-EF);
      rq = sqrt(aqs(jj));
      neq = neqf(ijj);
      
      for iw = 1:Nfreq;
	avg_cs(iw)  += neq*vc*trapz(wse,dIepsI./(ws(iw)-ekq - csign*wse + I*eta));
      endfor;
      %# treat delta-function analytically
      avg_c += neq*vc*C*rq./(ws-ekq-csign*alpha*rq + I*eta);
    endif;
  endfor;
  avg_cs(:) /= -lenS*pi;
  avg_c(:)  /= -lenS*pi;
  avg_t = avg_c + avg_cs;
  
endfunction;