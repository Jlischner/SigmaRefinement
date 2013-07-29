more off;
addpath("/auto/jlischner/Dropbox/RefinementCode/Sigma/Version4");

%# input variables:
%#-----------------
input;
dIepsI = mat(:,2);
wsef = mat(:,1);
[wse,epsmat] = getIepsI(qref_list,kirr,ref_dir,unref_dir);
EF = EFmid + EFdirac;
qcut = kcut;

%# derived variables:
%#-------------------
Sigma1 = zeros(size(ws));
Sigma2 = zeros(size(ws));
band = 5; %# which band in sum
global tlc  = 1/Nkx/5;
global tlcf = 1/Nkx/S(1)/4;
[p_irr,kfull] = fullbz(syms,kirr);
global partner_irr = p_irr;
global ks = kfull;

if(shift_flag == 0);
  global ks_shifted = ks;
  global partner_shifted = p_irr;
endif;

[BZks,neq_ks,indx_list,partner_BZk,partner_ks] = BZk(syms,ks,k,tlc);
Lz = 2*pi/Gmat(3,3);%# z-length of unit cell in bohr
vF =  0.38859; %# corresponding to 0.85*10^8 cm/s, only needed to calc eta
vol = (2*pi)^3/abs(det(Gmat));
Nkirr = size(kirr)(1);
global Nk = size(ks)(1);
global G  = norm(Gmat(1,:));
global zc = Lz/2;
global fact = 1/(vol*Nk);
global lenS = prod(S);
global d = 1/Nkx;
global Nfreq = size(ws)(1);
eta = G/Nkx/S(1) * vF * 27.21 * 2;
if(unref_flag == 1)
  eta = etaBGW;
endif;
[wc2,anchor] = min(abs(wse-wc_eps));
Nfreqe = size(wse)(1);
Ntrans = Nk;
mBZflag = 0;

%# set up fine grid
ms=[0:prod(S)-1];
ms=ms.';
m1=rem(ms,S(1));
m2=rem(floor(ms/S(1)),S(2));
global t = m1/S(1);
global u = m2/S(2);

Sf = [Nkx*S(1) Nkx*S(2)];
lenSf = prod(Sf);
ms=[0:prod(Sf)-1];
ms=ms.';
m1=rem(ms,Sf(1));
m2=rem(floor(ms/Sf(1)),Sf(2));
tf = m1/Sf(1);
uf = m2/Sf(2);
ksf = [tf uf zeros(lenSf,1)];

printf("symmetry reducing fine grid...\n");
if( load_bzkf != 1);
  [BZkf_o,neqf_o,indx_list_o,partner1_o,partner2_o] = BZk(syms,ksf,k,tlcf);  
  save BZkf_output BZkf_o neqf_o indx_list_o partner1_o partner2_o;

  if(unref_flag == 1)
    BZkf_o = BZks;
    neqf_o = neq_ks;
  endif;

else;
  load BZkf_output;
endif;
printf("done!\n");

global BZkf = BZkf_o;
global neqf = neqf_o
global NBZf = size(BZkf)(1);

avg_c  = zeros(size(ws));
indxv = [];
avg_n = zeros(size(ws),Ntrans);

count = 1;
for indx1 = 1:Nk;

  #printf("doing %d out of %d \n",indx1,Nk);
  if(indx1 == 1);
    printf("mBZ treatment for indx1=%d \n",indx1);
    mBZindx1 = indx1;
    mBZflag  = 1;
    indx1 = 2;
  endif;
  
  q1 = ks(indx1,:);
  qirr = kirr(partner_irr(indx1),:); %# find partner in irr wedge

  if( norm(qirr*Gmat) < qcut );

    %# find four corners of coarse grid square:
    q2 = bringBZ( q1 + [d 0 0] );
    if( norm(q2*Gmat) < 0.01 && unref_flag == 0);
      printf("mBZ treatment for indx1=%d \n",indx1);
      mBZindx1 = indx1;
      mBZflag  = 1;
      indx1 = 2;
    endif;
    
    q3 = bringBZ( q1 + [d d 0] );
    if( norm(q3*Gmat) < 0.01 && unref_flag == 0);
      printf("mBZ treatment for indx1=%d \n",indx1);
      mBZindx1 = indx1;
      mBZflag  = 1;
      indx1 = 2;
    endif;

    q4 = bringBZ( q1 + [0 d 0] );
    if( norm(q4*Gmat) < 0.01 && unref_flag == 0);
      printf("mBZ treatment for indx1=%d \n",indx1);
      mBZindx1 = indx1;
      mBZflag  = 1;
      indx1 = 2;
    endif;
    
    eps1 = epsmat(:,partner_irr(indx1) );
 
    [dis2,ind2] = min( sqrt( sum( (ks-ones(Nk,1)*q2).^2,2) ) ); 
    if( dis2(1) > tlc);
      printf("could not find k2 in list of kpoints! \n");
    else;
      indx2 = ind2(1);
      eps2 = epsmat(:,partner_irr(indx2) );
    endif;
    
    [dis3,ind3] = min( sqrt( sum( (ks-ones(Nk,1)*q3).^2,2) ) );
    if( dis3(1) > tlc);
      printf("could not find k3 in list of kpoints! \n");
    else;
      indx3 = ind3(1);
      eps3 = epsmat(:,partner_irr(indx3) );
    endif;
    
    [dis4,ind4] = min( sqrt( sum( (ks-ones(Nk,1)*q4).^2,2) ) );
    if( dis4(1) > tlc);
      printf("could not find k4 in list of kpoints! \n");
    else;
      indx4 = ind4(1);
      eps4 = epsmat(:,partner_irr(indx4) );
    endif;
    
    kq1 = bringBZ( k - q1 );
    kq2 = bringBZ( k - q2 );
    kq3 = bringBZ( k - q3 );
    kq4 = bringBZ( k - q4 );

    if( mBZflag == 1 && unref_flag == 1);
      kq1 = k;
    endif;

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
    
    
    #% interpolate energies:
    Ekqs= (1-t).*(1-u)*e1q + t.*(1-u)*e2q + t.*u*e3q + (1-t).*u*e4q;

    #% calculate Coulomb potential on fine grid:
    qcoul = q1;
    qcoul(qcoul>0.5) -= 1;
 
    qs_check = ones(lenS,1)*q1 + d*[t u zeros(lenS,1)]; 
    qs = ones(lenS,1)*qcoul + d*[t u zeros(lenS,1)];
    aqs = sqrt( sum( (qs*Gmat).^2,2));
    vcoul = 8*pi./aqs.^2 .* (1-exp(-aqs*zc)) * fact;

    avg_c = zeros(size(ws));
    for jj = 1:lenS;

      qjj = qs_check(jj,:);
      %# check if qjj is in symmetry reduced fine grid (BZkf):
      [djj,ijj] = min( sqrt( sum( (BZkf-ones(NBZf,1)*qjj).^2,2) ) );
      if( djj < tlcf );
	
	tc = t(jj);
	uc = u(jj);
	vc = vcoul(jj);
	ekq= Ekqs(jj);
	csign = sign(ekq-EF);
	neq = neqf(ijj);
	if(mBZflag == 1 && unref_flag == 1); %# for trivial refinement
	  neq = 1;
	endif;
	IepsI = inteps(tc,uc,eps1,eps2,eps3,eps4,wse,anchor)';
	
	for iw = 1:Nfreq;
	  avg_c(iw)  += neq*vc*trapz(wse,IepsI./(ws(iw)-ekq - csign*wse + I*eta));
	endfor;
	
      endif;
    endfor;
    avg_c *= -1/pi/lenS;

    if(mBZflag == 1);
      if(unref_flag == 0);
	avg_c = mBZfine(mBZindx1,ws,band,wsef,eta,dIepsI,w_eV,int_eV,EF);
      endif;
      indx1 = mBZindx1;
      Sigma2 += avg_c;
      mBZflag = 0;
    else;
      Sigma1 += avg_c;
    endif;

    neqfac = neq_ks( partner_BZk(indx1) );

    if( indx1 == 1);
      avg_n(:,count) = avg_c/neqfac;
      indxv = [indxv; indx1];
      count += 1;
    else;
      [diff,indxBZk] = min( abs( indxv - partner_ks(indx1) ));
      if( diff < 0.1);
	avg_n(:,indxBZk) += avg_c/neqfac;
      else;
	avg_n(:,count) = avg_c/neqfac;
	indxv = [indxv; indx1];
	count += 1;
      endif;
    endif;

  endif; %# qcut-loop
endfor; %# indx1-loop

%# write out stuff for BGW
Ntrans = count-1;
avg_r  = reshape(avg_n(:,1:Ntrans),Ntrans*Nfreq,1);

avg_r2 = [real(avg_r) imag(avg_r)];
fid = fopen("avg_file",'w');
fprintf(fid,"%20.12f %20.12f \n",avg_r2');
fclose(fid);
 
%# first entry in indxv is number of transitions
%# other entries are their indices
indxv = [Ntrans; indxv];
fid = fopen("indx_file",'w');
fprintf(fid,"%d \n",indxv');
fclose(fid);

printf("init_frequency_eval   %f \n",ws(1)-EF);
printf("delta_frequency_eval  %f \n",ws(2)-ws(1));
printf("number_frequency_eval %d \n",length(ws));

more on;