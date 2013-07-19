function [wse,IepsI] = getIepsI(qref_list,k_irr,ref_dir,unref_dir);

  global Nkx;
  global unref_flag;

  Nkref = size(qref_list)(1); %# number of refined epsilon-heads
  if(unref_flag == 1);
    Nkref = 0;
  endif;

  wse = load( strcat(ref_dir,int2str(2)))(:,1);
  Nfreqe = size(wse)(1);
  Nkirr = size(k_irr)(1);
  IepsI = zeros(Nfreqe,Nkirr);
  
  %# load unrefined heads:                 
  for n = 1:Nkirr-1;
    
    if( Nkx == 72);
      epsdat = load( strcat(unref_dir,int2str(n)));
      inteps = interp1(epsdat(:,1),epsdat(:,3),wse);
      inteps(isna(inteps))  = 0;
      inteps(isnan(inteps)) = 0;
      IepsI(:,n+1) = inteps;
    else;
      epsdat = load( strcat(unref_dir,int2str(n)));
    endif;
    
  endfor;
  
  %# load refined heads:      
  for n = 1:Nkref;
    
    if(Nkx == 72);
      epsdat = load(strcat(ref_dir,int2str(qref_list(n))));
      IepsI(:, qref_list(n))   = epsdat(:,3);
    else;
      epsdat = load(strcat(ref_dir,int2str(qref_list(n))));
      IepsI(:, qref_list(n)+1) = epsdat(:,3);
    endif;
  endfor
  
endfunction;