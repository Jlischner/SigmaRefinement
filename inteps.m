function out = inteps(t,u,f1,f2,f3,f4,ws,indx);

Nfreq = size(ws)(1);

f = zeros(Nfreq,4);
f = [abs(f1) abs(f2) abs(f3) abs(f4)];

%# find locations of maxima
xmax = zeros(4,2);
ymax = zeros(4,2);

for n = 1:4;
  
  [a,b] = max(f(1:indx,n ) );
  xmax(n,1) = ws(b);
  ymax(n,1) = a;
  [a,b] = max(f(indx:Nfreq,n ) );
  xmax(n,2) = ws(indx:Nfreq)(b);
  ymax(n,2) = a;
endfor;

%# find center of function corresponding to (t,u)
x1 = (1-t).*(1-u)*xmax(1,1)  + t.*(1-u)*xmax(2,1)  + t.*u*xmax(3,1)  + (1-t).*u*xmax(4,1);
x2 = (1-t).*(1-u)*xmax(1,2)  + t.*(1-u)*xmax(2,2)  + t.*u*xmax(3,2)  + (1-t).*u*xmax(4,2);
y1 = (1-t).*(1-u)*ymax(1,1)  + t.*(1-u)*ymax(2,1)  + t.*u*ymax(3,1)  + (1-t).*u*ymax(4,1);
y2 = (1-t).*(1-u)*ymax(1,2)  + t.*(1-u)*ymax(2,2)  + t.*u*ymax(3,2)  + (1-t).*u*ymax(4,2);

%# get interpolations from each corner
for n = 1:4;

  d0 = [x1-xmax(n,1) y1-ymax(n,1)];
  d1 = [x2-xmax(n,2) y2-ymax(n,2)];
  P  = [ws f(:,n)];
  A0 = d0 - (d0-d1)/(xmax(n,1)-xmax(n,2))*xmax(n,1);
  B0 = (d0-d1)/(xmax(n,1)-xmax(n,2));
  dis= ones(size(ws))*A0 + ws*B0;
  Ppt= P + dis;
  Pp = interp1(Ppt(:,1),Ppt(:,2),ws);
  Pp(isna(Pp)) = 0;
  Pp(isnan(Pp))= 0;
  fp(n,:) = Pp;
endfor;

fi = (1-t).*(1-u)*fp(1,:)  + t.*(1-u)*fp(2,:)  + t.*u*fp(3,:)  + (1-t).*u*fp(4,:);
out = -fi;
endfunction;