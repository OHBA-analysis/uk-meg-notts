function crosscorr_snugplot(a,b,c,pc,pcedge)
% snugplot(a,b,c,pc)
% like subplot(a,b,c)
% pc is the fraction of white space along each axis.
% pcedge is the extra fraction of white space at the edges.
% TB04

if(c>(a*b))
  error('ya twat');
end
if(nargin==3)
  pc=0.25;
  pcedge=0;
end

if(nargin==4)
  pcedge=0;
end

percy=(1-pc)/a;
percx=(1-pc)/b;

ygap=pc/(a+1);
xgap=pc/(b+1);

row=ceil(c/b);
row=a-row+1;
col=c-(ceil(c/b)-1)*b;

subplot('position',[(pcedge+xgap+(xgap+percx)*(col-1))/(1+pcedge) (pcedge+ygap+(ygap+percy)*(row-1))/(1+pcedge) percx/(1+pcedge) percy/(1+pcedge)]);

set(gcf,'PaperPositionMode','auto');

