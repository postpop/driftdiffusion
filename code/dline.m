function hhh=dline(y,in1,in2)
% lineHandle = dline(y, linetype, label)
%
% Draws a diagonal line on the current axes at the location specified by 'y'. 

% modified from vline/hline by
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if nargin==0
   lim = [get(gca,'YLim') get(gca,'XLim')];
   y = [min(lim), max(lim)];
   linetype='-k';
   label='';
end

if length(y)>2  % vector input
   for I=1:length(y)
      switch nargin
         case 1
            linetype='-k';
            label='';
         case 2
            if ~iscell(in1)
               in1={in1};
            end
            if I>length(in1)
               linetype=in1{end};
            else
               linetype=in1{I};
            end
            label='';
         case 3
            if ~iscell(in1)
               in1={in1};
            end
            if ~iscell(in2)
               in2={in2};
            end
            if I>length(in1)
               linetype=in1{end};
            else
               linetype=in1{I};
            end
            if I>length(in2)
               label=in2{end};
            else
               label=in2{I};
            end
      end
      h(I)=hline(y(I),linetype,label);
   end
else
   switch nargin
      case 1
         linetype='-k';
         label='';
      case 2
         linetype=in1;
         label='';
      case 3
         linetype=in1;
         label=in2;
   end
   
   
   
   
   g=ishold(gca);
   hold on
   
   x=get(gca,'xlim');
   h=plot(y,y,linetype);
   if ~isempty(label)
      yy=get(gca,'ylim');
      yrange=yy(2)-yy(1);
      yunit=(y-yy(1))/yrange;
      if yunit<0.2
         text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
      else
         text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
      end
   end
   
   if g==0
      hold off
   end
   set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else

if nargout
   hhh=h;
end
