function fig = UnitPlotTuning_ball( unit, protocol, graphinfo, fig, AnalysisType )
% UnitPlotTuning plots the tuning curves for a Unit
%
% UnitPlotTuning( unit )
%
% UnitPlotTuning( unit, protocol )
%
% UnitPlotTuning( unit, protocol, graphinfo ) 
% DEFAULT: [] (prompts the user)
%
% UnitPlotTuning( unit, protocol, graphinfo, fig )
% DEFAULT: [] (makes a new figure)
%
% UnitPlotTuning( unit, protocol, graphinfo, fig, AnalysisType )
% AnalysisType can be 'Mean', 'First Harmonic', or 'Second Harmonic'
% DEFAULT: '' (prompts the user)
%
% part of Spikes
%
% 09/00 MC
% 2001-10-10 VM use protocol.estfreqs if possible
% 2003-03 VM made it fit for traces
% 2007-04 MC added argument AnalysisType
% 2010-06 MC cosmetic changes
% 2010-06 MC minor changes so it can work with traces rather than spikes

if nargin < 5
    AnalysisType = '';
end

if nargin < 4
   fig = [];
end

if nargin < 3
    graphinfo = [];
end

if nargin < 2
    protocol = ProtocolLoad(unit);
end

%% 

if isempty(protocol)
    return;
end

if isempty(graphinfo)
   graphinfo = ProtocolGetGraphInfo(protocol);
   if isempty(graphinfo), return; end
end

if isempty(protocol.pfilefreqs)
   AnalysisType = 'Mean';
end

if isempty(fig)
    fig = FigTuning;
end

if isempty(AnalysisType)
   AnalysisType=questdlg('What response measure do you want to plot?', ...
      'Spike Sorter', ...
      'Mean','First Harmonic', 'Second Harmonic', 'Mean');
end

% errflag = 'std';
errflag = 'sem'; % MC reverted back to SEM 2010-03-08

switch AnalysisType
case 'Mean'
   [yy, ee] = UnitGetDC( unit, errflag );
case 'First Harmonic'
   if isempty(protocol.estfreqs)
      [yy, ee] = UnitGetHarm( 1, unit, protocol, errflag, 'broad' );
   else
      [yy, ee] = UnitGetHarm( 1, unit, protocol, errflag);
   end
   yy = abs(yy);
case 'Second Harmonic'
   if isempty(protocol.estfreqs)
      [yy, ee] = UnitGetHarm( 2, unit, protocol, errflag, 'broad' );
   else
      [yy, ee] = UnitGetHarm( 2, unit, protocol, errflag);
   end
   yy = abs(yy);
otherwise
   error('Huh?');
end

%% make the figure

figure(fig);

nrows = size(graphinfo,1);
ncols = size(graphinfo,2);
ax = zeros(nrows,ncols);

for irow = 1:nrows
   for icol = 1:ncols
      
      xx 	= graphinfo(irow,icol).xx;
      logflag = islogspaced(xx);
      uxx = unique(xx);
      if logflag && uxx(1)==0, xx(xx==0)=uxx(2)/2; end
      xticks = unique(xx);
      xticks = xticks(1:2:end);
      xticklabels = xticks;
      if logflag && uxx(1)==0, xticklabels(1)=0; end
      
      pos 	= graphinfo(irow,icol).pos;

      ax(irow,icol) = gridplot(nrows,ncols,irow,icol);
      
      % this used to be plottune( xx, yy, ee, pos );
      % xx does not include the blanks. pos does not either. yy and ee do.
      
      if isempty(xx), return; end
      
      nonblanks = setdiff(1:protocol.nstim,protocol.blankstims);
      
      if isempty(pos)
         % we should not plot any data in this ax
         x = sort(xx);
         y = NaN*[1:length(xx)];
         e = NaN*[1:length(xx)];
      else
         [ x, perm ]  = sort( xx(pos) );
         y = yy(nonblanks(pos)); y = y(perm);
         e = ee(nonblanks(pos)); e = e(perm);
      end
      
      if length(unique(x)) == length(x)
         
         pp = errorbar(x,y,e,'ks-');
         set(pp,'markerfacecolor','k');
         
      else
         pp = errorbar(x,y,e,'ks'); hold on
         set(pp,'markerfacecolor','k');
         xvalues = unique(x);
         yvalues = [];
         for xvalue = xvalues
            yvalues(end+1) = mean(y(x==xvalue));
         end
         plot(xvalues,yvalues,'k-'); 
      end
      hold on;
      if ~isempty(protocol.blankstims)
         plot([min(xx) max(xx)], mean(yy(protocol.blankstims))*[1 1],'k--')
      end
      
      xlabel( graphinfo(irow,icol).xlabel );
      if isfield(graphinfo,'title')
         title( graphinfo(irow,icol).title );
      end
      
      if logflag, set(gca,'xscale','log'); end
      
      set(gca,'xlim',[min(xx)-0.2*range(xx) max(xx)+0.2*range(xx)],'xtick',xticks,...
         'xticklabel',xticklabels);
      
      if isfield(graphinfo(irow,icol),'toprightcomment') && ~isempty(graphinfo(irow,icol).toprightcomment)
         
         text( 1, 1, graphinfo(irow,icol).toprightcomment,...
            'units','normalized',...
            'horizontalalignment','left',...
            'verticalalignment','top' );
      end
      if isfield(graphinfo(irow,icol),'midrightcomment')  && ~isempty(graphinfo(irow,icol).midrightcomment)
         text( 1, 0.5, graphinfo(irow,icol).midrightcomment,...
            'units','normalized',...
            'horizontalalignment','left',...
            'verticalalignment','middle' );
         
      end
   end
end

lims = [ min(yy(:)),max(yy(:))];
set(ax,'ylim', lims+[-0.2 0.2]*diff(lims) ); 

set(ax,'userdata','scaleme!');
% set(ax,'plotboxaspectratio',[1 1 1]); % new (7.11.00)

matchy(ax, 'bottom');

%set(fig,'menubar','none');

ylim = get(ax(1),'ylim');
set(findobj(fig,'tag','txtMaxResp'),'string',ylim(2));

%% write a description 

shortdesc = sprintf('Exp %d-%d Cell %s',unit.iseries, unit.iexp, unit.id );

% desc = [ desc ' - ' protocol.description ];

set(fig,'numbertitle','off','name',shortdesc)
supertitle(shortdesc,1);
