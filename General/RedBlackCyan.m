function cm = RedBlackCyan(n)

if nargin<1 n = 100; end

cm = [n:-1:0,zeros(1,n) ;
    zeros(1,n), zeros(1,n+1);
    zeros(1,n), 0:n]'/n;

cm(:,1) = cm(end:-1:1,1);
cm(:,2) = cm(end:-1:1,2);
cm(:,3) = cm(end:-1:1,3);

colormap(cm);

% cm = ([n*ones(1,n), n:-1:0 ; ...
%       0:n, n-1:-1:0; ...
%       0:n, ones(1,n)*n]' / n).^gamma;