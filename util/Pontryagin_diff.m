function [y, Y] = Pontryagin_diff(varargin)


global ellOptions;

if ~isstruct(ellOptions)
    evalin('base', 'ellipsoids_init;');
end

if nargin < 2
    error('MINKDIFF: first and second arguments must be single ellipsoids.');
end

E1 = varargin{1};
E2 = varargin{2};

if ~(isa(E1, 'ellipsoid')) | ~(isa(E2, 'ellipsoid'))
    error('MINKDIFF: first and second arguments must be single ellipsoids.');
end
  

if isbigger(E1, E2) == 0
    switch nargout
        case 0,
            fprintf('Geometric difference of these two ellipsoids is empty set.');
            return;
            
        case 1,
            y = [];
            return;
            
        otherwise,
            y = [];
            Y = [];
            return;
            
    end
end
[q1,Q1] = double(E1);
[q2,Q2] = double(E2);
% check the center fo ellipose

if norm(q1-q2)>1e-6
    error('differenet points in difference')
end
if nargin > 2
    if isstruct(varargin{3})
      Options = varargin{3};
    else
      Options = [];
    end
  else
    Options = [];
  end

  if ~isfield(Options, 'newfigure')
    Options.newfigure = 0;
  end

  if ~isfield(Options, 'fill')
    Options.fill = 0;
  end

  if ~isfield(Options, 'show_all')
    Options.show_all = 0;
  end

  if ~isfield(Options, 'color')
    Options.color = [1 0 0];
  end

  if ~isfield(Options, 'shade')
    Options.shade = 0.4;
  else
    Options.shade = Options.shade(1, 1);
  end


clr  = [1 0 0];
  m    = dimension(E1);
  n    = dimension(E2);
if m ~= n
    error('MINKDIFF: ellipsoids must be of the same dimension.');
end

if nargout == 0
    ih = ishold;
end
  if (Options.show_all ~= 0) & (nargout == 0)
    plot([E1 E2], 'b');
    hold on;
    if Options.newfigure ~= 0
      figure;
    else
      newplot;
    end
  end
    if ellOptions.verbose > 0
    if nargout == 0
      fprintf('Computing and plotting geometric difference of two ellipsoids...\n');
    else
      fprintf('Computing geometric difference of two ellipsoids...\n');
    end
  end
%%
switch n
    case 1,
        y       = q1 - q2;
        Y(1, 1) = q1 - q2 + sqrt(Q2) - sqrt(Q1);
        Y(1, 2) = q1 - q2 + sqrt(Q1) - sqrt(Q2);
        if nargout == 0
            h = ell_plot(Y);
            hold on;
            set(h, 'Color', clr, 'LineWidth', 2);
            h = ell_plot(y, '*');
            set(h, 'Color', clr);
        end
        
    otherwise
        y   = q1 - q2;
           
%         options.x0 = y;
        % random direction in each dimensionality
        K = 10*n^n;
        X = 2*rand(n,K)-1;
        l = X./repmat(sqrt(sum(X.^2,1)),n,1);
        l = rm_bad_directions(Q1, Q2, l);
%         options.num_samples = size(l,2);
%         Y = getLevelSetPdiff(E1.x,E1.getPoly(breaks),l,options);

        if size(l, 2) > 0
            [~, Y] = pdiffrho(E1, l);
            [~, X] = pdiffrho(E2, l);
            Y      = Y - X;
        else
            Y = y;
        end
        if nargout == 0
            vs   = size(Y, 2);
            if vs > 1
                chll = convhulln(Y');
                patch('Vertices', Y', 'Faces', chll, ...
                    'FaceVertexCData', clr(ones(1, vs), :), 'FaceColor', 'flat', ...
                    'FaceAlpha', Options.shade(1, 1));
            else
                h = ell_plot(y, '*');
                set(h, 'Color', clr);
            end
            hold on;
            shading interp;
            lighting phong;
            material('metal');
            view(3);
        end
        
end

if nargout == 0
    if ih == 0
        hold off;
    end
end

if nargout == 1
    y = Y;
end
if nargout == 0
    clear y Y;
end

return;
