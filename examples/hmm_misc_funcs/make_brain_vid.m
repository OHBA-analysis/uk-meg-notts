function make_brain_vid(fileName)

% run make_brain_graph first to get the figure up!

gcf(); % bring to front
view(274,22); % side-on view

% specify a set of views between which to interpolate
el = 10;
viewZ = [-90 el; ...
           0 el; ...
          90 el; ...
         180 el; ...
         270 el];
duration = 8; 
         
% Options = struct('FrameRate', 30, ...
%                  'Duration', duration, ...
%                  'Periodic', false, ...
%                  'Format', 'Motion JPEG AVI');
%              
% CaptureFigVid(viewZ, fileName, Options);

view(270, el);
nFrames = 120;
make_rotating_gif(gcf, fileName, nFrames, duration);


%% make still image
% view(270, el);
saveas(gcf, strrep(fileName, 'gif', 'png'), 'png');