function make_rotating_gif(handle, fileName, nFrames, oneLoopDuration, useAA)
%MAKE_ROTATING_GIF makes animated gif of rotated figure


if nargin < 5 || ~exist('useAA', 'var'),
    useAA = false;
end

% bring to front
figure(handle);



% pull out elements to rotate
imageContents = [findobj(gca, 'Type', 'surface'); ...
                 findobj(gca, 'Type', 'line');    ...
                 findobj(gca, 'Type', 'patch')];

% fix axes
set(gca, 'xlim', get(gca, 'xlim'));
set(gca, 'ylim', get(gca, 'ylim'));
set(gca, 'zlim', get(gca, 'zlim'));

% define parameters for movie
axis = [0 0 1];
angleIncrement = floor(360 ./ nFrames);
delayTime = oneLoopDuration ./ nFrames;

% save location
[saveDir, fileNameStem] = fileparts(fileName);
saveFileName = fullfile(saveDir, [fileNameStem '.gif']);
             
% write first frame
[A, map] = getImageMap(useAA);
imwrite(A, map, saveFileName, 'gif', 'DelayTime', delayTime, 'LoopCount', Inf);

% rotate and add to gif
for iFrame = 1:nFrames,
    rotate(imageContents, axis, angleIncrement);
    [A, map] = getImageMap(useAA);
    imwrite(A, map, saveFileName, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
end


function [A, map, im] = getImageMap(useAA)
if useAA
    h = myaa;
    cleanMe = onCleanup(@() close(h));
end
frame = getframe();
im = frame2im(frame);

[A, map] = rgb2ind(im, 256, 'nodither');
