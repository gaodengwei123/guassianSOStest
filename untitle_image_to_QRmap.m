clear all
close all
clc
% I=imread('C:\Users\Dengwei Gao\Pictures\Camera Roll\gaodengwei.jpg');
I=imread('C:\Users\Dengwei Gao\Desktop\aaaab.jpg');
figure
hold on
imshow(I,'DisplayRange',[0 10])
range = [0,0,500,500];
rectangle('Position',range)

if ndims(I) == 3
    I = rgb2gray(I);
end
I2=imresize(I,[500,500]);
[x,y]=size(I2);
X=1:x;
Y=1:y;
[xx,yy]=meshgrid(Y,X);
hh=im2double(I2);
hh(hh>0.8)=10;
hh(hh<0.8)=100;
hh(hh==100) = 1;
hh(hh==10) = 0;


figure;
BRMmap = round(hh);
save BRMmap
% mesh(xx,500-yy,hh);
% yy = 500-yy;
% map = robotics.BinaryOccupancyGrid(map,1/0.1);
% % prmSimple = robotics.PRM(map,100); % 100 PRM nodes in map
% show(map)

% colorbar
% figure;
% imshow(i)
x = reshape(xx,size(xx,1)*size(xx,2),1);
y = reshape(yy,size(yy,1)*size(yy,2),1);
h = 100*reshape(hh,size(hh,1)*size(hh,2),1);

% tic();
% [x, y, h, hm, xm, ym] = generate_terrain(7, 513, 0, 1, 0.05);
% toc();
% figure;
% trisurf(delaunay(x, y), x, y, h);
% colormap gray;                    % Default color map.
% set(gca, 'Position', [0 0 1 1]); % Fill the figure window.
% axis equal vis3d;                % Set aspect ratio.
% shading interp;                  % Interpolate color across faces.
% camlight left;                   % Add a light over to the left somewhere.
% lighting gouraud;                % Use decent lighting.
