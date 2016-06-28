
% 参数中name, 是gif文件的名字，frames就是抓取的帧，dt为每帧间的间隔。
function res=writegif(name,frames,dt)
  nframe = length(frames);
 
  for i=1:nframe
    [image,map] = frame2im(frames(i));
    [im,map2]          =  rgb2ind(image,128);
    if i==1
      imwrite(im,map2,name,'GIF','WriteMode','overwrite','DelayTime',dt,'LoopCount',inf);
    else
      imwrite(im,map2,name,'WriteMode','append','DelayTime',dt); %,'LoopCount',inf);       
    end
  end