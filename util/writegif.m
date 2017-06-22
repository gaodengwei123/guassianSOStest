function writegif(name,frames,dt)
nframe = length(frames);

for i=1:nframe
    image = frame2im(frames(i));
    [im,map2] =  rgb2ind(image,128);
    if i==1
        imwrite(im,map2,name,'gif','DelayTime',dt,'LoopCount',inf);
    else
        imwrite(im,map2,name,'gif','WriteMode','append','DelayTime',dt);
    end
end

end