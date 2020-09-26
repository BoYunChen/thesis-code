function genGIF(N,frames,dt,str)
    for i=1:N
    [image,map]=frame2im(frames(i));
    [im,map2]=rgb2ind(image,128);
    if i==1
        imwrite(im,map2,str,'gif','writeMode','overwrite','delaytime',dt,'loopcount',inf);
    else
        imwrite(im,map2,str,'gif','writeMode','append','delaytime',dt);
    end
    end
end