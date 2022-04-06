ORG_IM = imread('office.png'); % readt the original image
A = im2double(ORG_IM);
[m,n]=size(A); %size in pixels for A.
Anext=A; %use A=u(t) and Anext=u(t+h) (at the beginnig they are the same).

% setting the parameters and creaating the gradient of the image at any pixel(x,y)
landa=0.5;
[dx_im_smth, dy_im_smth] = imgradientxy(A);
gr_im_smth = (dx_im_smth.^2 + dy_im_smth.^2);
% element by element g() calculation
[mmm, nnn] = size(gr_im_smth);
g=zeros(mmm,nnn);
for i=1:mmm
    for j=1:nnn
        g(i,j) = 1./(1+(gr_im_smth(i,j)./(landa^.2)));
    end
end

imshow(g);




