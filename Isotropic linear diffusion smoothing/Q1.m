ORG_IM = imread('office_noisy.png'); % read the office noisy image
OFFICE_IM = imread('office.png'); % read the original image

A = im2double(ORG_IM);
AA = im2double(OFFICE_IM);

[m,n]=size(A); %size in pixels for A.
Anext=A; % use A=u(t) and Anext=u(t+h) (at the beginnig they are the same).

% print some initial objective quality evaluation MSE, SNR, and PSNR
err = immse(A,AA);
[peaksnr, snr] = psnr(A,AA);

fprintf('---------------------at First---------------------');
fprintf('\n The mean-squared error between office image and office noisy image is %0.4f\n', err);
fprintf('\n peak signal-to-noise ratio between office image and office noisy image is %0.4f\n',peaksnr);
fprintf(' \n The SNR value between office image and office noisy image is %0.4f\n', snr);

% setting the parameters
hx=1;
D=1;
ht=((hx^2)/(4*D))-0.005;
r=(ht*D)/hx^2;

for t=1:200 %time advance 
    for j=2:n-1 %go through the pixels, but avoiding the boundary ones
        for i=2:m-1
            Anext(i,j)=A(i,j)+r*(A(i,j+1)+A(i+1,j)+A(i,j-1)+A(i-1,j)-4*A(i,j)); % update the weights
        end
    end
    A=Anext; %set the updated weight as the current weight for the next round
    if((t==1)||(t==5)||(t==10)||(t==30) ||(t==100))
        
        % plotting two images
        figure; 
        subplot(1,2,1);
        imshow(im2double(ORG_IM));title('Office Noisy image');
        subplot(1,2,2)
        imshow(A);title(strcat('Office Noisy image after a diffusion time = ',int2str(t)));
        pause(3);
        
        err = immse(A,AA);
        [peaksnr, snr] = psnr(A,AA);
        
        fprintf('---------------------in the for---------------------');
        fprintf('\n time is t = %d', t);
        fprintf('\n The mean-squared error between office image and office noisy diffused image is %0.4f\n', err);
        fprintf('\n peak signal-to-noise ratio between office image and office noisy diffused image is %0.4f\n',peaksnr);
        fprintf('\n The SNR value between office image and office noisy diffused image is %0.4f\n', snr);
       
        K = imabsdiff(im2double(ORG_IM),A);
        figure
        imshow(K,[]);title(strcat('The imabsdiff difference between office noisy image and Smooth image at time t = ', int2str(t)));
        pause(3);
    end
end