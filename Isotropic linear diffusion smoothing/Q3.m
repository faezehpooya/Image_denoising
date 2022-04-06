ORG_IM = imread('office_noisy.png'); % readt the original image
A = im2double(ORG_IM);
[m,n]=size(A); %size in pixels for A.
Anext=A; %use A=u(t) and Anext=u(t+h) (at the beginnig they are the same).

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
    if((t==2)||(t==4)||(t==6)||(t==10))
        
        Sigma_value = sqrt(2*t);
        %h = fspecial('gaussian',[-3*ceil(Sigma_value) 3*ceil(Sigma_value)],Sigma_value); % create the [3*sigma 3*sigma] Gaussian kernel     
        %I_smooth = imfilter(im2double(ORG_IM),h); % perform the Gaussian smoothing on the office noisy image  
        I_smooth = imgaussfilt(im2double(ORG_IM),Sigma_value);
        
        % plotting two images
        figure; 
        subplot(1,2,1);
        imshow(I_smooth);title(strcat('Office Noisy image with Gaussian smoothing with Sigma value= ',num2str(Sigma_value)));
        subplot(1,2,2);
        imshow(A);title(strcat('Office Noisy image after a diffusion time = ',int2str(t)));
        pause(3);
       
       
        err = immse(A, I_smooth);
        [ssimval, ssimmap] = ssim(A,I_smooth);
        
        fprintf('---------------------in the for---------------------');
        fprintf('\n time is %d', t);
        fprintf('\n The mean-squared error is %0.4f\n', err);
        fprintf('\n Structural Similarity index is %0.4f\n', ssimval);
         
        figure, imshow(ssimmap,[]);
        title(strcat('ssim Index Map between diffused noisy image and Smoothed image at time t = ',int2str(t),' Sigma-value= ', num2str(Sigma_value), '- Mean ssim Value is: ', num2str(ssimval)));
        pause(3);
    end
end