function I = PoissonMixingGradients
    
% read images 
target= im2double(imread('target_2.jpg')); 
source= im2double(imread('source_2.jpg')); 
mask=imread('mask_2.bmp');

row_offset=130;
col_offset=10; 

source_scale=0.6;

source =imresize(source,source_scale);
mask =imresize(mask,source_scale);


N=sum(mask(:)); % N: Number of unknown pixels == variables


% enumerating pixels in the mask
mask_id = zeros(size(mask));
mask_id(mask) = 1:N;   
    
% neighborhood size for each pixel in the mask
[ir,ic] = find(mask);

Np = zeros(N,1); 

for ib=1:N
    
    i = ir(ib);
    j = ic(ib);
    
    Np(ib)=  double((row_offset+i> 1))+ ...
             double((col_offset+j> 1))+ ...
             double((row_offset+i< size(target,1))) + ...
             double((col_offset+j< size(target,2)));
end




% compute matrix A

% We use a sparse matrix since most of the data is zero
% and we only want to track non zero values (optimizing resources)

A = sparse(1:N,1:N,Np,N,N);
for curr_pixel=1:N
    % get row and column of current pixel in the mask
    row = ir(curr_pixel);
    col = ic(curr_pixel);
    A_row = A(curr_pixel,:);
    
    if(mask(row-1,col)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row-1,col));
        A_row(A_row_ind)=-1;
    end
    
    if(mask(row+1,col)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row+1,col));
        A_row(A_row_ind)=-1;
    end
    
    if(mask(row,col-1)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row,col-1));
        A_row(A_row_ind)=-1;
    end
    
    if(mask(row,col+1)~=0)
        A_row_ind = mask_id(sub2ind(size(mask),row,col+1));
        A_row(A_row_ind)=-1;
    end
    
    A(curr_pixel,:)=A_row;
    
end




% output intialization
seamless = target; 


for color=1:3 % solve for each colorchannel

    % compute b for each color
    b=zeros(N,1);
    
    for ib=1:N
    
    i = ir(ib);
    j = ic(ib);
    
      % All 4 if-statements represent a laplacian of the image as they 
      % substract the pixels of the neighbors, but also sum the current 
      % pixel. i.e. equivalent to laplacian kernel= [0,-1,0;-1,4,-1;0,-1,0]
      if (i>1)
          % Left laplacian neighbors
          target_laplacian_left = target(row_offset+i,col_offset+j,color) - ...
              target(row_offset+i-1,col_offset+j,color);
          source_laplacian_left = source(i,j,color)-source(i-1,j,color);
          mixed_laplacian_left = [target_laplacian_left,source_laplacian_left];
          [~,max_idx] = max(abs(mixed_laplacian_left));
          b(ib)=b(ib)+ target(row_offset+i-1,col_offset+j,color)*(1-mask(i-1,j))+...
                        mixed_laplacian_left(max_idx);
                          
      end

      if (i<size(mask,1))
          % Right laplacian neighbors
          target_laplacian_right = target(row_offset+i,col_offset+j,color) - ...
              target(row_offset+i+1,col_offset+j,color);
          source_laplacian_right = source(i,j,color)-source(i+1,j,color);
          mixed_laplacian_right = [target_laplacian_right,source_laplacian_right];
          [~,max_idx] = max(abs(mixed_laplacian_right));
          b(ib)=b(ib)+  target(row_offset+i+1,col_offset+j,color)*(1-mask(i+1,j))+ ...
                        mixed_laplacian_right(max_idx);
      end

      if (j>1)
          % Top laplacian neighbors
          target_laplacian_top = target(row_offset+i,col_offset+j,color) - ...
              target(row_offset+i,col_offset+j-1,color);
          source_laplacian_top = source(i,j,color)-source(i,j-1,color);
          mixed_laplacian_top = [target_laplacian_top,source_laplacian_top];
          [~,max_idx] = max(abs(mixed_laplacian_top));
          b(ib)= b(ib) +  target(row_offset+i,col_offset+j-1,color)*(1-mask(i,j-1))+...
                        mixed_laplacian_top(max_idx);
      end


      if (j<size(mask,2))
          % Bottom laplacian neighbors
          target_laplacian_bot = target(row_offset+i,col_offset+j,color) - ...
              target(row_offset+i,col_offset+j+1,color);
          source_laplacian_bot = source(i,j,color)-source(i,j+1,color);
          mixed_laplacian_bot = [target_laplacian_bot,source_laplacian_bot];
          [~,max_idx] = max(abs(mixed_laplacian_bot));
          b(ib)= b(ib)+ target(row_offset+i,col_offset+j+1,color)*(1-mask(i,j+1))+...
                        mixed_laplacian_bot(max_idx);
      end     


 

    end

    
     % solve linear system A*x = b;
    % your CODE begins here

    x = A\b;
    % your CODE ends here

   
    


    
    % impaint target image
    
     for ib=1:N
           seamless(row_offset+ir(ib),col_offset+ic(ib),color) = x(ib);
     end
     
     figure(1), imshow(target);
     title('Before gradient mixing blending');
     figure(2), imshow(seamless);
     title('After gradient mixing blending');
     I = seamless;
end


