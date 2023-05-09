function z = conv_fft(x,y)
%
% Performs a convolution of x and y using FFT methods
% JLG - written specifically for GNSS-R simulator. Assumes column vectors
%
maxlength  = max(size(x)) + max(size(y))-1; 

x0pad = zeros(maxlength,1); 
y0pad = zeros(maxlength,1); 

x0pad(1:max(size(x))) = x; 
y0pad(1:max(size(y))) = y; 

xft = fft(x0pad);
yft = fft(y0pad); 

z = ifft(xft.*yft); 

end

