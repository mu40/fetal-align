% Convolve with a Gaussian kernel along dimension. FWHM in pixels/voxels.
% Information about frequency convention and DFT representation:
% en.wikipedia.org/wiki/Fourier_transform#Square-integrable_functions
% gnu.org/software/gsl/manual/html_node/Overview-of-complex-data-FFTs.html
function im = gaussblur(im, fwhm, dim)

if all(fwhm == 0)
    return
end

if isscalar(fwhm)
	fwhm = fwhm * ones(1, ndims(im));
end

imsize = size(im);
if nargin() < 3
    dim = find(imsize > 1);
end

% Convolve with normalized Gaussian. FT of Gaussian is Gaussian, but
% standard deviation and height change. Construct kernel in frequency
% domain to save time. MATLAB uses unitary, ordinary frequency convention
% For normalized Gaussian in image space, height in k-space is 1:
% a = 1 ./ (sigma*sqrt(2*pi));
% ka = a .* sigma * sqrt(2*pi);
sigma = fwhm ./ (2*sqrt(2*log(2)));
ksigma = 1 ./ (sigma*2*pi);

extract = arrayfun(@(p){1:p}, imsize);
for i = dim
    % Pad with mirror image as zeros reduce intensity around edges.
    mirror = flip(im, i);
    padded = cat(i, im, mirror);
    dimsize = 2 * imsize(i);
    % Frequency step is df = 1/T.
    k = (-floor(dimsize/2) : ceil(dimsize/2)-1) / dimsize;
    kernel = exp(-0.5.*(k/ksigma(i)).^2);
    % Want same ordering as FFT output: f = (0, 1, 2, ..., -2, -1)/T.
    kernel = ifftshift(kernel);
    ind = ones(1, dimsize);
    ind(i) = dimsize;
    kernel = reshape(kernel, ind);
    kspace = fft(padded, [], i);
    kspace = kspace .* kernel;
    im = ifft(kspace, [], i);
    im = im(extract{:});    
end
im = real(im); % Should be real anyway.
