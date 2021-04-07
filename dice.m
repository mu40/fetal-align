% Compute Dice similarity metric between arrays of same size.
function out = dice(bw1, bw2)

bw1 = bw1 > 0;
bw2 = bw2 > 0;
out =  2 * nnz(bw1 & bw2) / (nnz(bw1) + nnz(bw2));
