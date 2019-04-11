% Return prior knowledge of fetal anatomy in mm for specific GA in wk.
% Snijders et al. (1994) and Robinson et al. (2008):
% https://www.pedrad.org/Portals/5/Fetal%20MRI/fetal%20eyes.pdf
function [ofd,bpd,odiam,odist] = anatomy(ga)

if ga < 14 || ga > 39
    error('ERROR: GA out of range 14-39 weeks.');
end

odiam = 0.47*ga - 0.7;
bod = 1.47*ga + 1.76;
odist = bod - odiam;
xga = 14:39;
bpd = [28 31 34 36 39 42 45 48 51 54 57 60 63 66 69 72 74 77 79 81 83 ...
    85 86 87 88 89];
ofd = [35 39 42 46 50 54 57 61 65 69 73 77 81 84 87 91 94 96 99 101 103 ...
    105 106 107 107 107];
bpd = interp1(xga, bpd, ga);
ofd = interp1(xga, ofd, ga);
