function [p, pblock] = ham_err(SNR)

% ========= SCRIPT  ============
% SNRdB  = 9.6   % this is EbNodB
pb = 0.5 * erfc(sqrt(10^(SNR/10)));
ecnodb = SNR + 10*log10(11/15);
p = 0.5 * erfc(sqrt(10^(ecnodb/10)));
 
% formula is nchoosek(15,t+1)*p^(t+1)*(1-p)^(15-(t+1)), where
% t = 1 is the correction capability of the code
pblock = 105 * p^2 * (1-p)^13;
pbit = (2/15)*pblock;
% =============================

% for a better approximation (for larger p values, lower SNR):
pblock = 0;
for e=2:15
    pblock = pblock + nchoosek(15,e)*p^e*(1-p)^(15-e);
end
% in script above, I just took the first term, the dominant term