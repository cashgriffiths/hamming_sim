function [p, pblock] = ham_err_11_15(SNR)

ecnodb = SNR + 10*log10(11/15);

p = 0.5 * erfc(sqrt(10^(ecnodb/10)));

pblock = 55 * p^2 * (1-p)^9;