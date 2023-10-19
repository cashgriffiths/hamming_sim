x = 3:7;
y_h = [1.27e-1 5.94e-2 2.17e-2 6.05e-3 1.07e-3];
y_s = [5.58e-2 1.78e-2 3.45e-3 5.46e-4 4.86e-5];
y_u = [3.97e-1 2.90e-1 1.58e-1 8.36e-2 4.81e-2];

y_ht = zeros(1,5);
y_st = zeros(1,5);

n=15;

N = zeros(1,n);
N(1) = 1;
N(2) = 15;

for snr = x
    [p, pb] = ham_err(snr);
    y_ht(snr-2) = pb;
end

A = zeros(1,n);
A(1) = 1;
A(2) = 0;
for idx = 1:(n-1)
    A(idx + 2) = (nchoosek(n, idx) - A(idx+1) - (n-idx+1)*A(idx)) / (idx+1);
end

for snr = x
    pblk = 0;
    for w = 3:n
        pblk = pblk + A(w+1)*0.5*erfc(sqrt(w * 11*10^(snr/10)/15));
    end
    y_st(snr-2) = pblk;
end

semilogy(x, y_h,'LineWidth',2), grid, hold on
semilogy(x, y_ht,'LineWidth',2)
semilogy(x, y_s,'LineWidth',2)
semilogy(x, y_st,'LineWidth',2)
semilogy(x, y_u,'LineWidth',2)
xlim([2 8])
xlabel('Eb/N0 (dB)')
ylabel('Block error probability')
legend('Hard', 'Hard (Theory)', 'Soft', 'Soft (Theory)', 'Uncoded')
axis([3 7 1e-5 1])
