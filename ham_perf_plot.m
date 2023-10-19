xs = 3:6;
xh = 6:9;
yh = [.013 .045 .242 2.911 72.782];
ys = [.069 .139 .633 4.278 44.518];
yh = [.013 .045 .242 2.911];
ys = [.069 .139 .633 4.278];

hold on
grid
semilogy(xs, ys, 'LineWidth',2)
semilogy(xh, yh, 'LineWidth',2)
% plot(xs, ys)
% plot(xh, yh)
xlabel('Eb/N0 (dB)')
ylabel('Execution time (s)')
legend('Soft', 'Hard')

% axis([3 10 0 80])