%stretched a little bit
[errordelta1(1),errornaive1(1),errorugt1(1),h1(1)] = irre_converge_test(0.0001,[0.1,0]);
[errordelta1(2),errornaive1(2),errorugt1(2),h1(2)] = irre_converge_test(0.0005,[0.1,0]);
[errordelta1(3),errornaive1(3),errorugt1(3),h1(3)] = irre_converge_test(0.001,[0.1,0]);
[errordelta1(4),errornaive1(4),errorugt1(4),h1(4)] = irre_converge_test(0.005,[0.1,0]);
[errordelta1(5),errornaive1(5),errorugt1(5),h1(5)] = irre_converge_test(0.01,[0.1,0]);
[errordelta1(6),errornaive1(6),errorugt1(6),h1(6)] = irre_converge_test(0.05,[0.1,0]);
[errordelta1(7),errornaive1(7),errorugt1(7),h1(7)] = irre_converge_test(0.1,[0.1,0]);
figure(1)
plot(log(h1),log(errordelta1),log(h1),log(errornaive1),log(h1),log(errorugt1))
axis equal
legend('DSC','OSC','GT')
title('Convergence for annulus test on irregular mesh')
xlabel('log avg edge length')
ylabel('log linf error')

%streched a lot
[errordelta(1),errornaive(1),errorugt(1),h(1)] = irre_converge_test(0.0001,[0.05,0]);
[errordelta(2),errornaive(2),errorugt(2),h(2)] = irre_converge_test(0.0005,[0.05,0]);
[errordelta(3),errornaive(3),errorugt(3),h(3)] = irre_converge_test(0.001,[0.05,0]);
[errordelta(4),errornaive(4),errorugt(4),h(4)] = irre_converge_test(0.005,[0.05,0]);
[errordelta(5),errornaive(5),errorugt(5),h(5)] = irre_converge_test(0.01,[0.05,0]);
[errordelta(6),errornaive(6),errorugt(6),h(6)] = irre_converge_test(0.05,[0.05,0]);
[errordelta(7),errornaive(7),errorugt(7),h(7)] = irre_converge_test(0.1,[0.05,0]);
figure(2)
plot(log(h),log(errordelta),log(h),log(errornaive),log(h),log(errorugt))
axis equal
legend('DSC','OSC','GT')
title('Convergence for annulus test on irregular mesh')
xlabel('log avg edge length')
ylabel('log linf error')