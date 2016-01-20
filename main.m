clear
clc

actN = 100;
alpha = 5;

distro = gem(actN, alpha);
data = data_generate(distro);

distro_mean = dp_post(data, alpha, actN);

plot(1:actN, distro, '-o', 1:actN, distro_mean, '--*')
title('Comparing the theoretical and the sampled distributions')
ylabel('probabilities')
xlabel('atoms')
legend('theoretical', 'sampled')