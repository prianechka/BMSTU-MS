pkg load statistics

EPS = 1e-6;


X = [-0.68, 0.71, 2.27, 0.38, 0.14, 0.06, 1.21, -0.59, 0.44, 1.98, 1.00, ...
       -0.88, -0.08, 1.87, -0.74, 0.83, -1.45, 0.58, 0.48, 3.26, 0.02, 0.26, ...
        2.96, 1.78, 0.58, 0.08, -1.60, 1.26, 1.28, -0.36, 0.15, -0.38, -1.04, ...
        0.95, -2.17, -0.30, 1.09, 0.39, 1.06, 0.98, -2.55, 2.62, -1.58, 3.75, ...
       -1.43, 0.92, 2.75, -0.55, 1.48, -0.96, 0.50, 2.67, -0.58, 0.41, -0.46, ...
       -0.48, 1.68, -0.08, 1.76, 0.08, -1.15, 0.66, 1.54, 0.17, -0.20, 1.34, ...
        1.08, 1.59, -0.05, 0.15, -0.35, 0.58, -0.87, 1.73, -0.27, 0.00, -0.67, ...
        0.13, 1.75, -0.59, 1.31, 1.20, 0.53, 0.14, -0.35, 1.00, -0.01, 0.21, ...
        1.58, -0.02, 1.28, 1.34, -1.66, 0.30, 0.08, 0.66, -0.26, 1.54, 1.22, ...
        1.24, 0.11, 0.79, -0.83, 1.41, 0.17, 0.55, 1.60, 1.26, 1.06, 0.39, ...
        -0.77, 1.49, 0.92, -1.58, 1.19, 0.13, 0.26, -2.14, 0.08, -1.75];

Xmin = min(X);
Xmax = max(X);

fprintf("----------------------------------------\n")
fprintf("Минимальное значение выборки: %d \n", Xmin)
fprintf("Максимальное значение выборки: %d \n", Xmax)
fprintf("----------------------------------------\n")

R = Xmax - Xmin;
printf("Размах выборки: %d \n", R);

N = length(X);
mu = mean(X);
sigma2 = var(X);
sigma = sqrt(sigma2);
correctedSigma2 = N / (N - 1) * sigma2;

fprintf("----------------------------------------\n")
fprintf("Оценка математического ожидания: %f \n", mu)
fprintf("Смещенная оценка дисперсии: %f \n", sigma2)
fprintf("Исправленная оценка дисперсии: %f \n", correctedSigma2)

m = floor(log2(N)) + 2;

intervalBounds = [];
tmp = Xmin;
intervalDelta = R / m;
for i = 1:(m + 1)
  intervalBounds(i) = tmp;
  tmp += intervalDelta;
end

intervalValuesNum = [];

for i = 1:(m - 1)
  tmpCount = 0;
  
  for j = 1:N
    if ((intervalBounds(i) < X(j)) || (abs(intervalBounds(i) - X(j)) < EPS)) ...
      && (X(j) < intervalBounds(i + 1))
      tmpCount += 1;
    endif
  endfor
  
  intervalValuesNum(i) = tmpCount;
endfor

tmpCount = 0;

for j = 1:N
  if (intervalBounds(m) < X(j) || abs(intervalBounds(m) - X(j)) < EPS) && ...
        (X(j) < intervalBounds(m + 1) || abs(intervalBounds(m + 1) - X(j)) < EPS)
    tmpCount += 1;
  endif
endfor

intervalValuesNum(m) = tmpCount;

fprintf("----------------------------------------\n");
fprintf("(г) группировка значений выборки в m = [log_2 n] + 2 интервала:\n");

for i = 1:(m - 1)
  fprintf("    [%f : %f) - %d значений\n",intervalBounds(i), ...
                              intervalBounds(i + 1), intervalValuesNum(i));
end

fprintf("    [%f : %f] - %d значений\n", intervalBounds(m), ...
                                intervalBounds(m + 1), intervalValuesNum(m));

fprintf("----------------------------------------\n");


fprintf("(д) построение гистограммы и графика функции плотности\n");
fprintf("    распределения вероятностей нормальной случайной величины\n");

figure('position',[100,100,1600,1200]);
title ("Гистограмма и график функции плотности нормальной случайной величины");
hold on;
grid on;
        
middleIntervalValues = zeros(1, m);
intervalHeight = zeros(1, m);

for i = 1:m
  intervalHeight(i) = intervalValuesNum(i) / (N * intervalDelta);
endfor

for i = 1:m
  middleIntervalValues(i) = intervalBounds(i + 1) - (intervalDelta / 2);
endfor

fprintf("    высоты столбцов гистограммы:\n");

for i = 1:m
  fprintf("    [%d] : %f\n", i, intervalHeight(i));
endfor

set(gca, "xtick", intervalBounds);
set(gca, "ytick", intervalHeight);
set(gca, "xlim", [min(intervalBounds) - 1, max(intervalBounds) + 1]);
set(gca, "fontsize", 16)
bar(middleIntervalValues, intervalHeight, 1, "facecolor", "blue", ... 
    "edgecolor", "w", "displayname", "Гистограмма");


rangeX = Xmin:(sigma / 100):Xmax;
normalPdf = normpdf(rangeX, mu, sigma);
plot(rangeX, normalPdf, "color", "r", "linewidth", 4, "displayname",  ...
     "Функция плотности нормальной случ.вел.");

myLegend = legend ("location", "northeastoutside");
legend(myLegend, "location", "northeastoutside");
xlabel('X')
ylabel('P')
print -djpg hist.jpg
hold off;

fprintf("----------------------------------------\n");
fprintf("Значения эмпирической функции распределения в точках:\n");

figure('position',[100,100,1600,1200]);
title ("График эмпирической функции распределения и функции распределения нормальной случайной величины");
hold on;
grid on;

m += 2;
intervalCumHeigth = zeros(1, m + 1);

intervalBounds = [(Xmin - intervalDelta) intervalBounds (Xmax + intervalDelta)];
intervalValuesNum = [0 0 intervalValuesNum 0];

curHeigth = 0;

for i = 2:m
  curHeigth += intervalValuesNum(i);
  intervalCumHeigth(i) = curHeigth / N;
end

intervalCumHeigth(m + 1) = 1;

rangeNormX = (Xmin - intervalDelta):(sigma / 100):(Xmax + intervalDelta);
normCdf = normcdf(rangeNormX, mu, sigma);
plot(rangeNormX, normCdf, "color", "r", "linewidth", 2, "displayname",  ...
     "Функция распределения нормальной случ.вел.");

for i = 2:m
  fprintf("x = %f : F(x) = %f\n", intervalBounds(i), intervalCumHeigth(i));
end

set(gca, "xtick", intervalBounds);
set(gca, "ylim", [0, 1.1]);
set(gca, "ytick", intervalCumHeigth);
set(gca, "fontsize", 16)
stairs(intervalBounds, intervalCumHeigth, "color", "blue", "linewidth", 4, ...
      "displayname", "График эмпирической функции распределения");

myLegend = legend ("location", "northeast");
legend(myLegend, "location", "northeast");
xlabel('X')
ylabel('F')
print -djpg cdf.jpg
hold off;