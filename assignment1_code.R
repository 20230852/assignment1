library(ggplot2)

# ---- 1(a) 데이터 생성 및 CSV 저장 ----
dir.create("data",  showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

set.seed(1)
n = 200
x = seq(0, 1, length.out = n)
y = sin(2*pi*x) + rnorm(n, sd = 0.15)

dat = data.frame(x = x, y = y)
write.csv(dat, "data/sim_data.csv", row.names = FALSE)
head(dat)

# ---- 1(b) 산점도와 회귀곡선 ----
dat = read.csv("data/sim_data.csv")

p = ggplot(dat, aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, colour = "red") +
  labs(title = "Scatter plot and fitted regression curve (LOESS)",
       x = "x", y = "y") +
  theme_minimal()
p

# ---- 1(c) 그림 저장 ----
ggsave("plots/regression_fit.png", plot = p, width = 7, height = 5, dpi = 300)
file.exists("plots/regression_fit.png")

# ---- 2(a) 버블 정렬 ----
bubble_sort = function(x, decreasing = FALSE) {
  n = length(x)
  if (n < 2) return(x)
  for (i in 1:(n - 1)) {
    swapped = FALSE
    for (j in 1:(n - i)) {
      # 순서가 어긋났는지 판단 (오름차순: 앞 > 뒤, 내림차순: 앞 < 뒤)
      wrong = if (decreasing) x[j] < x[j + 1] else x[j] > x[j + 1]
      if (wrong) {
        tmp = x[j]; x[j] = x[j + 1]; x[j + 1] = tmp
        swapped = TRUE
      }
    }
    if (!swapped) break   # 조기 종료
  }
  x
}

set.seed(1)
x = runif(10)
x
bs_asc  = bubble_sort(x, decreasing = FALSE)
bs_desc = bubble_sort(x, decreasing = TRUE)
bs_asc
bs_desc
# 내장 sort()와 비교
identical(bs_asc,  sort(x))
identical(bs_desc, sort(x, decreasing = TRUE))

# ---- 2(b) 퀵 정렬 ----
quick_sort = function(x, decreasing = FALSE) {
  if (length(x) <= 1) return(x)          # 종료 조건
  pivot = x[1]
  rest  = x[-1]
  if (decreasing) {
    left  = rest[rest >  pivot]          # pivot보다 큰 값이 왼쪽
    right = rest[rest <= pivot]
  } else {
    left  = rest[rest <  pivot]          # pivot보다 작은 값이 왼쪽
    right = rest[rest >= pivot]
  }
  c(quick_sort(left, decreasing), pivot, quick_sort(right, decreasing))
}

set.seed(1)
x = runif(10)
qs_asc  = quick_sort(x, decreasing = FALSE)
qs_desc = quick_sort(x, decreasing = TRUE)
qs_asc
qs_desc
identical(qs_asc,  sort(x))
identical(qs_desc, sort(x, decreasing = TRUE))

# ---- 3(a) 수치 미분 함수 ----
numerical_diff = function(f, x, h = 1e-6,
                          method = c("forward", "backward", "central")) {
  method = match.arg(method)
  switch(method,
    forward  = (f(x + h) - f(x)) / h,
    backward = (f(x) - f(x - h)) / h,
    central  = (f(x + h) - f(x - h)) / (2 * h)
  )
}

f      = function(x) cos(x) - x
fprime = function(x) -sin(x) - 1       # 해석적 도함수

xs = seq(0, 2*pi, length.out = 100)
d_exact    = fprime(xs)
d_forward  = numerical_diff(f, xs, method = "forward")
d_backward = numerical_diff(f, xs, method = "backward")
d_central  = numerical_diff(f, xs, method = "central")

# 최대 절대오차 비교
data.frame(method = c("forward", "backward", "central"),
           max_abs_error = c(max(abs(d_forward  - d_exact)),
                             max(abs(d_backward - d_exact)),
                             max(abs(d_central  - d_exact))))

# ---- 3(a) 수치 미분 함수 ----
df_d = data.frame(
  x = rep(xs, 4),
  value = c(d_exact, d_forward, d_backward, d_central),
  type = rep(c("Analytic", "Forward", "Backward", "Central"), each = length(xs))
)
ggplot(df_d, aes(x, value, colour = type, linetype = type)) +
  geom_line(linewidth = 0.8) +
  labs(title = "Numerical vs analytic derivative of cos(x) - x",
       x = "x", y = "f'(x)", colour = NULL, linetype = NULL) +
  theme_minimal()

df_e = data.frame(
  x = rep(xs, 3),
  err = pmax(c(abs(d_forward - d_exact), abs(d_backward - d_exact),
               abs(d_central - d_exact)), 1e-16),
  method = rep(c("Forward", "Backward", "Central"), each = length(xs))
)
ggplot(df_e, aes(x, err, colour = method)) +
  geom_line() + scale_y_log10() +
  labs(title = "Absolute error of numerical derivatives (h = 1e-6)",
       x = "x", y = "|error| (log scale)", colour = NULL) +
  theme_minimal()

# ---- 3(b) Newton-Raphson 함수 ----
newton_raphson = function(f, fprime = NULL, x0, maxiter = 100,
                          h = 1e-6, epsilon = 1e-10) {
  # 도함수가 없으면 중심차분 사용
  if (is.null(fprime)) {
    fprime = function(x) numerical_diff(f, x, h = h, method = "central")
  }
  x_old = x0
  for (t in 1:maxiter) {
    d = fprime(x_old)
    if (abs(d) < .Machine$double.eps) stop("도함수 값이 0에 가까워 진행할 수 없습니다.")
    x_new = x_old - f(x_old) / d
    if (abs(x_new - x_old) < epsilon) {
      return(list(root = x_new, f_root = f(x_new), iterations = t, converged = TRUE))
    }
    x_old = x_new
  }
  warning("maxiter 안에 수렴하지 않았습니다.")
  list(root = x_new, f_root = f(x_new), iterations = maxiter, converged = FALSE)
}

# ---- 3(c) cos(x) - x = 0 의 해 (x_0 = 0.5) ----
res_num   = newton_raphson(f, fprime = NULL,   x0 = 0.5)   # 수치미분 버전
res_exact = newton_raphson(f, fprime = fprime, x0 = 0.5)   # 도함수 제공 버전

print(res_num,   digits = 12)
print(res_exact, digits = 12)

# 두 해의 차이
abs(res_num$root - res_exact$root)

# ---- 4(a)-(c) 구현 ----
left_rect = function(f, a, b, n) {
  h = (b - a) / n
  x = a + (0:(n - 1)) * h            # x_0, ..., x_{n-1}
  h * sum(f(x))
}

trapezoid = function(f, a, b, n) {
  h = (b - a) / n
  x = a + (0:n) * h                   # x_0, ..., x_n
  fx = f(x)
  h / 2 * (fx[1] + 2 * sum(fx[2:n]) + fx[n + 1])
}

simpson = function(f, a, b, n) {
  if (n %% 2 != 0) stop("Simpson 방법은 n이 짝수여야 합니다.")
  h = (b - a) / n
  x = a + (0:n) * h
  fx = f(x)                           # fx[i+1] = f_i
  odd  = fx[seq(2, n, by = 2)]        # i = 1, 3, ..., n-1
  even = if (n > 2) fx[seq(3, n - 1, by = 2)] else 0   # i = 2, 4, ..., n-2
  h / 3 * (fx[1] + 4 * sum(odd) + 2 * sum(even) + fx[n + 1])
}

# ---- 4(d) sin x 를 [0,pi]에서 적분 (n=100) ----
true_val = 2
res_d = data.frame(
  method = c("Left Rectangle", "Trapezoid", "Simpson"),
  value  = c(left_rect(sin, 0, pi, 100),
             trapezoid(sin, 0, pi, 100),
             simpson(sin, 0, pi, 100))
)
res_d$abs_error = abs(res_d$value - true_val)
print(res_d, digits = 12)

# ---- 4(e) n에 따른 오차 비교 ----
ns = c(10, 30, 60, 100, 150, 200)
err = data.frame(
  n = ns,
  Left_Rectangle = sapply(ns, function(n) abs(left_rect(sin, 0, pi, n) - true_val)),
  Trapezoid      = sapply(ns, function(n) abs(trapezoid(sin, 0, pi, n) - true_val)),
  Simpson        = sapply(ns, function(n) abs(simpson(sin, 0, pi, n) - true_val))
)
knitr::kable(err, digits = 12, caption = "n별 절대오차 |적분근사 - 2|")

# ---- 4(e) n에 따른 오차 비교 ----
err_long = data.frame(
  n = rep(ns, 3),
  error = c(err$Left_Rectangle, err$Trapezoid, err$Simpson),
  method = rep(c("Left Rectangle", "Trapezoid", "Simpson"), each = length(ns))
)
ggplot(err_long, aes(n, error, colour = method)) +
  geom_line() + geom_point() +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Integration error of sin(x) on [0, pi]",
       x = "n (log scale)", y = "absolute error (log scale)", colour = NULL) +
  theme_minimal()

# ---- 5(a) Cholesky 분해 ----
A = matrix(c(4,2,2,2,5,1,2,1,3), 3)
A
U = chol(A)          # A = U'U (상삼각)
L = t(U)             # 하삼각
L
L %*% t(L)           # A와 같은지 확인
all.equal(L %*% t(L), A)

# ---- 5(b) 전진대입 forward ----
forward = function(L, b) {
  n = length(b)
  z = numeric(n)
  for (i in 1:n) {
    s = if (i > 1) sum(L[i, 1:(i - 1)] * z[1:(i - 1)]) else 0
    z[i] = (b[i] - s) / L[i, i]
  }
  z
}

b = c(1, -2, 3)
z = forward(L, b)
z
forwardsolve(L, b)                       # 내장 함수와 비교
all.equal(z, forwardsolve(L, b))

# ---- 5(c) 후진대입 backward ----
backward = function(L, z) {
  n = length(z)
  x = numeric(n)
  for (i in n:1) {
    s = if (i < n) sum(L[(i + 1):n, i] * x[(i + 1):n]) else 0
    x[i] = (z[i] - s) / L[i, i]
  }
  x
}

x_back = backward(L, z)
x_back
backsolve(t(L), z)                       # L'x = z (t(L)은 상삼각)
all.equal(x_back, backsolve(t(L), z))

# ---- 5(d) Ax=b 풀이 및 solve()와 비교 ----
x_chol  = backward(L, forward(L, b))
x_solve = solve(A, b)

x_chol
x_solve
all.equal(x_chol, x_solve)
max(abs(A %*% x_chol - b))               # 잔차 확인

# ---- 6(a) Gaussian kernel ----
gaussian_kernel = function(x, xprime = x, rho = 1) {
  x      = as.matrix(x)
  xprime = as.matrix(xprime)
  d2 = outer(rowSums(x^2), rowSums(xprime^2), "+") - 2 * x %*% t(xprime)
  d2 = pmax(d2, 0)                      # 수치오차로 인한 음수 방지
  exp(-rho * d2)
}

gaussian_kernel(c(0, 1, 2), rho = 1)    # 간단한 확인

# ---- 6(b) KRR 적합 함수 (class "krr") ----
krr = function(X, y, rho = 1, lambda = 1e-4) {
  X = as.matrix(X)
  n = nrow(X)
  K = gaussian_kernel(X, X, rho)
  alpha = solve(K + lambda * diag(n), y)
  structure(
    list(alpha = as.vector(alpha), X = X, y = y,
         rho = rho, lambda = lambda,
         fitted.values = as.vector(K %*% alpha),
         call = match.call()),
    class = "krr"
  )
}

# ---- 6(c) predict 메서드 ----
predict.krr = function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)   # newdata 없으면 훈련 적합값
  newdata = as.matrix(newdata)
  Kx = gaussian_kernel(newdata, object$X, object$rho)  # m x n
  as.vector(Kx %*% object$alpha)
}

# ---- 6(d) plot 메서드 ----
plot.krr = function(x, ...) {
  if (ncol(x$X) != 1) stop("plot.krr은 1차원 입력에 대해서만 지원합니다.")
  grid = matrix(seq(min(x$X), max(x$X), length.out = 400), ncol = 1)
  pred = predict(x, grid)
  plot(x$X[, 1], x$y, pch = 16, col = "grey50",
       xlab = "x", ylab = "y",
       main = sprintf("Kernel Ridge Regression (rho = %g, lambda = %g)",
                      x$rho, x$lambda), ...)
  lines(grid[, 1], pred, col = "red", lwd = 2)
  legend("topright", legend = c("data", expression(hat(f)(x))),
         pch = c(16, NA), lty = c(NA, 1), col = c("grey50", "red"), bty = "n")
  invisible(x)
}

# ---- 6(e) 시뮬레이션 및 적합 ----
set.seed(1)
n = 150
X = matrix(runif(n, -1, 1), ncol = 1)
ftrue = function(x) sin(2*pi*x) + 0.5*cos(4*pi*x)
y = ftrue(X[, 1]) + rnorm(n, sd = 0.1)

fit = krr(X, y, rho = 1, lambda = 1e-4)   # 기본값 사용
class(fit)

# 새 격자에서 예측
Xnew  = matrix(seq(-1, 1, length.out = 200), ncol = 1)
fnew  = predict(fit, Xnew)
head(fnew)

# 예측 정확도: 참 함수 대비 MSE
mean((fnew - ftrue(Xnew[, 1]))^2)

# ---- 6(e) 시뮬레이션 및 적합 ----
plot(fit)
curve(ftrue(x), from = -1, to = 1, add = TRUE, col = "blue", lty = 2, lwd = 2)
legend("bottomright", legend = "true f(x)", lty = 2, col = "blue", lwd = 2, bty = "n")

# ---- 6(e) 시뮬레이션 및 적합 ----
par(mfrow = c(1, 2))
for (r in c(1, 10)) {
  ft = krr(X, y, rho = r, lambda = 1e-4)
  plot(ft)
  curve(ftrue(x), from = -1, to = 1, add = TRUE, col = "blue", lty = 2)
}
par(mfrow = c(1, 1))
