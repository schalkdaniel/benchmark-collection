solve1 = "
arma::mat solveCholesky (const arma::mat& X, const arma::mat& y, const arma::mat U, const bool use_fast)
{
  arma::mat out;
  if (use_fast) {
    out =  arma::solve(arma::trimatl(U.t()), X.t() * y, arma::solve_opts::fast);
    out =  arma::solve(arma::trimatu(U), out, arma::solve_opts::fast);
  } else {
    out =  arma::solve(arma::trimatl(U.t()), X.t() * y);
    out =  arma::solve(arma::trimatu(U), out);
  }
  return out;
}
"

solve2 = "
arma::mat solveInverse (const arma::mat& X, const arma::mat& y, const arma::mat XtX_inv)
{
  return XtX_inv * X.t() * y;
}
"

Rcpp::cppFunction(code = constant1, depends = "RcppArmadillo")
Rcpp::cppFunction(code = constant2, depends = "RcppArmadillo")
Rcpp::cppFunction(code = src1, depends = "RcppArmadillo")
Rcpp::cppFunction(code = src2, depends = "RcppArmadillo")

solveCholesky(X, y, U, FALSE)


n_sim = 50000L
n_num = 50L

mbSolve = function (n_sim, n_numi, reps = 50L, unit = "ms", message = FALSE) {
  mydata = data.frame(matrix(rnorm(n_sim * n_num), ncol = n_num))
  target = runif(n = 1, min = 10, max = 30) + as.matrix(mydata) %*% rgamma(n_num, shape = 2)

  X = model.matrix(~ ., data = mydata)
  XtX_inv = invCpp(t(X) %*% X)
  U = cholCpp(t(X) %*% X)
  y = cbind(target)

  b1 = solveCholesky(X, y, U, TRUE)
  b2 = solveInverse(X, y, XtX_inv)

  if(message && all.equal(b1, b2)) {
    cat("Parameter are EQUAL!!! :)\n")
  }

  res = microbenchmark::microbenchmark(
    solveCholesky(X, y, U, TRUE),
    solveInverse(X, y, XtX_inv),
    times = reps,
    unit = unit
  )
  res_s = summary(res)
  res_df = data.frame(
    min = res_s$min,
    mean = res_s$mean,
    median = res_s$median,
    max = res_s$max,
    method = c("Cholesky", "Inverse"),
    nrow = nrow(X),
    ncol = ncol(X),
    unit = attr(res_s, "unit")
  )
  return (res_df)
}

mbSolve(n_sim, n_num)

n_sims = c(1000, 10000, 50000, 100000, 500000, 1000000)
n_nums = seq(from = 2L, to = 50L, by = 1L)

ll_sim = list()

pb = txtProgressBar(min = min(n_sims), max = max(n_sims), style = 3L)
for (n_sim in n_sims) {
  setTxtProgressBar(pb, n_sim)
  for (n_num in n_nums) {
    ll_sim[[paste0(n_sim, "_", n_num)]] = mbSolve(n_sim, n_num, reps = 20L, unit = "s")
  }
}



library(ggplot2)

df_plt = do.call(rbind, ll_sim)
df_plt$facet = as.factor(paste0("Number of observations: ", df_plt$nrow))
df_plt$facet = factor(df_plt$facet, levels(df_plt$facet)[c(1,2,5,3,6,4)])

gg = ggplot(data = df_plt, aes(x = ncol, y = median, color = method, fill = method)) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  xlab("Number of Parameter") +
  ylab("Runtime in Seconds\n(Median)") +
  labs(color = "", fill = "") +
  theme_light() +
  facet_wrap(~ facet, ncol = 2, scales = "free") +
  ggtitle("Runtime comparison of solving Ab = x by using\ngiven Cholesky decomposition and the inverse of A")

ggsave(filename = "figures/chol-vs-inverse.pdf", plot = gg, width = 8, height = 8)

