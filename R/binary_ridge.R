src = "
arma::vec binaryRidge(const arma::uvec& classes, const arma::vec& response, const double penalty)
{
  const unsigned int n_classes = arma::max(classes) + 1;
  const unsigned int n = response.size();

  arma::vec class_table(n_classes);
  class_table.fill(penalty);
  arma::vec xty(n_classes, arma::fill::zeros);

  unsigned int idx;

  for (unsigned int i = 0; i < n; i++) {
    idx = classes(i);
    class_table(idx) += 1;
    xty(idx) += response(i);
  }
  return (1 / class_table) % xty;
}
"
Rcpp::cppFunction(code = src, depends = "RcppArmadillo", rebuild = TRUE)

nsim = 1e8L

cls = sample(seq_len(5), nsim, replace = TRUE) - 1
y = rnorm(nsim)

library(dplyr)
library(data.table)

dft = data.frame(x = cls, y = y)
dtt = as.data.table(dft)

microbenchmark::microbenchmark(
  "my cpp" = { binaryRidge(cls, y, 0) },
#   "lm" = { lm(y ~ 0 + as.factor(x), data = dft) },
#   "aggregate" = { aggregate(y, by = list(cls), FUN = mean) },
#   "dplyr" = { dft %>% group_by(x) %>% summarize(mn = mean(y)) },
  "data.table" = { dtt[,.(Mean = mean(y)),.(x)] },
  times = 10L
)
