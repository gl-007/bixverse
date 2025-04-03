# correlation, covariance, equality tests --------------------------------------

## simple versions -------------------------------------------------------------

set.seed(123)
mat <- matrix(data = rnorm(100),
              nrow = 10,
              ncol = 10)
rownames(mat) <- sprintf("sample_%i", 1:10)
colnames(mat) <- sprintf("feature_%i", 1:10)

# Pearson
expect_equivalent(rs_cor(mat, spearman = FALSE), cor(mat))
# Spearman
expect_equivalent(rs_cor(mat, spearman = TRUE), cor(mat, method = "spearman"))
# Co-variance
expect_equivalent(rs_covariance(mat), cov(mat))

## upper triangle versions -----------------------------------------------------

# Check if the upper triangle class behaves as expected
cor_data <- rs_cor_upper_triangle(mat, spearman = FALSE, shift = 1L)

cor_class <- bixverse:::upper_triangular_cor_mat$new(cor_coef = cor_data,
                                                     features = colnames(mat),
                                                     shift = 1L)

expect_equal(cor_class$get_cor_matrix(.verbose = FALSE), cor(mat))
