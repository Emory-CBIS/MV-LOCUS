# MultiView.LOCUS A Low-Rank, Sparse Blind Source Separation for Multi-View Brain Connectome.

`MultiView.LOCUS` is an R package implementing a Low-Rank, Sparse Blind Source Separation framework designed for analyzing Multi-View brain connectivity data. The approach decompose the Multi-View brain connectivity data in to both shared and unique connectivity trait sources across views, addressing challenges inherent in brain connectome analyses such as high dimensionality and noise, providing interpretable and powerful insights into neurodevelopmental and neuropsychiatric studies.


-   I. Installation
-   II. Method
-   III. Detailed Descriptions of the Functions
-   IV. A Toy Example

## I. Package Installation
You can easily install `MultiView.LOCUS` from GitHub with:

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("Emory-CBIS/MultiView.LOCUS ")
library(MultiView.LOCUS)
```


## II. Method

###  MultiView.LOCUS Method
MultiView.LOCUS (** Multi-View Low-rank Decomposition of Brain Connectivity Matrices with Universal Sparsity by Canonical Correlation Analysis**) is a specialized form of Blind Source Separation designed explicitly for Multi-View brain connectivity data. The method uses a low-rank decomposition combined with sparsity-inducing penalties, enabling it to identify robust, interpretable connectivity patterns to each view.

Formally, Locus-CCA identifies $q_k$ latent source **$S^{(k)}$** of each view $k$:

<img src="Fig/decomp.png" width="650" align="center"/>
<img src="Fig/sources.png" width="650" align="center"/>


#### Method Highlights

- **Low-rank factorization**:  
  Locus-CCA employs low-rank decomposition to capture intrinsic, structured patterns within connectivity matrices. This efficiently reduces model complexity and enhances interpretability.

- **Universal Sparsity Regularization**:  
  Sparsity regularization (L1, Hardthreshold or SCAD penalties) is applied element-wise to canonical weights, ensuring robust and interpretable connectivity patterns that represent meaningful neural circuitry associated with clinical or behavioral phenotypes.

- **Common versus view-specific sources**:  
  The common and view-specific sources are considered simutaneously in our model.




### Functions Overview
The structure of the package is as follows, and detailed descriptions of the function arguments are provided in the section below:

-   **Main Function:**
    -   `multi_view_decomposition`: performs Multi-View LOCUS on Multi-View brain connectivity.
-   **Tuning Parameter Selection:**
    -   `calculate_bic`: selects the tuning parameters $\phi$ and $\psi$.
-   **Helper Functions:**
    -   `Ltrinv` and `Ltrans`: transform the brain connectivity to vectorized upper triangle and transform it back.
    -  `plot_conn`:  plots the canonical weights on brain connectivity in the form of heatmap for adjancency connectivity matrix.
-   **Function called**
    -  `Locus_preprocess`:  Preprocessing of connectivity data.
    -  `joint_update_approx`: The  function of fitting Multi-View LOCUS.
    -  `joint_initial`: Initialization for parameters in our model. 
## III. Detailed Descriptions of the Functions

### 1. Multi-View LOCUS function

```         
multi_view_decomposition(Y, q, q_common, V, MaxIteration=5000, penalty="SCAD", phi = 0.9, psi = 1, gamma =3,
                  espli1=5e-4, espli2=5e-4, rho=0.95, silent=FALSE)
```

Arguments
- Y: A list of length K (number of views), where each element is an 𝑛×𝑝 matrix of group-level brain connectivity data for a view.
Each row corresponds to a subject, and each column corresponds to an edge in the connectivity network.
To construct each matrix in Y from subject-level adjacency matrices (size V×V), use the Ltrans() function to vectorize the upper-triangular elements (excluding the diagonal).

- q: A vector of integers of length K. Specifies the number of view-specific connectivity traits to extract for each view.

- q_common: An integer.Specifies the number of shared connectivity traits across all views.

- V: An integer. The number of nodes in the brain network.The number of edges 𝑝 should satisfy $V(V-1)/2$


- MaxIteration: An integer (default = 5000).  The maximum number of iterations for the decomposition algorithm.

- penalty: A string specifying the sparsity regularization method for the source traits. Options include:

    -"NULL": No sparsity enforced.

    -"Hardthreshold": Hard-thresholding penalty.

    -"L1": Lasso penalty (elementwise l1-norm).

    -"SCAD": Smoothly Clipped Absolute Deviation penalty (default), introduced by Fan and Li (2001).

- phi: A numeric value (default = 0.9).
Regularization parameter controlling the strength of the sparsity penalty.

- psi: A numeric value (default = 1).
A coupling parameter enforcing synergy of shared components across different views.

- gamma: A numeric value (default = 3).
Used only when penalty = "SCAD".
Controls the concavity of the SCAD penalty. Must be greater than 2.

- espli1: A numeric value (default = 5e-4).
Tolerance for convergence based on the change in the mixing matrices.

- espli2: A numeric value (default = 5e-4).
Tolerance for convergence based on the change in the source traits.

- rho: A numeric value between 0 and 1 (default = 0.95).  A threshold parameter for determining the low-rank structure of the source traits. A higher rho encourages capturing more variance (leading to a higher rank).

- silent: Logical (TRUE or FALSE, default = FALSE). If FALSE, progress messages are printed during model fitting.


### 2. BIC_cal function

```
calculate_bic(Y, model)
```

- Y: 	A list of original input matrices for each view. Each element should be an 𝑁×𝑝 matrix, where N is the number of subjects and p is the number of connectivity edges (vectorized upper triangle).

-model: The result list returned by multi_view_decomposition(), containing fitted mixing matrices (A) and source matrices (S) for each view.

`BIC_cal` function serves as a valuable guide for tuning the parameters $\phi$ and $\psi$.  The function outputs a single BIC value.  A model with lower BIC value is prefered. However, it is worth noting that in certain datasets, the choice may not be straightforward solely based on BIC. Tuning parameters can also be selected based on visual inspection of the extracted connectivity traits to achieve the desired level of sparsity and appealing neuroscience interpretation.



## IV. A Simulation Example

In this section, we provide a toy example to demonstrate the implementation of the package. We generated toy example data **X**, **Y**, and **z** based on  estimated lantent connectivity traits from real brain connectivity and real clinical subscale dataset on cognition. 
Specifically, we generated connectivity matrices based on the real connectivity traits, using [Power's brain atlas](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3222858/). Each connectivity  trait is symmetric with dimensions of $node \times node$, where $node = 264$ is the number of nodes.   The input $X$ matrix would be of dimension $n \times p$, where $n = 300$ subjects and $p = V(V-1)/2$ edges. Suppose we have $n$ connectivity matrices from each of the $n$ subjects, where each matrix is a $node \times node$ symmetric matrix. To generate our input matrix $Y$, we use the `Ltrans()` function to extract the upper triangular elements of each  matrix and convert them into a row vector of length $p = \frac{(node-1)node}{2}$. We then concatenate these vectors across subjects to obtain the group connectivity data **X**. Similarly, **Y** is a matrix of subscale scores 

``` r
# library 
library(locusCCA.CVRtesting)
library(MASS)  # For ginv() function in data generating only
# generate the toy example data 
S_real_agg <- readRDS(system.file("data", "S_real_agg.rds", package = "locusCCA.CVRtesting"))
  original_Y <- readRDS(system.file("data", "original_Y.rds", package = "locusCCA.CVRtesting"))

# Define parameters
n <- 300
q <- 10
p <- ncol(S_real_agg)
m <- 6
node <- 264
# Simulate X and Y using known structures
# Generate synthetic signals based on the real dataset
U <- t(S_real_agg[ 1:m,]) / 1000
sample1 <- sample(2:13, q)
eigen_Y <- eigen(t(original_Y[, sample1]) %*% original_Y[, sample1])
V <- eigen_Y$vectors[ sample(1:10, q),sample(1:10, m)]

# Simulate X, Y, and z using known structures
set.seed(111)
fx <- matrix(rnorm(n * m, 0, 1), nrow = n)
fy <- fx + matrix(rnorm(n * m, 0, 0.6), nrow = n)
X <- 500 * fx[, 1:m] %*% (ginv(U)) + matrix(rnorm(n * p, 0, 0.01), nrow = n)
Y <- fy[, 1:m] %*% ginv(V) + matrix(rnorm(n * q, 0, 0.01), nrow = n)
weights = rnorm(2,1,0.1)
component = sample(1:6,2)
beta =2000*apply((U[,component] %*% diag(weights)), 1, sum)
z = X %*% beta + rnorm(n,sd = 0.1)

  
# check the dimension
dim(X)
dim(Y)
```

We propose to select the number of canonical correlation components  $m$  based on the  number of PCs needed to explain 95% variance of **Y**.

```r
determine_pca_components <- function(Y, variance_threshold = 0.95) {
  # Perform PCA
  pca_result <- prcomp(Y, center = TRUE, scale. = TRUE)
  
  # Calculate proportion of variance explained
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  
  # Compute cumulative variance explained
  cumulative_variance <- cumsum(variance_explained)
  
  # Find the number of components needed
  num_components <- min(which(cumulative_variance >= variance_threshold))
  
  return(list(num_components = num_components,
              cumulative_variance = cumulative_variance))
}
determine_pca_components(Y)
```

Next, we proceed to use the BIC-type criterion to select the hyperparameters `rho`. In this toy example, we  explore various values for $\rho$ to observe their impact on the BIC value. We recommend initially considering the range $seq(0, 0.05, 0.005)$ to evaluate the BIC.

``` r
# bic selection
rho_seq = seq(0, 0.05, 0.005)
BIC_list = c()
for (rho in rho_seq) {
    result_bic = Locus_CCA(X, Y, node = node, m = m, rho =rho,
                      penalt = "L1", proportion = 0.95,
                      silent = FALSE, tol = 1e-3)
  BIC_list = c(BIC_list,BIC_cal(X,Y,result_bic$U,result_bic$V))}
rho = rho_seq(which.min(BIC_list))
```


It is worth noting that the BIC criterion serves as a valuable guide in selecting the tuning the parameters $\rho$. However, the choice may not always be straightforward solely based on BIC in practice. Therefore, besides the BIC criterion, users can also employ supplementary selection strategies, such as specifying tuning parameters based on the desired sparsity level and the neuroscience interpretations they aim to achieve in the extracted connectivity traits.



Next, we perform the Locus-CCA using the parameters we have just selected.

``` r
## Run Locus-CCA 
result <- Locus_CCA(X, Y, node = node, m = m, rho = rho,
                      penalt = "L1", proportion = 0.95,
                      silent = FALSE, tol = 1e-3)
print(dim(result$U)) #p by m
print(dim(result$V)) #q by m
print(result$CC) #canonical correlation matrix
```
We visualize the canonical direction weights on brain connectivity  based on the Power's atlas. Please note that the visualization code is prepared based on the Power's atlas, and please modify as needed if other atlases are used. 

```r
plots <- lapply(1:m, function(j){
  conn_matrix <- Ltrinv(result$U[, j], node, FALSE)
  plot_conn(conn_matrix)
})

combined_plot <- grid.arrange(grobs = plots, ncol = 3, nrow = 2)

ggsave("combined_plot_with_margin.png", combined_plot,
       width = 15, height = 10,
       units = "in", dpi = 300, limitsize = FALSE)

```
<img src="fig/combined_plot_with_margin.png" width="650" align="center"/>


The CVR testing procedure is then implemented to evaluate the significance of each canonical variants in characterizing the overall response **z**, which gives T_stats of m canonical components. 

## Run CVR testing 
```r
T_stats <- CVR_testing(result$U, X, z)
p.adjust(2*(1-pnorm(abs(T_stats))),,method='fdr') #adjusted p-values using FDR correction to correct multiple (m) testing.
```
