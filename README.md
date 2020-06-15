
Let $A \subset \mathbb{R}^d$ and let $\{Y(\textbf{x}): \textbf{x} \in A\}$ be a Gaussian process with mean function $\mu_{\boldsymbol{\beta}}(\textbf{x}) = \textbf{m}(\textbf{x})^T\boldsymbol{\beta},$ where $\textbf{m}:\mathbb{R}^d \to \mathbb{R}^p$ and $\boldsymbol{\beta} \in \mathbb{R}^p,$ and covariance function $K_{\boldsymbol{\theta}}$. Let $\textbf{x}_1,\dots, \textbf{x}_n \in A$ and suppose we observe $Y_1 = Y(\textbf{x}_1),\dots, Y_n = Y(\textbf{x}_n)$. Let $\textbf{Y} = (Y_1,\dots, Y_n)^T,$ $\textbf{X} = (\textbf{m}(\textbf{x}_1), \dots,  \textbf{m}(\textbf{x}_n) )^T$ and let $\boldsymbol{\Sigma}(\boldsymbol{\theta})$ be an $n\times n$ matrix with $i,j$th entry given by $K_{\boldsymbol{\theta}}(\textbf{x}_i,\textbf{x}_j)$. Then we have
\begin{align}
    \textbf{Y} \sim \mathcal{N}(\textbf{X}\boldsymbol{\beta} , \boldsymbol{\Sigma}(\boldsymbol{\theta})).
\end{align}

Suppose that the $n \times p$ matrix $\textbf{X}$ is full rank. To carry out REML estimation, we need to first write down the joint density of a set of contrasts $\textbf{KY}$ where $\textbf{K}$ any $n-p\times n$ full rank matrix such that $\mathbb{E}[\textbf{KY}] = 0$. Suppose that the first $p$ rows of $\textbf{X}$ are linearly independent, and let $\textbf{Y}_{[j]}$ denote $(Y_1,\dots, Y_j)^T.$ Then the BLUP of $Y_{p+j}$ given $\textbf{Y}_{[p+j-1]}$ exists for $j = 1,\dots, n-p$.
\begin{enumerate}
    \item 
    Let $\boldsymbol{\Sigma}_{[p+j-1]}$ denote the covariance matrix for $\textbf{Y}_{[p+j-1]}$.
    \item 
    Let $\textbf{X}_{[p+j-1]}$ denote the design matrix for for $\textbf{Y}_{[p+j-1]}$.
    \item 
    Let $\textbf{k}_j$ denote the vector of covariances between $\textbf{Y}_{[p+j-1]}$ and $Y_{p+j}$.
    \item 
    Let $\textbf{X}_j$ denote the $j$th row of $\textbf{X}$.
\end{enumerate}
For $j = 1,\dots, n-p,$ let $\boldsymbol{\lambda}_j$ be the first $p+j-1$ entries of the vector 
\begin{align*}
    \begin{pmatrix}
    \boldsymbol{\Sigma}_{[p+j-1]} & \textbf{X}_{[p+j-1]}\\
    \textbf{X}_{[p+j-1]}^T & 0
    \end{pmatrix}^{-1}\begin{pmatrix}\textbf{k}_j\\ \textbf{X}_j\end{pmatrix}
\end{align*}
Let $\textbf{K}$ be an $n-p\times n$ where the $j$th row is given by $(-\boldsymbol{\lambda}^T_j,1,0,\dots,0)$. Then $\textbf{K}$  is full rank and $\mathbb{E}[\textbf{KY}] = 0,$ so $\textbf{W}$ is a suitable set of contrasts. The $j$th entry of $\textbf{W} = \textbf{KY}$ is just the error of the BLUP of $Y_{p+j}$ based on $\textbf{Y}_{[p+j-1]}$.
Consequently, the entries of $\textbf{W}$ are uncorrelated with each other. Since they are also jointly normal, they are independent. Let $\textbf{V} = \textbf{K}\boldsymbol{\Sigma} \textbf{K}^T$.  If $r(\textbf{w};\boldsymbol{\beta},\boldsymbol{\theta})$ denotes the joint density of $\textbf{W},$ then
\begin{align}
    \log r(\textbf{w};\boldsymbol{\beta},\boldsymbol{\theta})&= -\frac{n-p}{2}\log (2\pi) - \frac{1}{2}\log |\textbf{V}| - \frac{1}{2}\textbf{W}^T\textbf{V}^{-1}\textbf{W}\\
    &=\sum_{j = 1}^{n-p}\frac{1}{2}\Big(-\log(2\pi) -\log { \textbf{V}_{jj}} - \textbf{V}_{jj}^{-1}\textbf{W}_{j}^2\Big)
\end{align}
Note that $\textbf{V}_{jj}$ is just the variance of the error of the BLUP of $Y_{p+j}$ based on $\textbf{Y}_{[p+j-1]},$ or equivalently, the mse of the BLUP . Now let $\textbf{S}_{[p+j-1]}\subset \textbf{Y}_{[p+j-1]}$ have $ \text{b} = \min(p+j-1,m)$ entries for $j=1,\dots,n-p$ where $m<< n-p$ ($b$ corresponds to bsize-1 in the code).
Vecchia's approximation of (6) is 
\begin{align}
    \log r(\textbf{w};\boldsymbol{\beta},\boldsymbol{\theta})\approx \sum_{j = 1}^{n-p}\frac{1}{2}\Big(-\log(2\pi) -\log { \tilde{\textbf{V}}_{jj}} - \tilde{\textbf{V}}_{jj}^{-1}\tilde{\textbf{W}}_{j}^2\Big)
\end{align}
where $\tilde{\textbf{W}}_j$ is the error of the BLUP of $Y_{p+j}$ based on $\textbf{S}_{[p+j-1]}$ and $\tilde{\textbf{V}}_{jj}$ is the variance of this error. We can obtain $\tilde{\textbf{W}}_j$ and $\tilde{\textbf{V}}_{jj}$ as follows:
\begin{enumerate}
    \item 
    Let $\tilde{\boldsymbol{\Sigma}}_{[p+j-1]}$ denote the covariance matrix for $\textbf{S}_{[p+j-1]}$.
    \item 
    Let $\tilde{\textbf{X}}_{[p+j-1]}$ denote the design matrix for for $\textbf{S}_{[p+j-1]}$.
    \item 
    Let $\tilde{\textbf{k}}_j$ denote the vector of covariances between $\textbf{S}_{[p+j-1]}$ and $Y_{p+j}$.
    \item 
    Let $\textbf{X}_j$ denote the $j$th row of $\textbf{X}$.
\end{enumerate}
For $j = 1,\dots, n-p,$ let $\tilde{\boldsymbol{\lambda}}_j$ be the first $b$ entries of the vector 
\begin{align*}
    \begin{pmatrix}
    \tilde{\boldsymbol{\Sigma}}_{[p+j-1]} & \tilde{\textbf{X}}_{[p+j-1]}\\
    \tilde{\textbf{X}}_{[p+j-1]}^T & 0
    \end{pmatrix}^{-1}\begin{pmatrix}\tilde{\textbf{k}}_j\\ \textbf{X}_j\end{pmatrix}
\end{align*}
Then $\tilde{\textbf{W}}_j = (-\tilde{\boldsymbol{\lambda}}_j^T,1)\textbf{S}_{[p+j]}$ and $\Tilde{\textbf{V}}_{jj} = (-\tilde{\boldsymbol{\lambda}}_j^T,1)\boldsymbol{\Sigma}_{[p+j]} (-\tilde{\boldsymbol{\lambda}}_j^T,1)^T $. 

We can embed $(-\tilde{\boldsymbol{\lambda}}_j^T,1)$ in an $n$-row-vector of zeros $\textbf{C}_j$ to make $\tilde{\textbf{W}}_j = \textbf{C}_j\textbf{Y}$. Let $\textbf{C}$ be an $n \times (n-p)$ matrix with rows $\textbf{C}_j$.  Then $\tilde{\textbf{W}} :=(\tilde{\textbf{W}}_{11},\dots, \tilde{\textbf{W}}_{n-p,n-p})^T = \textbf{CY} \sim \mathcal{N}(0,\textbf{C}\boldsymbol{\Sigma}\textbf{C}^T)$  where $\tilde{\textbf{V}} :=\textbf{C}\boldsymbol{\Sigma}\textbf{C}^T$ is a diagonal matrix with diagonal entries $\tilde{\textbf{V}}_{jj}$. There are formulas for obtaining $\tilde{\textbf{V}}_{jj}$ directly.


The  formula for the gradient is contained in the Stein et al. paper :
\begin{align*}
   \frac{ \partial} {\partial \theta_k }rl(\theta, w) = -\frac{1}{2} \Big[\Big(\tilde{\textbf{V}}_{jj}^{-1} \cdot  \frac{ \partial} {\partial \theta_k}  \tilde{\textbf{V}}_{jj}\Big) + \Big(2\tilde{\textbf{W}}_j \cdot \tilde{\textbf{V}}_{jj}^{-1} \cdot  \frac{ \partial} {\partial \theta_k} \tilde{\textbf{W}}_j\Big) - \Big( \tilde{\textbf{W}}_j^2 \cdot \tilde{\textbf{V}}_j^{-2} \cdot  \frac{ \partial} {\partial \theta_k} \tilde{\textbf{V}}_j) \Big]
\end{align*}
The Fisher Information matrix can be obtained using the fact that $\tilde{\textbf{W}} \sim \mathcal{N}(0,\textbf{C}\boldsymbol{\Sigma}\textbf{C}^T)$:
\begin{align*}
    \boldsymbol{\mathcal{I}}_{kl}=\frac{1}{2}\sum_{j=1}^{n-p}\Big[\tilde{\textbf{V}}_{jj}^{-1}\cdot \frac{ \partial} {\partial \theta_k}  \tilde{\textbf{V}}_{jj}\cdot\tilde{\textbf{V}}_{jj}^{-1}\cdot  \frac{ \partial} {\partial \theta_l}  \tilde{\textbf{V}}_{jj}\Big].
\end{align*}

The restricted likelihood, gradient and Fisher Information can be computed in one pass through the data, if the formulas above are correct.
