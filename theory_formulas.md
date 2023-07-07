
**To render file, press Ctrl+Shift+V in VS Code!**

# Coordinate transformation
## 2D integrals
The theory below is adapted from Section 5.1 in [this document](/FEM_pdfs/Gagandeep%20Singh%2C%20Short%20Introduction%20to%20Finite%20Element%20Method.pdf).
Integrals are computed by performing coordinate transformation from arbitrary triangle $T$ to reference triangle $\bar{T}$. The corners of $T$ are $(r_i,z_i)$ with $i=1,2,3$.

$\int_T f(r,z) dr dz = \int_{\bar{T}} f(r(\xi,\eta),z(\xi,\eta)) \; |J| \; d\xi d\eta$

with $J = \begin{pmatrix} \partial_\xi r & \partial_\eta r \\ \partial_\xi z & \partial_\eta z \end{pmatrix} = \begin{pmatrix} r_2-r_1 & r_3-r_1 \\ z_2-z_1 & z_3-z_1 \end{pmatrix}$

Each triangle has 3 non-zero basis functions $\varphi_i(r,z)$ that are equal to 1 in one corner and 0 in the other corners. The equivalent functions in the reference triangle are
- $\bar{N_{1}}(\xi,\eta) = 1-\xi-\eta$
- $\bar{N_{2}}(\xi,\eta) = \xi$
- $\bar{N_{3}}(\xi,\eta) = \eta$

### **Integral 1**
$\int_T r \; \begin{pmatrix}\sigma_r & 0\\ 0 & \sigma_z\end{pmatrix} \nabla\varphi_i \cdot \nabla\varphi_j dr dz\\$
$ = \int_{\bar{T}} (r_1\bar{N_1}+r_2\bar{N_2}+r_3\bar{N_3}) \; (\sigma_r \partial_r\varphi_i \partial_r\varphi_j + \sigma_z \partial_z\varphi_i \partial_z\varphi_j) \; |J| \; d\xi d\eta\\$
$ = \frac{1}{6} (r_1+r_2+r_3) \; (\sigma_r \partial_r\varphi_i \partial_r\varphi_j + \sigma_z\partial_z\varphi_i \partial_z\varphi_j) \; |J| $

Where
- $\partial_r\varphi_1 = \frac{1}{|J|} (z_2-z_3) \;\;\; \partial_z\varphi_1 = \frac{1}{|J|} (r_3-r_2)$
- $\partial_r\varphi_2 = \frac{1}{|J|} (z_3-z_1) \;\;\; \partial_z\varphi_2 = -\frac{1}{|J|} (r_3-r_1)$
- $\partial_r\varphi_3 = -\frac{1}{|J|} (z_2-z_1) \;\;\; \partial_z\varphi_3 = \frac{1}{|J|} (r_2-r_1)$

### **Integral 2**
$\int_T r \; \varphi_i \; \varphi_j dr dz\\$
$= \int_{\bar{T}} (r_1\bar{N_1}+r_2\bar{N_2}+r_3\bar{N_3}) \; \bar{N_{i}} \; \bar{N_{j}} \; |J| \; d\xi d\eta\\$

Integrals $\int_{\bar{T}} \bar{N_{k}} \bar{N_{i}} \bar{N_{j}} d\xi d\eta$ are easy:
- Equal to $\frac{1}{20}$ for $i=k=j$
- Equal to $\frac{1}{60}$ for pair of equals
- Equal to $\frac{1}{120}$ for three different $i,k,j$

### **Integral 3**
$\int_T r \; \varphi_j dr dz\\$
$= \int_{\bar{T}} (r_1\bar{N_1}+r_2\bar{N_2}+r_3\bar{N_3}) \; \bar{N_{j}} \; |J| \; d\xi d\eta\\$
$= \frac{1}{24} \; |J| \; (r_1+r_2+r_3+r_j)$

## 1D integrals
Integrals are computed by performing coordinate transformation from arbitrary edge $E$ to the interval $\xi \in [0,1]$. The endpoints of $E$ are $(r_i,z_i)$ with $i=1,2$.

$\int_E f(r,z) dE = \int_0^1 f(r(\xi),z(\xi)) \; ||E||_2 \; d\xi$

with $||E||_2 = \sqrt{(r_2-r_1)^2 + (z_2-z_1)^2}$

Each triangle has 2 non-zero basis functions $\varphi_i(r,z)$ that are equal to 1 in one endpoint and 0 in the other endpoint. The equivalent functions in the reference triangle are
- $\bar{N_{1}}(\xi,\eta) = 1-\xi$
- $\bar{N_{2}}(\xi,\eta) = \xi$

### **Integral 1**
$\int_E r \; \varphi_j \; dE\\$
$=\int_0^1 (r_1\bar{N_1}+r_2\bar{N_2}) \; \bar{N_j} \; ||E||_2 \; d\xi\\$
$=\frac{1}{6}(r_1+r_2+r_j) \; ||E||_2$

### **Integral 2**
$\int_E r \; \varphi_i \; \varphi_j \; dE\\$
$ = \int_0^1 (r_1\bar{N_1}+r_2\bar{N_2}) \; \bar{N_i} \; \bar{N_j} \; ||E||_2 \; d\xi\\$

Integrals $\int_0^1 \bar{N_i}\bar{N_k}\bar{N_j} d\xi$ are easy:
- Equal to $\frac{1}{4}$ if $i=k=j$
- Else equal to $\frac{1}{12}$

# Newton iterations
### Linear approximation of $R_u$ and $R_v$:
$$
\begin{align*}
R_u(C_u,C_v) = \frac{V_{mu}C_{u}}{(K_{mu}+C_u)(1 + \frac{C_v}{K_{mv}})} &\approx \frac{V_{mu}}{K_{mu}}C_u \\
R_v(C_u,C_v) = r_q R_u(C_u,C_v) + \frac{V_{mfv}}{1+ \frac{C_u}{K_{mfu}}} &\approx r_q \frac{V_{mu}}{K_{mu}} C_u + V_{mfv} \\
\end{align*}
$$
### Linearization of the system:

We linearilize the non-linearity around a point $\begin{pmatrix}c^k_u \\ c^k_v\end{pmatrix}$:

$\begin{pmatrix}
H_u(c_u, c_v) \\ H_v(c_u, c_v)
\end{pmatrix} \approx \begin{pmatrix}
H_u(c^k_u, c^k_v) \\ H_v(c^k_u, c^k_v)
\end{pmatrix} + J\begin{pmatrix}\Delta c_u \\ \Delta c_v\end{pmatrix}$

With $J = \begin{pmatrix}J_1 & J_2 \\ J_3 & J_4\end{pmatrix}
    = \begin{pmatrix}\frac{\delta H_u}{\delta c_u} & \frac{\delta H_u}{\delta c_v} \\
    \frac{\delta H_v}{\delta c_u} &  \frac{\delta H_v}{\delta c_v} \end{pmatrix}$


- $J_1(m, i) = \dfrac{\delta}{\delta c_i} \int_{\Omega}rR_u(C_u, C_v)\varphi_m(r, z)d\Omega$
- $J_2(m, i) = \dfrac{\delta}{\delta c_{M+i}} \int_{\Omega}rR_u(C_u, C_v)\varphi_m(r, z)d\Omega$
- $J_3(m, i) = \dfrac{\delta}{\delta c_i} \int_{\Omega}-rR_v(C_u, C_v)\varphi_m(r, z)d\Omega$
- $J_4(m, i) = \dfrac{\delta}{\delta c_{M+i}} \int_{\Omega}-rR_v(C_u, C_v)\varphi_m(r, z)d\Omega$

We look at the contribution to these integrals by each of the elements. For one triangle with active basis functions $\varphi_1, \varphi_2, \varphi_3$ and corresponding coefficients $c_1, c_2, c_3, c_{M+1}, c_{M+2}, c_{M+3}$, we get the following contributions for each pair ($\varphi_m$,$\varphi_i$):

- $\alpha = V_{mu}K_{mu}|J|\int_{\bar{T}} \dfrac{(r_1\hat{N}_1 + r_2\hat{N}_2 + r_3\hat{N}_3)\hat{N}_i\hat{N}_m}{(1 + \frac{1}{K_{mv}}\sum_{j=1}^{3}c_{M+j}\hat{N}_j)(K_{mu} + \sum_{j=1}^{3}c_{j}\hat{N}_j)^2}d\xi d\eta$

- $\beta = -\dfrac{V_{mu}}{K_{mv}}|J|\int_{\bar{T}} \dfrac{(r_1\hat{N}_1 + r_2\hat{N}_2 + r_3\hat{N}_3)(\sum_{j=1}^{3}c_{j}\hat{N}_j)\hat{N}_i\hat{N}_m}{(1 + \frac{1}{K_{mv}}\sum_{j=1}^{3}c_{M+j}\hat{N}_j)^2(K_{mu} + \sum_{j=1}^{3}c_{j}\hat{N}_j)}d\xi d\eta$

- $\gamma = \dfrac{V_{mfv}}{K_{mfu}}|J|\int_{\bar{T}} \dfrac{(r_1\hat{N}_1 + r_2\hat{N}_2 + r_3\hat{N}_3)\hat{N}_i\hat{N}_m}{(1 + \frac{1}{K_{mfu}}\sum_{j=1}^{3}c_{j}\hat{N}_j)^2}d\xi d\eta$

- $J_1(m,i) = J_1(m,i) + \alpha$

- $J_2(m,i) = J_2(m,i) + \beta$

- $J_3(m,i) = J_3(m, i) -r_q\alpha + \gamma$

- $J_4(m,i) = J_4(m,i) - r_q\beta$

All calculations can be found in the file [theory_jacobian.pdf](theory_jacobian.pdf). Similarly, for evaluating the nonlinear function itself given a point $(c_u,c_v)$, we get the following contributions from a basis function $\varphi_m$ in a triangle:

- $\delta = V_{mu}|J|\int_{\bar{T}} \dfrac{(r_1\hat{N}_1 + r_2\hat{N}_2 + r_3\hat{N}_3)(\sum_{j=1}^{3}c_{j}\hat{N}_j)\hat{N}_m}{(1 + \frac{1}{K_{mv}}\sum_{j=1}^{3}c_{M+j}\hat{N}_j)(K_{mu} + \sum_{j=1}^{3}c_{j}\hat{N}_j)}d\xi d\eta$
- $\varepsilon = V_{mfv}|J|\int_{\bar{T}} \dfrac{(r_1\hat{N}_1 + r_2\hat{N}_2 + r_3\hat{N}_3)\hat{N}_m}{(1 + \frac{1}{K_{mfu}}\sum_{j=1}^{3}c_{j}\hat{N}_j)}d\xi d\eta$
- $H_u(m) = H_u(m) + \delta$
- $H_v(m) = H_v(m) -r_q\delta - \varepsilon$