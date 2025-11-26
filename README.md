# KF

Double minimisation of $\chi^2=\phi+\lambda c$
- For each $\lambda$, the minimisation will find a place where $\phi$ and $c$ are tangent---inner minimisation;
- Then among these $\lambda$, the one that minimises $|c|$ (note the abs $|\cdot|$ here, not in inner min.) is chosen---external minimisation;
- By definition, $c\rightarrow 0$ is achieved in the external min.
- The residual $|c|\sim\epsilon$, multiplied by the corresponding $\lambda$ moves $\phi$ away from a pure $\chi^2=\phi$ minimisation. But the final location is still where $\phi$ and $c$ touch tangentially.
- The external convergence is judged by MINUIT flag and $\epsilon<\varepsilon$, where $\varepsilon$ is chosen to accommodate the numerical accuracy. The choice of $\varepsilon$ will affect the convergence rate.
