# KF

Double minimisation of $\chi^2=\phi+\lambda c$
- Inner min.: for each $\lambda$, the minimisation will find a place where $\phi$ and $c$ are tangent;
- Outer min.: then among these $\lambda$, the one that minimises $|c|$ is chosen. Note the critical abs $|\cdot|$ here and only here;
- By definition, $c\rightarrow 0$ is achieved in the outer min.
- The residual $|c|\sim\delta$, multiplied by the corresponding $\lambda$, moves $\phi$ away from a pure $\chi^2=\phi$ minimisation. But the final location is still where $\phi$ and $c$ touch tangentially.
- The outer convergence is judged by MINUIT flag and $\delta<\varepsilon$, where $\varepsilon$ is chosen to accommodate the numerical accuracy. The choice of $\varepsilon$ will affect the convergence rate.
