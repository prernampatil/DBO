#DBO1D_Codes

We present codes to solve stochastic partial differential equations using the dynamically/bi-orthonormal method. For more details about the results please refer to this [paper](https://arxiv.org/abs/1910.04299)

**Case I: Stochastic linear advection equation**

$$ \frac{\partial u}{\partial t} + V(\omega) \frac{\partial u}{\partial x} = 0 \quad \qquad x \in [0, 2\pi] \quad \mbox{and} \quad  t\in[0,t_f],\\$$
$$          u(x,0) &= \sin(x), \quad \qquad x \in [0, 2\pi]$$
![MeanError](LinearAdvection/ErrorPlots/MeanError.pdf)

**Case II: Stochastic Burgers' equation with manufactured solution** 

$$ \frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} &= \nu \frac{\partial^2 u}{\partial x^2} + f(x,t; \omega), \quad \qquad x \in [0, 2\pi] \quad \mbox{and} \quad  t\in[0,t_f] \\$$
$$     u(x,0;\omega) &= g(x),    \quad \qquad x \in [0, 2\pi]. $$
![VarianceError](https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/LinearAdvection/ErrorPlots/VarError.pdf)

**Case III: Burgers' equation with stochastic forcing**

$$\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} &= \nu \frac{\partial^2 u}{\partial x^2} + \frac{(1+\xi)}{2} \sin(2\pi t),   \quad \qquad x \in [0, 2\pi] \quad \mbox{and} \quad  t\in[0,t_f]\\$$
$$    u(x,0;\omega) &= g(x) \quad \qquad x \in [0, 2\pi]$$
