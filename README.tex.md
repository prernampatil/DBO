#DBO1D_Codes

We present codes to solve stochastic partial differential equations using the dynamically/bi-orthonormal method. For more details about the results please refer to this [paper](https://arxiv.org/abs/1910.04299).

**Case I: Stochastic linear advection equation**

$$ \frac{\partial u}{\partial t} + V(\omega) \frac{\partial u}{\partial x} = 0 \quad \qquad x \in [0, 2\pi] \quad \mbox{and} \quad  t\in[0,t_f],\\$$
$$          u(x,0) &= \sin(x), \quad \qquad x \in [0, 2\pi]$$

<table style="width:100%">
  <tr>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/LinearAdvection/ErrorPlots/MeanError.png" alt="Mean Error" width="430"/></th>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/LinearAdvection/ErrorPlots/VarError.png" alt="Variance Error" width="430"/></th> 
  </tr>
  <tr>
    <th>Mean Error</th>
    <th>Variance Error</th> 
  </tr>
  </table>


**Case II: Stochastic Burgers' equation with manufactured solution** 

$$ \frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} &= \nu \frac{\partial^2 u}{\partial x^2} + f(x,t; \omega), \quad \qquad x \in [0, 2\pi] \quad \mbox{and} \quad  t\in[0,t_f] \\$$
$$     u(x,0;\omega) &= g(x),    \quad \qquad x \in [0, 2\pi]. $$

<table style="width:100%">
  <tr>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/ErrorPlots/eps_3/Eigenvalues.png" alt="Mean Error" width="430"/></th>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/ErrorPlots/eps_5/Eigenvalues2PI.png" alt="Variance Error" width="430"/></th> 
  </tr>
  <tr>
    <th>Eigenvalues <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-3}" border="0"/>  </th>
    <th>Eigenvalues <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-5}" border="0"/>  </th> 
  </tr>
  </table>
  
  <table style="width:100%">
  <tr>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/ErrorPlots/eps_3/MeanError.png" alt="Mean Error" width="430"/></th>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/ErrorPlots/eps_5/MeanError2PI.png" alt="Variance Error" width="430"/></th> 
  </tr>
  <tr>
    <th>Mean Error <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-3}" border="0"/> </th>
    <th>Mean Error <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-5}" border="0"/> </th> 
  </tr>
  </table>
  
  <table style="width:100%">
  <tr>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/ErrorPlots/eps_3/VarError.png" alt="Mean Error" width="430"/></th>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/ErrorPlots/eps_5/VarError2PI.png" alt="Variance Error" width="430"/></th> 
  </tr>
  <tr>
    <th>Variance Error <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-3}" border="0"/> </th>
    <th>Variance Error <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-5}" border="0"/> </th> 
  </tr>
  </table>
  
  <table style="width:100%">
  <tr>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/Basisplots/epsilon10_3/SpatialModes.gif" alt="Mean Error" width="430"/></th>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/Basisplots/epsilon10_5/SpatialModes.gif" alt="Variance Error" width="430"/></th> 
  </tr>
  <tr>
    <th>Spatial Modes <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-3}" border="0"/> </th>
    <th>Spatial Modes <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-5}" border="0"/> </th> 
  </tr>
  </table> 

<table style="width:100%">
  <tr>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/Phasespace/epsilon10_3/Phasespace.gif" alt="Mean Error" width="430"/></th>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/ManufacturedSolutionPI/Phasespace/epsilon10_5/Phasespace.gif" alt="Variance Error" width="430"/></th> 
  </tr>
  <tr>
    <th>Spatial Modes <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-3}" border="0"/> </th>
    <th>Spatial Modes <img src="http://latex.codecogs.com/svg.latex?\epsilon = 10^{-5}" border="0"/> </th> 
  </tr>
  </table> 



**Case III: Burgers' equation with stochastic forcing**

$$\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} &= \nu \frac{\partial^2 u}{\partial x^2} + \frac{(1+\xi)}{2} \sin(2\pi t),   \quad \qquad x \in [0, 2\pi] \quad \mbox{and} \quad  t\in[0,t_f]\\$$
$$    u(x,0;\omega) &= g(x) \quad \qquad x \in [0, 2\pi]$$

<p align="center">
  <img width="430" src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/StochasticForcing/ErrorPlots/N9LongTime_nu_004_Eigenvalues.png"> </p>

<table style="width:100%">
  <tr>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/StochasticForcing/ErrorPlots/N9LongTime_nu_004_ErrorMean.png" alt="Mean Error" width="430"/></th>
    <th><img src="https://github.com/ppatil1708/DBO/blob/master/DBO_1DCodes/StochasticForcing/ErrorPlots/N9LongTime_nu_004_ErrorVariance.png" alt="Variance Error" width="430"/></th> 
  </tr>
  <tr>
    <th>Mean Error</th>
    <th>Variance Error</th> 
  </tr>
  </table>
