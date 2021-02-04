<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.6/MathJax.js?config=TeX-AMS_CHTML"></script>

# Basic equation in Lorenz (1996)
The basic equation in Lorenz (1996) is expressed as: 
\begin{equation}
\dfrac{dX_k}{dt} =\left(X_{k+1}-X_{k-2} \right) X_{k-1}-X_k+F,\quad (k=1,\cdots ,N), \tag{L.1} \label{eq:exam-3-2-3-1}
\end{equation}
where the prediction variable \\(X_k\\) has the periodic condition: 
\begin{equation}
X_{k-N}=X_k=X_{k+N}. \tag{L.2} \label{eq:exam-3-2-3-2}
\end{equation}


# Time integration for the forecast
The Lorenz96 equation (\ref{eq:exam-3-2-3-1}): 
\begin{equation}
\dfrac{d\textbf{x}}{dt} =f(t,\textbf{x}),\quad \textbf{x} \equiv (X_1,\cdots ,X_k,\cdots ,X_N)^T \tag{F.1} \label{eq:app-1-1}
\end{equation}
is integrated by the standard 4th-order Runge-Kutta scheme: 
\begin{equation}
\textbf{x} _{i+1} =\textbf{x} _i+\dfrac{\Delta t}{6} \left(\textbf{k} _1+2\textbf{k} _2+2\textbf{k} _3+\textbf{k} _4 \right) =M(\textbf{x} _i) , \nonumber
\end{equation}
\begin{equation}
\textbf{k} _1\equiv f(t_i,\textbf{x} _i) , \nonumber
\end{equation}
\begin{equation}
\textbf{k} _2\equiv f(t_i+\Delta t/2,\textbf{x} _i+(\Delta t/2)\textbf{k} _1) , \nonumber
\end{equation}
\begin{equation}
\textbf{k} _3\equiv f(t_i+\Delta t/2,\textbf{x} _i+(\Delta t/2)\textbf{k} _2) , \nonumber
\end{equation}
\begin{equation}
\textbf{k} _4\equiv f(t_i+\Delta t,\textbf{x} _i+\Delta t\textbf{k} _3) . \nonumber
\end{equation}


# Initialization and configuration of an observation system simulation experiment (OSSE)
* Spin-up experiment (required by making initial conditions for the perfect and forecast models)
  * Initial condition: \\(\textbf{x}^{\mathrm{s}}_0=\left[1.1,\; 1.0,\; \cdots ,\; 1.0 \right] ^T \\)
  * Time integration: \\(N_s\\) steps (for sufficiently long period)
* Perfect model experiment (\\(\textbf{x}^{\mathrm{t}}_i\\))
  * Initial condition: \\(\textbf{x}^{\mathrm{t}}_0=\textbf{x}^{\mathrm{s}}_N \\)
  * Time integration: \\(N_t\\) steps
* Forecast model experiment (\\(\textbf{x}^{\mathrm{f}}_i\\))
  * Initial condition: Temporal average of \\(\textbf{x}^{\mathrm{s}} \\) for a period of 0 to \\(N_s\\)
  * Time integration: \\(N_t\\) steps
* Pseudo-observation (\\(\textbf{y}^{\mathrm{o}}_i\\))
  * Observation grids are located at model grids (without any interpolations)
  * \\(\textbf{y}^{\mathrm{o}}_i=\textbf{x}^{\mathrm{t}}_i+\textbf{r} \\)
    * \\(\textbf{r} \\) is composed of random perturbations with a normal distribution (\\(N(0, \sigma _R)\\))
    * \\(\sigma _R \\) is standard deviations of observation


# The data assimilation-forecast cycles
## (Extended) Kalman Filter
* Forecast equations:
  \begin{equation}
  \textbf{x}^{\mathrm{f}} _{i+1}=M(\textbf{x}^{\mathrm{a}} _i), \tag{KF.1} \label{eq:KF-1}
  \end{equation}
  \begin{equation}
  \textbf{P} ^{\mathrm{f}} _{i+1}\approx \textbf{M} \textbf{P} ^{\mathrm{a}} _i \textbf{M} ^T. \tag{KF.2} \label{eq:KF-2}
  \end{equation}

* Kalman gain (\\(\textbf{K} _i\\)) equation:
  \begin{equation}
  \textbf{K} _i=\textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i\left(\textbf{H} _i \textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i+\textbf{R} _i \right) ^{-1}. \tag{KF.3} \label{eq:KF-3}
  \end{equation}

* Analysis equations: 
  \begin{equation}
  \textbf{x}^{\mathrm{a}} _i=\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \left[\textbf{y}^{\mathrm{o}} _i -H_i(\textbf{x}^{\mathrm{f}} _i) \right] . \tag{KF.4} \label{eq:KF-4}
  \end{equation}
  \begin{equation}
  \textbf{P} ^{\mathrm{a}} _i=\left(\textbf{I} -\textbf{K} _i\textbf{H} _i \right) \textbf{P} ^{\mathrm{f}} _i. \tag{KF.5} \label{eq:KF-5}
  \end{equation}

* Covariance inflation: 
  \begin{equation}
  \textbf{P}^{\mathrm{a}} _i\leftarrow (1+\Delta )\textbf{P}^{\mathrm{a}} _i,\quad (0<\Delta ). \tag{KF.6} \label{eq:KF-6}
  \end{equation}

* Symbols
  * \\(\textbf{x}^{\mathrm{f}}\\): Forecast (i.e., first guess) variables (\\(N\\)-dimension vector), 
  * \\(\textbf{x}^{\mathrm{a}}\\): Analysis variables (\\(N\\)-dimension vector), 
  * \\(\textbf{y}^{\mathrm{o}}\\): Observation variables (\\(p\\)-dimension vector), 
  * \\(\textbf{P} ^{\mathrm{f}} \\): Background covariance matrix (\\(N\times N\\)), 
  * \\(\textbf{P} ^{\mathrm{a}} \\): Analysis covariance matrix (\\(N\times N\\)), 
  * \\(\textbf{R} \\): Observation covariance matrix (\\(p\times p\\)), 
  * \\(M\\): Model operator for time integration (linear or non-linear), 
  * \\(\textbf{M} \equiv \partial M/\partial \textbf{x} \\): Tangent linear operator corresponding to the Model operator (\\(N\times N\\)), 
  * \\(H\\): Observation operator (linear or non-linear), 
  * \\(\textbf{H} \equiv \partial H/\partial \textbf{x} \\): Tangent linear operator corresponding to the Observation operator (\\(p\times N\\)), 


## Singular Evolutive Extended Kalman (SEEK) Filter
* Forecast equations:
  \begin{equation}
  \textbf{x}^{\mathrm{f}} _{i+1}=M(\textbf{x}^{\mathrm{a}} _i), \tag{SEEKF.1} \label{eq:SEEKF-1}
  \end{equation}
  \begin{equation}
  \hat{\textbf{U}}' _{i+1}\approx \textbf{M} \hat{\textbf{U}} _i. \tag{SEEKF.2} \label{eq:SEEKF-2}
  \end{equation}

* Kalman gain (\\(\textbf{K} _i\\)) equation:
  \begin{equation}
  \textbf{K} _i=\hat{\textbf{U}}' _i\hat{\textbf{D}}' _i(\hat{\textbf{U}}' _i)^T \textbf{H} ^T_i\textbf{R} ^{-1}_i. \tag{SEEKF.3} \label{eq:SEEKF-3}
  \end{equation}

* Analysis equations: 
  \begin{equation}
  \textbf{x}^{\mathrm{a}} _i=\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \left[\textbf{y}^{\mathrm{o}} _i -H_i(\textbf{x}^{\mathrm{f}} _i) \right] , \tag{SEEKF.4} \label{eq:SEEKF-4}
  \end{equation}
  \begin{equation}
  \hat{\textbf{D}}' _i\approx \hat{\textbf{D}} _{i-1}-\hat{\textbf{D}} _{i-1}(\hat{\textbf{U}}' _i)^T\textbf{H} ^T_i\left[\textbf{H} _i\hat{\textbf{U}}' _i\hat{\textbf{D}} _{i-1}(\hat{\textbf{U}}' _i)^T\textbf{H} ^T_i+\textbf{R} _i \right] ^{-1}\textbf{H} _i\hat{\textbf{U}}' _i\hat{\textbf{D}} _{i-1}, \tag{SEEKF.5} \label{eq:SEEKF-5}
  \end{equation}
  \begin{equation}
  \hat{\textbf{D}}' _i=\textbf{L} \textbf{L} ^T, \quad (\mathrm{Cholesky\; decomposition}), \tag{SEEKF.6} \label{eq:SEEKF-6}
  \end{equation}
  \begin{equation}
  \textbf{L} ^T(\hat{\textbf{U}}' _i)^T\hat{\textbf{U}}' _i\textbf{L} =\textbf{V} \textbf{E} \textbf{V} ^T, \quad (\mathrm{Eigenvalue\; decomposition}), \tag{SEEKF.7} \label{eq:SEEKF-7}
  \end{equation}
  \begin{equation}
  \hat{\textbf{U}} _i=\hat{\textbf{U}}' _i\textbf{L} \textbf{V} \textbf{E} ^{-1/2},\quad \hat{\textbf{D}} _i=\textbf{E} . \tag{SEEKF.8} \label{eq:SEEKF-8}
  \end{equation}

* Symbols
  * \\(\textbf{x}^{\mathrm{f}}\\): Forecast (i.e., first guess) variables (\\(N\\)-dimension vector), 
  * \\(\textbf{x}^{\mathrm{a}}\\): Analysis variables (\\(N\\)-dimension vector), 
  * \\(\textbf{y}^{\mathrm{o}}\\): Observation variables (\\(p\\)-dimension vector), 
  * \\(\textbf{P} ^{\mathrm{f}} \\): Background covariance matrix (\\(N\times N\\)), 
  * \\(\textbf{P} ^{\mathrm{a}} \\): Analysis covariance matrix (\\(N\times N\\)), 
  * \\(\textbf{R} \\): Observation covariance matrix (\\(p\times p\\)), 
  * \\(M\\): Model operator for time integration (linear or non-linear), 
  * \\(\textbf{M} \equiv \partial M/\partial \textbf{x} \\): Tangent linear operator corresponding to the Model operator (\\(N\times N\\)), 
  * \\(H\\): Observation operator (linear or non-linear), 
  * \\(\textbf{H} \equiv \partial H/\partial \textbf{x} \\): Tangent linear operator corresponding to the Observation operator (\\(p\times N\\)), 


## Local Ensemble Transform Kalman Filter (LETKF)
* Forecast equations (for each ensemble member, m):
  \begin{equation}
  \textbf{X}^{\mathrm{f}} _{i+1}=M(\textbf{X}^{\mathrm{a}} _i), \tag{LETKF.1} \label{eq:LETKF-1}
  \end{equation}
  \begin{equation}
  \textbf{X}\equiv \left[\textbf{x}^{(1)}|\cdots |\textbf{x}^{(m)} \right] =\overline{\textbf{X}} +\delta \textbf{X} , \tag{LETKF.2} \label{eq:LETKF-2}
  \end{equation}
  where \\(\overline{(\; )}\\) means ensemble mean. 

* Analysis equations: 
  \begin{equation}
  \textbf{X}^{\mathrm{a}} _i=\overline{\textbf{X}} ^{\mathrm{f}} _i+\delta \textbf{X} ^{\mathrm{f}} _i\left[\textbf{U} \textbf{D} ^{-1}\textbf{U} ^T(\textbf{H} _i\delta \textbf{X} ^{\mathrm{f}} _i)^T(\textbf{R} _i)^{-1}(\textbf{Y} ^{\mathrm{o}}_i-\overline{H_i(\textbf{X} ^{\mathrm{f}} _i)} )+\; \sqrt[]{m-1} \textbf{U} \textbf{D} ^{1/2}\textbf{U} ^T \right] , \tag{LETKF.3} \label{eq:LETKF-3}
  \end{equation}
  \begin{equation}
  (m-1)\textbf{I}+(\textbf{H} _i\delta \textbf{X} ^{\mathrm{f}} _i)^T(\textbf{R} _i)^{-1}\textbf{H} _i\delta \textbf{X} ^{\mathrm{f}} _i=\textbf{U} \textbf{D} \textbf{U} ^T, \qquad (\mathrm{Eigenvalue\; decomposition}). \tag{LETKF.4} \label{eq:LETKF-4}
  \end{equation}

* Sub equations (Not required in the analysis procedure):
  \begin{equation}
  \textbf{K} _i=\delta \textbf{X} ^{\mathrm{f}} _i\textbf{U} \textbf{D} ^{-1}\textbf{U} ^T(\textbf{H} _i\delta \textbf{X} ^{\mathrm{f}} _i)^T(\textbf{R} _i)^{-1} . \tag{LETKF.5} \label{eq:LETKF-5}
  \end{equation}
  \begin{equation}
  \textbf{P} \equiv \dfrac{1}{m-1} (\delta \textbf{X} )(\delta \textbf{X} )^T. \tag{LETKF.6} \label{eq:LETKF-6}
  \end{equation}

* Covariance inflation (Multiplicative inflation): 
  \begin{equation}
  \delta \textbf{X}^{\mathrm{a}} _i\leftarrow (1+\Delta )\delta \textbf{X}^{\mathrm{a}} _i,\quad (0<\Delta ). \tag{LETKF.7} \label{eq:LETKF-7}
  \end{equation}

* Symbols
  * \\(\textbf{x}^{\mathrm{f}(k)}\\): Forecast (i.e., first guess) variables in the \\(k\\)-th ensemble member (\\(N\\)-dimension vector), 
  * \\(\textbf{x}^{\mathrm{a}(k)}\\): Analysis variables in the \\(k\\)-th ensemble member (\\(N\\)-dimension vector), 
  * \\(\textbf{Y}^{\mathrm{o}}\\): Observation variables (\\(p\times m\\)-dimension vector), 
  * \\(\textbf{R} \\): Observation covariance matrix (\\(p\times p\\)), 
  * \\(M\\): Model operator for time integration (linear or non-linear), 
  * \\(H\\): Observation operator (linear or non-linear), 
  * \\(\textbf{H} \equiv \partial H/\partial \textbf{x} \\): Tangent linear operator corresponding to the Observation operator (\\(p\times N\\)), 
  * \\(m\\): Total ensemble member.


## A hybrid Ensemble Kalman Filter (EnKF)

(Under construction)
* Background covariance
  \begin{equation}
  \textbf{P} ^{\mathrm{f}} _i=\beta \textbf{P} _{\mathrm{stat}}+(1-\beta )\textbf{P} ^{\mathrm{f}} _{\mathrm{flow},i}, \quad (0\leq \beta \leq 1). \tag{HEnKF.1} \label{eq:HEnKF-1}
  \end{equation}

* Symbols
  * \\(\textbf{x}^{\mathrm{f}(k)}\\): Forecast (i.e., first guess) variables in the \\(k\\)-th ensemble member (\\(N\\)-dimension vector), 
  * \\(\textbf{x}^{\mathrm{a}(k)}\\): Analysis variables in the \\(k\\)-th ensemble member (\\(N\\)-dimension vector), 
  * \\(\textbf{Y}^{\mathrm{o}}\\): Observation variables (\\(p\times m\\)-dimension vector), 
  * \\(\textbf{P}_{\mathrm{stat}} \\): Statistical or climatorogical background covariance matrix (\\(N\times N\\)), 
  * \\(\textbf{P} ^{\mathrm{f}} _{\mathrm{flow},i} \\): Flow-dependent background covariance matrix (\\(N\times N\\)), 
  * \\(\textbf{R} \\): Observation covariance matrix (\\(p\times p\\)), 
  * \\(M\\): Model operator for time integration (linear or non-linear), 
  * \\(H\\): Observation operator (linear or non-linear), 
  * \\(\textbf{H} \equiv \partial H/\partial \textbf{x} \\): Tangent linear operator corresponding to the Observation operator (\\(p\times N\\)), 
  * \\(m\\): Total ensemble member.

