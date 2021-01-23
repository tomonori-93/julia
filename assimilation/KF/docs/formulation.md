<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.6/MathJax.js?config=TeX-AMS_CHTML"></script>

(Under construction...)

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
\begin{split}
&\textbf{x} _{i+1} =\textbf{x} _i+\dfrac{\Delta t}{6} \left(\textbf{k} _1+2\textbf{k} _2+2\textbf{k} _3+\textbf{k} _4 \right) =M(\textbf{x} _i) ,\\
&\textbf{k} _1\equiv f(t_i,\textbf{x} _i) , \\
&\textbf{k} _2\equiv f(t_i+\Delta t/2,\textbf{x} _i+(\Delta t/2)\textbf{k} _1) =f(t^{(2)}_i,\textbf{x} ^{(2)}_i) , \\
&\textbf{k} _3\equiv f(t_i+\Delta t/2,\textbf{x} _i+(\Delta t/2)\textbf{k} _2) =f(t^{(3)}_i,\textbf{x} ^{(3)}_i) , \\
&\textbf{k} _4\equiv f(t_i+\Delta t,\textbf{x} _i+\Delta t\textbf{k} _3) =f(t^{(4)}_i,\textbf{x} ^{(4)}_i) .
\end{split}
\tag{F.2} \label{eq:app-1-2}
\end{equation}


# The data assimilation-prediction cycle in the Kalman Filter
- Prediction equations:
\begin{equation}
\textbf{x}^{\mathrm{f}} _{i+1}=M(\textbf{x}^{\mathrm{a}} _i), \tag{KF.1} \label{eq:KF-1}
\end{equation}
\begin{equation}
\textbf{P} ^{\mathrm{f}} _{i+1}\approx \textbf{M} \textbf{P} ^{\mathrm{a}} _i \textbf{M} ^T. \tag{KF.2} \label{eq:KF-2}
\end{equation}

- Kalman gain (\\(\textbf{K} _i\\)) equation:
\begin{equation}
\textbf{K} _i=\textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i\left(\textbf{H} _i \textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i+\textbf{R} _i \right) ^{-1}. \tag{KF.3} \label{eq:KF-3}
\end{equation}

- Analysis equations: 
\begin{equation}
\textbf{x}^{\mathrm{a}} _i=\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \left[\textbf{y}^{\mathrm{o}} _i -H_i(\textbf{x}^{\mathrm{f}} _i) \right] . \tag{KF.4} \label{eq:KF-4}
\end{equation}
\begin{equation}
\textbf{P} ^{\mathrm{a}} _i=\left(\textbf{I} -\textbf{K} _i\textbf{H} _i \right) \textbf{P} ^{\mathrm{f}} _i. \tag{KF.5} \label{eq:KF-5}
\end{equation}

- Symbols
  - \\(\textbf{x}^{\mathrm{f}}\\): Forecast (i.e., first guess) variables (\\(N\\)-dimension vector), 
  - \\(\textbf{x}^{\mathrm{a}}\\): Analysis variables (\\(N\\)-dimension vector), 
  - \\(\textbf{y}^{\mathrm{o}}\\): Observation variables (\\(p\\)-dimension vector), 
  - \\(\textbf{P} ^{\mathrm{f}} \\): Background covariance matrix (\\(N\times N\\)), 
  - \\(\textbf{P} ^{\mathrm{a}} \\): Analysis covariance matrix (\\(N\times N\\)), 
  - \\(\textbf{R} \\): Observation covariance matrix (\\(p\times p\\)), 
  - \\(M\\): Model operator for time integration (linear or non-linear), 
  - \\(\textbf{M} \equiv \partial M/\partial \textbf{x} \\): Tangent linear operator corresponding to the Model operator (\\(N\times N\\)), 
  - \\(H\\): Observation operator (linear or non-linear), 
  - \\(\textbf{H} \equiv \partial H/\partial \textbf{x} \\): Tangent linear operator corresponding to the Observation operator (\\(p\times N\\)), 

# The data assimilation-prediction cycle in the Kalman Filter
- Prediction equations:
\begin{equation}
\textbf{x}^{\mathrm{f}} _{i+1}=M(\textbf{x}^{\mathrm{a}} _i), \tag{KF.1} \label{eq:KF-1}
\end{equation}
\begin{equation}
\textbf{P} ^{\mathrm{f}} _{i+1}\approx \textbf{M} \textbf{P} ^{\mathrm{a}} _i \textbf{M} ^T. \tag{KF.2} \label{eq:KF-2}
\end{equation}

- Kalman gain (\\(\textbf{K} _i\\)) equation:
\begin{equation}
\textbf{K} _i=\textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i\left(\textbf{H} _i \textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i+\textbf{R} _i \right) ^{-1}. \tag{KF.3} \label{eq:KF-3}
\end{equation}

- Analysis equations: 
\begin{equation}
\textbf{x}^{\mathrm{a}} _i=\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \left[\textbf{y}^{\mathrm{o}} _i -H_i(\textbf{x}^{\mathrm{f}} _i) \right] . \tag{KF.4} \label{eq:KF-4}
\end{equation}
\begin{equation}
\textbf{P} ^{\mathrm{a}} _i=\left(\textbf{I} -\textbf{K} _i\textbf{H} _i \right) \textbf{P} ^{\mathrm{f}} _i. \tag{KF.5} \label{eq:KF-5}
\end{equation}

- Symbols
  - \\(\textbf{x}^{\mathrm{f}}\\): Forecast (i.e., first guess) variables (\\(N\\)-dimension vector), 
  - \\(\textbf{x}^{\mathrm{a}}\\): Analysis variables (\\(N\\)-dimension vector), 
  - \\(\textbf{y}^{\mathrm{o}}\\): Observation variables (\\(p\\)-dimension vector), 
  - \\(\textbf{P} ^{\mathrm{f}} \\): Background covariance matrix (\\(N\times N\\)), 
  - \\(\textbf{P} ^{\mathrm{a}} \\): Analysis covariance matrix (\\(N\times N\\)), 
  - \\(\textbf{R} \\): Observation covariance matrix (\\(p\times p\\)), 
  - \\(M\\): Model operator for time integration (linear or non-linear), 
  - \\(\textbf{M} \equiv \partial M/\partial \textbf{x} \\): Tangent linear operator corresponding to the Model operator (\\(N\times N\\)), 
  - \\(H\\): Observation operator (linear or non-linear), 
  - \\(\textbf{H} \equiv \partial H/\partial \textbf{x} \\): Tangent linear operator corresponding to the Observation operator (\\(p\times N\\)), 

# The data assimilation-prediction cycle in the Kalman Filter
- Prediction equations:
\begin{equation}
\textbf{x}^{\mathrm{f}} _{i+1}=M(\textbf{x}^{\mathrm{a}} _i), \tag{KF.1} \label{eq:KF-1}
\end{equation}
\begin{equation}
\textbf{P} ^{\mathrm{f}} _{i+1}\approx \textbf{M} \textbf{P} ^{\mathrm{a}} _i \textbf{M} ^T. \tag{KF.2} \label{eq:KF-2}
\end{equation}

- Kalman gain (\\(\textbf{K} _i\\)) equation:
\begin{equation}
\textbf{K} _i=\textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i\left(\textbf{H} _i \textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i+\textbf{R} _i \right) ^{-1}. \tag{KF.3} \label{eq:KF-3}
\end{equation}

- Analysis equations: 
\begin{equation}
\textbf{x}^{\mathrm{a}} _i=\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \left[\textbf{y}^{\mathrm{o}} _i -H_i(\textbf{x}^{\mathrm{f}} _i) \right] . \tag{KF.4} \label{eq:KF-4}
\end{equation}
\begin{equation}
\textbf{P} ^{\mathrm{a}} _i=\left(\textbf{I} -\textbf{K} _i\textbf{H} _i \right) \textbf{P} ^{\mathrm{f}} _i. \tag{KF.5} \label{eq:KF-5}
\end{equation}

- Symbols
  - \\(\textbf{x}^{\mathrm{f}}\\): Forecast (i.e., first guess) variables (\\(N\\)-dimension vector), 
  - \\(\textbf{x}^{\mathrm{a}}\\): Analysis variables (\\(N\\)-dimension vector), 
  - \\(\textbf{y}^{\mathrm{o}}\\): Observation variables (\\(p\\)-dimension vector), 
  - \\(\textbf{P} ^{\mathrm{f}} \\): Background covariance matrix (\\(N\times N\\)), 
  - \\(\textbf{P} ^{\mathrm{a}} \\): Analysis covariance matrix (\\(N\times N\\)), 
  - \\(\textbf{R} \\): Observation covariance matrix (\\(p\times p\\)), 
  - \\(M\\): Model operator for time integration (linear or non-linear), 
  - \\(\textbf{M} \equiv \partial M/\partial \textbf{x} \\): Tangent linear operator corresponding to the Model operator (\\(N\times N\\)), 
  - \\(H\\): Observation operator (linear or non-linear), 
  - \\(\textbf{H} \equiv \partial H/\partial \textbf{x} \\): Tangent linear operator corresponding to the Observation operator (\\(p\times N\\)), 

