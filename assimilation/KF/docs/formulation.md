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


# The data assimilation-prediction cycle in the Kalman Filter
- Prediction equations:
\begin{equation}
\textbf{x}^{\mathrm{f}} _{i+1}=M(\textbf{x}^{\mathrm{a}} _i), \tag{KF.1} \label{eq:3-2-1}
\end{equation}
\begin{equation}
\textbf{P} ^{\mathrm{f}} _{i+1}\approx \textbf{M} \textbf{P} ^{\mathrm{a}} _i \textbf{M} ^T. \tag{KF.2} \label{eq:3-2-12}
\end{equation}

  - \\(\textbf{x}^{\mathrm{f}}\\): Forecast (i.e., first guess) variables (\\(N\\)-dimension vector), 
  - \\(\textbf{x}^{\mathrm{a}}\\): Analysis variables (\\(N\\)-dimension vector), 
  - \\(\textbf{P} ^{\mathrm{f}} \\): Background covariance matrix (\\(N\times N\\)), 
  - \\(\textbf{P} ^{\mathrm{a}} \\): Analysis covariance matrix (\\(N\times N\\)), 
  - \\(M\\): Model operator (linear or non-linear), 
  - \\(\textbf{M} \\): Tangent linear operator corresponding to the Model operator, 

- Kalman gain (\\(\textbf{K} _i\\)) equation:
\begin{equation}
\textbf{K} _i=\textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i\left(\textbf{H} _i \textbf{P} ^{\mathrm{f}} _i \textbf{H} ^T_i+\textbf{R} _i \right) ^{-1}. \tag{KF.3} \label{eq:3-2-9}
\end{equation}

- Kalman filter equation: 
\begin{equation}
\textbf{x}^{\mathrm{a}} _i=\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \left[\textbf{y}^{\mathrm{o}} _i -H_i(\textbf{x}^{\mathrm{f}} _i) \right] , \tag{3.2.3} \label{eq:3-2-3}
\end{equation}


- Analysis equation
