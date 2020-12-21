<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.6/MathJax.js?config=TeX-AMS_CHTML"></script>

(Under construction...)

# Basic equation in Lorenz (1996)
The basic equation in Lorenz (1996) is expressed as: 
\begin{equation}
\dfrac{dX_k}{dt} =\left(X_{k+1}-X_{k-2} \right) X_{k-1}-X_k+F,\quad (k=1,\cdots ,N), \label{eq:exam-3-2-3-1}
\end{equation}
where the prediction variable \\(X_k\\) has the periodic condition: 
\begin{equation}
X_{k-N}=X_k=X_{k+N}. \label{eq:exam-3-2-3-2}
\end{equation}


# The data assimilation-prediction cycle in the Kalman Filter
- Prediction equations:
\begin{equation}
\textbf{x}^{\mathrm{f}} _{i+1}=M(\textbf{x}^{\mathrm{a}} _i), \tag{3.2.1} \label{eq:3-2-1}
\end{equation}

where \\(\textbf{x}^{\mathrm{f}}\\) and \\(\textbf{x}^{\mathrm{a}}\\) are forecast (i.e., first guess) and analysis variables, respectively.

- Kalman gain equation

- Kalman filter equation: 
\begin{equation}
\textbf{x}^{\mathrm{a}} _i=\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \left[\textbf{y}^{\mathrm{o}} _i -H_i(\textbf{x}^{\mathrm{f}} _i) \right] =\textbf{x}^{\mathrm{f}} _i+\textbf{K} _i \textbf{d} _i, \tag{3.2.3} \label{eq:3-2-3}
\end{equation}


- Analysis equation
