# Numerical integration functions for EnKF
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/12/17
# Modification: 2021/01/18
# License: LGPL2.1
#############################
# functions for calculation #
#############################

module EnKF_functions

using LinearAlgebra
using Statistics

export En_mean, Get_ensemble_covariance

function En_mean(x,N,M)  # calculation of ensemble mean (M members) for x[1:N]
    # x: 状態変数 (N x M)
    m_res = fill(0.0,N,M)
    
    for i in 1:N
        m_res[i,1:M] = fill(mean(x[i,1:M]),1,M)
    end
    
    return m_res
end

function Get_ensemble_covariance(x,N,Nsamp)
    # x(N,Nsamp): 状態変数
    # N: x の個数
    # Nsamp: 共分散計算のサンプル数
    B_res = reshape(zeros(N,N),N,N)
    xmean = reshape(zeros(N,Nsamp),N,Nsamp)
    delx = reshape(zeros(N,Nsamp),N,Nsamp)
    
    xmean = En_mean(x,N,Nsamp)
    delx = x - xmean

    B_res = delx * delx'
    B_res = B_res ./ (Nsamp - 1)

    return B_res
end

end
