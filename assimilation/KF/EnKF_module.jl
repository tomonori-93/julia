# Numerical integration functions for EnKF
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/12/17
# License: LGPL2.1
#############################
# functions for calculation #
#############################

module EnKF_functions

using LinearAlgebra
using Statistics

export En_mean

function En_mean(x,N,M)  # calculation of ensemble mean (M members) for x[1:N]
    # x: 状態変数 (N x M)
    m_res = fill(0.0,N,M)
    
    for i in 1:N
        m_res[i,1:M] = fill(mean(x[i,1:M]),1,M)
    end
    
    return m_res
end

end
