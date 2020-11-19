# Numerical integration functions for Lorenz (1996)
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/11/18
# License: LGPL2.1
#############################
# functions for calculation #
#############################

module Lorenz96_functions

using LinearAlgebra
using Statistics

export L96_EU1, L96_RK4, L96_EU1_tangent, L96_RK4_tangent, L96_TL_check, L96_get_background_stat

function L96_tendency(x,F,N)  # calculation of tendency in Lorenz (1996)
    # x: 状態変数
    # F: 強制項の値 (スカラー)
    # N: x の個数
    m_res = fill(0.0,N)
    
    for i in 3:N-1
        m_res[i] = - x[i] + (x[i+1] - x[i-2]) * x[i-1] + F
    end
    # i=1
    m_res[1] = - x[1] + (x[2] - x[N-1]) * x[N] + F
    # i=2
    m_res[2] = - x[2] + (x[3] - x[N]) * x[1] + F
    # i=N
    m_res[N] = - x[N] + (x[1] - x[N-2]) * x[N-1] + F
    
    return m_res
end

function L96_EU1(x,dt,F,N)  # Lorenz (1996) with explicit Euler scheme
    # x: 状態変数
    # dt: 時間ステップ
    # F: 強制項の値 (スカラー)
    # N: x の個数
    m_res = fill(0.0,N)
    
    m_res = x + dt * L96_tendency(x,F,N)
    
    return m_res
end

function L96_RK4(x,dt,F,N)  # Lorenz (1996) with the 4th-order Runge-Kutta scheme
    # x: 状態変数
    # dt: 時間ステップ
    # F: 強制項の値 (スカラー)
    # N: x の個数
    m_res = fill(0.0,N)
    
    k1 = L96_tendency(x,F,N)
    k2 = L96_tendency(x+0.5*dt*k1,F,N)
    k3 = L96_tendency(x+0.5*dt*k2,F,N)
    k4 = L96_tendency(x+dt*k3,F,N)
    
    m_res = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    
    return m_res
end

function L96_EU1_tangent(x,dt,N)  # The tangential operator for L96_EU1
    # x: 状態変数
    # dt: 時間ステップ
    # N: x の個数
    mt_res = fill(0.0,N,N)
    
    # 対角成分
    for i in 1:N
        mt_res[i,i] = (1.0 - dt)
    end
    # k-2
    for i in 3:N
        mt_res[i,i-2] = - dt * x[i-1]
    end
    # k-1
    for i in 3:N-1
        mt_res[i,i-1] = dt * (x[i+1] - x[i-2])
    end
    # k+1
    for i in 2:N-1
        mt_res[i,i+1] = dt * x[i-1]
    end
    # other terms
    mt_res[1,2] = dt * x[N]
    mt_res[N,1] = dt * x[N-1]
    mt_res[1,N] = dt * (x[2] - x[N-1])
    mt_res[2,1] = dt * (x[3] - x[N])
    mt_res[N,N-1] = dt * (x[1] - x[N-2])
    mt_res[1,N-1] = - dt * x[N]
    mt_res[2,N] = - dt * x[1]
    
    return mt_res
end

function L96_RK4_Fdx(x,N)  # The gradient for the forcing term (i.e., dx/dt = F <- this) in L96_RK4
    # x: 状態変数
    # dt: 時間ステップ
    # N: x の個数
    mt_res = fill(0.0,N,N)
    
    # 対角成分
    for i in 1:N
        mt_res[i,i] = -1.0
    end
    # k-2
    for i in 3:N
        mt_res[i,i-2] = -x[i-1]
    end
    # k-1
    for i in 3:N-1
        mt_res[i,i-1] = x[i+1] - x[i-2]
    end
    # k+1
    for i in 2:N-1
        mt_res[i,i+1] = x[i-1]
    end
    # other terms
    mt_res[1,2] = x[N]
    mt_res[N,1] = x[N-1]
    mt_res[1,N] = x[2] - x[N-1]
    mt_res[2,1] = x[3] - x[N]
    mt_res[N,N-1] = x[1] - x[N-2]
    mt_res[1,N-1] = -x[N]
    mt_res[2,N] = -x[1]
    
    return mt_res
end

function L96_RK4_tangent(x,dt,F,N)  # The tangential operator for L96_RK4
    # x: 状態変数
    # dt: 時間ステップ
    # F: 強制項の値 (スカラー)
    # N: x の個数
    mt_res = fill(0.0,N,N)
    dt6 = dt / 6.0
    dt2 = 0.5 * dt
    
    I_mat = Matrix{Float64}(I,N,N)  # 単位行列の利用

    k1 = L96_tendency(x,F,N)
    x2 = x + dt2 .* k1
    k2 = L96_tendency(x2,F,N)
    x3 = x + dt2 .* k2
    k3 = L96_tendency(x3,F,N)
    x4 = x + dt .* k3
    k4 = L96_tendency(x4,F,N)
    
    F1 = L96_RK4_Fdx(x,N)
    F2 = L96_RK4_Fdx(x2,N)
    F3 = L96_RK4_Fdx(x3,N)
    F4 = L96_RK4_Fdx(x4,N)
    
    F12 = (I_mat + dt2 .* F1) * F2
    F123 = (I_mat + dt2 .* F12) * F3
    F1234 = (I_mat + dt .* F123) * F4
    
    mt_res = I_mat + dt6 .* (F1 + 2.0 .* F12 + 2.0 .* F123 + F1234)

    return mt_res
end

function L96_TL_check(x,deltax,dt,nt,F,N,scheme)  # The TL check for L96_{EU1,RK4}_tangent
    # x: 状態変数
    # deltax: x に対する微小変数
    # dt: 時間ステップ
    # nt: 検証する時間ステップ数
    # F: 強制項の値 (スカラー)
    # N: x の個数
    # scheme: 検証するスキームの種類, "EU1" or "RK4"
    xt = x
    delx = deltax
    println("delx ",delx)
    
    if scheme == "EU1"
        println("Start TL check for EU1")
    elseif scheme == "RK4"
        println("Start TL check for RK4")
    end
    
    for i in 1:nt
        if scheme == "EU1"
            Lxt = L96_EU1(xt,dt,F,N)  # Time integration for x
            Lxtd = L96_EU1(xt+delx,dt,F,N)  # Time integration for x + deltax
            Mop = L96_EU1_tangent(xt,dt,N)  # TLM for x
        elseif scheme == "RK4"
            Lxt = L96_RK4(xt,dt,F,N)  # Time integration for x
            Lxtd = L96_RK4(xt+delx,dt,F,N)  # Time integration for x + deltax
            Mop = L96_RK4_tangent(xt,dt,F,N)  # TLM for x
        end
        Ldx = Lxtd - Lxt
        Mdx = Mop * delx
        println("TL check:i=",i,", (dL/dx)δx=",Ldx'*Ldx,", Mδx=",Mdx'*Mdx)
        
        # Updating x and deltax
        xt = Lxt
        delx = Mdx
        println("dx check ",delx' * delx)
    end
    
end

function L96_get_background_stat(x,dt,del_x,F,N,Nsamp,Ninteg)
    # x(N,Nsamp): 状態変数
    # dt: 時間ステップ
    # del_x: 誤差の振幅
    # F: 強制項の値 (スカラー)
    # N: x の個数
    # Nsamp: 共分散計算のサンプル数
    # Nnteg: 誤差の時間発展積分ステップ数
    Mop_linear = reshape(zeros(N,N),N,N)
    gb_res = reshape(zeros(N,N),N,N)
    lambda_max = fill(0.0,Nsamp)
    eigvec_max = reshape(zeros(N,Nsamp),N,Nsamp)
    xt = reshape(zeros(N,1),N,1)
    xe = reshape(zeros(N,1),N,1)
    gb = reshape(zeros(N,N,Nsamp),N,N,Nsamp)
    
    for i in 1:Nsamp
        Mop_linear = L96_RK4_tangent(x[1:N,i],dt,F,N)
        lambda = eigen(Symmetric(Mop_linear' * Mop_linear),N:N)
        lambda_max = lambda.values[1]  # 最大固有値
        eigvec_max[1:N,i] = lambda.vectors[1:N,1]  # 最大固有値の固有ベクトル
        xt = x[1:N,i]
        xe = xt + del_x .* eigvec_max[1:N,i]
        
        for j in 1:Ninteg
            xt = L96_RK4(xt,dt,F,N)
            xe = L96_RK4(xe,dt,F,N)
        end
        
        gb[1:N,1:N,i] = (xe - xt) * (xe - xt)'
    end

    for i in 1:N
        for j in 1:N
            gb_res[i,j] = mean(gb[i,j,1:Nsamp])
        end
    end

    return gb_res
end

end
