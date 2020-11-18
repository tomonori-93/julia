# Numerical integration functions for Lorenz (1963)
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/11/18
# License: LGPL2.1
#############################
# functions for calculation #
#############################

module Lorenz63_functions

using LinearAlgebra

export L63_EU1, L63_RK4, L63_EU1_tangent, L63_TL_check #L63_RK4_tangent, 

function L63_tendency(x,cx,cy,cz)  # calculation of tendency in Lorenz (1963)
    # x: 状態変数
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    m_res = fill(0.0,3)
    
    # i=1
    m_res[1] = - cx * (x[1] - x[2])
    # i=2
    m_res[2] = - x[2] - x[1] * x[3] + cy * x[1]
    # i=3
    m_res[3] = - cz * x[3] + x[1] * x[2]
    
    return m_res
end

function L63_EU1(x,dt,cx,cy,cz)  # Lorenz (1963) with explicit Euler scheme
    # x: 状態変数
    # dt: 時間ステップ
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    m_res = fill(0.0,3)
    
    m_res = x + dt * L63_tendency(x,cx,cy,cz)
    
    return m_res
end

function L63_RK4(x,dt,cx,cy,cz)  # Lorenz (1963) with the 4th-order Runge-Kutta scheme
    # x: 状態変数
    # dt: 時間ステップ
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    m_res = fill(0.0,3)
    
    k1 = L63_tendency(x,cx,cy,cz)
    k2 = L63_tendency(x+0.5*dt*k1,cx,cy,cz)
    k3 = L63_tendency(x+0.5*dt*k2,cx,cy,cz)
    k4 = L63_tendency(x+dt*k3,cx,cy,cz)
    
    m_res = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    
    return m_res
end

function L63_EU1_tangent(x,dt,cx,cy,cz)  # The tangential operator for L63_full
    # x: 状態変数
    # dt: 時間ステップ
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    mt_res = fill(0.0,3,3)
    
    mt_res[1,1] = 1.0 - cx * dt
    mt_res[2,2] = 1.0 - dt
    mt_res[3,3] = 1.0 - cz * dt
    mt_res[1,2] = cx * dt
    mt_res[2,1] = (cy - x[3]) * dt
    mt_res[2,3] = - x[1] * dt
    mt_res[3,1] = x[2] * dt
    mt_res[3,2] = x[1] * dt
    
    return mt_res
end

function L63_TL_check(x,deltax,dt,nt,cx,cy,cz)  # The TL check for L63_tangent
    # x: 状態変数
    # deltax: x に対する微小変数
    # dt: 時間ステップ
    # nt: 検証する時間ステップ数
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    xt = x
    delx = deltax
    println("delx ",delx)
    
    for i in 1:nt
        Lxt = L63_EU1(xt,dt,cx,cy,cz)  # Time integration for x
        Lxtd = L63_EU1(xt+delx,dt,cx,cy,cz)  # Time integration for x + deltax
        Ldx = Lxtd - Lxt
        Mop = L63_EU1_tangent(xt,dt,cx,cy,cz)  # TLM for x
        Mdx = Mop * delx
        println("TL check:i=",i,", (dL/dx)δx=",Ldx'*Ldx,", Mδx=",Mdx'*Mdx)
        
        # Updating x and deltax
        xt = Lxt
        delx = Mdx
        println("dx check ",delx' * delx)
    end
    
end

end
