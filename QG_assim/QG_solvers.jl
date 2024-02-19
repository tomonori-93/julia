module QG_solvers

using LinearAlgebra

export zeta_to_psi1, zeta_to_psi2, grad1, laplace1, zeta_tend, integrate_RK4

function zeta_to_psi1( dx, dy, nx, ny, zeta )
# Diagnosis of zeta to psi

    # intent(in) zeta = reshape(zeros(nx,ny),nx,ny)
    psi = reshape(zeros(nx,ny),nx,ny)
    na = (nx-2) * (ny-2)
    A = reshape(zeros(na,na),na,na)
    b = reshape(zeros(na,1),na,1)
    xb = reshape(zeros(na,1),na,1)
    dx2_inv = 1.0 / (dx * dx)
    dy2_inv = 1.0 / (dy * dy)

    # inner domain
    for j in 3:ny-2
        for i in 3:nx-2
            # Mapping of i and j to A[ix,jy]
            # ix: zeta(i,j) of b(ix), jy: psi(i,j) of x(jy)
            ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
            b[ix] = zeta[i,j]
            
            # (1) psi(i,j)
            jy1 = ix
            A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)
            
            # (2) psi(i-1,j)
            jy2 = (i-2) + (j-2) * (nx-2)
            A[ix,jy2] = dx2_inv
            
            # (3) psi(i+1,j)
            jy3 = (i) + (j-2) * (nx-2)
            A[ix,jy3] = dx2_inv
            
            # (4) psi(i,j-1)
            jy4 = (i-1) + (j-3) * (nx-2)
            A[ix,jy4] = dy2_inv
            
            # (5) psi(i,j+1)
            jy5 = (i-1) + (j-1) * (nx-2)
            A[ix,jy5] = dy2_inv
        end
    end

    # boundaries
    i = 2  # west boundary
    for j in 3:ny-2
        ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
        b[ix] = zeta[i,j]

        # (1) psi(2,j)
        jy1 = ix
        A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)
    
        # (2) psi(nx-1,j)
        jy2 = (nx-2) + (j-2) * (nx-2)
        A[ix,jy2] = dx2_inv
    
        # (3) psi(3,j)
        jy3 = (i) + (j-2) * (nx-2)
        A[ix,jy3] = dx2_inv
    
        # (4) psi(2,j-1)
        jy4 = (i-1) + (j-3) * (nx-2)
        A[ix,jy4] = dy2_inv
    
        # (5) psi(2,j+1)
        jy5 = (i-1) + (j-1) * (nx-2)
        A[ix,jy5] = dy2_inv
    end

    i = nx-1  # east boundary
    for j in 3:ny-2
        ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
        b[ix] = zeta[i,j]
    
        # (1) psi(nx-1,j)
        jy1 = ix
        A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)
    
        # (2) psi(nx-2,j)
        jy2 = (i-2) + (j-2) * (nx-2)
        A[ix,jy2] = dx2_inv
    
        # (3) psi(2,j)
        jy3 = 1 + (j-2) * (nx-2)
        A[ix,jy3] = dx2_inv
    
        # (4) psi(nx-1,j-1)
        jy4 = (i-1) + (j-3) * (nx-2)
        A[ix,jy4] = dy2_inv
    
        # (5) psi(nx-1,j+1)
        jy5 = (i-1) + (j-1) * (nx-2)
        A[ix,jy5] = dy2_inv
    end

    j = 2  # south boundary
    for i in 3:nx-2
        ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
        b[ix] = zeta[i,j]
    
        # (1) psi(i,2)
        jy1 = ix
        A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)
    
        # (2) psi(i-1,2)
        jy2 = (i-2) + (j-2) * (nx-2)
        A[ix,jy2] = dx2_inv
    
        # (3) psi(i+1,2)
        jy3 = (i) + (j-2) * (nx-2)
        A[ix,jy3] = dx2_inv
    
        # (4) psi(i,1)
        b[ix] = b[ix] - psi[i,1] * dy2_inv
    
        # (5) psi(i,3)
        jy5 = (i-1) + (j-1) * (nx-2)
        A[ix,jy5] = dy2_inv
    end

    j = ny-1  # north boundary
    for i in 3:nx-2
        ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
        b[ix] = zeta[i,j]
    
        # (1) psi(i,ny-1)
        jy1 = ix
        A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)
    
        # (2) psi(i-1,ny-1)
        jy2 = (i-2) + (j-2) * (nx-2)
        A[ix,jy2] = dx2_inv
    
        # (3) psi(i+1,ny-1)
        jy3 = (i) + (j-2) * (nx-2)
        A[ix,jy3] = dx2_inv
    
        # (4) psi(i,ny-2)
        jy4 = (i-1) + (j-3) * (nx-2)
        A[ix,jy4] = dy2_inv
    
        # (5) psi(i,ny)
        b[ix] = b[ix] - psi[i,ny] * dy2_inv
    end

    # corner boundaries
    i = 2
    j = 2  # SW-corner
    ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
    b[ix] = zeta[i,j]

    # (1) psi(2,2)
    jy1 = ix
    A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)

    # (2) psi(nx-1,2)
    jy2 = (nx-2) + (j-2) * (nx-2)
    A[ix,jy2] = dx2_inv

    # (3) psi(3,2)
    jy3 = (i) + (j-2) * (nx-2)
    A[ix,jy3] = dx2_inv

    # (4) psi(2,1)
    b[ix] = b[ix] - psi[i,1] * dy2_inv

    # (5) psi(2,3)
    jy5 = (i-1) + (j-1) * (nx-2)
    A[ix,jy5] = dy2_inv

    i = nx-1
    j = 2  # SE-corner
    ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
    b[ix] = zeta[i,j]

    # (1) psi(nx-1,2)
    jy1 = ix
    A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)

    # (2) psi(nx-2,2)
    jy2 = (i-2) + (j-2) * (nx-2)
    A[ix,jy2] = dx2_inv

    # (3) psi(2,2)
    jy3 = 1 + (j-2) * (nx-2)
    A[ix,jy3] = dx2_inv

    # (4) psi(nx-1,1)
    b[ix] = b[ix] - psi[i,1] * dy2_inv

    # (5) psi(nx-1,3)
    jy5 = (i-1) + (j-1) * (nx-2)
    A[ix,jy5] = dy2_inv

    i = 2
    j = ny-1  # NW-corner
    ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
    b[ix] = zeta[i,j]

    # (1) psi(2,ny-1)
    jy1 = ix
    A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)

    # (2) psi(nx-1,ny-1)
    jy2 = (nx-2) + (j-2) * (nx-2)
    A[ix,jy2] = dx2_inv

    # (3) psi(3,ny-1)
    jy3 = (i) + (j-2) * (nx-2)
    A[ix,jy3] = dx2_inv

    # (4) psi(2,ny-2)
    jy4 = (i-1) + (j-3) * (nx-2)
    A[ix,jy4] = dy2_inv

    # (5) psi(2,ny)
    b[ix] = b[ix] - psi[i,ny] * dy2_inv

    i = nx-1
    j = ny-1  # SW-corner
    ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
    b[ix] = zeta[i,j]

    # (1) psi(nx-1,ny-1)
    jy1 = ix
    A[ix,jy1] = - 2.0 * (dx2_inv + dy2_inv)

    # (2) psi(nx-2,ny-1)
    jy2 = (i-2) + (j-2) * (nx-2)
    A[ix,jy2] = dx2_inv

    # (3) psi(2,ny-1)
    jy3 = 1 + (j-2) * (nx-2)
    A[ix,jy3] = dx2_inv

    # (4) psi(nx-1,ny-2)
    jy4 = (i-1) + (j-3) * (nx-2)
    A[ix,jy4] = dy2_inv

    # (5) psi(nx-1,ny)
    b[ix] = b[ix] - psi[i,ny] * dy2_inv
    xb = inv(A) * b


    for j in 2:ny-1
        for i in 2:nx-1
            ix = (i-1) + (j-2) * (nx-2)  # fixed at a (i,j)
            psi[i,j] = xb[ix]
        end
    end
    
    return psi
end

function zeta_to_psi2( n_itr, errmax, dx, dy, nx, ny, zeta, psi_i )
# Diagnosis of zeta to psi

    dx2_inv = 1.0 / (dx * dx)
    dy2_inv = 1.0 / (dy * dy)
    dx2 = dx * dx
    dy2 = dy * dy
    dx2dy2 = 1.0 / (dx2_inv + dy2_inv)

    psi_now = reshape(zeros(nx,ny),nx,ny)
    psi_tmp = reshape(zeros(nx,ny),nx,ny)
    psi_err = reshape(zeros(nx,ny),nx,ny)
    psi_ixm = reshape(zeros(nx-2,ny-2),nx-2,ny-2)
    psi_ixp = reshape(zeros(nx-2,ny-2),nx-2,ny-2)
    psi_jym = reshape(zeros(nx-2,ny-2),nx-2,ny-2)
    psi_jyp = reshape(zeros(nx-2,ny-2),nx-2,ny-2)
    psi = reshape(zeros(nx,ny),nx,ny)
    psi .= psi_i
    psi_tmp .= psi_i
    diffmax = 0.0

    psi_now .= psi_i
    psi_ixm .= psi_i[1:nx-2,2:ny-1]
    psi_ixp .= psi_i[3:nx,2:ny-1]
    psi_jym .= psi_i[2:nx-1,1:ny-2]
    psi_jyp .= psi_i[2:nx-1,3:ny]

    for i in 1:n_itr
        psi_now[2:nx-1,2:ny-1] .= 0.5 * ( ( psi_ixm + psi_ixp ) * dx2_inv + ( psi_jym + psi_jyp ) * dy2_inv - zeta[2:nx-1,2:ny-1] ) * dx2dy2
        
        psi_now[1:nx,1] .= psi_i[1:nx,1]
        psi_now[1:nx,ny] .= psi_i[1:nx,ny]
        psi_now[1,1:ny] .= psi_now[nx-1,1:ny]
        psi_now[nx,1:ny] .= psi_now[2,1:ny]
        
        psi_ixm .= psi_now[1:nx-2,2:ny-1]
        psi_ixp .= psi_now[3:nx,2:ny-1]
        psi_jym .= psi_now[2:nx-1,1:ny-2]
        psi_jyp .= psi_now[2:nx-1,3:ny]
        #println(i, psi_now[nx÷2,ny÷2],"\n")

        psi_err .= abs.(psi_now - psi_tmp)
        diffmax = psi_err[findall(isequal(maximum(psi_err)),psi_err)]
        if diffmax[1] <= errmax
           break
        end
        psi_tmp .= psi_now
    end

    psi[2:nx-1,2:ny-1] .= 0.5 * ( ( psi_ixm + psi_ixp ) * dx2_inv + ( psi_jym + psi_jyp ) * dy2_inv - zeta[2:nx-1,2:ny-1] ) * dx2dy2
    psi[1:nx,1] .= psi_i[1:nx,1]
    psi[1:nx,ny] .= psi_i[1:nx,ny]
    psi[1,1:ny] .= psi_i[1,1:ny]
    psi[nx,1:ny] .= psi_i[nx,1:ny]

    return psi
end


# functions for time integration
function grad1( dirord, dl, nx, ny, nz, val )  # calculation of 1d gradient
    # dirord: the order of the coordinate for the gradient (1: x, 2: y, 3: z)
    # dl: line element for dirord
    # nx,ny,nz: x-y-z grid numbers
    # val(nx,ny,nz): a variable

    dl_inv = 1.0 / dl

    grad1_res = fill(0.0,nx,ny,nz)

    valm = fill(0.0,nx,ny,nz)
    valp = fill(0.0,nx,ny,nz)

    if dirord == 1
        valm[2:nx,1:ny,1:nz] .= val[1:nx-1,1:ny,1:nz]
        valp[1:nx-1,1:ny,1:nz] .= val[2:nx,1:ny,1:nz]
        valm[1,1:ny,1:nz] .= val[nx-1,1:ny,1:nz]
        valp[nx,1:ny,1:nz] .= val[2,1:ny,1:nz]
    elseif dirord == 2
        valm[1:nx,2:ny,1:nz] .= val[1:nx,1:ny-1,1:nz]
        valp[1:nx,1:ny-1,1:nz] .= val[1:nx,2:ny,1:nz]
        valm[1:nx,1,1:nz] .= val[1:nx,1,1:nz]
        valp[1:nx,ny,1:nz] .= val[1:nx,ny,1:nz]
    elseif dirord == 3
        valm[1:nx,1:ny,2:nz] .= val[1:nx,1:ny,1:nz-1]
        valp[1:nx,1:ny,1:nz-1] .= val[1:nx,1:ny,2:nz]
        valm[1:nx,1:ny,1] .= val[1:nx,1:ny,1]
        valp[1:nx,1:ny,nz] .= val[1:nx,1:ny,nz]
    end

    grad1_res .= 0.5 * ( valp - valm ) * dl_inv
    
end


function laplace1( dirord, dl, nx, ny, nz, val )  # calculation of 1d gradient
    # dirord: the order of the coordinate for the gradient (1: x, 2: y, 3: z)
    # dl: line element for dirord
    # nx,ny,nz: x-y-z grid numbers
    # val(nx,ny,nz): a variable

    dl2_inv = 1.0 / (dl * dl)

    grad1_res = fill(0.0,nx,ny,nz)

    valm = fill(0.0,nx,ny,nz)
    valp = fill(0.0,nx,ny,nz)

    if dirord == 1
        valm[2:nx,1:ny,1:nz] .= val[1:nx-1,1:ny,1:nz]
        valp[1:nx-1,1:ny,1:nz] .= val[2:nx,1:ny,1:nz]
        valm[1,1:ny,1:nz] .= val[nx-1,1:ny,1:nz]
        valp[nx,1:ny,1:nz] .= val[2,1:ny,1:nz]
    elseif dirord == 2
        valm[1:nx,2:ny,1:nz] .= val[1:nx,1:ny-1,1:nz]
        valp[1:nx,1:ny-1,1:nz] .= val[1:nx,2:ny,1:nz]
        valm[1:nx,1,1:nz] .= val[1:nx,1,1:nz]
        valp[1:nx,ny,1:nz] .= val[1:nx,ny,1:nz]
    elseif dirord == 3
        valm[1:nx,1:ny,2:nz] .= val[1:nx,1:ny,1:nz-1]
        valp[1:nx,1:ny,1:nz-1] .= val[1:nx,1:ny,2:nz]
        valm[1:nx,1:ny,1] .= val[1:nx,1:ny,1]
        valp[1:nx,1:ny,nz] .= val[1:nx,1:ny,nz]
    end

    grad1_res .= 0.25 * ( valp + valm - 2.0 * val ) * dl2_inv
    
end


function zeta_tend( nu, dx, dy, dz, nx, ny, nz, zeta, psi )  # calculation of the forcing terms in the governing equations
    # nx,ny,nz: x-y-z grid numbers
    # dx,dy,dz: x-y-z grid intervals
    # zeta(nx,ny,nz): vorticity
    # psi(nx,ny,nz): streamfunction

    dzeta_res = fill(0.0,nx,ny,nz)
    dzdx = fill(0.0,nx,ny,nz)
    dzdy = fill(0.0,nx,ny,nz)
    dpdx = fill(0.0,nx,ny,nz)
    dpdy = fill(0.0,nx,ny,nz)
    pdzdx = fill(0.0,nx,ny,nz)
    pdzdy = fill(0.0,nx,ny,nz)
    zdpdx = fill(0.0,nx,ny,nz)
    zdpdy = fill(0.0,nx,ny,nz)
    dpdzxdy = fill(0.0,nx,ny,nz)
    dpdzydx = fill(0.0,nx,ny,nz)
    dzdpxdy = fill(0.0,nx,ny,nz)
    dzdpydx = fill(0.0,nx,ny,nz)
    Δzeta = fill(0.0,nx,ny,nz)
    Δ2zeta = fill(0.0,nx,ny,nz)
    Δ3zeta = fill(0.0,nx,ny,nz)
    adv = fill(0.0,nx,ny,nz)
    diff = fill(0.0,nx,ny,nz)

    dzdx .= grad1( 1, dx, nx, ny, nz, zeta )
    dzdy .= grad1( 2, dy, nx, ny, nz, zeta )
    dpdx .= grad1( 1, dx, nx, ny, nz, psi )
    dpdy .= grad1( 2, dy, nx, ny, nz, psi )

    pdzdx .= psi .* dzdx
    pdzdy .= psi .* dzdy
    zdpdx .= zeta .* dpdx
    zdpdy .= zeta .* dpdy

    dpdzxdy .= grad1( 2, dy, nx, ny, nz, pdzdx )
    dpdzydx .= grad1( 1, dx, nx, ny, nz, pdzdy )
    dzdpxdy .= grad1( 2, dy, nx, ny, nz, zdpdx )
    dzdpydx .= grad1( 1, dx, nx, ny, nz, zdpdy )

    Δzeta .= laplace1( 1, dx, nx, ny, nz, zeta ) .+ laplace1( 2, dy, nx, ny, nz, zeta )
    Δ2zeta .= laplace1( 1, dx, nx, ny, nz, Δzeta ) .+ laplace1( 2, dy, nx, ny, nz, Δzeta )
    Δ3zeta .= laplace1( 1, dx, nx, ny, nz, Δ2zeta ) .+ laplace1( 2, dy, nx, ny, nz, Δ2zeta )

    adv .= - (
              (dpdx .* dzdy .- dpdy .* dzdx)  # J1
              .+(dpdzydx .- dpdzxdy)  # J2
              .+(dzdpxdy .- dzdpydx)  # J3
             ) ./ 3.0

    diff .= nu .* Δ3zeta  # Diffusion

    dzeta_res .= adv .+ diff

    return dzeta_res
end

function integrate_RK4( nu, ntr, emax, dx, dy, dz, nx, ny, nz, zeta, psi, dt )  # integrating each prognostic equation with the RK4
    # nu: diffusion coefficient
    # ntr: the iteration times
    # emax: maximum difference for the iteration of zeta to psi
    # nx,ny,nz: x-y-z grid numbers
    # dx,dy,dz: x-y-z grid intervals
    # zeta(nx,ny,nz): vorticity
    # psi(nx,ny,nz): streamfunction
    # dt: time step for integration

    zeta_res = fill(0.0,nx,ny,nz)
    psi_res = fill(0.0,nx,ny,nz)
    zeta_work = fill(0.0,nx,ny,nz)
    psi_work = fill(0.0,nx,ny,nz)
    df1 = fill(0.0,nx,ny,nz)
    df2 = fill(0.0,nx,ny,nz)
    df3 = fill(0.0,nx,ny,nz)
    df4 = fill(0.0,nx,ny,nz)

    zeta_work .= zeta
    psi_work .= psi
    for i in 1:nz
        psi_work[1:nx,1:ny,i] .= zeta_to_psi2( ntr, emax, dx, dy, nx, ny, zeta_work[1:nx,1:ny,i], psi_work[1:nx,1:ny,i] )
    end
    df1 .= zeta_tend( nu, dx, dy, dz, nx, ny, nz, zeta_work, psi_work )

    zeta_work .= zeta + 0.5 * df1
    for i in 1:nz
        psi_work[1:nx,1:ny,i] .= zeta_to_psi2( ntr, emax, dx, dy, nx, ny, zeta_work[1:nx,1:ny,i], psi_work[1:nx,1:ny,i] )
    end
    df2 .= zeta_tend( nu, dx, dy, dz, nx, ny, nz, zeta_work, psi_work )

    zeta_work .= zeta + 0.5 * df2
    for i in 1:nz
        psi_work[1:nx,1:ny,i] .= zeta_to_psi2( ntr, emax, dx, dy, nx, ny, zeta_work[1:nx,1:ny,i], psi_work[1:nx,1:ny,i] )
    end
    df3 .= zeta_tend( nu, dx, dy, dz, nx, ny, nz, zeta_work, psi_work )

    zeta_work .= zeta + df3
    for i in 1:nz
        psi_work[1:nx,1:ny,i] .= zeta_to_psi2( ntr, emax, dx, dy, nx, ny, zeta_work[1:nx,1:ny,i], psi_work[1:nx,1:ny,i] )
    end
    df4 .= zeta_tend( nu, dx, dy, dz, nx, ny, nz, zeta_work, psi_work )

    psi_res .= psi_work
    zeta_res .= zeta + (dt/6.0) * (df1 + 2.0 * df2 + 2.0 * df3 + df4)
    
    return zeta_res, psi_res
end


end


