using SparseArrays: SparseMatrixCSC, sparse
using ArnoldiMethod: SR
using Base: OneTo
using LinearAlgebra: eigen

function fmrg_layout(g::AbstractGraph,
                       locs_x::Array{Float64,1},
                       locs_y::Array{Float64,1},
                       MAXITER=500;
                       C=2.0,
                       INITTEMP=2.0)

    nvg = nv(g)
    adj_matrix = adjacency_matrix(g)
    work_x = locs_x
    work_y = locs_y
    # The optimal distance bewteen vertices
    k = C * sqrt(4.0 / nvg)
    k² = k * k

    # Store forces and apply at end of iteration all at once
    force_x = zeros(nvg)
    force_y = zeros(nvg)

    # Iterate MAXITER times
    @inbounds for iter = 1:MAXITER
        # Calculate forces
        for i = 1:nvg
            force_vec_x = 0.0
            force_vec_y = 0.0
            for j = 1:nvg
                i == j && continue
                d_x = locs_x[j] - locs_x[i]
                d_y = locs_y[j] - locs_y[i]
                dist²  = (d_x * d_x) + (d_y * d_y)
                dist = sqrt(dist²)

                if !( iszero(adj_matrix[i,j]) && iszero(adj_matrix[j,i]) )
                    # Attractive + repulsive force
                    F_d = dist² / k - k² / dist # original FR algorithm
                    # F_d = dist / k - k² / dist²
                else
                    # Just repulsive
                    F_d = -k² / dist  # original FR algorithm
                    # F_d = -k² / dist²
                end
                force_vec_x += F_d*d_x
                force_vec_y += F_d*d_y
            end
            force_x[i] = force_vec_x
            force_y[i] = force_vec_y
        end
        # Cool down
        temp = INITTEMP / iter
        # Now apply them, but limit to temperature
        for i = 1:nvg
            fx = force_x[i]
            fy = force_y[i]
            force_mag  = sqrt((fx * fx) + (fy * fy))
            scale      = min(force_mag, temp) / force_mag
            work_x[i] += force_x[i] * scale
            work_y[i] += force_y[i] * scale
        end
    end

    # Scale to unit square
    min_x, max_x = minimum(work_x), maximum(work_x)
    min_y, max_y = minimum(work_y), maximum(work_y)
    function scaler(z, a, b)
        2.0*((z - a)/(b - a)) - 1.0
    end
    map!(z -> scaler(z, min_x, max_x), work_x, work_x)
    map!(z -> scaler(z, min_y, max_y), work_y, work_y)

    return work_x, work_y
end
