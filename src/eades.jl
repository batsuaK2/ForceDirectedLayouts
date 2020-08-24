using SparseArrays: SparseMatrixCSC, sparse
using ArnoldiMethod: SR
using Base: OneTo
using LinearAlgebra: eigen


function eades_layout(g::AbstractGraph,
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

           c1 = 2
           c2 = 1
           c3 = 1
           c4 = 0.1

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
                           F_x = -1 * c1 * (abs(d_x) / c2)
                           F_y = -1 * c1 * (abs(d_y) / c2)
                       else
                           F_x = c3 / (d_x * d_x)
                           F_y = c3 / (d_y * d_y)
                       end
                       force_vec_x += F_x
                       force_vec_y += F_y
                   end
                   force_x[i] += c4 * force_vec_x
                   force_y[i] += c4 * force_vec_y
               end
               # Cool down
           end

           for i = 1:nvg
               work_x[i] += force_x[i]
               work_y[i] += force_y[i]
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


using Random: MersenneTwister

function eades_layout(g::AbstractGraph, seed::Integer, kws...)
    rng = MersenneTwister(seed)
    eades_layout(g, 2 .* rand(rng, nv(g)) .- 1.0, 2 .* rand(rng,nv(g)) .- 1.0; kws...)
end
