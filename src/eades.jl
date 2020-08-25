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

           c1 = 1
           c2 = 100
           c3 = 2
           c4 = 0.1

           # Iterate MAXITER times
           @inbounds for iter = 1:MAXITER
               # Calculate forces
               for i = 1:nvg
                   force_vec_x = 0.0
                   force_vec_y = 0.0
                   for j = 1:nvg
                       i == j && continue
                       d_x = work_x[j] - work_x[i]
                       d_y = work_y[j] - work_y[i]
                       if !( iszero(adj_matrix[i,j]) && iszero(adj_matrix[j,i]) )
                           # F_x = abs(c1 * log(abs(d_x) / c2))
                           # F_y = abs(c1 * log(abs(d_y) / c2))
                           F_x = abs( c1 * d_x)
                           F_y = abs( c1 * d_y)

                       else
                           F_x = c3 / (d_x * d_x)
                           F_y = c3 / (d_y * d_y)

                           # F_x = abs(c3 / log(abs(d_x * d_x)))
                           # F_y = abs(c3 / log(abs(d_y * d_y)))
                       end
                       if (work_x[i] < work_x[j])
                           force_vec_x += F_x
                       else
                           force_vec_x -= F_x
                       end
                       if (work_y[i] < work_y[j])
                           force_vec_y += F_y
                       else
                           force_vec_y -= F_y
                       end
                   end
                   force_x[i] += c4 * force_vec_x
                   force_y[i] += c4 * force_vec_y
               end
               for i = 1:nvg
                   work_x[i] += force_x[i]
                   work_y[i] += force_y[i]
               end
               # Cool down
           end


           # Scale to unit square
           min_x, max_x = minimum(work_x), maximum(work_x)
           min_y, max_y = minimum(work_y), maximum(work_y)
           function scaler(z, a, b)
               3.0*((z - a)/(b - a)) - 1.0
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
