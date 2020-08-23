using SparseArrays: SparseMatrixCSC, sparse
using ArnoldiMethod: SR
using Base: OneTo
using LinearAlgebra: eigen

function kamkaw_layout(g::AbstractGraph,
                       locs_x::Array{Float64,1},
                       locs_y::Array{Float64,1},
                       MAXITER=500;
                       C=2.0,
                       INITTEMP=2.0)
       nvg = nv(g)
       adj_matrix = adjacency_matrix(g)

       work_x = locs_x
       work_y = locs_y

       d = [(i,j) for i=1:nvg for j=1:nvg]


       iter = 1
       temp = INITTEMP / iter

       for i = 1:nvg
           fx = force_x[i]
           fy = force_y[i]
           force_mag  = sqrt((fx * fx) + (fy * fy))
           scale      = min(force_mag, temp) / force_mag
           work_x[i] += force_x[i] * scale
           work_y[i] += force_y[i] * scale
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
