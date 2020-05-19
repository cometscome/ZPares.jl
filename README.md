# ZPares

[![Build Status](https://travis-ci.com/cometscome/ZPares.jl.svg?branch=master)](https://travis-ci.com/cometscome/ZPares.jl)
[![Codecov](https://codecov.io/gh/cometscome/ZPares.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cometscome/ZPares.jl)
[![Coveralls](https://coveralls.io/repos/github/cometscome/ZPares.jl/badge.svg?branch=master)](https://coveralls.io/github/cometscome/ZPares.jl?branch=master)

Julia wrapper for [z-Pares](https://zpares.cs.tsukuba.ac.jp).

I confirmed that this works on MacOS 10.14.6.

Julia 1.4 or higher is needed. 

```
z-Pares: Parallel Eigenvalue Solver

z-Pares is a package for solving generalized eigenvalue problems. z-Pares is designed to compute a few eigenvalues and eigenvectors of sparse matrices. The symmetries and definitenesses of the matrices can be exploited suitably. z-Pares implements a complex moment based contour integral eigensolver. z-Pares computes eigenvalues inside a user-specified contour path and corresponding eigenvectors. The most important feature of z-Pares is two-level Message Passing Interface (MPI) distributed parallelism.
```

Sample

```julia
using ZPares
using Test
using SparseArrays

function test()

    mat_size = 500
    A = zeros(ComplexF64,mat_size,mat_size)
    B = zeros(ComplexF64,mat_size,mat_size)
    for i=1:mat_size
        for j=1:mat_size
            if i==j
                A[i,j] = 2
            elseif abs(i-j) == 1
                A[i,j] = 1
            end 
        end
    end

    for i=1:mat_size
        B[i,i] = 1
    end

#    @time ZPares.eigsolve(A,B)
    left = 3.5
    right = 3.6
    ρ = (right-left)/2
    γ = (right+left)/2
    L=8
    N=32
    M=16
    Lmax=32
    γ = 0
    ρ = 0.2
    left = γ-ρ
    right = γ+ρ
    @time eigval,X,num_ev,res = ZPares.eigensolve(sparse(A),left,right)
    println("index:   eigenvalues : residuals")
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i])
    end


    @time eigval,X,num_ev,res = ZPares.eigensolve(sparse(A),left,right,ishermitian=true)

    println("index:   eigenvalues : residuals")
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i])
    end
    #exit()
   
    @time eigval,X,num_ev,res = ZPares.eigensolve(sparse(A),sparse(B),left,right)
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i])
    end    
    @time eigval,X,num_ev,res = ZPares.eigensolve(sparse(A),sparse(B),left,right,ishermitian=true)
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i])
    end 
    #e,v = eigen(A)
    #println(e)

    
end

test()

```
