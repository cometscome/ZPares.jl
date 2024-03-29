
# ZPares
[![CI](https://github.com/cometscome/ZPares.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/cometscome/ZPares.jl/actions/workflows/CI.yml)

Julia wrapper for [z-Pares](https://zpares.cs.tsukuba.ac.jp).

# What is this?

We can get eigenvalues with a given domain Gamma. Now Gamma is a elipsoidal shape. 
```left``` means the left side and ```right``` means the right side. 
If you use ```ishermitian = true```, left is the minimum eigenvalue in the domain Gamma and right is the maximum in this. 

Julia 1.6 or higher is needed. 

```
z-Pares: Parallel Eigenvalue Solver

z-Pares is a package for solving generalized eigenvalue problems. z-Pares is designed to compute a few eigenvalues and eigenvectors of sparse matrices. The symmetries and definitenesses of the matrices can be exploited suitably. z-Pares implements a complex moment based contour integral eigensolver. z-Pares computes eigenvalues inside a user-specified contour path and corresponding eigenvectors. The most important feature of z-Pares is two-level Message Passing Interface (MPI) distributed parallelism.
```
# How to install

```
add https://github.com/cometscome/ZPares.jl
```

# Sample

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
