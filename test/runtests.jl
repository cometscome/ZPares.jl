using ZPares
using Test
using LinearAlgebra
using SparseArrays

function test()
    eps = 1e-6

    mat_size = 500
    A = spzeros(ComplexF64,mat_size,mat_size)
    B = spzeros(ComplexF64,mat_size,mat_size)

    for i=1:mat_size
        j=i+1
        if 1 <= j <= mat_size
            A[i,j] = -1
        end

        j=i-1
        if 1 <= j <= mat_size
            A[i,j] = -1
        end

    end
    #=
    for i=1:mat_size
        for j=1:mat_size
            if i==j
                A[i,j] = 2
            elseif abs(i-j) == 1
                A[i,j] = 1
            end 
        end
    end
    =#

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

    @time  e,v = eigen(Matrix(A))
    println("exact results")
    println("index:   eigenvalues")
    exacte = []
    for i=1:length(e)
        if γ - ρ <= e[i] <= γ +  ρ
            println(i,"\t",e[i])
            push!(exacte,e[i])
        end
    end



    @time eigval,X,num_ev,res = ZPares.eigensolve(A,left,right)
    println("index:   eigenvalues : residuals : exact diff")
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i],"\t",eigval[i]-exacte[i])
    end

    @test sum(abs.(eigval .- exacte))/sum(abs.(exacte)) < eps


    @time eigval,X,num_ev,res = ZPares.eigensolve(A,left,right,ishermitian=true)



    println("index:   eigenvalues : residuals : exact diff")
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i],"\t",eigval[i]-exacte[i])
    end

    @test sum(abs.(eigval .- exacte))/sum(abs.(exacte)) < eps

    #exit()
   
    @time eigval,X,num_ev,res = ZPares.eigensolve(A,B,left,right)
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i],"\t",eigval[i]-exacte[i])
    end    

    @test sum(abs.(eigval .- exacte))/sum(abs.(exacte)) < eps


    @time eigval,X,num_ev,res = ZPares.eigensolve(sparse(A),sparse(B),left,right,ishermitian=true)
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i],"\t",eigval[i]-exacte[i])
    end 
    #e,v = eigen(A)
    #println(e)

    @test sum(abs.(eigval .- exacte))/sum(abs.(exacte)) < eps


    
end



@testset "ZPares.jl" begin
    # Write your own tests here.
    test()

end
