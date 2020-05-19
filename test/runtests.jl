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
   
    @time eigval,X,num_ev = ZPares.eigensolve(sparse(A),sparse(B),left,right)
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i])
    end    
    @time eigval,X,num_ev = ZPares.eigensolve(sparse(A),sparse(B),left,right,ishermitian=true)
    for i=1:num_ev
        println(i,"\t",eigval[i],"\t",res[i])
    end 
    #e,v = eigen(A)
    #println(e)

    
end



@testset "ZPares.jl" begin
    # Write your own tests here.
    test()

end
