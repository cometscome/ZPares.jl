module ZPares
    using IterativeSolvers
    using LinearAlgebra
    using LinearMaps

    const currentdir = pwd()
    const insdir = homedir()*"/.julia/packages/ZPares"


    function eigensolve(A,emin,emax;L=8,N=32,M=16,Lmax=32,ishermitian=false)
        function mulC(z)
            function mulCy!(y::AbstractVector, x::AbstractVector)
                @. y[:] = x[:]
                mul!(y, A, x, -1, z) #A B α + C β y = (z - H)x
                return y
            end
            (y,x) -> mulCy!(y,x)
        end 
        function mulCdag(z)
            function mulCy!(y::AbstractVector, x::AbstractVector)
                @. y[:] = x[:]
                if ishermitian
                    mul!(y, A, x, -1, -z) #A B α + C β y = (z - H)x
                else
                    mul!(y, A', x, -1, -z) #A B α + C β y = (z - H)x
                end
                return y
            end
            (y,x) -> mulCy!(y,x)
        end


        if ishermitian == false
            typeforeig = ComplexF64
            emin += 0im
        else
            typeforeig = Float64
        end


        #initialize_zpares_prm(L,N,M,Lmax)
        cd(insdir)
        ccall((:initialize_zpares_prm,"zpares_wrapper.so"),Nothing, 
            (
            Ref{Int64}, #L
            Ref{Int64}, #N
            Ref{Int64}, #M
            Ref{Int64}), #LMAX
            L,N,M,Lmax)
        

        tasks = zeros(Int32,7)

        #subroutine get_ZPARES_TASK(tasks) bind(c,name='get_ZPARES_TASK')
        ccall((:get_ZPARES_TASK,"zpares_wrapper.so"),Nothing, 
            (
            Ref{Int32},), #tasks
            tasks)

        ZPARES_TASK_FINISH = tasks[1]
        ZPARES_TASK_FACTO = tasks[2]
        ZPARES_TASK_SOLVE = tasks[3]
        ZPARES_TASK_SOLVE_H = tasks[4]
        ZPARES_TASK_MULT_A = tasks[5]
        ZPARES_TASK_NONE = tasks[6]
        ZPARES_TASK_MULT_B = tasks[7]

        #println(tasks)

        mat_size = size(A)[1]
        
        #rwork(mat_size, Lmax), cwork(mat_size, Lmax)
        rwork = zeros(ComplexF64,mat_size,Lmax)
        cwork = zeros(ComplexF64,mat_size,Lmax)

        #integer function zpares_get_ncv_wrapper() bind(c, name = 'zpares_get_ncv_wrapper')
        ncv = ccall((:zpares_get_ncv_wrapper,"zpares_wrapper.so"),Int64, 
            ()
            )
        eigval = zeros(typeforeig,ncv)
        res = zeros(Float64,ncv)
        X = zeros(ComplexF64,mat_size,ncv)
        z = zeros(ComplexF64,1)
        num_ev =zeros(Int32,1)
        info =zeros(Int32,1)
        itask =zeros(Int32,1)
        C = similar(A)
        #C = zeros(ComplexF64,size(A)[1],size(A)[2])
        xs =zeros(Int32,1)
        ws =zeros(Int32,1)
        nc =zeros(Int32,1)

        itask[1] = ZPARES_TASK_NONE
        while itask[1] != ZPARES_TASK_FINISH

            if ishermitian == true
                            #=
                    subroutine wrapper_zpares_zrciheev &
                    (ncv,Lmax,mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info) bind(c,name = 'wrapper_zpares_zrciheev')
                    =#
                ccall((:wrapper_zpares_zrciheev,"zpares_wrapper.so"),Nothing, 
                (Ref{Int64},#ncv
                Ref{Int64},#Lmax
                Ref{Int64},#mat_size
                Ref{ComplexF64}, #z
                Ref{ComplexF64}, #rwork
                Ref{ComplexF64}, #cwork
                Ref{Float64}, #emin
                Ref{Float64}, #emax
                Ref{Int32}, #num_ev
                Ref{Float64}, #eigval
                Ref{ComplexF64}, #X
                Ref{Float64},#res
                Ref{Int32},#itask
                Ref{Int32},#xs
                Ref{Int32},#ws
                Ref{Int32},#nc
                Ref{Int32}), #info
                ncv,Lmax,mat_size, z, rwork, cwork, emin, 
                    emax, num_ev, eigval, X, res, itask,xs,ws,nc,info)
            else
                ccall((:wrapper_zpares_zrcigeev,"zpares_wrapper.so"),Nothing, 
                (Ref{Int64},#ncv
                Ref{Int64},#Lmax
                Ref{Int64},#mat_size
                Ref{ComplexF64}, #z
                Ref{ComplexF64}, #rwork
                Ref{ComplexF64}, #cwork
                Ref{ComplexF64}, #emin
                Ref{Float64}, #emax
                Ref{Int32}, #num_ev
                Ref{ComplexF64}, #eigval
                Ref{ComplexF64}, #X
                Ref{Float64},#res
                Ref{Int32},#itask
                Ref{Int32},#xs
                Ref{Int32},#ws
                Ref{Int32},#nc
                Ref{Int32}), #info
                ncv,Lmax,mat_size, z, rwork, cwork, emin, 
                    emax, num_ev, eigval, X, res, itask,xs,ws,nc,info)
            end
            #println(z)
            if itask[1] ==ZPARES_TASK_FACTO
                #println("ZPARES_TASK_FACTO")
                #=
                @. C[:,:] = -A[:,:]
                for i=1:mat_size
                    C[i,i] += z[1]
                end
                =#
                
            elseif   itask[1] == ZPARES_TASK_SOLVE
                C = LinearMap(mulC(z[1]) , mat_size; ismutating=true)

                #println("ZPARES_TASK_SOLVE")
                for jj=ws[1]:ws[1]+nc[1]-1
                    cwork[:,jj] = gmres(C,cwork[:,jj],tol = 1e-16) 
                end
            elseif itask[1] == ZPARES_TASK_SOLVE_H
                C = LinearMap(mulCdag(z[1]) , mat_size; ismutating=true)
                #println("ZPARES_TASK_SOLVE_H")
                for jj=ws[1]:ws[1]+nc[1]-1
                    cwork[:,jj] = gmres(C,cwork[:,jj],tol = 1e-16) 
#                    cwork[:,jj] = gmres(C',cwork[:,jj],tol = 1e-16) 
                end
            elseif itask[1]  ==ZPARES_TASK_MULT_A
                rwork[:,ws[1]:ws[1]+nc[1]-1] = A*X[:,xs[1]:xs[1]+nc[1]-1]

            end
            
        end

        for i=1:num_ev[1]
            println(i,"\t",eigval[i],"\t",res[i])
        end
#        println(eigval[1:num_ev[1]])
        cd(currentdir)

        return eigval[1:num_ev[1]],X[:,1:num_ev[1]],num_ev[1]

    end

    function eigensolve(A,B,left,right;L=8,N=32,M=16,Lmax=32,ishermitian=false)

        function mulC(z)
            function mulCy!(y::AbstractVector, x::AbstractVector)
                mul!(y,B,x)
#                @. y[:] = x[:]
                mul!(y, A, x, -1, z) #A B α + C β y = (zB - H)x
                return y
            end
            (y,x) -> mulCy!(y,x)
        end 

        function mulCdag(z)
            function mulCy!(y::AbstractVector, x::AbstractVector)
                mul!(y,B',x)
#                @. y[:] = x[:]
                if ishermitian
                    mul!(y, A, x, -1, -z) #A B α + C β y = (zB - H)x
                else
                    mul!(y, A', x, -1, -z) #A B α + C β y = (zB - H)x
                end
                return y
            end
            (y,x) -> mulCy!(y,x)
        end





        cd(insdir)
        #initialize_zpares_prm(L,N,M,Lmax)
        ccall((:initialize_zpares_prm,"zpares_wrapper.so"),Nothing, 
            (
            Ref{Int64}, #L
            Ref{Int64}, #N
            Ref{Int64}, #M
            Ref{Int64}), #LMAX
            L,N,M,Lmax)

        tasks = zeros(Int32,7)

        #subroutine get_ZPARES_TASK(tasks) bind(c,name='get_ZPARES_TASK')
        ccall((:get_ZPARES_TASK,"zpares_wrapper.so"),Nothing, 
            (
            Ref{Int32},), #tasks
            tasks)

        ZPARES_TASK_FINISH = tasks[1]
        ZPARES_TASK_FACTO = tasks[2]
        ZPARES_TASK_SOLVE = tasks[3]
        ZPARES_TASK_SOLVE_H = tasks[4]
        ZPARES_TASK_MULT_A = tasks[5]
        ZPARES_TASK_NONE = tasks[6]
        ZPARES_TASK_MULT_B = tasks[7]

        #println(tasks)

        mat_size = size(A)[1]
        
        #rwork(mat_size, Lmax), cwork(mat_size, Lmax)
        rwork = zeros(ComplexF64,mat_size,Lmax)
        cwork = zeros(ComplexF64,mat_size,Lmax)

        #integer function zpares_get_ncv_wrapper() bind(c, name = 'zpares_get_ncv_wrapper')
        ncv = ccall((:zpares_get_ncv_wrapper,"zpares_wrapper.so"),Int64, 
            ()
            )
        eigval = zeros(ComplexF64,ncv)
        res = zeros(Float64,ncv)
        X = zeros(ComplexF64,mat_size,ncv)
        z = zeros(ComplexF64,1)
        num_ev =zeros(Int32,1)
        info =zeros(Int32,1)
        itask =zeros(Int32,1)
        C = similar(A)
        #C = zeros(ComplexF64,size(A)[1],size(A)[2])
        xs =zeros(Int32,1)
        ws =zeros(Int32,1)
        nc =zeros(Int32,1)

        itask[1] = ZPARES_TASK_NONE
        while itask[1] != ZPARES_TASK_FINISH
            #=
            subroutine wrapper_zpares_zrciheev &
            (ncv,Lmax,mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info) bind(c,name = 'wrapper_zpares_zrciheev')
            =#
            ccall((:wrapper_zpares_zrcigegv,"zpares_wrapper.so"),Nothing, 
            (Ref{Int64},#ncv
            Ref{Int64},#Lmax
            Ref{Int64},#mat_size
            Ref{ComplexF64}, #z
            Ref{ComplexF64}, #rwork
            Ref{ComplexF64}, #cwork
            Ref{ComplexF64}, #emin
            Ref{ComplexF64}, #emax
            Ref{Int32}, #num_ev
            Ref{ComplexF64}, #eigval
            Ref{ComplexF64}, #X
            Ref{Float64},#res
            Ref{Int32},#itask
            Ref{Int32},#xs
            Ref{Int32},#ws
            Ref{Int32},#nc
            Ref{Int32}), #info
            ncv,Lmax,mat_size, z, rwork, cwork, left, 
                right, num_ev, eigval, X, res, itask,xs,ws,nc,info)
            #println(z)
            if itask[1] ==ZPARES_TASK_FACTO
                #println("ZPARES_TASK_FACTO")
                #@. C[:,:] = -A[:,:] + z[1]*B[:,:]
#                for i=1:mat_size
#                    C[i,i] += z[1]
#                end
            elseif   itask[1] == ZPARES_TASK_SOLVE
                C = LinearMap(mulC(z[1]) , mat_size; ismutating=true)

                #println("ZPARES_TASK_SOLVE")
                for jj=ws[1]:ws[1]+nc[1]-1
                    #x = C \ cwork[:,jj]
                    cwork[:,jj] = gmres(C,cwork[:,jj],tol = 1e-16) 
                    #x = bicgstabl(C, cwork[:,jj] , 2)
                    #cwork[:,jj] =  x[:] #C \ cwork[:,jj] #bicgstabl(C, cwork[:,jj] , 2)
                end
            elseif itask[1] == ZPARES_TASK_SOLVE_H
                C = LinearMap(mulCdag(z[1]) , mat_size; ismutating=true)
                #println("ZPARES_TASK_SOLVE_H")
                for jj=ws[1]:ws[1]+nc[1]-1
                    #x = C' \ cwork[:,jj]
                    cwork[:,jj] = gmres(C',cwork[:,jj],tol = 1e-16) 
#                    x = bicgstabl(C', cwork[:,jj] , 2)
                    #println(norm(x-x2))
                    #cwork[:,jj] =  x[:]#C' \ cwork[:,jj]#bicgstabl(C', cwork[:,jj] , 2)
                end
            elseif itask[1]  ==ZPARES_TASK_MULT_A
                rwork[:,ws[1]:ws[1]+nc[1]-1] = A*X[:,xs[1]:xs[1]+nc[1]-1]
            elseif itask[1]  ==ZPARES_TASK_MULT_B
                rwork[:,ws[1]:ws[1]+nc[1]-1] = B*X[:,xs[1]:xs[1]+nc[1]-1]
            end
            
        end

        for i=1:num_ev[1]
            println(i,"\t",eigval[i],"\t",res[i])
        end
#        println(eigval[1:num_ev[1]])
        cd(currentdir)

        return eigval[1:num_ev[1]],X[:,1:num_ev[1]],num_ev[1]

    end

end # module
