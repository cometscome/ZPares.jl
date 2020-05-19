const ZParesversion = "0.9.6a"
const ZParespath = "libzpares.a"

const ZParestar = "zpares_"*ZParesversion*".tar.gz"
const ZParesurl = "https://zpares.cs.tsukuba.ac.jp/?download=242"

if isfile(ZParespath) == false
    if isfile("zpares_"*ZParesversion) == false
        if isfile(ZParestar) == false
            println("Downloading zpares...")
            download( ZParesurl,ZParestar)
            println("done.")
        end
        run(`tar zxf  $ZParestar`) 
    end
    cd("zpares_"*ZParesversion)
    run(`cp Makefile.inc/make.inc.gfortran.seq make.inc`)
    run(`ls`)
    #run(`pwd`)
    run(`make`) 
    run(`cp ./lib/libzpares.a ../`)
    run(`cp ./src/zpares.mod ../`)
    cd("../")

#        run(`cp ./lib/libzpares.a ../../src/`)
#        run(`cp ./src/zpares.mod ../../src/`)
    #cd(pacpass)


end

if isfile("zpares_wrapper.so") == false
    #cd(pacpass*"/src")
    #run(`ls`)
    run(`pwd`)
    run(`gfortran -L./ -lzpares zpares_wrapper.f90 -o zpares_wrapper.so -shared -fPIC -llapack -lblas`)
    insdir = homedir()*"/.julia/packages/ZPares"
    run(`cp zpares_wrapper.so $insdir/src/`)
    run(`cp libzpares.a $insdir/src/`)
    run(`cp zpares.mod $insdir/src/`)
    run(`cp zpares_wrapper.so $(insdir)/test/`)
    run(`cp libzpares.a $(insdir)/test/`)
    run(`cp zpares.mod $(insdir)/test/`)
end