"""
    mantis(design::DataFrame,filepath::String,destination::ContainerName)

Create Mantis dipsense instructions for SBS plate dispensing

  ## Arguments 
  * `design`: a (# of runs) x (# of stocks) dataframe containing the volume of each stock for each run in µL.
  * `filepath`: the ouput path of the dispense file. The function automatically strips any incorrect file extension and adds the appropriate one.
  * `destination`: The destination plate type. See `keys(containers)` for available options. 
"""
function mantis(design::DataFrame,filepath::String,destination::ContainerName)
        n=nrow(design)
        R,C=destination.container.shape
        colnames=names(design)
        if n != R*C
            error(ArgumentError("Number of wells in the design ($n) does not match the size of the destination plate $(R*C). Fill the design with 0's to designate empty wells."))
        end
        directory=dirname(filepath)
        expt_name=basename(filepath)
        expt_name,ext=splitext(expt_name)
        ext=".dl.txt"
        name=string(expt_name,ext)
        if ~isdir(directory)
            mkdir(directory)
          end 
        filename=joinpath(directory,name)
        n_stocks=ncol(design)
        delay_header=vcat(n_stocks,repeat([0,""],n_stocks))
        outfile=open(filename,"w")
        platefilename=destination.mantis
        if platefilename==""
            error(ArgumentError("Selected plate type $(destination.container.name) is not supported by Mantis."))
        end 
        print(outfile,join(["[ Version: 5 ]"],'\t'),"\r\n")
        print(outfile,platefilename,"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")
        print(outfile,join([1],'\t'),"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")

        for i in 1:n_stocks 
            vols=Vector(design[:,i])
            vols=round.(reshape(vols,R,C),digits=1)
            print(outfile,join([colnames[i],"","Normal"],'\t'),"\r\n")
            print(outfile,join(["Well",1],'\t'),"\r\n")
            for r in 1:R
                print(outfile,join(vols[r,:],'\t'),"\r\n")
            end 
        end 
        close(outfile)
end 


