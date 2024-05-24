


#= Adam Dama Python Code 
        with open(path, "w", newline="", encoding="utf-8") as f:
            # Header
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["Version            :", 6])
            writer.writerow(
                ["Plate type name    :", plate_type.formulatrix_plate_template_name]
            )
            writer.writerow(["Priority Delays    :", 0])

            # Pipette volumes
            for idx, r in enumerate(reagents):
                writer.writerow(["Reagent Name    :", r])
                writer.writerow(["Barcode            :"])
                writer.writerow(["Priority           :", 1])
                vols = worklist.iloc[:, idx].values
                vols = vols.reshape(plate_type.shape)
                writer.writerows(vols)
=# 
"""
    tempest(design::DataFrame,filepath::String,destination::ContainerName)

Create Tempest dipsense instructions for SBS plate dispensing

  ## Arguments 
  * `design`: a (# of runs) x (# of stocks) dataframe containing the volume of each stock for each run in µL.
  * `filepath`: the ouput path of the dispense file. The function automatically strips any incorrect file extension and adds the appropriate one.
  * `destination`: The destination plate type. See `keys(containers)` for available options. 
"""
function tempest(design::DataFrame,filepath::String,destination::ContainerName)
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
    platefilename=destination.tempest
    if platefilename==""
        error(ArgumentError("Selected plate type $(destination.container.name) is not supported by Tempest."))
    end 
    print(outfile,join(["Version            :", 6],'\t'),"\r\n")
    print(outfile,join(["Plate type name    :", platefilename],'\t'),"\r\n")
    print(outfile,join(["Priority Delays    :", 0],'\t'),"\r\n")


    for i in 1:n_stocks 
        vols=Vector(design[:,i])
        vols=reshape(vols,R,C)
        print(outfile,join(["Reagent Name    :", colnames[i]],'\t'),"\r\n")
        print(outfile,join(["Barcode            :"],'\t'),"\r\n")
        print(outfile,join(["Priority           :", 1],'\t'),"\r\n")
        for r in 1:R
            print(outfile,join(vols[r,:],'\t'),"\r\n")
        end 
    end 
    close(outfile)
end 

"""
    multi_tempest(designs::Vector{DataFrame},dlnames::Vector{String},directory::String,mdlname::String,destination::ContainerName)

Create Tempest multidipsense instructions for SBS plate dispensing

  ## Arguments 
  * `design`: a vector of k dataframes where each dataframe is a (# of runs) x (# of stocks) dataframe containing the volume of each stock for each run in µL.
  * `dlnames`" The name of each dispense list. There should be k names, one for each design.
  * `directory`: the ouput directory of the dispense file. The function automatically strips any incorrect file extension and adds the appropriate one.
  * `mdlname`: The name of the multidispense list file. The funcion outputs the multidispense list in directory
  * `destination`: The destination plate type. See `keys(containers)` for available options. 
"""
function multi_tempest(designs::Vector{DataFrame},dlnames::Vector{String},directory::String,mdlname::String,destination::ContainerName)
    n=length(designs)

    filepaths=joinpath.((directory,),dlnames)

    ext=".dl.txt"
    basenames=string.(dlnames,(ext,))

    mdlbasename=string(mdlname,".mdl.txt")
    outpath=joinpath(directory,mdlbasename)

    TempestDispense.(designs,filepaths,(destination,))

    outfile=open(outpath,"w")

    print(outfile,join(["LP:$(destination.tempest).pd.txt"],'\t'),"\r\n")
    for i in 1:n
        print(outfile,join(["P:$i","TP:1","DL:$(basenames[i])"],'\t'),"\r\n")
    end
    close(outfile)
end 

