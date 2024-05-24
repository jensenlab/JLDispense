# JLDispense.jl
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jensenlab/JensenLabDispense/blob/main/LICENSE)

 

# Contents 
[Description](#description) \
[Installing JLDispense](#installing-jldispense) \
[Example Usage](#example-usage)


# Description 
Drivers for Jensen Lab Liquid Handlers 

# Installing JLDispense
 Requires installation of  [Julia](https://julialang.org/downloads/). Once Julia is installed. Install JensenLabDispense by navigating to package mode:  

```julia 
add https://github.com/jensenlab/JLDispense
```

# Example Usage 
```julia 
    using JLDispense, DataFrames

```
## Cobra  
 ``` cobra``` generates a directory populated with the appropriate dispense files to control the ARI/Hudson Cobra liquid handler. use ```?cobra``` to see documentation

```julia
    design = DataFrame(ones(384,96),:auto) # dispense 1 ul from each well of a 96 wp into each well of a 384 wp
    lq_classes= ["Water" for _ in 1:96] # each source has a liquid class. 

    cobra(design, "<PATH TO DIRECTORY>", containers[:WP96], containers[:WP384],lq_classes)
```


## Mantis 
```mantis``` generates a .dl.txt file for the Fomulatrix Mantis liquid handler. use ```?mantis``` for documentation
```julia
    design = DataFrame(ones(96,6),:auto) # dispense 1 ul from 6 sources to each well of a 96 wp

    mantis(design, "<PATH TO DIRECTORY>", containers[:WP96]) 
```

## Tempest 
```tempest``` generates a .dl.txt file for the Fomulatrix Mantis liquid handler. use ```?tempest``` for documentation
```julia
    design = DataFrame(ones(96,6),:auto) # dispense 1 ul from 6 sources to each well of a 96 wp

    tempest(design, "<PATH TO DIRECTORY>", containers[:WP96]) 
```

```multi_tempest``` generates a .mdl.txt and associated .dl.txt files to run multiple unique dispense lists using the plate changer on the Formulatrix Tempest liquid handler. use ```?multi_tempest``` for documentation. 
```julia
    design1 = DataFrame(ones(96,6),:auto) # dispense 1 ul from 6 sources to each well of a 96 wp
    design2 = DataFrame(2*ones(96,6),:auto) # dispense 2 ul from 6 sources to each well of a 96 wp
    design3 = DataFrame(3*ones(96,6),:auto) # dispense 3 ul from 6 sources to each well of a 96 wp
    designs=[design1,design2,design3]
    names=["design1","design2","design3"]
    multi_tempest(designs,names, "<PATH TO DIRECTORY>","multi_design", containers[:WP96]) 
```


## Nimbus

`nimbus` generates a .csv file that the Hamilton Nimbus Software parses to control the machine. Use `?nimbus` for documentation. 

```julia 
    design = DataFrame(ones(96,6),:auto) # dispense 1 ul from 6 sources to each well of a 96 wp
    disp_list=NimbusDispenseList(design,"TubeRack50ML_0001","Cos_96_DW_2mL_0001")
    
    nimbus(disp_list,"<PATH TO .CSV FILE>")
```