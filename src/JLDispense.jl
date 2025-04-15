module JLDispense
import Base: showerror,length
using 
    CSV,
    DataFrames,
    Gurobi,
    JLIMS,
    JLConstants,
    JuMP,
    Random,
    Unitful,
    UnitfulParsableString,
    JSON,
    UUIDs,
    Plots,
    TextWrap

import Plots.plot 

include("./RandomProtocolNames/random_protocol_name.jl")
include("./Configuration/Configuration.jl")
include("./utils.jl")
#include("./DispenseSolver/dispense_solver.jl")
include("./Cobra/Cobra.jl")
include("./Mantis/Mantis.jl")
include("./Tempest/Tempest.jl")
include("./Nimbus/Nimbus.jl")
include("./GilsonPlateMaster/PlateMaster.jl")
include("./Micropipettors/SingleChannel.jl")
include("./Micropipettors/EightChannel.jl")



#types
export random_protocol_name

#export MixingError,OverdraftError,InsufficientIngredientError,ContainerError,StockCompatibilityError
#export strain_array, ingredient_array
#export mixer

#export dispense_solver, minimize_transfers!,minimize_labware!,minimize_overdrafts!,minimize_sources!,enforce_maxShots!,minimize_overshots!,minimize_labware_crossover!,minimize_robots!


end # module JLDispense
