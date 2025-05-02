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
    StatsBase,
    LinearAlgebra,
    TextWrap

import Plots.plot 

include("./RandomProtocolNames/random_protocol_name.jl")
include("./Equipment/Configuration.jl")
include("./utils.jl")
include("./Slotting/slotting.jl")
#include("./DispenseSolver/dispense_solver.jl")
include("./Equipment/Cobra.jl")
include("./Equipment/Mantis.jl")
include("./Equipment/Tempest.jl")
include("./Equipment/Nimbus.jl")
include("./Equipment/PlateMaster.jl")
include("./Equipment/SingleChannel.jl")
include("./Equipment/EightChannel.jl")
include("./Equipment/NullRobot.jl")
include("./DispenseSolver/utils.jl")
include("./DispenseSolver/objectives.jl")
include("./DispenseSolver/dispense_solver.jl")


#types
export random_protocol_name

export min_operations!,min_sources!,min_labware_crossover!

#export MixingError,OverdraftError,InsufficientIngredientError,ContainerError,StockCompatibilityError
#export strain_array, ingredient_array
#export mixer

#export dispense_solver, minimize_transfers!,minimize_labware!,minimize_overdrafts!,minimize_sources!,enforce_maxShots!,minimize_overshots!,minimize_labware_crossover!,minimize_robots!


end # module JLDispense
