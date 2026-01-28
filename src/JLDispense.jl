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

#include("./DispenseSolver/dispense_solver.jl")
include("./Equipment/Cobra.jl")
include("./Equipment/Mantis.jl")
include("./Equipment/Tempest.jl")
include("./Equipment/Nimbus.jl")
include("./Equipment/PlateMaster.jl")
include("./Equipment/SingleChannel.jl")
include("./Equipment/EightChannel.jl")
include("./Equipment/NullRobot.jl")
include("./DispenseSolver/interface.jl")
include("./DispenseSolver/utils.jl")
include("./DispenseSolver/objectives.jl")
include("./DispenseSolver/dispense_solver.jl")
include("./Slotting/slotting.jl")
include("./Compilier/protocol_compilier.jl")
include("./Scheduler/scheduler.jl")


#types
export random_protocol_name

export min_operations!,min_sources!,min_labware_crossover!

export dispense_solver 
export protocol_compiler 
export scheduler 

export cobra
export nimbus 
export platemaster 
export p1000, p200, p20, p2
export  mantis_lv, mantis_hv 
export  tempest_lv , tempest_hv
export  multichannel_p1000_h,multichannel_p1000_v, multichannel_p100_h ,multichannel_p100_v , multichannel_p10_h , multichannel_p10_v 

end # module JLDispense
