module JLDispense
import Base.showerror
using 
    CSV,
    DataFrames,
    Gurobi,
    JLIMS,
    JuMP,
    Random,
    Unitful,
    UnitfulParsableString,
    JSON
include("./RandomProtocolNames/random_protocol_name.jl")
include("./Robot/Robot.jl")
include("./DispenseSolver/dispense_solver.jl")
include("./Cobra/Cobra.jl")
include("./Mantis/Mantis.jl")
include("./Tempest/Tempest.jl")
include("./Nimbus/Nimbus.jl")
include("./Human/Human.jl")




#types
export random_protocol_name
export showerror
export Robot,RobotProperties,RobotConfiguration,DeckPosition,Human,HumanProperties,HumanConfiguration,Cobra,CobraProperties,CobraConfiguration,Mantis,MantisProperties,MantisConfiguration
export feasibility, MixingError,OverdraftError,InsufficientIngredientError,ContainerError,StockCompatibilityError
export mixer
export cobra, cobra_default
export human,human_default, omnipotent_robot
export mantis, mantis_default
export tempest,multi_tempest,tempest_default
export nimbus,nimbus_default
export dispense_solver, minimize_transfers!,minimize_labware!,minimize_overdrafts!,minimize_sources!,enforce_maxShots!,minimize_overshots!,minimize_labware_crossover!


end # module JLDispense
