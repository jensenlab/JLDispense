module JLDispense
import Base.showerror
using 
    CSV,
    DataFrames,
    Gurobi,
    JLIMS,
    JuMP,
    RandomProtocolNames,
    Unitful,
    JSON

include("./Robot/Robot.jl")
include("./Mixer/Mixer.jl")
include("./Feasibility/Feasibility.jl")
include("./Cobra/Cobra.jl")
#include("./Mantis/Mantis.jl")
#include("./Tempest/Tempest.jl")
#include("./Nimbus/Nimbus.jl")
include("./Human/Human.jl")




#types
export showerror
export Robot,RobotProperties,RobotConfiguration,DeckPosition,Human,HumanProperties,HumanConfiguration,Cobra,CobraProperties,CobraConfiguration
export feasibility, MixingError,OverdraftError,InsufficientIngredientError,ContainerError,StockCompatibilityError
export mixer
export cobra, cobra_default
export human,human_default


end # module JLDispense
