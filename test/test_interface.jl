import JLDispense: volconc_to_stock, volconc_to_quant, stock_to_quant, stock_to_volconc, quant_to_stock, quant_to_volconc


n = 20
vc_names = ["volume","water","alanine","glucose","sodium_chloride","a","b"]

v = fill(200,n) 
w = fill(100,n) 
chems = 10* rand(n,5) 

m = hcat(v,w,chems)

df = DataFrame(m, vc_names)
units = DataFrame(["µL" "%" "g/l"  "mg/mL" "mM" "g/l" "ng/µl"],vc_names)

stocks = volconc_to_stock(df,units;chem_context=[JLIMS,JLConstants]) # parse into stocks 

vc, un  = stock_to_volconc(stocks;chem_context=[JLIMS,JLConstants]) # parse back into dataframes 

stocks2 = volconc_to_stock(vc,un;chem_context=[JLIMS,JLConstants]) # parse back into stocks 

q, q_un = stock_to_quant(stocks2;chem_context=[JLIMS,JLConstants]) 


stocks3 = quant_to_stock(q,q_un;chem_context=[JLIMS,JLConstants])


@testset "Interface" begin 
    @test all(isapprox.(stocks,stocks2))
    @test all(isapprox.(stocks2,stocks3)) 
end 



