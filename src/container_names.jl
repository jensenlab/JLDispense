struct ContainerName 
    name::String
    shape::Tuple{Int64,Int64}
    cobra::String
    mantis::String
    tempest::String
    nimbus::String
end 




containers=Dict{Symbol,ContainerName}()

containers[:WP96]=ContainerName("WP96",(8,12),"96 Costar","PT3-96-Assay.pd.txt","PT3-96-Assay.pd.txt","")
containers[:WP384]=ContainerName("WP384",(16,24),  "384 Well p/n 3575 3576", "PT9-384-Assay.pd.txt","PT9-384-Assay.pd.txt","")
containers[:brPCR96]=ContainerName("brPCR96",(8,12),    "","breakaway_pcr_96.pd.txt","breakaway_pcr_96.pd.txt","")
containers[:dWP96_2ml]=ContainerName("dWP96_2ml",(8,12),"Deep Well 2 ml","PT3-96-Assay-DW.pd.txt","","")

containers[:conical_15ml]=ContainerName("concial_15ml",(1,1),"","","","")
containers[:conical_50ml]=ContainerName("conical_50ml",(1,1),"","","","")

