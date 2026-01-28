

Ma,Md,Sa,Sd = JLDispense.masks(cobra.head,JLConstants.WP96(1,"test"))


# the 

C = 4 
P = 60 
W = 96 

a_counts = zeros(Int64, 8,12)

for p in 1: P 
    for c in 1:C
        for w in 1:W 

            if Ma(w,p,c) 
                a_counts[w] +=1 
            end 
        end 
    end 
end 

C = 4 
P = 132 
W = 96 

d_counts = zeros(Int64, 8,12)

for p in 1: P 
    for c in 1:C
        for w in 1:W 

            if Md(w,p,c) 
                d_counts[w] +=1 
            end 
        end 
    end 
end 