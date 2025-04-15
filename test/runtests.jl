using JLDispense, JLIMS, JLConstants , Unitful, Plots 


x=[Conical50(1,"conical1"),Conical50(1,"conical2"),Conical50(1,"conical3"),WP96(1,"plate1")]


JLDispense.can_place(x[1],JLDispense.tuberack50mL_0001)

y=JLDispense.nimbus_slotting_greedy(x,JLDispense.nimbus)

plot(y)