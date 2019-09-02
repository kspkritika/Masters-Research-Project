using DelimitedFiles
using DifferentialEquations
using Optim

include("Dynamic_Cybernetic.jl")
include("Dynamic_multi_ode.jl")
include("Objective.jl")
include("loadData.jl")

#DefineTime
tStart = 0.3;
tStop = 9.5;
tStep = 0.1;


DataDict = Dynamic_Cybernetic(tStart,tStep,tStop)
x0 = DataDict["InitialConditions"]
Z = DataDict["ModeMatrix"]
(num_reations,num_modes) = size(Z)
# IC_Ensemble = Data_dict["IC_Ensemble"]
S = DataDict["Stoich"]
tspan = (0.3,9.5)
# P = DataDict["Parameters"]

Param = DataDict["Parameter_Ensemble"]

P_opt = zeros(10,10)
#
# time = zeros(93,100)
# sol_bio = zeros(93,100)
# sol_glc = zeros(93,100)
# sol_lac = zeros(93,100)
# sol_ser = zeros(93,100)
# sol_asn = zeros(93,100)
# sol_gln = zeros(93,100)
# sol_glu = zeros(93,100)
# sol_ala = zeros(93,100)
# sol_asp = zeros(93,100)
# sol_gly = zeros(93,100)
# sol_nh3 = zeros(93,100)
#
for i = 1:10
    P = Param[:,i]
    cost(P) = Objective(P)
    res = optimize(cost,P,NelderMead())
    P_new = Optim.minimizer(res)
    P_opt[:,i] = P_new
    # @show P_opt
end
    # prob = ODEProblem(odefunc,x0,tspan,P_new)
#     sol = solve(prob, saveat=0.1)
#     time[:,i] = sol.t
#     sol_bio[:,i] = sol[14,:]
#     sol_glc[:,i] = sol[4,:]
#     sol_lac[:,i] = sol[5,:]
#     sol_ser[:,i] = sol[6,:]
#     sol_asn[:,i] = sol[7,:]
#     sol_gln[:,i] = sol[8,:]
#     sol_glu[:,i] = sol[9,:]
#     sol_ala[:,i] = sol[10,:]
#     sol_asp[:,i] = sol[11,:]
#     sol_gly[:,i] = sol[12,:]
#     sol_nh3[:,i] = sol[13,:]
# end
#
# writedlm( "Solution_Ensemble/biomass.csv",  sol_bio, ',')
# writedlm( "Solution_Ensemble/glucose.csv",  sol_glc, ',')
# writedlm( "Solution_Ensemble/lactate.csv",  sol_lac, ',')
# writedlm( "Solution_Ensemble/serine.csv",  sol_ser, ',')
# writedlm( "Solution_Ensemble/asparagine.csv",  sol_asn, ',')
# writedlm( "Solution_Ensemble/glutamine.csv",  sol_gln, ',')
# writedlm( "Solution_Ensemble/glutamate.csv",  sol_glu, ',')
# writedlm( "Solution_Ensemble/alanine.csv",  sol_ala, ',')
# writedlm( "Solution_Ensemble/asparate.csv",  sol_asp, ',')
# writedlm( "Solution_Ensemble/glycine.csv",  sol_gly, ',')
# writedlm( "Solution_Ensemble/nh3.csv",  sol_nh3, ',')
# P_opt = zeros(9,100)

# for i = 1:100
# IC = IC_Ensemble[:,i]

# cost(P) = Objective(P)
# res = optimize(cost,P,NelderMead())
# P_new = Optim.minimizer(res)

# P_opt[:,i] = P_new

# prob = ODEProblem(odefunc,x0,tspan,P_new)
# sol = solve(prob)

# enz1 = sol[1,:];
# enz2 = sol[2,:];
# enz3 = sol[3,:];
# xglc = sol[4,:];
# xlac = sol[5,:];
# xser = sol[6,:];
# xasn = sol[7,:];
# xgln = sol[8,:];
# xglu = sol[9,:];
# xala = sol[10,:];
# xasp = sol[11,:];
# xgly = sol[12,:];
# xnh3 = sol[13,:];
# xbio = sol[14,:];
#
# t = sol.t;
#
# using PyPlot
# figure(1)
# plot(t,xbio,color="black")
# scatter(tBio,Bio,color="black")
# #
# figure(2)
# plot(t,xglc,color="black")
# scatter(tGlc,Glc,color="black")
# # #
# figure(3)
# plot(t,xlac,color="black")
# scatter(tLac,Lac,color="black")
# # #
# figure(4)
# plot(t,xser,color="black")
# scatter(tSer,Ser,color="black")
# #
# figure(5)
# plot(t,xasn,color="black")
# scatter(tAsn,Asn,color="black")
# #
# figure(6)
# plot(t,xgln,color="black")
# scatter(tGln,Gln,color="black")
# # #
# figure(7)
# plot(t,xglu,color="black")
# scatter(tGlu,Glu,color="black")
# #
# figure(8)
# plot(t,xala,color="black")
# scatter(tAla,Ala,color="black")
# # #
# figure(9)
# plot(t,xasp,color="black")
# scatter(tAsp,Asp,color="black")
# # #
# figure(10)
# plot(t,xgly,color="black")
# scatter(tGly,Gly,color="black")
# # #
# figure(11)
# plot(t,xnh3,color="black")
# scatter(tnh3,nh3,color="black")
#
# figure(12)
# plot(t,enz1,color="black")
# plot(t,enz2,color="orange")
# plot(t,enz3,color="green")










# include("Dynamic_Cybernetic.jl")
# # # #
# # # #DefineTime
# tStart = 0.3;
# tStop = 9.5;
# tStep = 0.1;
# # # # # tSim = collect(tStart:tStep:tStop)
# # # #
# DataDict = Dynamic_Cybernetic(tStart,tStep,tStop)
# x0 = DataDict["InitialConditions"]
# Z = DataDict["ModeMatrix"]
# (num_reations,num_modes) = size(Z)
# S = DataDict["Stoich"]
# tspan = (0.3,9.5)
# Param = DataDict["Parameter_Ensemble"]
# # #
time = zeros(93,10)
sol_bio = zeros(93,10)
sol_glc = zeros(93,10)
sol_lac = zeros(93,10)
sol_ser = zeros(93,10)
sol_asn = zeros(93,10)
sol_gln = zeros(93,10)
sol_glu = zeros(93,10)
sol_ala = zeros(93,10)
sol_asp = zeros(93,10)
sol_gly = zeros(93,10)
sol_nh3 = zeros(93,10)
# #
for i = 1:10
    P = P_opt[:,i]
    prob = ODEProblem(odefunc,x0,tspan,P)
    sol = solve(prob, saveat=0.1)
    time[:,i] = sol.t
    sol_bio[:,i] = sol[14,:]
    sol_glc[:,i] = sol[4,:]
    sol_lac[:,i] = sol[5,:]
    sol_ser[:,i] = sol[6,:]
    sol_asn[:,i] = sol[7,:]
    sol_gln[:,i] = sol[8,:]
    sol_glu[:,i] = sol[9,:]
    sol_ala[:,i] = sol[10,:]
    sol_asp[:,i] = sol[11,:]
    sol_gly[:,i] = sol[12,:]
    sol_nh3[:,i] = sol[13,:]
end
#
# # writedlm( "Solution_Ensemble/biomass.csv",  sol_bio, ',')
# # writedlm( "Solution_Ensemble/glucose.csv",  sol_glc, ',')
# # writedlm( "Solution_Ensemble/lactate.csv",  sol_lac, ',')_
# # writedlm( "Solution_Ensemble/serine.csv",  sol_ser, ',')
# # writedlm( "Solution_Ensemble/asparagine.csv",  sol_asn, ',')
# # writedlm( "Solution_Ensemble/glutamine.csv",  sol_gln, ',')
# # writedlm( "Solution_Ensemble/glutamate.csv",  sol_glu, ',')
# # writedlm( "Solution_Ensemble/alanine.csv",  sol_ala, ',')
# # writedlm( "Solution_Ensemble/asparate.csv",  sol_asp, ',')
# # writedlm( "Solution_Ensemble/glycine.csv",  sol_gly, ',')
# # writedlm( "Solution_Ensemble/nh3.csv",  sol_nh3, ',')
#
# # Lets calculate the mean of the solution space obtained
bio_mean = mean(sol_bio',dims=1)
glc_mean = mean(sol_glc',dims=1)
lac_mean = mean(sol_lac',dims=1)
ser_mean = mean(sol_ser',dims=1)
asn_mean = mean(sol_asn',dims=1)
gln_mean = mean(sol_gln',dims=1)
glu_mean = mean(sol_glu',dims=1)
ala_mean = mean(sol_ala',dims=1)
asp_mean = mean(sol_asp',dims=1)
gly_mean = mean(sol_gly',dims=1)
nh3_mean = mean(sol_nh3',dims=1)
#
# # Lets Calculate the standard deviation of the given system and get values for 95% confidence interval
bio_std = zeros(93,1)
#Biomass
for i = 1:93
    bio_std[i] = std(sol_bio'[:,i]).*0.8158            #0.2576
end

glc_std = zeros(93,1)
#Glucose
for i = 1:93
    glc_std[i] = std(sol_glc'[:,i]).*0.8158
end

#Lactate
lac_std = zeros(93,1)
for i = 1:93
    lac_std[i] = std(sol_lac'[:,i]).*0.8158
end

#Serine
ser_std = zeros(93,1)
for i = 1:93
    ser_std[i] = std(sol_ser'[:,i]).*0.8158
end

#Asparagine
asn_std = zeros(93,1)
for i = 1:93
    asn_std[i] = std(sol_asn'[:,i]).*0.8158
end

#Glutamine
gln_std = zeros(93,1)
for i = 1:93
    gln_std[i] = std(sol_gln'[:,i]).*0.8158
end

#Glutamate
glu_std = zeros(93,1)
for i = 1:93
    glu_std[i] = std(sol_glu'[:,i]).*0.8158
end

#Alanine
ala_std = zeros(93,1)
for i = 1:93
    ala_std[i] = std(sol_ala'[:,i]).*0.8158
end

#Asparate
asp_std = zeros(93,1)
for i = 1:93
    asp_std[i] = std(sol_asp'[:,i]).*0.8158
end

#Glycine
gly_std = zeros(93,1)
for i = 1:93
    gly_std[i] = std(sol_gly'[:,i]).*0.8158
end

#nh3
nh3_std = zeros(93,1)
for i = 1:93
    nh3_std[i] = std(sol_nh3'[:,i]).*0.8158
end

time = time[:,1]

# Lets plot the mean solutions and experimental Data
using PyPlot
figure(1)
upper_bio =  bio_mean' + bio_std
lower_bio =  bio_mean' - bio_std
ub_bio = upper_bio[:,1]
lb_bio = lower_bio[:,1]
fill_between(time, ub_bio, lb_bio,color="lightblue", alpha=1)
plot(time,bio_mean',color="black", label = "Simulated")
scatter(tBio,Bio,color="black", label = "Experimental")
errorbar(tBio,Bio,yerr = 0.2*Bio, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Biomass Profile")
tight_layout()
legend()
savefig("bio.svg")
savefig("bio.pdf")
savefig("bio.png")

figure(2)
upper_glc =  glc_mean' + glc_std
lower_glc =  glc_mean' - glc_std
ub_glc = upper_glc[:,1]
lb_glc = lower_glc[:,1]
fill_between(time, ub_glc, lb_glc,color="lightblue", alpha=1)
plot(time,glc_mean',color="black", label = "Simulated")
scatter(tGlc,Glc,color="black", label = "Experimental")
errorbar(tGlc,Glc,yerr = 0.2*Glc, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Glucose Profile")
tight_layout()
legend()
savefig("glc.svg")
savefig("glc.pdf")
savefig("glc.png")

figure(3)
upper_lac =  lac_mean' + lac_std
lower_lac =  lac_mean' - lac_std
ub_lac = upper_lac[:,1]
lb_lac = lower_lac[:,1]
fill_between(time, ub_lac, lb_lac,color="lightblue", alpha=1)
plot(time,lac_mean',color="black", label = "Simulated")
scatter(tLac,Lac,color="black", label = "Experimental")
errorbar(tLac,Lac,yerr = 0.2*Lac, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Lactate Profile")
tight_layout()
legend()
savefig("lac.svg")
savefig("lac.pdf")
savefig("lac.png")

figure(4)
upper_ser =  ser_mean' + ser_std
lower_ser =  ser_mean' - ser_std
ub_ser = upper_ser[:,1]
lb_ser = lower_ser[:,1]
fill_between(time, ub_ser, lb_ser,color="lightblue", alpha=1)
plot(time,ser_mean',color="black", label = "Simulated")
scatter(tSer,Ser,color="black", label = "Experimental")
errorbar(tSer,Ser,yerr = 0.2*Ser, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Serine Profile")
tight_layout()
legend()
savefig("ser.svg")
savefig("ser.pdf")
savefig("ser.png")

figure(5)
upper_asn =  asn_mean' + asn_std
lower_asn =  asn_mean' - asn_std
ub_asn = upper_asn[:,1]
lb_asn = lower_asn[:,1]
fill_between(time, ub_asn, lb_asn,color="lightblue", alpha=1)
plot(time,asn_mean',color="black", label = "Simulated")
scatter(tAsn,Asn,color="black", label = "Experimental")
errorbar(tAsn,Asn,yerr = 0.2*Asn, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Asparagine Profile")
tight_layout()
legend()
savefig("asn.svg")
savefig("asn.pdf")
savefig("asn.png")

figure(6)
upper_gln =  gln_mean' + gln_std
lower_gln =  gln_mean' - gln_std
ub_gln = upper_gln[:,1]
lb_gln = lower_gln[:,1]
fill_between(time, ub_gln, lb_gln,color="lightblue", alpha=1)
plot(time,gln_mean',color="black", label = "Simulated")
scatter(tGln,Gln,color="black", label = "Experimental")
errorbar(tGln,Gln,yerr = 0.2*Gln, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Glutamine Profile")
tight_layout()
legend()
savefig("gln.svg")
savefig("gln.pdf")
savefig("gln.png")

figure(7)
upper_glu =  glu_mean' + glu_std
lower_glu =  glu_mean' - glu_std
ub_glu = upper_glu[:,1]
lb_glu = lower_glu[:,1]
fill_between(time, ub_glu, lb_glu,color="lightblue", alpha=1)
plot(time,glu_mean',color="black", label = "Simulated")
scatter(tGlu,Glu,color="black", label = "Experimental")
errorbar(tGlu,Glu,yerr = 0.2*Glu, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Glutamate Profile")
tight_layout()
legend()
savefig("glu.svg")
savefig("glu.pdf")
savefig("glu.png")

figure(8)
upper_ala =  ala_mean' + ala_std
lower_ala =  ala_mean' - ala_std
ub_ala = upper_ala[:,1]
lb_ala = lower_ala[:,1]
fill_between(time, ub_ala, lb_ala,color="lightblue", alpha=1)
plot(time,ala_mean',color="black", label = "Simulated")
scatter(tAla,Ala,color="black", label = "Experimental")
errorbar(tAla,Ala,yerr = 0.2*Ala, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Alanine Profile")
tight_layout()
legend()
savefig("ala.pdf")
savefig("ala.png")
# #
figure(9)
upper_asp =  asp_mean' + asp_std
lower_asp =  asp_mean' - asp_std
ub_asp = upper_asp[:,1]
lb_asp = lower_asp[:,1]
fill_between(time, ub_asp, lb_asp,color="lightblue", alpha=1)
plot(time,asp_mean',color="black", label = "Simulated")
scatter(tAsp,Asp,color="black", label = "Experimental")
errorbar(tAsp,Asp,yerr = 0.2*Asp, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Asparate Profile")
tight_layout()
legend()
savefig("asp.pdf")
savefig("asp.png")
# #
figure(10)
upper_gly =  gly_mean' + gly_std
lower_gly =  gly_mean' - gly_std
ub_gly = upper_gly[:,1]
lb_gly = lower_gly[:,1]
fill_between(time, ub_gly, lb_gly,color="lightblue", alpha=1)
plot(time,gly_mean',color="black", label = "Simulated")
scatter(tGly,Gly,color="black", label = "Experimental")
errorbar(tGly,Gly,yerr = 0.2*Gly, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Glycine Profile")
tight_layout()
legend()
savefig("gly.pdf")
savefig("gly.png")
# #
figure(11)
upper_nh3 =  nh3_mean' + nh3_std
lower_nh3 =  nh3_mean' - nh3_std
ub_nh3 = upper_nh3[:,1]
lb_nh3 = lower_nh3[:,1]
fill_between(time, ub_nh3, lb_nh3,color="lightblue", alpha=1)
plot(time,nh3_mean',color="black", label = "Simulated")
scatter(tnh3,nh3,color="black", label = "Experimental")
errorbar(tnh3,nh3,yerr = 0.2*nh3, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Ammonia Profile")
tight_layout()
legend()
savefig("amm.pdf")
savefig("amm.png")
