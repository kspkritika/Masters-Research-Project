using DelimitedFiles
using DifferentialEquations
using Statistics
# using Sundials
# using ODE
include("loadData.jl")

function odefunc(dxdt,x,P,t)

    # mode_matrix = float(readdlm("Data/dynamic_mode1.csv",','))[2649:3194,:];
    #
    # Z = deepcopy(mode_matrix)
    #
    # (num_reations,num_modes) = size(Z)
    #
    # # normalisation on the basis of biomass flux
    # a = 6.3108538979; # the values after 4690 are not normalized since it just becomes 1. make this correction in ecoli too
    # b = 2.3943091593;
    # c = 1.7519335312;
    #
    # for j=1:546
    #     Z[j,1] = Z[j,1]/a
    # end
    #
    # for j=1:546
    #     Z[j,2] = Z[j,2]/b
    # end
    #
    # for j=1:546
    #     Z[j,3] = Z[j,3]/c
    # end

    # # blocking product based exchange reactions for mode 1
    # Z[2922,1] = 0*Z[2922,1] #GLUTAMATE EXCHANGE
    # Z[3004,1] = 0*Z[3004,1] #LACTATE EXCHANGE
    # Z[2737,1] = 0*Z[2737,1] #ALANINE EXCHANGE
    # Z[2756,1] = 0*Z[2756,1] #ASPARATE ECHANGE
    # Z[2924,1] = 0*Z[2924,1] #GLYCINE EXCHANGE
    # Z[3044,1] = 0*Z[3044,1] #NH3 EXCHANGE
    #
    # # blocking product based exchange reactions for mode 2
    # Z[2922,1] = 0*Z[2922,2] #GLUTAMATE EXCHANGE
    # Z[2737,1] = 0*Z[2737,2] #ALANINE EXCHANGE
    # Z[2756,1] = 0*Z[2756,2] #ASPARATE ECHANGE
    # Z[2919,1] = 0*Z[2919,2] #GLUTAMINE ECHANGE
    #
    # # blocking product based exchange reactions for mode 3
    # Z[2922,1] = 0*Z[2922,3] #GLUTAMATE EXCHANGE
    # Z[2919,1] = 0*Z[2919,3] #GLUTAMINE ECHANGE
    # Z[2924,1] = 0*Z[2924,3] #GLYCINE EXCHANGE
    # Z[3044,1] = 0*Z[3044,3] #NH3 EXCHANGE

    #Load Stoichiometric Matrix
    # stm_matrix = float(readdlm("Data/Network.csv",','));
    # S_idx = zeros(10,546)
    # S_idx[1,:] = stm_matrix[293,2649:3194]#glc
    # S_idx[2,:] = stm_matrix[320,2649:3194]  #*-1 #lac is produced :please fix this after fixing the modes. shouldn't require the negative sign
    # S_idx[3,:] = stm_matrix[199,2649:3194] # ser
    # S_idx[4,:] = stm_matrix[193,2649:3194] # asparagine
    # S_idx[5,:] = stm_matrix[197,2649:3194] # glutamine
    # S_idx[6,:] = stm_matrix[249,2649:3194] # glutamate
    # S_idx[7,:] = stm_matrix[170,2649:3194] # alanine
    # S_idx[8,:] = stm_matrix[211,2649:3194] # asparate
    # S_idx[9,:] = stm_matrix[247,2649:3194] # glycine
    # S_idx[10,:] = stm_matrix[319,2649:3194] # nh3
    # S = abs.(S_idx)

    # Alias the species vector -
    M_glc = x[4];
    M_lac = x[5];
    M_ser = x[6];
    M_asn = x[7];
    M_gln = x[8];
    M_glu = x[9];
    M_ala = x[10];
    M_asp = x[11];
    M_gly = x[12];
    M_nh3 = x[13];
    M_bio = x[14];

    # Declare Parameters
    ke1 = P[1]
    alpha = 0.005; #P[2]
    beta1 = P[2]; #2.772;
    kmax1 = P[3]
    Kser1 = 1; #P[5]
    Kglc1 = 5; #P[6]
    Kasn1 = 0.3; #P[7]
    Kgln1 = 0.4; #P[8]
    kmax2 = P[4]
    Kser2 = 1; #P[10]
    Kglc2 = 20; #P[11]
    Kasn2 = 0.3; #P[12]
    Klac2 = 4; #P[13]
    Kgly2 = 3; #P[14]
    Knh2 = 0.6; #P[15]
    kmax3 = P[5]
    Kser3 = 1; #P[17]
    Kglc3 = 20; #P[18]
    Kasn3 = 0.3; #P[19]
    Klac3 = 4; #P[20]
    Kala3 = 2.5; #P[21]
    Kasp3 = 1; #P[22]
    Kgly3 = 3; #P[23]
    ke2 = P[6]
    ke3 = P[7]
    beta2 = P[8]; #2.772
    beta3 = P[9]; #2.772
    s_o = P[10]

    # Initialize rate_vectors
    rM = Float64[]
    rE = Float64[]
    rG = Float64[]

    #Metabolite Reaction
    fill!(rM,0.0)
    push!(rM,kmax1*x[1]*(M_glc/(Kglc1+M_glc))*(M_ser/(Kser1+M_ser))*(M_asn/(Kasn1+M_asn))*(M_gln/(Kgln1+M_gln)))
    push!(rM,kmax2*x[2]*(M_glc/(Kglc2+M_glc))*(M_ser/(Kser2+M_ser))*(M_asn/(Kasn2+M_asn))*(M_lac/(Klac2+M_lac))*(M_gly/(Kgly2+M_gly))*(M_nh3/(Knh2+M_nh3)))
    push!(rM,kmax3*x[3]*(M_glc/(Kglc3+M_glc))*(M_ser/(Kser3+M_ser))*(M_asn/(Kasn3+M_asn))*(M_lac/(Klac3+M_lac))*(M_ala/(Kala3+M_ala))*(M_asp/(Kasp3+M_asp))*(M_gly/(Kgly3+M_gly)))

    # push!(rM,kmax1*x[1]*(M_glc/(Kglc1+M_glc))*(M_ser/(Kser1+M_ser))*(M_asn/(Kasn1+M_asn))*(M_gln/(Kgln1+M_gln)))
    # push!(rM,kmax2*x[2]*(M_glc/(Kglc2+M_glc))*(M_ser/(Kser2+M_ser))*(M_asn/(Kasn2+M_asn)))
    # push!(rM,kmax3*x[3]*(M_glc/(Kglc3+M_glc))*(M_ser/(Kser3+M_ser))*(M_asn/(Kasn3+M_asn))*(M_lac/(Klac3+M_lac)))


    #EnzymeReactionRate
    fill!(rE,0.0)
    push!(rE,ke1*(M_glc/(Kglc1+M_glc))*(M_ser/(Kser1+M_ser))*(M_asn/(Kasn1+M_asn))*(M_gln/(Kgln1+M_gln)))
    push!(rE,ke2*(M_glc/(Kglc2+M_glc))*(M_ser/(Kser2+M_ser))*(M_asn/(Kasn2+M_asn))*(M_lac/(Klac2+M_lac))*(M_gly/(Kgly2+M_gly))*(M_nh3/(Knh2+M_nh3)))
    push!(rE,ke3*(M_glc/(Kglc3+M_glc))*(M_ser/(Kser3+M_ser))*(M_asn/(Kasn3+M_asn))*(M_lac/(Klac3+M_lac))*(M_ala/(Kala3+M_ala))*(M_asp/(Kasp3+M_asp))*(M_gly/(Kgly3+M_gly)))

    # push!(rE,ke*(M_glc/(Kglc1+M_glc))*(M_ser/(Kser1+M_ser))*(M_asn/(Kasn1+M_asn))*(M_gln/(Kgln1+M_gln)))
    # push!(rE,ke*(M_glc/(Kglc2+M_glc))*(M_ser/(Kser2+M_ser))*(M_asn/(Kasn2+M_asn)))
    # push!(rE,ke*(M_glc/(Kglc3+M_glc))*(M_ser/(Kser3+M_ser))*(M_asn/(Kasn3+M_asn))*(M_lac/(Klac3+M_lac)))


    # Growthrate
    fill!(rG,0.0)
    push!(rG,(5.3108538979*rM[1]))
    push!(rG,(4.3943091593*rM[2]))
    push!(rG,(2.7519335312*rM[3]))

    # for i = 1:num_modes
    #     push!(rG,(Z[4690,i]*rM[i])) ##biomass for a producing cell line
    # end

    #Initialize Control Vector
    cybernetic_var = Float64[]
    fill!(cybernetic_var,0.0)
    # push!(cybernetic_var,(abs.(Z[475,1]+Z[268,1]+Z[106,1]+Z[271,1]))*rM[1])
    # push!(cybernetic_var,(abs.(Z[475,2]+Z[268,2]+Z[106,2]+Z[356,2]+Z[276,2]+Z[396,2]))*rM[2])
    # push!(cybernetic_var,(abs.(Z[475,3]+Z[268,3]+Z[106,3]+Z[356,3]+Z[89,3]+Z[108,3]+Z[276,3]))*rM[3])

    push!(cybernetic_var,(abs.(Z[271,1]))*rM[1])
    push!(cybernetic_var,(abs.(Z[106,2]))*rM[2])
    push!(cybernetic_var,(abs.(Z[106,3]+Z[356,3]))*rM[3])



    # push!(cybernetic_var,(abs.(Z[475,1]+Z[268,1]+Z[106,1]+Z[271,1]))*rM[1])
    # push!(cybernetic_var,(abs.(Z[475,2]+Z[268,2]+Z[106,2]))*rM[2])
    # push!(cybernetic_var,(abs.(Z[475,3]+Z[268,3]+Z[106,3]+Z[356,3]))*rM[3])




    #Initialize cybernetic variables u and v
    # u = zeros(num_modes)
    # v = zeros(num_modes)

    u = zeros(3)
    v = zeros(3)


    for i = 1:3
        u[i] = cybernetic_var[i]/sum(cybernetic_var)
        v[i] = cybernetic_var[i]/maximum(cybernetic_var)
    end

    # @show u
    # @show v

    #Redefine rates (rate*control)
    Enzyme_rate = rE.*u;
    Growth_rate = rG.*v;
    Reaction_rate = rM.*v;
    mu = sum(Growth_rate)

# This conditional statement takes care of the fact that the dilution starts at day 3.
# there is no dilution before day 3. Vo = 1.7. dV/dt = F = 0.1. D = F/(Vo + Ft)
    if t < 3
        D = 0
    else
        D = (0.1/(1.7 + (0.1*t)))
    end

    # dxdt = similar(x)
    dxdt[1] = alpha .+ (rE[1]*u[1]) .- (beta1+mu)*x[1]
    dxdt[2] = alpha .+ (rE[2]*u[2]) .- (beta2+mu)*x[2]
    dxdt[3] = alpha .+ (rE[3]*u[3]) .- (beta3+mu)*x[3]
    # dxdt[1:3] = alpha .+ Enzyme_rate .- (beta+mu)*x[1:3]
    # dxdt[4:13] =  0.05555*s_o .+ S[1:10,:]*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[4] =  D*(s_o - x[4]) .+ S[1,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[5] =  D*(s_o - x[5]) .+ S[2,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[6] =  D*(s_o - x[6]) .+ S[3,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[7] =  D*(s_o - x[7]) .+ S[4,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[8] =  D*(s_o - x[8]) .+ S[5,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[9] =  D*(s_o - x[9]) .+ S[6,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[10] =  D*(s_o - x[10]) .+ S[7,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[11] =  D*(s_o - x[11]) .+ S[8,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[12] =  D*(s_o - x[12]) .+ S[9,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[13] =  D*(s_o - x[13]) .+ S[10,:]'*Z*Reaction_rate*M_bio  # 0.05555*s_o .-
    dxdt[14] = (mu - D)*M_bio
    # dxdt
    return dxdt
end

include("Dynamic_Cybernetic.jl")
# # # #DefineTime
tStart = 0.3;
tStop = 9.5;
tStep = 0.1;
# # # # # tSim = collect(tStart:tStep:tStop)
DataDict = Dynamic_Cybernetic(tStart,tStep,tStop)
x0 = DataDict["InitialConditions"]
Z = DataDict["ModeMatrix"]
(num_reations,num_modes) = size(Z)
S = DataDict["Stoich"]
tspan = (0.3,9.5)
P = DataDict["Parameters"]
prob = ODEProblem(odefunc,x0,tspan,P)
sol = solve(prob, saveat=0.1)
#
# solutions for three mode case
enz1 = sol[1,:];
enz2 = sol[2,:];
enz3 = sol[3,:];
xglc = sol[4,:];
xlac = sol[5,:];
xser = sol[6,:];
xasn = sol[7,:];
xgln = sol[8,:];
xglu = sol[9,:];
xala = sol[10,:];
xasp = sol[11,:];
xgly = sol[12,:];
xnh3 = sol[13,:];
xbio = sol[14,:];

t = sol.t;

using PyPlot
figure(1)
plot(t,xbio,color="black")
scatter(tBio,Bio,color="black")
errorbar(tBio,Bio,yerr = 0.2*Bio, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Biomass Profile")
tight_layout()
legend()
savefig("bio.pdf")
savefig("bio.png")
# #
figure(2)
plot(t,xglc,color="black")
scatter(tGlc,Glc,color="black")
errorbar(tGlc,Glc,yerr = 0.2*Glc, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Glucose Profile")
tight_layout()
legend()
savefig("glc.pdf")
savefig("glc.png")
# # #
figure(3)
plot(t,xlac,color="black")
scatter(tLac,Lac,color="black")
errorbar(tLac,Lac,yerr = 0.2*Lac, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Lactate Profile")
tight_layout()
legend()
savefig("lac.pdf")
savefig("lac.png")
# #
figure(4)
plot(t,xser,color="black")
scatter(tSer,Ser,color="black")
errorbar(tSer,Ser,yerr = 0.2*Ser, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Serine Profile")
tight_layout()
legend()
savefig("ser.pdf")
savefig("ser.png")
#
figure(5)
plot(t,xasn,color="black")
scatter(tAsn,Asn,color="black")
errorbar(tAsn,Asn,yerr = 0.2*Asn, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Asparagine Profile")
tight_layout()
legend()
savefig("asn.pdf")
savefig("asn.png")
#
figure(6)
plot(t,xgln,color="black")
scatter(tGln,Gln,color="black")
errorbar(tGln,Gln,yerr = 0.2*Gln, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Glutamine Profile")
tight_layout()
legend()
savefig("gln.pdf")
savefig("gln.png")
# #
figure(7)
plot(t,xglu,color="black")
scatter(tGlu,Glu,color="black")
errorbar(tGlu,Glu,yerr = 0.2*Glu, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Glutamate Profile")
tight_layout()
legend()
savefig("glu.pdf")
savefig("glu.png")
# #
figure(8)
plot(t,xala,color="black")
scatter(tAla,Ala,color="black")
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
plot(t,xasp,color="black")
scatter(tAsp,Asp,color="black")
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
plot(t,xgly,color="black")
scatter(tGly,Gly,color="black")
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
plot(t,xnh3,color="black")
scatter(tnh3,nh3,color="black")
errorbar(tnh3,nh3,yerr = 0.2*nh3, fmt = "o", color="black")
xlabel("Time (days)")
ylabel("Abundance (mM)")
title("Ammonia Profile")
tight_layout()
legend()
savefig("amm.pdf")
savefig("amm.png")

figure(12)
plot(t,enz1,color="black")
plot(t,enz2,color="orange")
plot(t,enz3,color="green")







# Param = DataDict["Parameter_Ensemble"]
# # #
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
# # #
# for i = 1:100
#     P = Param[:,i]
#     prob = ODEProblem(odefunc,x0,tspan,P)
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
# #
# # # writedlm( "Solution_Ensemble/biomass.csv",  sol_bio, ',')
# # # writedlm( "Solution_Ensemble/glucose.csv",  sol_glc, ',')
# # # writedlm( "Solution_Ensemble/lactate.csv",  sol_lac, ',')_
# # # writedlm( "Solution_Ensemble/serine.csv",  sol_ser, ',')
# # # writedlm( "Solution_Ensemble/asparagine.csv",  sol_asn, ',')
# # # writedlm( "Solution_Ensemble/glutamine.csv",  sol_gln, ',')
# # # writedlm( "Solution_Ensemble/glutamate.csv",  sol_glu, ',')
# # # writedlm( "Solution_Ensemble/alanine.csv",  sol_ala, ',')
# # # writedlm( "Solution_Ensemble/asparate.csv",  sol_asp, ',')
# # # writedlm( "Solution_Ensemble/glycine.csv",  sol_gly, ',')
# # # writedlm( "Solution_Ensemble/nh3.csv",  sol_nh3, ',')
# #
# # # Lets calculate the mean of the solution space obtained
# bio_mean = mean(sol_bio',dims=1)
# glc_mean = mean(sol_glc',dims=1)
# lac_mean = mean(sol_lac',dims=1)
# ser_mean = mean(sol_ser',dims=1)
# asn_mean = mean(sol_asn',dims=1)
# gln_mean = mean(sol_gln',dims=1)
# glu_mean = mean(sol_glu',dims=1)
# ala_mean = mean(sol_ala',dims=1)
# asp_mean = mean(sol_asp',dims=1)
# gly_mean = mean(sol_gly',dims=1)
# nh3_mean = mean(sol_nh3',dims=1)
# #
# # # Lets Calculate the standard deviation of the given system and get values for 95% confidence interval
# bio_std = zeros(93,1)
# #Biomass
# for i = 1:93
#     bio_std[i] = std(sol_bio'[:,i]).*0.2576
# end
#
# glc_std = zeros(93,1)
# #Glucose
# for i = 1:93
#     glc_std[i] = std(sol_glc'[:,i]).*0.2576
# end
#
# #Lactate
# lac_std = zeros(93,1)
# for i = 1:93
#     lac_std[i] = std(sol_lac'[:,i]).*0.2576
# end
#
# #Serine
# ser_std = zeros(93,1)
# for i = 1:93
#     ser_std[i] = std(sol_ser'[:,i]).*0.2576
# end
#
# #Asparagine
# asn_std = zeros(93,1)
# for i = 1:93
#     asn_std[i] = std(sol_asn'[:,i]).*0.2576
# end
#
# #Glutamine
# gln_std = zeros(93,1)
# for i = 1:93
#     gln_std[i] = std(sol_gln'[:,i]).*0.2576
# end
#
# #Glutamate
# glu_std = zeros(93,1)
# for i = 1:93
#     glu_std[i] = std(sol_glu'[:,i]).*0.2576
# end
#
# #Alanine
# ala_std = zeros(93,1)
# for i = 1:93
#     ala_std[i] = std(sol_ala'[:,i]).*0.2576
# end
#
# #Asparate
# asp_std = zeros(93,1)
# for i = 1:93
#     asp_std[i] = std(sol_asp'[:,i]).*0.2576
# end
#
# #Glycine
# gly_std = zeros(93,1)
# for i = 1:93
#     gly_std[i] = std(sol_gly'[:,i]).*0.2576
# end
#
# #nh3
# nh3_std = zeros(93,1)
# for i = 1:93
#     nh3_std[i] = std(sol_nh3'[:,i]).*0.2576
# end
#
# time = time[:,1]
#
# # Lets plot the mean solutions and experimental Data
# using PyPlot
# figure(1)
# upper_bio =  bio_mean' + bio_std
# lower_bio =  bio_mean' - bio_std
# ub_bio = upper_bio[:,1]
# lb_bio = lower_bio[:,1]
# fill_between(time, ub_bio, lb_bio,color="lightblue", alpha=1)
# plot(time,bio_mean',color="black", label = "Simulated")
# # scatter(tBio,Bio,color="black", label = "Experimental")
# errorbar(tBio,Bio,yerr = 0.2*Bio, fmt = "o", color="black", label = "Experimental")
# xlabel("Time (hours)")
# ylabel("Abundance (mM)")
# title("Biomass Profile")
# tight_layout()
# legend()
# # savefig("bio.svg")
# # savefig("bio.pdf")
# # savefig("bio.png")
#
# figure(2)
# upper_glc =  glc_mean' + glc_std
# lower_glc =  glc_mean' - glc_std
# ub_glc = upper_glc[:,1]
# lb_glc = lower_glc[:,1]
# fill_between(time, ub_glc, lb_glc,color="lightblue", alpha=1)
# plot(time,glc_mean',color="black", label = "Simulated")
# # scatter(tGlc,Glc,color="black", label = "Experimental")
# errorbar(tGlc,Glc,yerr = 0.2*Glc, fmt = "o", color="black", label = "Experimental")
# xlabel("Time (hours)")
# ylabel("Abundance (mM)")
# title("Glucose Profile")
# tight_layout()
# legend()
# # savefig("glc.svg")
# # savefig("glc.pdf")
# # savefig("glc.png")
#
# figure(3)
# upper_lac =  lac_mean' + lac_std
# lower_lac =  lac_mean' - lac_std
# ub_lac = upper_lac[:,1]
# lb_lac = lower_lac[:,1]
# fill_between(time, ub_lac, lb_lac,color="lightblue", alpha=1)
# plot(time,lac_mean',color="black", label = "Simulated")
# # scatter(tLac,Lac,color="black", label = "Experimental")
# errorbar(tLac,Lac,yerr = 0.2*Lac, fmt = "o", color="black", label = "Experimental")
# xlabel("Time (hours)")
# ylabel("Abundance (mM)")
# title("Lactate Profile")
# tight_layout()
# legend()
# # savefig("lac.svg")
# # savefig("lac.pdf")
# # savefig("lac.png")
#
# figure(4)
# upper_ser =  ser_mean' + ser_std
# lower_ser =  ser_mean' - ser_std
# ub_ser = upper_ser[:,1]
# lb_ser = lower_ser[:,1]
# fill_between(time, ub_ser, lb_ser,color="lightblue", alpha=1)
# plot(time,ser_mean',color="black", label = "Simulated")
# # scatter(tSer,Ser,color="black", label = "Experimental")
# errorbar(tSer,Ser,yerr = 0.2*Ser, fmt = "o", color="black", label = "Experimental")
# xlabel("Time (hours)")
# ylabel("Abundance (mM)")
# title("Serine Profile")
# tight_layout()
# legend()
# # savefig("ser.svg")
# # savefig("ser.pdf")
# # savefig("ser.png")
#
# figure(5)
# upper_asn =  asn_mean' + asn_std
# lower_asn =  asn_mean' - asn_std
# ub_asn = upper_asn[:,1]
# lb_asn = lower_asn[:,1]
# fill_between(time, ub_asn, lb_asn,color="lightblue", alpha=1)
# plot(time,asn_mean',color="black", label = "Simulated")
# # scatter(tAsn,Asn,color="black", label = "Experimental")
# errorbar(tAsn,Asn,yerr = 0.2*Asn, fmt = "o", color="black", label = "Experimental")
# xlabel("Time (hours)")
# ylabel("Abundance (mM)")
# title("Asparagine Profile")
# tight_layout()
# legend()
# # savefig("asn.svg")
# # savefig("asn.pdf")
# # savefig("asn.png")
#
# figure(6)
# upper_gln =  gln_mean' + gln_std
# lower_gln =  gln_mean' - gln_std
# ub_gln = upper_gln[:,1]
# lb_gln = lower_gln[:,1]
# fill_between(time, ub_gln, lb_gln,color="lightblue", alpha=1)
# plot(time,gln_mean',color="black", label = "Simulated")
# # scatter(tGln,Gln,color="black", label = "Experimental")
# errorbar(tGln,Gln,yerr = 0.2*Gln, fmt = "o", color="black", label = "Experimental")
# xlabel("Time (hours)")
# ylabel("Abundance (mM)")
# title("Glutamine Profile")
# tight_layout()
# legend()
# # savefig("gln.svg")
# # savefig("gln.pdf")
# # savefig("gln.png")
#
# figure(7)
# upper_glu =  glu_mean' + glu_std
# lower_glu =  glu_mean' - glu_std
# ub_glu = upper_glu[:,1]
# lb_glu = lower_glu[:,1]
# fill_between(time, ub_glu, lb_glu,color="lightblue", alpha=1)
# plot(time,glu_mean',color="black", label = "Simulated")
# # scatter(tGlu,Glu,color="black", label = "Experimental")
# errorbar(tGlu,Glu,yerr = 0.2*Glu, fmt = "o", color="black", label = "Experimental")
# xlabel("Time (hours)")
# ylabel("Abundance (mM)")
# title("Glutamate Profile")
# tight_layout()
# legend()
# # savefig("glu.svg")
# # savefig("glu.pdf")
# # savefig("glu.png")



# prob = ODEProblem(odefunc,x0,tspan,P)
# #
# # # sol = solve(prob,CVODE_BDF())
# #
# sol = solve(prob, saveat=0.1)
# #
# # #f(t,x) = odefunc(x,t,P)
# # #t,X = ode45(f,x0,tSim; points=:specified)
# #
# # # , abstol = 1e-3, reltol = 1e-3, dt = 0.1)


# P = DataDict["Parameters"]
#
# prob = ODEProblem(odefunc,x0,tspan,P)
# sol = solve(prob)
#
# # # solutions for three mode case
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
# #
# #
# using PyPlot
# figure(1)
# plot(t,xbio,color="black")
# scatter(tBio,Bio,color="black")
# # #
# figure(2)
# plot(t,xglc,color="black")
# scatter(tGlc,Glc,color="black")
# # # #
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
# # #
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
