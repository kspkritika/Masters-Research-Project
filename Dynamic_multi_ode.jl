using DelimitedFiles
using DifferentialEquations
using Statistics
# using Sundials
# using ODE
include("loadData.jl")

function odefunc(dxdt,x,P,t)

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
    alpha = 0.005; 
    beta1 = P[2]; 
    kmax1 = P[3]
    Kser1 = 1; 
    Kglc1 = 5; 
    Kasn1 = 0.3;
    Kgln1 = 0.4; 
    kmax2 = P[4]
    Kser2 = 1; 
    Kglc2 = 20; 
    Kasn2 = 0.3; 
    Klac2 = 4; 
    Kgly2 = 3; 
    Knh2 = 0.6; 
    kmax3 = P[5]
    Kser3 = 1; 
    Kglc3 = 20; 
    Kasn3 = 0.3;
    Klac3 = 4; 
    Kala3 = 2.5;
    Kasp3 = 1; 
    Kgly3 = 3; 
    ke2 = P[6]
    ke3 = P[7]
    beta2 = P[8]; 
    beta3 = P[9]; 
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

    #EnzymeReactionRate
    fill!(rE,0.0)
    push!(rE,ke1*(M_glc/(Kglc1+M_glc))*(M_ser/(Kser1+M_ser))*(M_asn/(Kasn1+M_asn))*(M_gln/(Kgln1+M_gln)))
    push!(rE,ke2*(M_glc/(Kglc2+M_glc))*(M_ser/(Kser2+M_ser))*(M_asn/(Kasn2+M_asn))*(M_lac/(Klac2+M_lac))*(M_gly/(Kgly2+M_gly))*(M_nh3/(Knh2+M_nh3)))
    push!(rE,ke3*(M_glc/(Kglc3+M_glc))*(M_ser/(Kser3+M_ser))*(M_asn/(Kasn3+M_asn))*(M_lac/(Klac3+M_lac))*(M_ala/(Kala3+M_ala))*(M_asp/(Kasp3+M_asp))*(M_gly/(Kgly3+M_gly)))

    # Growthrate
    fill!(rG,0.0)
    push!(rG,(5.3108538979*rM[1]))
    push!(rG,(4.3943091593*rM[2]))
    push!(rG,(2.7519335312*rM[3]))

    #Initialize Control Vector
    cybernetic_var = Float64[]
    fill!(cybernetic_var,0.0)

    push!(cybernetic_var,(abs.(Z[271,1]))*rM[1])
    push!(cybernetic_var,(abs.(Z[106,2]))*rM[2])
    push!(cybernetic_var,(abs.(Z[106,3]+Z[356,3]))*rM[3])

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







