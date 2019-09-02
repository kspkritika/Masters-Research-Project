using PyCall
@pyimport numpy as np

function Objective(P)

    prob = ODEProblem(odefunc,x0,tspan,P)
    # prob = ODEProblem(odefunc,IC,tspan,P)
    sol = solve(prob)

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

    glc_sim = np.interp(tGlc,t,xglc)
    lac_sim = np.interp(tLac,t,xlac)
    ser_sim = np.interp(tSer,t,xser)
    asn_sim = np.interp(tAsn,t,xasn)
    gln_sim = np.interp(tGln,t,xgln)
    glu_sim = np.interp(tGlu,t,xglu)
    ala_sim = np.interp(tAla,t,xala)
    asp_sim = np.interp(tAsp,t,xasp)
    gly_sim = np.interp(tGly,t,xgly)
    nh3_sim = np.interp(tnh3,t,xnh3)
    bio_sim = np.interp(tBio,t,xbio)

    z1 = sum((glc_sim-Glc).^2)
    z2 = sum((lac_sim-Lac).^2)
    z3 = sum((ser_sim-Ser).^2)
    z4 = sum((asn_sim-Asn).^2)
    z5 = sum((gln_sim-Gln).^2)
    z6 = sum((glu_sim-Glu).^2)
    z7 = sum((ala_sim-Ala).^2)
    z8 = sum((asp_sim-Asp).^2)
    z9 = sum((gly_sim-Gly).^2)
    z10 = sum((nh3_sim-nh3).^2)
    z11 = sum((bio_sim-Bio).^2)

        cost = z1 + z2 + z3 + z4 + z5 +  z11 + z6 + z7 + z8 + z9 + z10 
    return cost
end
