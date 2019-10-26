using DelimitedFiles

function Dynamic_Cybernetic(tStart,tStop,tStep)

#FED BATCH PARAMETERS. These have been estimated by hand, trial and error and with the help of optim
    Parameters = [
    1;      # 1 ke
    0.05*5;  # 2 beta
    0.528; # 3 kmax1
    0.288; # 4 kmax2
    0.72; # 5 kmax3
    10; # 6 ke2
    30; # 7 ke3
    0.05*5; # 8 beta2
    0.05*5; # 9 BETA3
    9.5;                 # 10 s_o feed substrate, this value has been taken from a previous example on fed batch cybernetic model for hybridoma cells
    ]


    Param_Ensemble = zeros(10,10)
    for i=1:10
        # Param_Ensemble[:,i] = abs.(Parameters .+ 0.1*randn())
        Param_Ensemble[1,i] = abs.(Parameters[1] .+ 0.99*randn());
        Param_Ensemble[2,i] = abs.(Parameters[2] .+ 0.99*randn());
        Param_Ensemble[3,i] = abs.(Parameters[3] .+ 0.99*randn()); #abs.(Parameters[3] .+ 0.001*randn())
        Param_Ensemble[4,i] = abs.(Parameters[4] + 0.99*randn());  #abs.(Parameters[4] .+ 0.001*randn())
        Param_Ensemble[5,i] = abs.(Parameters[5] + 0.99*randn());    #abs.(Parameters[5] .+ 0.001*randn())
        Param_Ensemble[6,i] = abs.(Parameters[6] + 0.99*randn());
        Param_Ensemble[7,i] = abs.(Parameters[7] + 0.99*randn());
        Param_Ensemble[8,i] = abs.(Parameters[8] + 0.99*randn());
        Param_Ensemble[9,i] = abs.(Parameters[9] + 0.99*randn());
        Param_Ensemble[10,i] = abs.(Parameters[10] + 0.99*randn());
    end

    #============================#
    # Lets load Z
    mode_matrix = float(readdlm("Data/FBA_MODES_old_new.csv",','));
    # mode_matrix = float(readdlm("Data/MCMC_MODES_old_new.csv",','));

    Z = deepcopy(mode_matrix)


# FED BATCH
    # normalisation on the basis of biomass flux
    a = 5.3108538979; # the values after 4690 are not normalized since it just becomes 1. make this correction in ecoli too
    b = 4.3943091593;
    c = 2.7519335312;

    for j=1:4725
        Z[j,1] = Z[j,1]/a
    end

    for j=1:4725
        Z[j,2] = Z[j,2]/b
    end

    for j=1:4725
        Z[j,3] = Z[j,3]/c
    end

    #=============================#
    #Load Stoichiometric Matrix
    #antibody reaction has been added to a reaction network of CHO-K1 cell line given by Hefzi et al. This reaction was taken from Nolan et al
    stm_matrix = float(readdlm("Data/Network.csv",','));
    S_idx = zeros(11,4725)
    S_idx[1,:] = stm_matrix[293,:]#glc
    S_idx[2,:] = stm_matrix[320,:]  #*-1 #lac is produced :please fix this after fixing the modes. shouldn't require the negative sign
    S_idx[3,:] = stm_matrix[199,:] # ser
    S_idx[4,:] = stm_matrix[193,:] # asparagine
    S_idx[5,:] = stm_matrix[197,:] # glutamine
    S_idx[6,:] = stm_matrix[249,:] # glutamate
    S_idx[7,:] = stm_matrix[170,:] # alanine
    S_idx[8,:] = stm_matrix[211,:] # asparate
    S_idx[9,:] = stm_matrix[247,:] # glycine
    S_idx[10,:] = stm_matrix[319,:] # nh3
    S_idx[11,:] = stm_matrix[2775,:] # antibody
    S = abs.(S_idx)


    #Initial Conditions for 3 modes
    x0 = [
    0.8; # e1
    0.8; # e2
    0.8; # e3
    70.0410228873; #4 glc
    11.5482629792; #5 lac
    10.7144927914; # 6 serine
    19.7833200136; # 7 asparagine
    3.3026417615; # 8 glutamine
    0.8328368955; # 9 glutamate
    0.8250474968; # 10 alanine
    2.1248193993; # 11 asparate
    3.401370313; # 12 glycine
    1.1260873627; # 13 nh3
    0.6173367204; # 14 antibody
    6.170588591; # 15 vcd
    5.7144726262; # 16 m_BIO overall
    ];


#=======================================#
#=======================================#
    # Parameters and initial conditions from Data_dict
    Data_dict = Dict()

    Data_dict["Parameters"] = Parameters
    Data_dict["InitialConditions"] = x0
    Data_dict["ModeMatrix"] = Z
    Data_dict["Stoich"] = S
    Data_dict["Parameter_Ensemble"] = Param_Ensemble

#=========================================#
  return Data_dict
end
