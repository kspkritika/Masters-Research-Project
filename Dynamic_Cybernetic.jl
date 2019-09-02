using DelimitedFiles

function Dynamic_Cybernetic(tStart,tStop,tStep)

    Parameters = [
    1;       
    0.05*5;  
    0.528; 
    0.288; 
    0.72; 
    10;
    30;
    0.05*5;
    0.05*5; 
    9.5;                 
    ]

    Param_Ensemble = zeros(10,10)
    for i=1:10
        # Param_Ensemble[:,i] = abs.(Parameters .+ 0.1*randn())
        Param_Ensemble[1,i] = Parameters[1] + randn();
        Param_Ensemble[2,i] = Parameters[2] + randn();
        Param_Ensemble[3,i] = Parameters[3] + randn(); 
        Param_Ensemble[4,i] = Parameters[4] + randn();  
        Param_Ensemble[5,i] = Parameters[5] + randn();    
        Param_Ensemble[6,i] = Parameters[6] + randn();
        Param_Ensemble[7,i] = Parameters[7] + randn();
        Param_Ensemble[8,i] = Parameters[8] + randn();
        Param_Ensemble[9,i] = Parameters[9] + randn();
        Param_Ensemble[10,i] = Parameters[10] + randn();
    end

    #============================#

    mode_matrix = float(readdlm("Data/Mode_Analysis_FBA.csv",','))[2649:3194,:];

    Z = deepcopy(mode_matrix)


    # normalisation on the basis of biomass flux
    a = 5.3108538979; 
    b = 4.3943091593;
    c = 2.7519335312;

    for j=1:546
        Z[j,1] = Z[j,1]/a
    end

    for j=1:546
        Z[j,2] = Z[j,2]/b
    end

    for j=1:546
        Z[j,3] = Z[j,3]/c
    end

    #=============================#
    #Load Stoichiometric Matrix
    stm_matrix = float(readdlm("Data/Network.csv",','));
    S_idx = zeros(10,546)
    S_idx[1,:] = stm_matrix[293,2649:3194]#glc
    S_idx[2,:] = stm_matrix[320,2649:3194]  #lac
    S_idx[3,:] = stm_matrix[199,2649:3194] # ser
    S_idx[4,:] = stm_matrix[193,2649:3194] # asparagine
    S_idx[5,:] = stm_matrix[197,2649:3194] # glutamine
    S_idx[6,:] = stm_matrix[249,2649:3194] # glutamate
    S_idx[7,:] = stm_matrix[170,2649:3194] # alanine
    S_idx[8,:] = stm_matrix[211,2649:3194] # asparate
    S_idx[9,:] = stm_matrix[247,2649:3194] # glycine
    S_idx[10,:] = stm_matrix[319,2649:3194] # nh3
    S = abs.(S_idx)

#=======================================#
#=======================================#
    # Parameters and initial conditions from Data_dict
    Data_dict = Dict()

    Data_dict["Parameters"] = Parameters
    Data_dict["InitialConditions"] = x0
    Data_dict["ModeMatrix"] = Z
    Data_dict["Stoich"] = S
    Data_dict["Parameter_Ensemble"] = Param_Ensemble
    # Data_dict["IC_Ensemble"] = IC_Ensemble

#=========================================#
  return Data_dict
end
