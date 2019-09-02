using DelimitedFiles

function Dynamic_Cybernetic(tStart,tStop,tStep)
    #============================#
    #EnzymeRateParameters :PARA with orig ksats
    # Parameters = [
    # 10;       #4.3;  # 1 ke :8
    # 0.005;     #0.02852; # 2 alpha:10
    # 0.5;  # 3 beta 2.772:1
    # 0.1; # 4 kmax 1
    # 1; # 5 serine Ksat1 5.622495495, 9.3410709178
    # 5; # 6 glucose Ksat2 41.626579913, 60.1520473825
    # 0.3; # 7 asparagine Ksat3 7.5179583177,18.1220565684
    # 0.4; # 8 glutamine Ksat4 0.4730679356, 1.9095941421
    # 0.07; # 9 kmax 2
    # 1; # 10 SERINE Ksat 5, 5.6224295495
    # 20; # 11 glc Ksat6, 42.626579913
    # 0.3; # 12 asn Ksat7, 7.5179583177
    # 4; # 13 lac Ksat8, 32.5306831612
    # 3; # 14 gly Ksat9, 3.5989039649
    # 0.6; # 15 nh3 Ksat10, 2.4196843244
    # 0.08; # 16 kmax 3
    # 1; # 17 ser Ksat11, 3.5166488962
    # 20; # 18 glc Ksat12, 33.2811185771
    # 0.3; # 19 asn Ksat13, 1.0632485079
    # 4; # 20 lac Ksat14, 2.0068776236,4
    # 2.5; # 21 ala Ksat15, 7.5419168128
    # 1; # 22 asp Ksat16, 3.9821819505
    # 3; #23GLY Ksat 17
    # 10; # 24 ke2
    # 100; # 25 ke3
    # 0.5; # beta2
    # 0.2; # BETA3
    # ]

# with  ; bio,gln,glu,asn,ser
    # Parameters = [
    # 1.8338544888611676;       #4.3;  # 1 ke :8
    # 1.6498055294467309;  # 2 beta 2.772:1
    # 0.4783811481891599; # 3 kmax 1
    # 0.005853514302477346; # 4 kmax 2
    # 0.15898948279771452; # 5 kmax 3
    # 114.11006222271486; # 6 ke2
    # 188.75020093588333; # 7 ke3
    # 5.282352057003303; # 8 beta2
    # 0.056542859863806895; # 9 BETA3
    # ]

    # # with  ; bio,lac, gln,glu,asn,ser,glc (VALUES FROM PRIOR)
    # Parameters = [
    # 1.49409946237286;       #4.3;  # 1 ke :8
    # 0.3029936637171954;  # 2 beta 2.772:1
    # 0.2808213926665067; # 3 kmax 1
    # 0.0037740788330917607; # 4 kmax 2
    # 0.08844195630181774; # 5 kmax 3
    # 353.5130954653908; # 6 ke2
    # 235.51381500111873; # 7 ke3
    # 8.003753673610975; # 8 beta2
    # 0.08349054054595537; # 9 BETA3
    # 95;                 # 10 s_o feed substrate
    # ]

    # # with  ; bio,lac, gln,glu,asn,ser,glc (TRIAL)
    # Parameters = [
    # 1;       #4.3;  # 1 ke :8
    # 0.05;  # 2 beta 2.772:1
    # 0.85; # 3 kmax 1
    # 0.75; # 4 kmax 2
    # 0.9; # 5 kmax 3
    # 1; # 6 ke2
    # 1; # 7 ke3
    # 0.05; # 8 beta2
    # 0.05; # 9 BETA3
    # 0.95;                 # 10 s_o feed substrate
    # ]

    # # ( I am getting a good range now)
    # Parameters = [
    # 1;       #4.3;  # 1 ke :8
    # 0.05;  # 2 beta 2.772:1
    # 0.85; # 3 kmax 1
    # 0.75; # 4 kmax 2
    # 0.9; # 5 kmax 3
    # 1; # 6 ke2
    # 1; # 7 ke3
    # 0.05; # 8 beta2
    # 0.05; # 9 BETA3
    # 9.5;                 # 10 s_o feed substrate
    # ]
# another okayish range set
    # Parameters = [
    # 1;       #4.3;  # 1 ke :8
    # 0.05;  # 2 beta 2.772:1
    # 0.031; # 3 kmax 1
    # 0.035; # 4 kmax 2
    # 0.038; # 5 kmax 3
    # 1; # 6 ke2
    # 1; # 7 ke3
    # 0.05; # 8 beta2
    # 0.05; # 9 BETA3
    # 9.5;                 # 10 s_o feed substrate
    # ]

    Parameters = [
    1;       #4.3;  # 1 ke :8
    0.05*5;  # 2 beta 2.772:1
    0.528; # 3 kmax 1 0.085
    0.288; # 4 kmax 2 0.075
    0.72; # 5 kmax 3 0.09
    10; # 6 ke2
    30; # 7 ke3
    0.05*5; # 8 beta2
    0.05*5; # 9 BETA3
    9.5;                 # 10 s_o feed substrate
    ]

    # Parameters = [
    # rand()*10;       #4.3;  # 1 ke :8
    # rand();  # 2 beta 2.772:1
    # rand(); # 3 kmax 1
    # rand()*0.001; # 4 kmax 2
    # rand()*0.01; # 5 kmax 3
    # rand()*1000; # 6 ke2
    # rand()*1000; # 7 ke3
    # rand()*10; # 8 beta2
    # rand()*0.1; # 9 BETA3
    # ]

    # # trying to get a physiologically relevant zone 1
    # Parameters = [
    # 0.21574664056400064;   #*(0.1);       #4.3;  # 1 ke :8
    # 0.1599796972464626;  # 2 beta 2.772:1
    # 0.8680410016070501; # 3 kmax 1
    # 0.47031784188947884;  #*(10); # 4 kmax 2
    # 3.0628958429107875;   #*(10); # 5 kmax 3
    # 4.866375755781506; # 6 ke2
    # 3.304751403666941; # 7 ke3
    # 8.96356361382608; # 8 beta2
    # 0.09385989453007569; # 9 BETA3
    # ]

    # # trying to get a physiologically relevant zone 2
    # Parameters = [
    # 0.006852986817764412;   #*(0.1);       #4.3;  # 1 ke :8
    # 0.04617018967805757;  # 2 beta 2.772:1
    # 1.2441132837016362; # 3 kmax 1
    # 0.10239008413352171;  #*(10); # 4 kmax 2
    # 10.6373352594878;   #*(10); # 5 kmax 3
    # 0.04621805330484753; # 6 ke2
    # 0.015886769983009662; # 7 ke3
    # 0.02144945464935608; # 8 beta2
    # 0.041219275514860265; # 9 BETA3
    # ]

    Param_Ensemble = zeros(10,10)
    for i=1:10
        # Param_Ensemble[:,i] = abs.(Parameters .+ 0.1*randn())
        Param_Ensemble[1,i] = Parameters[1] + randn();
        Param_Ensemble[2,i] = Parameters[2] + randn();
        Param_Ensemble[3,i] = Parameters[3] + randn(); #abs.(Parameters[3] .+ 0.001*randn())
        Param_Ensemble[4,i] = Parameters[4] + randn();  #abs.(Parameters[4] .+ 0.001*randn())
        Param_Ensemble[5,i] = Parameters[5] + randn();    #abs.(Parameters[5] .+ 0.001*randn())
        Param_Ensemble[6,i] = Parameters[6] + randn();
        Param_Ensemble[7,i] = Parameters[7] + randn();
        Param_Ensemble[8,i] = Parameters[8] + randn();
        Param_Ensemble[9,i] = Parameters[9] + randn();
        Param_Ensemble[10,i] = Parameters[10] + randn();
    end

    # Param_Ensemble = zeros(9,100)
    # for i=1:100
    #     # Param_Ensemble[:,i] = abs.(Parameters .+ 0.1*randn())
    #     Param_Ensemble[1,i] = 0.55*Parameters[1] .+ 0.9*Parameters[1].*rand()
    #     Param_Ensemble[2,i] = 0.55*Parameters[2] .+ 0.9*Parameters[2].*rand()
    #     Param_Ensemble[3,i] = 0.55*Parameters[3] .+ 0.9*Parameters[3].*rand() #abs.(Parameters[3] .+ 0.001*randn())
    #     Param_Ensemble[4,i] = 0.55*Parameters[4] .+ 0.9*Parameters[4].*rand()  #abs.(Parameters[4] .+ 0.001*randn())
    #     Param_Ensemble[5,i] = 0.55*Parameters[5] .+ 0.9*Parameters[5].*rand()    #abs.(Parameters[5] .+ 0.001*randn())
    #     # Param_Ensemble[3,i] = Parameters[3] .+ Parameters[3].*rand() #abs.(Parameters[3] .+ 0.001*randn())
    #     # Param_Ensemble[4,i] = Parameters[4] .+ Parameters[4].*rand()  #abs.(Parameters[4] .+ 0.001*randn())
    #     # Param_Ensemble[5,i] = Parameters[5] .+ Parameters[5].*rand()    #abs.(Parameters[5] .+ 0.001*randn())
    #     Param_Ensemble[6,i] = 0.55*Parameters[6] .+ 0.9*Parameters[6].*rand()
    #     Param_Ensemble[7,i] = 0.55*Parameters[7] .+ 0.9*Parameters[7].*rand()
    #     Param_Ensemble[8,i] = 0.55*Parameters[8] .+ 0.9*Parameters[8].*rand()
    #     Param_Ensemble[9,i] = 0.55*Parameters[9] .+ 0.9*Parameters[9].*rand()
    # end

    # with  ; TO make physiological sense
    # Parameters = [
    # 3.486439616717923*(0.1);       #4.3;  # 1 ke :8
    # 27;  # 2 beta 2.772:1s
    # 0.2552327014409511*(100); # 3 kmax 1
    # 12.657023044757516*(100); # 4 kmax 2
    # 1.0222349390305727*(100); # 5 kmax 3
    # 0.4250644012335958*(0.01); # 6 ke2
    # 0.010656867663803474*(0.01); # 7 ke3
    # 27; # 8 beta2
    # 27; # 9 BETA3
    # ]

    #EnzymeRateParameters : back to square one ::TRIED JULY 8
    # Parameters = [
    # 10.9;       #4.3;  # 1 ke
    # 0.005;     #0.02852; # 2 alpha
    # 0.05;  # 3 beta 2.772
    # 0.05; # 4 kmax 1
    # 1; # 5 serine Ksat1 5.622495495, 9.3410709178
    # 5; # 6 glucose Ksat2 41.626579913, 60.1520473825
    # 1.3; # 7 asparagine Ksat3 7.5179583177,18.1220565684
    # 0.1; # 8 glutamine Ksat4 0.4730679356, 1.9095941421
    # 0.04; # 9 kmax 2
    # 1; # 10 SERINE Ksat 5, 5.6224295495
    # 5; # 11 glc Ksat6, 42.626579913
    # 0.3; # 12 asn Ksat7, 7.5179583177
    # 4; # 13 lac Ksat8, 32.5306831612
    # 3; # 14 gly Ksat9, 3.5989039649
    # 0.7; # 15 nh3 Ksat10, 2.4196843244
    # 0.01; # 16 kmax 3
    # 1; # 17 ser Ksat11, 3.5166488962
    # 5; # 18 glc Ksat12, 33.2811185771
    # 0.3; # 19 asn Ksat13, 1.0632485079
    # 4; # 20 lac Ksat14, 2.0068776236,4
    # 2.5; # 21 ala Ksat15, 7.5419168128
    # 1; # 22 asp Ksat16, 3.9821819505
    # 3; #GLY Ksat 17
    # ]
    # ke = 4.3; #2.852
    # # ke = 2.852
    # alpha = 0.02852
    # beta = 2.772
    # # kmax = 1209600
    # kmax = 0.099; #1.209600
    # K = 0.13
    #============================#
    # Lets load Z
    # mode_matrix = float(readdlm("Data/dynamic_mode1.csv",','))[2649:3194,:];

    mode_matrix = float(readdlm("Data/Mode_Analysis_FBA.csv",','))[2649:3194,:];

    # mode_matrix = float(readdlm("Data_july8/multi_mode.csv",','))[2649:3194,:];

    Z = deepcopy(mode_matrix)


    # normalisation on the basis of biomass flux
    a = 5.3108538979; # the values after 4690 are not normalized since it just becomes 1. make this correction in ecoli too
    b = 4.3943091593;
    c = 2.7519335312;

    # # gln
    # # Z[271,1] = -4.2
    # Z[271,1] = -5   #-9 #-29
    # Z[271,2] = 0.08  #1  #50
    # Z[271,3] = 0.01   #1    #70
    #
    #
    # # nh3
    # Z[396,1] = 4    #6   #8
    # Z[396,2] = -1.1   #-5   #-4   #-67
    # Z[396,3] = -0.4   #-0.5   #0.001  #0.1
    #
    # # GLUTAMATE
    # Z[274,1] = 2   #3    #0.5
    # Z[274,2] = 0.5    #0.5
    # Z[274,3] = 0.097  #0.1 #0.097
    #
    # # asparate
    # Z[108,1] = 1    #0.8     #2.5    #2.5
    # Z[108,2] = 2    #2.5
    # Z[108,3] = -0.2  #-0.009
    #
    # #asparagine
    # Z[106,1] = -10   #-15   #-42.5
    # Z[106,2] = -4.5
    # Z[106,3] = -1.1    #-0.97
    #
    # #glucose
    # Z[268,1] = -45
    # Z[268,2] = -0.01  #-19
    # Z[268,3] = -0.2   #-1.5
    #
    #
    # #serine
    # Z[475,1] = -5     #-9  #-5   #-15
    # Z[475,2] = -1   #-0.1    #-17
    # Z[475,3] = -0.5    #1     #10
    #
    # # #glycine
    # Z[276,1] = 0.7
    # Z[276,2] = -0.15
    # Z[276,3] = -0.2
    #
    # # # #alanine
    # # Z[89,1] = 60
    # # Z[89,1] = 4
    # # Z[89,3] = 1
    #
    # #LACTATE
    # Z[356,1] =  40   #60     #150   #80    #90   #200
    # Z[356,2] = -0.1 #-1   #-100    #150      #890  #-9
    # Z[356,3] = -5    #-9   #-1#3



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
    S_idx[2,:] = stm_matrix[320,2649:3194]  #*-1 #lac is produced :please fix this after fixing the modes. shouldn't require the negative sign
    S_idx[3,:] = stm_matrix[199,2649:3194] # ser
    S_idx[4,:] = stm_matrix[193,2649:3194] # asparagine
    S_idx[5,:] = stm_matrix[197,2649:3194] # glutamine
    S_idx[6,:] = stm_matrix[249,2649:3194] # glutamate
    S_idx[7,:] = stm_matrix[170,2649:3194] # alanine
    S_idx[8,:] = stm_matrix[211,2649:3194] # asparate
    S_idx[9,:] = stm_matrix[247,2649:3194] # glycine
    S_idx[10,:] = stm_matrix[319,2649:3194] # nh3
    S = abs.(S_idx)


    #Initial Conditions for 3 modes
    x0 = [
    0.8; # e1
    0.8; # e2
    0.8; # e3
    70.0410228873; #1 glc
    11.5482629792; #2 lac
    10.7144927914; # 3 serine
    19.7833200136; # 4 asparagine
    3.3026417615; # 5 glutamine
    0.8328368955; # 6 glutamate
    0.8250474968; # 7 alanine
    2.1248193993; # 8 asparate
    3.401370313; # 9 glycine
    1.1260873627; # 10 nh3
    5.7144726262; # 11 m_BIO
    ];

    # IC_Ensemble = zeros(14,100)
    # for i=1:100
    #     IC_Ensemble[1,i] = x0[1]*rand()
    #     IC_Ensemble[2,i] = x0[2]*rand()
    #     IC_Ensemble[3,i] = x0[3]*rand()
    #     IC_Ensemble[4,i] = x0[4]
    #     IC_Ensemble[5,i] = x0[5]
    #     IC_Ensemble[6,i] = x0[6]
    #     IC_Ensemble[7,i] = x0[7]
    #     IC_Ensemble[8,i] = x0[8]
    #     IC_Ensemble[9,i] = x0[9]
    #     IC_Ensemble[10,i] = x0[10]
    #     IC_Ensemble[11,i] = x0[11]
    #     IC_Ensemble[12,i] = x0[12]
    #     IC_Ensemble[13,i] = x0[13]
    #     IC_Ensemble[14,i] = x0[14]
    # end

    #Initial Conditions for mode 1
    # x0 = [
    # 2; # e1
    # # 2; # e2
    # # 2; # e3
    # 70.0410228873; #1 glc
    # 11.5482629792; #2 lac
    # 10.7144927914; # 3 serine
    # 19.7833200136; # 4 asparagine
    # 3.3026417615; # 5 glutamine
    # 0.8328368955; # 6 glutamate
    # 0.8250474968; # 7 alanine
    # 2.1248193993; # 8 asparate
    # 3.401370313; # 9 glycine
    # 1.1260873627; # 10 nh3
    # 5.7144726262; # 11 m_BIO
    # ];

    # #Initial Conditions for mode 2
    # x0 = [
    # 2; # e1
    # # 2; # e2
    # # 2; # e3
    # 42.9112388389; #1 glc
    # 36.9591961968; #2 lac
    # 6.356446242; # 3 serine
    # 11.357176734; # 4 asparagine
    # 0.218996417; # 5 glutamine
    # 2.8158305383; # 6 glutamate
    # 3.9531585163; # 7 alanine
    # 3.2370324214; # 8 asparate
    # 3.8802531206; # 9 glycine
    # 3.019299547; # 10 nh3
    # 38.9878032741; # 11 m_BIO
    # ];

    # #Initial Conditions for mode 3
    # x0 = [
    # 2; # e1
    # # 2; # e2
    # # 2; # e3
    # 39.7661983476; #1 glc
    # 13.1866142429; #2 lac
    # 4.2923371416; # 3 serine
    # 1.9436857338; # 4 asparagine
    # 0.8914759815; # 5 glutamine
    # 3.8338963836; # 6 glutamate
    # 9.4996072089; # 7 alanine
    # 4.5508469918; # 8 asparate
    # 3.3343678262; # 9 glycine
    # 1.5948100864; # 10 nh3
    # 68.5797854082; # 11 m_BIO
    # ];

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
