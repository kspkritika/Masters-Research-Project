using DelimitedFiles

AntiData = readdlm("Experimental_Data/anti.csv",',')
tAnti = AntiData[:,1]
Anti = AntiData[:,2]

vcdData = readdlm("Experimental_Data/vcd.csv",',')
tvcd = vcdData[:,1]
vcd = vcdData[:,2]
VCD = vcdData[:,3]  # IN mM

BioData = readdlm("Experimental_Data/Biomass.csv",',')
tBio = BioData[:,1]
Bio = BioData[:,2]

SerData = readdlm("Experimental_Data/ser.csv",',')
tSer = SerData[:,1]
Ser = SerData[:,2]

GluData = readdlm("Experimental_Data/glu.csv",',')
tGlu = GluData[:,1]
Glu = GluData[:,2]

GlcData = readdlm("Experimental_Data/gluc.csv",',')
tGlc = GlcData[:,1]
Glc = GlcData[:,2]

AsnData = readdlm("Experimental_Data/asn.csv",',')
tAsn = AsnData[:,1]
Asn = AsnData[:,2]

GlnData = readdlm("Experimental_Data/gln.csv",',')
tGln = GlnData[:,1]
Gln = GlnData[:,2]

LacData = readdlm("Experimental_Data/lac.csv",',')
tLac = LacData[:,1]
Lac = LacData[:,2]

AlaData = readdlm("Experimental_Data/ala.csv",',')
tAla = AlaData[:,1]
Ala = AlaData[:,2]

AspData = readdlm("Experimental_Data/asp.csv",',')
tAsp = AspData[:,1]
Asp = AspData[:,2]

GlyData = readdlm("Experimental_Data/gly.csv",',')
tGly = GlyData[:,1]
Gly = GlyData[:,2]

nh3Data = readdlm("Experimental_Data/nh3.csv",',')
tnh3 = nh3Data[:,1]
nh3 = nh3Data[:,2]
