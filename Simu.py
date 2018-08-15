#!/usr/bin/env python
# -*-coding:utf-8 -*-
# ** author:Oasis
# *****************************
import msprime, sys
from math import log
from math import exp

seed = 1
#nhaps = 20002000
#nhaps = 8000
nhaps = int(sys.argv[1]) * 2
scale = 100
nbp = 50000000
mu=1.25e-8/scale # mutation rate per bp
rho=1e-8/scale   # recombination rate per bp
genlen=25
N0=7310
Thum=5920
Naf=14474

# years per generation
# initial population size, and reference effective size
# =148000/genlen; generations back to initial pop growth (advent of modern humans) # size of african population from t0 onwards
Tooa=2040
Nb=1861
mafb=1.5e-4 # migration rate between Africa and Out-of-Africa
#51000/genlen # number of generations back to Out of Africa # size of out of Africa population
Teu=920 #23000/genlen # number of generations back to Asia-Europe split
Neu=1032
Nas=554 # bottleneck population sizes after the split
mafeu=2.5e-5; mafas=7.8e-6; meuas=3.11e-5 # migration rates between Africa, Europe and Asia ras=0.0048 # growth rates per generation in Asia
reu=0.00292
Tex=150 # generations back to accelerated population growth rex=0.03 # accelerated growth rate
Tmod=15
rmod=0.13 # growth rate in most recent generations
# pop0 is Africa, pop1 is eur (uk), pop2 is east asn
# will generate nhaps from eur pop
samplesize = nhaps
othersize = 0
rex=0.03
ras=0.0048

pop_config = [
    msprime.PopulationConfiguration(sample_size=othersize,initial_size=Naf*exp(rex*Tex),\
                                growth_rate= rex),
    msprime.PopulationConfiguration(sample_size=samplesize,initial_size= Neu*exp(reu*(Teu-Tex))* exp(rex*(Tex-Tmod))* exp(rmod*Tmod), \
                                    growth_rate=rmod),
    msprime.PopulationConfiguration(sample_size=othersize,\
                                initial_size= Nas*exp((Teu-Tex)*ras)* exp(rex*Tex),growth_rate=rex)
]

mig_mat = [[0,mafeu,mafas],[mafeu,0,meuas],[mafas,meuas,0]]

# recent change in growth rate
recent_event = [ msprime.PopulationParametersChange(time=Tmod,growth_rate=rex,population_id=1) ]

# populations stop having accelerated growth (advent of agriculture)
ag_event = [
    msprime.PopulationParametersChange(time=Tex,growth_rate=ras,population_id=2),\
    msprime.PopulationParametersChange(time=Tex,growth_rate=reu,population_id=1), \
    msprime.PopulationParametersChange(time=Tex,growth_rate=0.0,population_id=0)
]

# Asia and Europe merge, migration changes, population size changes, growth stops
eu_event = [
    msprime.MigrationRateChange(time=Teu,rate=0.0), msprime.PopulationParametersChange(time=Teu,growth_rate=0.0, population_id=2),\
    msprime.MassMigration(time=Teu+0.0001,source=2,destination=1,proportion=1.0), \
    msprime.PopulationParametersChange(time=Teu+0.0001,initial_size=Nb,growth_rate=0.0, population_id=1), \
    msprime.MigrationRateChange(time=Teu+0.0001,rate=mafb,matrix_index=(0,1)), \
    msprime.MigrationRateChange(time=Teu+0.0001,rate=mafb,matrix_index=(1,0))
]

# Out of Africa event (looking back, Africa and Europe merge)
ooa_event = [
    msprime.MigrationRateChange(time=Tooa,rate=0.0),\
    msprime.MassMigration(time=Tooa+0.0001,source=1,destination=0,proportion=1.0)
]

# initial population size
init_event = [
    msprime.PopulationParametersChange(time=Thum,initial_size=N0,population_id=0)
]

# cat all the events together
events = recent_event + ag_event + eu_event + ooa_event + init_event

# run the simulation
treeseq = msprime.simulate(population_configurations=pop_config, migration_matrix=mig_mat, \
                           demographic_events=events, length=nbp, recombination_rate=rho, \
                           mutation_rate=mu, random_seed=seed)

# # debug log
# dd = msprime.DemographyDebugger(
#         population_configurations=pop_config,
#         migration_matrix=mig_mat,
#         demographic_events=events)
# dd.print_history()

# #vcffile = '/Users/yeweijian/Documents/MsprimeSimul/t.vcf'
with sys.stdout as vcffile:
    treeseq.write_vcf(vcffile, 2)