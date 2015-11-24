#ARTIFICIAL SELECTION

# #
# sim16b_4allele_mono_1of2

# 2 Chromosomes, each wih 30 loci.  Selection at loci 05, 21 and/or 37.

# Most loci (initially) monomorphic.
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural or artificial selection and sampling.  Selection occuring at different periods.  Equal events spacing.

# 2 Initial subpopns, expanding to 6

# Selection at locus 5 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 5, this occurs before the split also
# 
# phylip, fasta, fstat formatted output
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM 	  
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 8.85e-7 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 8.85e-8# Backward Mutation Rate
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
## Settings ###################################################################################

for iter in range(10):
	z1 = 500   # Number of Steps (Generations) before population splitting
	j1 = 400   # Number of Steps (Generations) before selection event 1
	j2 = 550   # Number of Steps (Generations) before selection event 2
	j3 = 700   # Number of Steps (Generations) before selection event 3
	
	d1 = 50 # duration of selection event 1 (in generations)
	d2 = 50 # duration of selection event 2 (in generations)
	d3 = 50 # duration of selection event 3 (in generations)
       
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_000000_sample_%d.dat' (iter, pop.dvars().gen), gui=False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_000000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_000000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_000000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    # 		sim.DiscardIf(True, subPops=[
    # 	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )


    
for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_050000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_050000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_050000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_050000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_0021000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_002100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_002100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_002100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    # 		sim.DiscardIf(True, subPops=[
    # 	        (0, 'Affected')], at = j1),
          sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_000037_sample_%d.dat', gui = False),        
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_000037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_000037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_000037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    # 		sim.DiscardIf(True, subPops=[
    # 	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
          sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_052100_sample_%d.dat', gui = False),        
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_052100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_052100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_052100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_052137_sample_%d.dat', gui = False),        
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_052137_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_052137_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_052137_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_050037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_050037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_050037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_050037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_210500_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_210500_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_210500_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_210500_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250    # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_370005_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_370005_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_370005_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_370005_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 100   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_372105_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_372105_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_372105_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_372105_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 700   # Number of Steps (Generations) before selection event 2
    j3 = 550   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_053721_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_053721_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_053721_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_053721_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_370521_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_370521_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_370521_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_370521_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_realmu_%s_210537_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_realmu_%s_210537_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_realmu_%s_210537_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_210537_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_realmu_%s_210537_sample_%d.phy    sim16b_4allele_mono_1of2_realmu_%s_210537_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



#ARTIFICIAL SELECTION

# #
# sim16b_4allele_mono
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural or artificial selection and sampling.  Selection occuring at different periods.  Equal events spacing.

# 2 Initial subpopns, expanding to 6

# Selection at locus 5 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 5, this occurs before the split also
# 
# phylip, fasta, fstat formatted output
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM 	  
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 8.85e-7 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 8.85e-8# Backward Mutation Rate
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
## Settings ###################################################################################
 

for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
        
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_000000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_000000_sample_%d.dat' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_realmu_%s_000000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_realmu_%s_000000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    # 		sim.DiscardIf(True, subPops=[
    # 	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )


    
for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_050000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_050000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_realmu_%s_050000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_realmu_%s_050000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_002100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_002100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_realmu_%s_002100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_realmu_%s_002100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    # 		sim.DiscardIf(True, subPops=[
    # 	        (0, 'Affected')], at = j1),
          sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_000037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_000037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_realmu_%s_000037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_realmu_%s_000037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    # 		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    # 		sim.DiscardIf(True, subPops=[
    # 	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
          sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_052100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_052100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_realmu_%s_052100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_realmu_%s_052100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_052137_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_052137_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_realmu_%s_052137_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_realmu_%s_052137_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_050037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_050037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_realmu_%s_050037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_realmu_%s_050037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_210500_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_210500_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_realmu_%s_210500_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_realmu_%s_210500_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250    # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_370005_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_370005_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_realmu_%s_370005_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_realmu_%s_370005_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 100   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_372105_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_372105_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_realmu_%s_372105_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_realmu_%s_372105_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 700   # Number of Steps (Generations) before selection event 2
    j3 = 550   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_053721_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_053721_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_realmu_%s_053721_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_realmu_%s_053721_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
    	global q, g
    	rate=2**15
    	# if .... stop growing
    	if gen >= g and all([x < q  for x in pop.subPopSizes()]):
    		return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
    	elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
    		return ([q for x in pop.subPopSizes(ancGen=-1)])
    	else:
    		return ([x for x in pop.subPopSizes(ancGen=-1)])
    		
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_realmu_%s_370521_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_realmu_%s_370521_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_realmu_%s_370521_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_realmu_%s_370521_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
    	    sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
    	    sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
    	    sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
    	    sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
    	    sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
    	    sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
    	],
    	preOps=[
    		sim.SNPMutator(u = u, v = v),
    		sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
    		sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
    	    sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    		sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    		sim.DiscardIf(True, subPops=[
    	        (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
    	
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
    	                sim.MendelianGenoTransmitter(),
    	                sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
    	postOps=[
    		simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
    	gen=(t+1),
    )




#ARTIFICIAL SELECTION

# #
# sim16b_4allele_diphyletic
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural or artificial selection and sampling.  Selection occuring at different periods.  Equal events spacing.

# 2 Initial subpopns, expanding to 6

# Selection at locus 5 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 5, this occurs before the split also
# 
# phylip, fasta, fstat formatted output
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM     
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 8.85e-7 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 8.85e-8# Backward Mutation Rate
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
## Settings ###################################################################################
 

for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
        
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_000000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_000000_sample_%d.dat' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_000000_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_000000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_000000_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_000000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )


    
for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_050000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_050000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_050000_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_050000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_050000_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_050000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
                        sim.DiscardIf(True, subPops=[
                (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_002100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_002100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_002100_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_002100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_002100_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_002100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1),
          sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_000037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_000037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_000037_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_000037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_000037_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_000037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
          sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_052100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_052100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_052100_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_052100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_052100_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_052100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_052137_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_052137_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_052137_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_052137_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_052137_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_052137_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_050037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_050037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_050037_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_050037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_050037_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_050037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_210500_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_210500_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_210500_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_210500_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_210500_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_210500_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250    # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_370005_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_370005_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_370005_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_370005_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_370005_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_370005_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 100   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_372105_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_372105_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_372105_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_372105_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_372105_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_372105_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 700   # Number of Steps (Generations) before selection event 2
    j3 = 550   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_053721_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_053721_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_053721_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_053721_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_053721_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_053721_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_370521_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_370521_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_370521_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_370521_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_370521_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_370521_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_realmu_%s_210537_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_realmu_%s_210537_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_realmu_%s_210537_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_210537_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_realmu_%s_210537_sample_%d.phy    sim16b_4allele_diphyletic_realmu_%s_210537_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )




#ARTIFICIAL SELECTION

# #
# sim16b_4allele_diphyletic_sep
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural or artificial selection and sampling.  Selection occuring at different periods.  Equal events spacing.

# 2 Initial subpopns, expanding to 6

# Selection at locus 5 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 5, this occurs before the split also
# 
# phylip, fasta, fstat formatted output
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM     
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 8.85e-7 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 8.85e-8# Backward Mutation Rate
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
## Settings ###################################################################################
 

for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
        
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_000000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_000000_sample_%d.dat' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_000000_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_000000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_000000_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_000000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1+20),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )


    
for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_050000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_050000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_050000_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_050000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_050000_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_050000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
                        sim.DiscardIf(True, subPops=[
                (1, 'Affected')], at = j1+20),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_002100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_002100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_002100_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_002100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_002100_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_002100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1+20),
          sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_000037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_000037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_000037_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_000037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_000037_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_000037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1+20),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
          sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_052100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_052100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_052100_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_052100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_052100_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_052100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_052137_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_052137_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_052137_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_052137_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_052137_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_052137_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_050037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_050037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_050037_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_050037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_050037_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_050037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_210500_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_210500_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_210500_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_210500_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_210500_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_210500_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250    # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_370005_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_370005_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_370005_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_370005_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_370005_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_370005_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 100   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_372105_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_372105_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_372105_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_372105_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_372105_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_372105_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 700   # Number of Steps (Generations) before selection event 2
    j3 = 550   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_053721_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_053721_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_053721_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_053721_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_053721_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_053721_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_370521_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_370521_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_370521_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_370521_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_370521_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_370521_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_210537_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_sep_realmu_%s_210537_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_sep_realmu_%s_210537_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_210537_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_sep_realmu_%s_210537_sample_%d.phy    sim16b_4allele_diphyletic_sep_realmu_%s_210537_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1+20),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )




#NATURAL SELECTION

# #
# sim16b_4allele_mono_1of2
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural or artificial selection and sampling.  Selection occuring at different periods.  Equal events spacing.

# 2 Initial subpopns, expanding to 6

# Selection at locus 5 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 5, this occurs before the split also
# 
# phylip, fasta, fstat formatted output
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM     
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 8.85e-7 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 8.85e-8# Backward Mutation Rate
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
## Settings ###################################################################################
 

for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
        
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_000000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_000000_sample_%d.phy' (iter, pop.dvars().gen), gui=False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_000000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_000000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_000000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)),# Positive Selection for allele 1 at locus given
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )


    
for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_050000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_050000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_050000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_050000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaSelector(loci=5, fitness = [1,1-s1,1-s1], wildtype = 1, begin = j1, end = (j1+d1)), # Positive Selection for allele 1 at locus given
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_002100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_002100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_002100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_002100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
          sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_000037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_000037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_000037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_000037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
          sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_052100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_052100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_052100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_052100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_052137_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_052137_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_052137_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_052137_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_050037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_050037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_050037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_050037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_210500_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_210500_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_210500_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_210500_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250    # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_370005_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_370005_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_370005_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_370005_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 100   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_372105_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_372105_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_372105_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_372105_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 700   # Number of Steps (Generations) before selection event 2
    j3 = 550   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_053721_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_053721_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_053721_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_053721_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_370521_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_370521_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_370521_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_370521_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_210537_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_1of2_nat_realmu_%s_210537_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_1of2_nat_realmu_%s_210537_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_210537_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_1of2_nat_realmu_%s_210537_sample_%d.phy    sim16b_4allele_mono_1of2_nat_realmu_%s_210537_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )




#NATURAL SELECTION

# #
# sim16b_4allele_mono
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural or artificial selection and sampling.  Selection occuring at different periods.  Equal events spacing.

# 2 Initial subpopns, expanding to 6

# Selection at locus 5 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 5, this occurs before the split also
# 
# phylip, fasta, fstat formatted output
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM     
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 1   # Initial Number of Subpopns
d = 6 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 8.85e-7 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 8.85e-8# Backward Mutation Rate
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
## Settings ###################################################################################
 

for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
        
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_000000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_000000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_000000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_000000_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_000000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )


    
for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_050000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_050000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_050000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_050000_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_050000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_002100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_002100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_002100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_002100_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_002100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
          sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_000037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_000037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_000037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_000037_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_000037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
          sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_052100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_052100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_052100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_052100_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_052100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_052137_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_052137_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_052137_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_052137_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_052137_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_050037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_050037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_050037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_050037_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_050037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_210500_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_210500_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_210500_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_210500_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_210500_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250    # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_370005_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_370005_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_370005_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_370005_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_370005_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 100   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_372105_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_372105_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_372105_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_372105_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_372105_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 700   # Number of Steps (Generations) before selection event 2
    j3 = 550   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_053721_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_053721_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_053721_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_053721_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_053721_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_370521_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_370521_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_370521_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_370521_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_370521_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_mono_nat_realmu_%s_210537_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_mono_nat_realmu_%s_210537_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_mono_nat_realmu_%s_210537_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_210537_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_mono_nat_realmu_%s_210537_sample_%d.phy    sim16b_4allele_mono_nat_realmu_%s_210537_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )





#NATURAL SELECTION

# #
# sim16b_4allele_diphyletic
# 
# Script to simulate a population with n subpops of i individuals, with population rebound, 
# natural or artificial selection and sampling.  Selection occuring at different periods.  Equal events spacing.

# 2 Initial subpopns, expanding to 6

# Selection at locus 5 by exclusion of minor allele (bottleneck) before splitting
# If other selection events before that at locus 5, this occurs before the split also
# 
# phylip, fasta, fstat formatted output
# 
# #
# Author: Richard Stephens
# Created: July 30, 2013 09:04:33 AM 
# Modified: July 30, 2013 09:04:47 AM     
# #


import simuOpt
import simuPOP as sim
import math, os
from simuPOP.utils import export, saveCSV, Exporter, migrIslandRates, importPopulation, viewVars
from simuPOP.sampling import drawRandomSample

       
## Settings ###################################################################################
n = 2   # Initial Number of Subpopns
d = 3 # Divisor for subpopn splitting
i = 100000 # Number of Indivs/Subpopn
l = 30    # Number of Loci per Chromosome
c = 2 # Number of Chromosomes
g = 10 # Number of Steps (Generations) before expansion
t = 2500 # Total Number of Steps (Generations)
u = 8.85e-7 # Forward Mutation Rate # estimate from Kahler et al. 1984
v = 8.85e-8# Backward Mutation Rate
m = 0.0001  # Migration rate - Based on Morrell et al 2003 - historical coalescent estimated average per-locus value in wild barley
e = 98 # Proportion of selfing
s1 = 0.008 # selection coefficient
q = 100000 # Maximum Population Size
r1 = 2.5e-4 #mean selection intensity as on chr 3H (based on a rate of 0.4Mb/cM) ie.(1/4/100)
## Settings ###################################################################################
 

for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
        
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_000000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_000000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_000000_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_000000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_000000_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_000000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )


    
for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_050000_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_050000_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_050000_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_050000_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_050000_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_050000_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
            sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
            sim.DiscardIf(True, subPops=[
                (0, 'Affected')], at = j1),
                        sim.DiscardIf(True, subPops=[
                (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_002100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_002100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_002100_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_002100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_002100_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_002100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1),
          sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_000037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_000037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_000037_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_000037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_000037_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_000037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
    #       sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
    #       sim.DiscardIf(True, subPops=[
    #           (0, 'Affected')], at = j1),
    #       sim.DiscardIf(True, subPops=[
    #           (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
          sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_052100_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_052100_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_052100_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_052100_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_052100_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_052100_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_052137_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_052137_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_052137_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_052137_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_052137_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_052137_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_050037_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_050037_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_050037_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_050037_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_050037_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_050037_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_210500_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_210500_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_210500_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_210500_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_210500_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_210500_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
    #       sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250    # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_370005_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_370005_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_370005_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_370005_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_370005_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_370005_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
    #       sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 100   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_372105_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_372105_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_372105_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_372105_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_372105_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_372105_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 700   # Number of Steps (Generations) before selection event 2
    j3 = 550   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_053721_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_053721_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_053721_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_053721_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_053721_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_053721_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 550   # Number of Steps (Generations) before selection event 2
    j3 = 250   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_370521_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_370521_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_370521_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_370521_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_370521_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_370521_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



for iter in range(10):

    z1 = 500   # Number of Steps (Generations) before population splitting
    j1 = 400   # Number of Steps (Generations) before selection event 1
    j2 = 250   # Number of Steps (Generations) before selection event 2
    j3 = 700   # Number of Steps (Generations) before selection event 3
    
    d1 = 50 # duration of selection event 1 (in generations)
    d2 = 50 # duration of selection event 2 (in generations)
    d3 = 50 # duration of selection event 3 (in generations)
    
    
    
    pop = sim.Population(size=[i]*n, loci=[l]*c, lociPos = list(range(0,30,1))*2, infoFields=['migrate_to','fitness'])
    pop.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter()
        ])
    )
    z=20
    def demo(pop,gen):
        global q, g
        rate=2**15
        # if .... stop growing
        if gen >= g and all([x < q  for x in pop.subPopSizes()]):
            return ([(x+rate) for x in pop.subPopSizes(ancGen=-1)])
        elif gen >= g and all([x > q  for x in pop.subPopSizes()]):
            return ([q for x in pop.subPopSizes(ancGen=-1)])
        else:
            return ([x for x in pop.subPopSizes(ancGen=-1)])
            
    def sampleAndExport(pop):
        sz = pop.subPopSizes()
        new_sz = [x//2000 for x in sz]
        sample = drawRandomSample(pop, new_sz)
        export(sample, format = 'fstat', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_210537_sample_%d.dat', gui = False),
        export(sample, format = 'phylip', output = 'sim16b_4allele_diphyletic_nat_realmu_%s_210537_sample_%d.phy' (iter, pop.dvars().gen), alleleNames = ('A','C','G','T'), gui=False),
        os.system('perl convert_diploid.pl   N   sim16b_4allele_diphyletic_nat_realmu_%s_210537_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_210537_merged_sample_%d.phy' % (iter, pop.dvars().gen)),
        os.system('perl phylip_to_fasta.pl   sim16b_4allele_diphyletic_nat_realmu_%s_210537_sample_%d.phy    sim16b_4allele_diphyletic_nat_realmu_%s_210537_sample_%d.fas' % (iter, pop.dvars().gen)),
        return True
    
    simu = sim.Simulator(pop, rep=1)
    simu.evolve(
        initOps=[
            sim.InitGenotype(loci = [3,7,21,36,50,54], freq=(0.0, 0.80, 0.20, 0.0)),
            sim.InitGenotype(loci = [4,6,17,20,22,35,37,39,51,53,55], freq=(0.2, 0.0, 0.8, 0.0)),
            sim.InitGenotype(loci = [5,19,23,38,52,58], freq=(0.0, 0.2, 0.0, 0.8)),
            sim.InitGenotype(loci = [0,12,24,29,30,42,43,57,59], freq=(1.0, 0.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [1,9,14,16,26,33,44,46,47,49], freq=(0.0, 1.0, 0.0, 0.0)),
            sim.InitGenotype(loci = [8,13,25,31,32,40,41,45], freq=(0.0, 0.0, 1.0, 0.0)),
            sim.InitGenotype(loci = [2,10,11,15,18,27,28,34,48,56], freq=(0.0, 0.0, 0.0, 1.0)),
        ],
        preOps=[
            sim.SNPMutator(u = u, v = v),
            sim.Migrator(rate=migrIslandRates(m, n), end = z1-1),
            sim.SplitSubPops(subPops = [k for k in (list(range(n)))], proportions=[(1.0/d)]*d, randomize=True, at=z1),
            sim.Migrator(rate=migrIslandRates(m, (n*d)), begin = z1+1),
        sim.MaPenetrance(loci=5, penetrance=[1, 1, 0.0], wildtype=1),
        sim.DiscardIf(True, subPops=[
            (0, 'Affected')], at = j1),
        sim.DiscardIf(True, subPops=[
            (1, 'Affected')], at = j1),
            sim.MaSelector(loci=21, fitness = [1,1-s1,1-s1], wildtype = 2, begin = j2, end = (j2+d2)),
            sim.MaSelector(loci=37, fitness = [1,1-s1,1-s1], wildtype = 0, begin = j3, end = (j3+d3)),
        ],
        
        #mating scheme structure for sexless mating with variable proportion of selfing, with recombination
        matingScheme=sim.HeteroMating([
            sim.SelfMating(ops=sim.Recombinator(intensity=r1), weight=e),
            sim.HomoMating(
                chooser=sim.CombinedParentsChooser(
                    sim.RandomParentChooser(),
                    sim.RandomParentChooser()),
                generator=sim.OffspringGenerator(
                    ops= [
                        sim.MendelianGenoTransmitter(),
                        sim.Recombinator(intensity=r1)
                    ],
                ),
                weight=100-e)
            ],
        subPopSize=demo
        ),
        postOps=[
            simOperator(func=sampleAndExport, at = [0,550,750,t])
        ],
        gen=(t+1),
    )



print('End of file')