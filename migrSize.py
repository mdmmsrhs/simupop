import simuPOP as sim
pop = sim.Population(size=[500, 1000], infoFields='migrate_to')
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.Migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(popSize=True),
        sim.PyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)