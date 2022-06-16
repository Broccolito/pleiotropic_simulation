"
In this simulation, we simulated 1000 subjects and generated two phenotypes, 
each effected by 60 genotypes, where 40 of them are unique, 
and 20 of them are pleitropic.

All genotypes' effect sizes are either 0 (not relevant), or 1 (relevant).
The allele frequencies all of genotypes are in between 0.1 and 0.9 (randomly distributed).

The inheribilities (h2) of both phenotypes are 0.3.
"

source("simulate.R")
set.seed(8314)

# Simulate Genotypes
n_genotype = 100
n = 1000
eafs = runif(100, min = 0.1, max = 0.9)
genotypes = simulate_genotypes(n_genotype = n_genotype, n = n, eafs = eafs)

# Simulate Phenotype1
effect_sizes = c(rep(1,60),rep(0,40))
heritability = 0.3
phenotype1 = simulate_phenotype(genotypes, effect_sizes = effect_sizes, heritability = heritability, phenotype_name = "PHENOTYPE1")

# Simulate Phenotype2
effect_sizes = c(rep(0,40),rep(1,60))
heritability = 0.3
phenotype2 = simulate_phenotype(genotypes, effect_sizes = effect_sizes, heritability = heritability, phenotype_name = "PHENOTYPE2")