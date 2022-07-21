source("simulate.R")
set.seed(8314)

# Simulate Genotypes
n_genotype = 185
n = 4000 # Total sample size
eafs = rep(0.5, n_genotype) # Effect allele frequency
genotypes = simulate_genotypes(n_genotype = n_genotype, n = n, eafs = eafs)

# Simulate Phenotypes
effect_size_matrix = matrix(nrow = 2, ncol = n_genotype)

# Define pleitropic effects
effect_size_matrix[1,1:5] = 5
effect_size_matrix[2,1:5] = 5

effect_size_matrix[1,6:95] = 1
effect_size_matrix[2,96:185] = 1

effect_size_matrix[is.na(effect_size_matrix)] = 0
effect_size_matrix = effect_size_matrix[,sample(1:n_genotype,
                                                size = n_genotype,
                                                replace = FALSE)]

# Simulate Phenotype1
phenotype1 = simulate_phenotype(genotypes, effect_sizes = effect_size_matrix[1,], 
                                heritability = 0.5, 
                                environment_effect = rep(0, n),
                                phenotype_name = "PHENOTYPE1")

# Simulate Phenotype2
phenotype2 = simulate_phenotype(genotypes, effect_sizes = effect_size_matrix[2,], 
                                heritability = 0.5, 
                                environment_effect = rep(0, n),
                                phenotype_name = "PHENOTYPE2")

# Confirm convariance between the phenotypes
cor(x = phenotype1$PHENOTYPE1, y = phenotype2$PHENOTYPE2)

