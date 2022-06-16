# Pleiotropic Simulation

This project aims to use simulation techniques to study pleiotropy effects in genetics.



#### Simulating genotypes

Genotypes can be simulated by specifying effective allele frequency (eaf), total sample size (n) and name of the SNP (snp_name)

```R
simulate_genotype = function(eaf = 0.2, n = 1000, snp_name = "GENOTYPE1"){
  aaf = 1-eaf
  
  f_AA = aaf^2
  f_Aa = 2*aaf*eaf
  f_aa = eaf^2
  
  n_aa = round(n*f_aa, 0)
  n_Aa = round(n*f_Aa, 0)
  n_AA = n - n_aa - n_Aa
  
  subject_name = paste0("X", seq(1,n))
  subject_genotype = sample(c(rep(2, n_aa), rep(1, n_Aa), rep(0, n_AA)), n, replace = FALSE)
  
  genotype = data.frame(SUBJECT = subject_name, GENOTYPE = subject_genotype)
  names(genotype)[2] = snp_name
  
  return(genotype)
}
```



#### Simulating phenotypes

Phenotypes are partially dictated by genotypes (heritable portion) and an error term. A phenotype can be simulated by specifying the genotype matrix, effect size of each variant, heritability and name of the trait..

```R
simulate_phenotype = function(genotypes, effect_sizes, heritability, phenotype_name = "PHENOTYPE1"){
  if(any(heritability<0, heritability>1)){
    return("Please re-specify heritability")
  }
  
  genotype_effect = simulate_genotype_effect(genotypes = genotypes, effect_sizes = effect_sizes)
  error = simulate_error(genotype_effect = genotype_effect, heritability = heritability)
  phenotype = genotype_effect + error
  
  subject_name = names(phenotype)
  names(phenotype) = NULL
  
  phenotype = data.frame(SUBJECT = subject_name, PHENOTYPE = phenotype)
  names(phenotype)[2] = phenotype_name
  
  return(phenotype)
}
```

where the heritable portion can be simulated by

```R
simulate_genotype_effect = function(genotypes, effect_sizes){
  subject_name = genotypes[,1]
  
  genotype_mat = as.matrix(genotypes[,-1])
  genotype_mat = t(genotype_mat)
  freq = (2 * rowSums(genotype_mat == 0, na.rm=TRUE) + rowSums(genotype_mat == 1, na.rm=TRUE)) / (2 * ncol(genotype_mat))
  allele_count = 2 * (genotype_mat == 0) + (genotype_mat == 1)
  weight = (allele_count - 2 * freq) / sqrt(2 * freq * (1 - freq))
  
  genotype_effect = colSums(weight * effect_sizes, na.rm = TRUE)
  names(genotype_effect) = subject_name
  
  return(genotype_effect)
}
```

and the error term can be simulated based on heritability

```R
simulate_error = function(genotype_effect, heritability){
  residual_var = var(genotype_effect) * (1/heritability - 1)
  residual = rnorm(length(genotype_effect), sd = sqrt(residual_var))
  
  return(residual)
}
```



#### Specification

In the current run, we simulated 1000 subjects and generated two phenotypes, each effected by 60 genotypes, where 40 of them are unique, and 20 of them are pleiotropic. All genotypes' effect sizes are either 0 (not relevant), or 1 (relevant). The allele frequencies all of genotypes are in between 0.1 and 0.9 (randomly distributed). The inheritabilities (h2) of both phenotypes are 0.3.

```R
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
```



### Reference

Climente-González, H. & Azencott, C.-A. martini: an R package for genome-wide association studies using SNP networks. *bioRxiv* 2021.01.25.428047 (2021) doi:10.1101/2021.01.25.428047