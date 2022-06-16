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

simulate_genotypes = function(n_genotype = 100, n = 1000, eafs = rep(0.2, 100)){
  if(length(eafs)!=n_genotype){
    return("Please re-specify eafs...\n")
  }
  
  if(any(c(eafs>1,eafs<0))){
    return("Please re-specify eafs...\n")
  }
  
  snp_names = paste0("GENOTYPE",seq(1,n_genotype))
  subject_name = paste0("X", seq(1,n))
  
  genotypes = as.data.frame(sapply(eafs, function(x){
    simulate_genotype(eaf = x, n = n)[,2]
  }))
  
  names(genotypes) = snp_names
  
  genotypes = cbind.data.frame(data.frame(SUBJECT = subject_name), genotypes)
  
  return(genotypes)
}


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

simulate_error = function(genotype_effect, heritability){
  residual_var = var(genotype_effect) * (1/heritability - 1)
  residual = rnorm(length(genotype_effect), sd = sqrt(residual_var))
  
  return(residual)
}


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
