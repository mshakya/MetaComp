#
# load Kraken assignment #1
#
dat <- data.frame(load_kraken_assignment(
  "../test_data/SSputum-dil-DNase-cDNA/kraken_mini/allReads-kraken_mini.list.txt"))

# columns import test
expect_that( dim(dat)[1], equals(2977) )

#rows import test
expect_that( dim(dat)[2], equals(4) )

# col names test
expect_that(colnames(dat), equals( c("LEVEL", "TAXA", "COUNT", "ABUNDANCE") ) )

# [0.0] subset family
family_table <- dplyr::filter( dat, LEVEL == "family")

# [0.1] subset taxa
ent_row <- dplyr::filter( family_table, TAXA == "Propionibacteriaceae")

# [0.2] test
# family Propionibacteriaceae 4123 8 31957
expect_that(ent_row$COUNT, equals(4123))
expect_that(ent_row$ABUNDANCE, equals(0.2913297))

#
# test the failure
#
expect_that(load_metaphlan_assignment("../test_data/nonexistantfile.txt"), throws_error())
