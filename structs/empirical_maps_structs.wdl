version 1.0

struct Dataset {
  String name
  String parent1
  String parent2
  String cross
  Array[String]? chromosomes
  Boolean multiallelics
}

# Describes one mapping population and how to identify its samples in key files.
#
# Fields:
#   dataset          – population metadata (name, parents, cross, etc.) passed to Maps
#   pedigree         – canonical Pedigree value in GBS key files (e.g. "CQxMQ")
#   pedigree_aliases – reciprocal / alternate notations for the same cross
#                      (e.g. ["MQxCQ"]); use [] if there are none
#
# Note: dataset.parent1 and dataset.parent2 must match the SeedLot values as
# they appear in the VCF sample header (spaces replaced by underscores).
struct PopulationSpec {
  Dataset       dataset
  String        pedigree
  Array[String] pedigree_aliases
}
