#ifndef __PARAMETERS_H
#define __PARAMETERS_H

//**************************************************************************
// HOUSE KEEPING PARAMETERS
// Number of available processors (threads to create); _WIDTH _must_ be divisible by _CPUS
// Note that this version runs only with 1 CPU, not optimized for multiple CPUs
#define _CPUS                           1

// Seed for the random number generator, 0 -> seed based on time
#define _RAND_SEED                      0

// Log() only logs to graph file when (Gen%_SYS) == 0. _SYS == 0 -> disable
#define _SYS                            1

// Log() only logs to graph file when (Gen%_DEME) == 0. _DEME == 0 -> disable
#define _DEME                           1

// Log() only logs to graph file when (Gen%_DISTANCE) == 0. _DISTANCE == 0 -> disable
#define _DISTANCE                       1

// Log() only logs to state file when (Gen%_IND) == 0.  _IND == 0  -> disable
#define _IND                            0
#define _IND_SAMPLE                     2

// Log() only logs to graph file when (Gen%_HISTO) == 0. _HISTO == 0 -> disable
#define _HISTO                          1

//**************************************************************************
// INITIALIZATION PARAMETERS
// Width and height of environment
// DO NOT FORGET TO ADJUST _MAXLEVEL ABOVE (32->5; 64->6)
#define _WIDTH                          32
#define _HEIGHT                         32

// Number of different characteristics (for _*e and _*p only)
#define _k                              4

// Maximum number of individuals per patch: this has to do with memory allocation.
// 128 hits the upper max regardless of _KO (even when it is 40, so make it 256!!!)
#define _MAX_INDIVIDUALS                256

// Initial number of individuals to initialize per patch
// watch out for _MAX_INDIVIDUALS for memory and _K0 for dynamics (10 for cats on KI)
#define _INIT_INDIVIDUALS               2

// this is the 'dispersal range'
// should be an ODD number
// previous model used 3, 9, and 15
#define _DRANGE                         11

//**************************************************************************
//LIFE HISTORY PARAMETERS
// maximum age of an individual (think in mating cycles!) if generations are overlapping
// (mouse: 12 (2); cats: 10 (5); foxes: 6 (3); rabbits: 18 (3); black rat: 12(2); goats: 18(9))
#define _MAXAGE                         12

// number of mating cycles/year if generations are overlapping
// has no implications whatsoever, to convert old gens to years
// ***change it in view.h as well!
// (mouse: 6; cats: 2; foxes: 2; rabbits: 6; black rat: 6; goats: 2)
#define _MATINGCYCLES                   6

// Birth rate
// (mouse: 6; cats: 5; foxes: 3; rabbits: 4; black rat: 5; goats: 2.1)
#define _b                              6

// (Boolean) NEW flag for converting _Ln to sex chromosomes without deleting the old code
// should be ON
#define _SEXCHROMOSOME                  1

//**************************************************************************
// PLAGUE PARAMETERS
#define _PlagueInitGen                  30
#define _PlagueEndGen                   36

#define _PlagueSeInitGen                600//60
#define _PlagueSeEndGen                 600//66

#define _PlagueThInitGen                600//90
#define _PlagueThEndGen                 600//96

#define _PlagueFoInitGen                600//120
#define _PlagueFoEndGen                 600//126

#define _PlagueFiInitGen                600//150
#define _PlagueFiEndGen                 600//156

#define _PlagueIntesity                 20 //10K in beverton-holt is reduced this much during nonplague years (for memmoory issues, can I do reverse?).
//i.e. if 2000 to 200000 is different from 10,000 to 200,000 say. the number of inoculants will change? density of inoculants?
//i.e. initial density initial inoculation effort etc.

//**************************************************************************
// Immigration PARAMETERS (immigration is allowed ONLY to the edge patches, so max no. of immigrants is: 2x(width+height))
#define _IMMIGRATION                    1
// BOOLEAN, turn on one of them: either per patch immigration probability is determined by the overall density, or is constant
#define _DENSITYDEPENDENTIMMIGRATION    0
#define _CONSTANTIMMIGRATION            1

// if constant immigration rate: not used if it is density-dependent immigration
#define _ImmigrationRate                0.05 //0.05, 0.1

//**************************************************************************
// GENE DRIVE PARAMETERS
// STRATEGY tCRISPR: gene drive is t-allele (m) and prolactin knockout (f)
#define _PROLACTINLOCUS                 1     // (Boolean)
#define _PROLACTINKNOCKOUT              1     // (Boolean) knocks out and skips infertile females
#define _BOTHSEXKNOCKOUT                0     // (Boolean) 1. default prolactin knockout in males (0=OFF), in both sexes if ON(=1)
#define _XKNOCKOUT                      0     // (Boolean) 2. x shredder xknockout is in males
#define _SHREDPROB                      0.    // prob of knockout in _XKNOCKOUT
#define _RECOMRATE                      0.0   // prob. of recombination between CRISPR and t-haplotype

#define _MALEDRIVE                      1     // (Boolean) inoculated individuals are male, otherwise random sex allocation

#define _DRIVETPROB                     0.95  // probability of segregation distortion
#define _PROBCUT                        0.8   // probability of successful cut by Cas9/gRNA
#define _LOFPROB                        1.0   // probability of loss of function after a cut

#define _TDISPERSAL                     0     // (Boolean) if OFF, _DISTANCEPENALTY is not used, everyone has the same _DRANGE
#define _DISTANCEPENALTY                0     // 1, or 2; max distance that a wildtype ind moves compared to a t-allele carrying ind

//**************************************************************************
// INOCULATION PARAMETERS
// Boolean() only inoculates drive carrying individuals in THREAD() in thread.c if ON
#define _INOC                           1

// NEW VERSION
// gene-drive individuals inoculation generations init and end
// temporal inoculation effort
#define _InocGenInit                    12
#define _InocGenEnd                     12//32//18

// spatial incoulation effort
#define _InocMod                        8
// Number of individuals to inoculate per patch (always same number of across multiple inoculations)
// watch out for _MAX_INDIVIDUALS + _K0 for and dynamics
// separated by comma if _NoInoc>1, e.g. if _NoInoc=3, then {2,1,1}
#define _InocIndividuals                1

// Boolean, second inoculation
#define _SecondInoculation              0

//number of cycles after first inoc ends (expected time for tcrispr to be lost is ~16 gens)
#define _SecondInoculationGap           6 //after the first one ends
#define _SecondInoculationDuration      6

// genotypes of inoculated individuals for the corresponding generations above
// 1: TCRISPR; 2: Tw2; separated by comma if _NoInoc>1, e.g. if _NoInoc=3, then {2,1,1}
#define _InocType                       1

//**************************************************************************
// island invasion version where the initial population has functional resistant prolactin alleles present
#define _MainlandInvasion               0

// initial frequency of island specific (nonresistant) prolactin on mainland
#define _FreqProlactin                  0.1

//**************************************************************************
// OTHER LIFE HISTORY PARAMETERS

// (0.53; K=19 )
#define _SURVIVALPROBABILITY            0.53
#define _DRIVEFITNESS                   1.0 // males only? or both sexes (currently males only)

// Maximum carrying capacity (#K=round(-50*x + 45.75)) 10 for cats on KI if w=0.1
// (cats: 12; rabbits: 0.6; )
#define _K0                             30//30//15

// (Boolean) set one of them to 1, rest to zero
// _RAND_DISPERSAL: random dispersal
// _DISTANCEDENSITYDISPERSAL: density dependent dispersal based on distance, define density threshold as well
#define _RAND_DISPERSAL                 0
#define _DISTANCEDENSITYDISPERSAL       1
#define _DISPERSALCOEF                  1.
#define _DENSITYCOEF                    1.

// (Boolean) dispersal mode: (1) Child dispersal, or (0) Parent dispersal
// (PARENT DISPERSAL IS NOT USED CURRENTLY)
// again for gene-drive models it should be 1, for parent dispersal
// check the code before making it 0
#define _CHILD_DISPERSAL                1

// prob. of BEHAVIORAL POLYANDRY per female. This is also the actual GENETIC POLYANDRY, since for simplicity
// BEHAVIORAL POLYANDRY = GENETIC POLYANDRY, ALL multiple matings results in sires "multiple paternity opportunity"
// multiple mating = two males only (default 0.46)
#define _POLYANDRYPROB                  0.63

// in multiple mating, the advantage of the first male's sperm (early sperm advantage, 1ST come 1ST serve, copulatory plug)
// if both the males have the same genotype (w,w or dr,dr), then m1 = 0.7, m2 = 1-m1 = 0.3
#define _FIRSTSPERMADVANTAGE            0.5

// this is the 'offset' if the males are (w,dr), if (w,dr)=(0.7+0.2, (1-0.7)-0.2); if if (dr,w)=(0.7-0.2, (1-0.7)+0.2);
// THE TWO ABOVE CANNOT BE GREATER THAN 1 (default 0.17)
#define _SPERMCOMPETITIONCOEF           0.17

//**************************************************************************
// DO NOT CHANGE THESE
// Number of loci per characteristic
#define _Le                             1
#define _Lp                             1
#define _Lm                             1
#define _Lk                             1
#define _Lf                             1
#define _Ln                             1

// Number of bits per allele; _must_ be a member of {1,2,4,8}
#define _Ae                             4
#define _Ap                             4
#define _Am                             4
#define _Ak                             4
#define _Af                             4
#define _An                             1

// Crossover rate: original code treats loci independently, there is free recombination between the loci/trait
// but not used here since there is 1 Loci for each trait
#define _Xi                             {0.,0.,0.,0.,0.,0.}

// When enabled, this will cause niche to use Micheal's faster version of exp()
#define  _FAST_EXP                      1

// Will cause some extra printing to be done if == 1 (2 for more)
#define _VERBOSE                        1

// Causes some timing information to be computed and printed
#define _BENCHMARK                      1

// Prefix to use for output file names (before .ind, etc).  See FormatName() for more details
#define _FILE_NAME                      "%g"
//**************************************************************************

#endif
