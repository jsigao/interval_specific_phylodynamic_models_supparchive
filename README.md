# Supplementary Archive for: Phylodynamic insights on the early spread of the COVID-19 pandemic and the efficacy of intervention measures
This archive contains the necessary materials to replicate the analyses and reproduce the results of our study.
These materials are located within four directories: the [`data`](#data) directory contains the genomic, geographic and sampling data used in our phylodynamic analyses, and also includes data on human-travel volume and intervention measures; the [`analyses`](#analyses) directory contains the XML scripts specifying the phylodynamic analyses we performed using the `BEAST` software package; the [`scripts`](#scripts) directory contains the `R` scripts that we used to curate data and to process the output from our analyses, and; the [`program`](#program) directory contain a modified version of the `BEAST` software package that implements various extensions required to perform some of the phylodynamic analyses in our study. Each of these directories is detailed below.
Note that we also provide this archive as a [Dryad repository](link here).

## <a name="data"></a>Data
The `data` directory is divided into three subdirectories, one for each primary type of data used in our study:
* SARS-CoV-2 genomic sequence data ([`data/sequence_data`](#sequence_data)),
* human-travel data ([`data/sequence_data`](#travel_data)), and
* intervention-measure data ([`data/measure_data`](#measure_data)).

### <a name="sequence_data"></a>SARS-CoV-2 Sequence Data
Our study involves analyses of two SARS-CoV-2 genomic datasets: the *reduced dataset*, with 1271 sequences, and the *entire dataset*, with 2598 sequences.
The reduced dataset is based on all SARS-CoV-2 genomic sequences&mdash;collected during the early phase of the COVID-19 pandemic&mdash;that were deposited in [GISAID](https://www.gisaid.org) as of April 19, 2020.
The entire dataset is based on all SARS-CoV-2 sequences collected during the early phase of the COVID-19 pandemic that were deposited in GISAID as of September 22, 2020.

GISAID policy prohibits us from directly sharing the `FASTA` files for the SARS-CoV-2 genomic sequences used in our study. Instead, we provide the corresponding GISAID accession numbers of those SARS-CoV-2 sequences: sequences comprising the reduced dataset are listed in `data/sequence_data/gisaid_acknowledgement_table_041920_030820.tsv`, and sequences comprising the entire dataset are listed in `data/sequence_data/gisaid_acknowledgement_table_092220_030820.tsv`.
The metadata for each sequence&mdash;including the sampling date and geographic location&mdash;are listed in `data/sequence_data/metainfo_041920.tsv` for the reduced dataset, and `data/sequence_data/metadata_092220.tsv` for the entire dataset.
The accession numbers of sequences that we excluded during the data-curation stage (because they are duplicates) are listed in `data/sequence_data/exclude_seqs.tsv` (modified from [Nextstrain ncov exclude list(https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt)].

### <a name="travel_data"></a>Human-Travel Data
Our study also involved three types of data regarding human-travel volume:
* data on the volume of global air travel ([`data/travel_data/global_flights`](#global_flights_data)),
* data on the volume of domestic travel in China ([`data/travel_data/china_domestic`](#china_domestic_data)), and
* data on the volume of domestic travel in the United States ([`data/travel_data/us_domestic`](#us_domestic_data)).

#### <a name="global_flights_data"></a>Global air-travel volume
Data on the daily volume of global air travel&mdash;sourced from FlightAware&mdash;is located in `data/travel_data/global_flights/nflights_daily_byaircraft.csv`; these records include the number of all commercial passenger flights (with the number of flights for each type of aircraft) on each day between December 30, 2019 to March 8, 2020.
We summarize data on the passenger occupancy (number of seats) in each type of aircraft in  `data/travel_data/global_flights/aircraft_nseats.csv`.
We then estimated the daily air-travel volume (assuming maximum capacity for all flights) as follows: first, we compute the product of the number of flights for each type of aircraft times the passenger occupancy for that type of aircraft, then we sum over all aircraft types.

#### <a name="china_domestic_data"></a>Volume of domestic travel in China
Data on daily domestic travel within China&mdash;sourced from [the Baidu Migration platform](https://qianxi.baidu.com/2020) for location-based services&mdash;is located in `data/travel_data/china_domestic/baidu_mobility_index.csv`; these records include the domestic mobility index for each day between January 1 to March 8, 2020.

#### <a name="us_domestic_data"></a>Volume of domestic travel in the United States
We extracted data on domestic travel within the US from two publicly available datasets, [Apple COVID-19 mobility trends reports](https://covid19.apple.com/mobility) and [Google COVID-19 community mobility reports](https://www.google.com/covid19/mobility/).
The Apple location-based services provide data on three types of mobility indices (walking, driving, and transit), which are summarized in `data/travel_data/us_domestic/apple_mobility_indices.csv`. These records are based on user direction requests in Apple Maps, where each daily mobility index is the relative number of direction requests made on that day compared to the corresponding baseline volume on January 13, 2020 (these indices are only made available by Apple from this date onward to help study mobility trends during the COVID‑19 pandemic).
Google location-based services provide data for six types of mobility indices (transit, retail and recreation, groceries and pharmacy, parks, workplaces, and residential), which are summarized in `data/travel_data/us_domestic/google_mobility_indices.csv`. These mobility indices are relative values, where the mobility index for a given day is computed as the number of visits on that day divided by the median number of visits on that day measured from January 3 to February 6, 2020 (these indices are only made available by Google from February 6 onward onward to help study mobility trends during the COVID‑19 pandemic).

### <a name="measure_data"></a>Public-Health Measures
We focussed on two types of public-health measures involving China during the early phase of COVID-19 pandemic: targeted-containment measures (*i.e.*, international bans on air-travel with China), and the domestic-mitigation measures within China.
The table `data/measure_data/international_airtravelban_withchina.csv` contains a collection of countries or territories that enacted international travel bans with China, as well as the associated initiation date and information source.
The table `data/measure_data/china_domestic.csv` contains a collection of provinces or cities of China that enacted domestic-mitigation measures, as well as the associated implementation period, measure type, and information source.

## <a name="analyses"></a> Phylodynamic Analyses
The `analyses` directory contains the XML scripts that we used to preform the phylodyamic analyses in the `BEAST` software package.
We include a subdirectory for each stage of analyses that we performed in our study, including:
* analyses to infer a dated phylogeny for the reduced SARS-CoV-2 dataset ([`analyses/step1_dated_phylogeny_inference`](#step1_dated_phylogeny_inference_analyses));
* analyses to evaluate the fit of candidate biogeographic models to the reduced SARS-CoV-2 dataset ([`analyses/step2_model_evaluation`](#step2_model_evaluation_analyses));
* analyses to infer the phylogeny, divergence times and biogeographic history for the entire SARS-CoV-2 dataset  ([`analyses/step3_joint_analyses`](#step3_joint_analyses)), under the preferred model identified by the model-evaluation analyses;
* ancillary analyses to estimate daily global dispersal rates for the entire SARS-CoV-2 dataset ([`analyses/ancillary1_daily_rate`](#ancillary1_daily_rate)); and
* ancillary analyses to assess the sensitivity of our results to incomplete and non-random sampling of SARS-CoV-2 sequences ([`analyses/ancillary2_reduced_dataset`](#ancillary2_reduced_dataset)).

We performed all analyses using the `BEAST` software package, although the specific version varied among analyses (see details below), with the `BEAGLE` library (compiled from [the `hmc-clock` branch, commit `dd36bf5`](https://github.com/beagle-dev/beagle-lib/tree/dd36bf5b8d88348c77a93eeeef0917d90df71a4f)) enabled to accelerate both CPU and GPU computation.
Note that the XML scripts that we incude in this repository include `?` as placeholders for each sequence (because GISAID policy prohibits us from publishing the sequence data). In order to reproduce our analyses, it is necessary to: (1) first, download the sequences from GISAID using the GISAID accession numbers provided in the [`data/sequence_data`](#sequence_data) subdirectory; (2) then, run the `R` scripts provided in the [`scripts/sequence_data_curation`](#sequence_data_curation_scripts) subdirectory to recreate our curated sequence alignments, and; (3) finally, replace the `?`s with the corresponding sequences to replicate our XML scripts.
To improve clarity, the name of each sequence in the XML scripts comprises three parts, denoting the GISAID accession ID, sampling geographic area, and sampling date, respectively. For example, the sequence listed under GISAID ID `EPI_ISL_402124` was sampled in China on December 30, 2019, so its name in the XML scripts is thus `402124_China_123019`.

### <a name="step1_dated_phylogeny_inference_analyses"></a>Step 1: Estimating the dated phylogeny for the reduced SARS-CoV-2 dataset
We inferred the dated phylogeny for the reduced SARS-CoV-2 dataset using the XML script `data/step1_dated_phylogeny_inference/coalExp_ucln_posterior_run1.xml`.
The resulting Maximum Clade Credibility (MCC) summary tree inferred from these analyses is located in `data/step1_dated_phylogeny_inference/coalExp_ucln_posterior_MCC.tre`.
We performed these analyses using [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1).

### <a name="step2_model_evaluation_analyses"></a>Step 2: Evaluating candidate geographic models using the reduced SARS-CoV-2 dataset
We explored the relative and absolute fit of the reduced SARS-CoV-2 dataset to 12 candidate biogeographic models corresponding to combinations of:
* symmetric and asymmetric Q matrices;
* default and alternative priors on the total number of dispersal routes;
* default and alternative priors on the average dispersal rate, and;
* models where the Q matrices and average dispersal rates are piecewise constant over 1, 2, or 4 pre-specified intervals.

We assessed the *relative fit* of these candidate biogeographic models to our reduced SARS-CoV-2 dataset by computing Bayes factors from the marginal likelihoods of competing models.
We estimated the marginal likelihood for each candidate biogeographic model by performing power-posterior analyses; these XML scrips are located in the subdirectory `analyses/step2_model_evaluation/marginal_likelihood`.
We also assessed the *absolute fit* of each candidate biogeographic model to our reduced SARS-CoV-2 dataset using posterior-predictive simulation.
We estimated the posterior probability distribution for each candidate biogeographic model by performing MCMC simulation; these XML scrips are located in the subdirectory `analyses/step2_model_evaluation/constant/posterior_predictive`.

Our evaluation of candidate biogeographic models conditioned on the MCC summary phylogeny inferred in [step 1](#step1_dated_phylogeny_inference_analyses).
Our analyses under the time-constant biogeographic models were performed using [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1), whereas those under the piece-wise constant biogeographic models were performed using the modified version of `BEAST` (see the [program section](#program) for details).

### <a name="step3_joint_analyses"></a>Step 3: Jointly inferring the phylogeny, divergence times and geographic history for the entire SARS-CoV-2 dataset
We inferred the phylogeny, divergence times and geographic history for the entire SARS-CoV-2 dataset. The phylodynamic model that we used for these analyses included the previously selected geographic (sub)model (see [step 2](#step2_model_evaluation_analyses)).
We performed these joint phylodynamic analyses using the modified version of `BEAST`; the corresponding XML files are located in `analyses/step3_joint_analyses`.
We provide the inferred posterior distribution of phylodynamic model parameters in our Dryad repository (file-size constraints prevent us from including it in this repository).

### <a name="ancillary1_daily_rate"></a>Ancillary analyses 1: Estimating daily global viral dispersal rates for the entire SARS-CoV-2 dataset
We estimated the daily global viral dispersal rate for the entire SARS-CoV-2 dataset under a more granular phylodynamic model that allows the average dispersal rate to vary from day to day (while keeping the other biogeographic model components identical to the joint inference) to assess the correlation between these estimates with the daily volume of global air travel during the early phase of the COVID-19 pandemic.
We performed these analyses using the modified version of `BEAST`; the corresponding XML scripts are located in `analyses/ancillary1_daily_rate`.
To accommodate phylogenetic uncertainty, we performed these analyses while averaging over (a subsample of) the marginal posterior distribution of dated phylogenies inferred in the joint analyses; we include this tree file in our Dryad repository (file-size constraints prevent us from including it in this repository).

### <a name="ancillary2_reduced_dataset"></a>Ancillary analyses 2: Assessing the impact of incomplete viral sampling
We also performed additional series of analyses to assess the sensitivity of our results to incomplete and non-random sampling sampling of SARS-CoV-2 sequences.
Specifically, we replicated the entire series of analyses that we performed on the entire SARS-CoV-2 dataset&mdash;including joint inference of phylogeny, divergence times and biogeographic history, as well as analyses to infer daily global viral dispersal rates&mdash;for the reduced SARS-CoV-2 dataset, which has fewer (1271 vs 2598) sequences and different temporal and spatial sampling intensities.
`analyses/ancillary2_reduced_dataset` contains the XML scripts for these sensitivity analyses.
These analyses were performed with the modified version of `BEAST`.

## <a name="scripts"></a>Scripts
The `scripts` subdirectory contains the `R` scripts we used in this study; it is further divided into three subdirectories, corresponding to the three major types of processing we did (either before or after running the `BEAST` analyses), including:
* the scripts that curate SARS-CoV-2 genomic sequences acquired from GISAID to generate the sequence alignment, sampling time, and sampling geographic area data to be used in the analyses ([`scripts/sequence_data_curation`](#sequence_data_curation_scripts)),
* the scripts that perform both stochastic mappings and posterior-predictive simulations using the inferred geographic model parameters and dated phylogenies ([`scripts/history_simulation`](#history_simulation_scripts)), and
* the scripts that process the `BEAST` output files to produce summaries of geographic model parameters ([`scripts/parameter_summary`](#parameter_summary_scripts)).

### <a name="sequence_data_curation_scripts"></a>Curating SARS-CoV-2 Sequence Data
The first step of sequence data curation is to filter the raw sequences downloaded from GISAID those that undermine the analyses due to various reasons (*e.g.*, a sequence that is incomplete or lacks the associated sampling time or location information).
The specific procedures differ slightly between the reduced dataset and the entire dataset, since we downloaded unaligned sequences from GISAID to generate the reduced dataset, but sequence alignment to generate the entire dataset.
We present both procedures here so that either way can be replicated.
For the reduced dataset (when the initial input are raw sequences), `scripts/sequence_data_curation/sequence_preprocess_beforealign.R` is run first as a preliminary filtration of the sequences, and then (after the alignment is inferred using the output of the preliminary filtration), `scripts/sequence_data_curation/sequence_process_aftealign.R` is run to filter the inferred alignment again (as some filtration conditions in this step require that site homology has been established).
For the entire dataset (when the initial input is sequence alignment), `scripts/sequence_data_curation/sequence_process_gisaidalignment.R` is run to filter it in one installment.
These scripts also trim the alignment after the filtration to generate nucleotide alignments that only keep the coding regions.

After the curation step is done, `scripts/sequence_data_curation/sequence_process_alignment_foranalysis.R` can be run to partition the entire alignment according to various partition schemes (*e.g.*, by codon positions and/or by gene regions), so that different substitution model among the partitions can be easily specified in the analyses.
`scripts/sequence_data_curation/discretetrait_writer.R` generates a data table containing the sampling age (in unit of days) and sampling geographic area (according to the area-grouping scheme defined in `scripts/sequence_data_curation/main_config.R`) information of each virus sample included in the processed sequence alignment.
Finally, `scripts/sequence_data_curation/generate_xml_data_head.R` reads in the partitioned sequence alignment (generated by `sequence_process_alignment_foranalysis.R`) and the sampling-information table (generated by `discretetrait_writer.R`) to generate the `taxa` part of a XML script that can replace the corresponding part of XML scripts (whose sequence data were represented by `?`s due to the terms of use of GISAID) provided in the `analyses` subdirectory to reproduce our analyses.

### <a name="history_simulation_scripts"></a>Simulating Geographic Histories
We simulated full dispersal history over the sampled dated phylogeny for two major purposes: inferring the posterior distribution of histories (through stochastic mapping) and assessing the adequacy of geographic models (through posterior-predictive simulation).
`scripts/history_simulation/history_simulator_functions.R` is the interface script that a user can execute the `history_simulator` routine within it to either perform stochastic mapping or posterior-predictive simulation (by setting the argument `conditional` to true or false).
Other scripts in this subdirectory contains subroutines that will be evoked by the interface routine.
One exception is `scripts/history_simulation/posterior_predictive_teststatistics_functions`, which provides functions to calculate posterior-predictive summary statistics (including the tipwise multinomial statistic and the parsimony statistic).

### <a name="parameter_summary_scripts"></a>Summarizing Geographic-Model Parameters
`scripts/parameter_summary/get_BFs_counts.R` processes the `BEAST` parameter log file (and simulated history log file if available) to summarize the geographic model parameter estimates for both constant and piecewise constant model, including (interval-specific) pairwise dispersal routes, (interval-specific) global dispersal rate, and (interval-specific) number of pairwise dispersal events.
`scripts/parameter_summary/get_daily_events_nlieanges.R` processes the `BEAST` simulated history log file (or the tree log file) to summarize the number of dispersal events occurred on each dispersal route (or the number of active viral lineages) in each arbitrarily specified time interval.

## <a name="program"></a>Modified BEAST Program
We included an executable `BEAST` program compiled from the `BEAST` source code we modified, which extends the machinery (that was designed for the constant model) to work for the piecewise constant model (including both phylogenetic likelihood computation and stochastic mapping).
Our modification of the source code was based on [the `master` branch, commit `c0e5f8d3a`](https://github.com/beast-dev/beast-mcmc/commit/c0e5f8d3af1f943bed5754d8e6466aab4eb776f4); the modified source code (from which the executable that we used in our analyses in this study was compiled) can be found as [the `epochal_model_extension` branch, commit `49e476e41`](https://github.com/jsigao/beast-mcmc/commit/49e476e4198d6328933d3d8fcaf3923374eaf5e2) of a forked `BEAST` source code repository.
As a quick example, this executable can be run by invoking:
```
java -jar beast.jar -working ./analyses/parameter_estimation/interval_specific/092220_030820/4interval_exphyper_asymmetric_poissonintermediate_orf1restcp_posterior_run1.xml
```
See [`BEAST` official tutorial website](http://beast.community/index.html) for further details about running `BEAST` analyses using a `Java` executable.
