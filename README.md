# Supplementary Archive for: Phylodynamic insights on the early spread of the COVID-19 pandemic and the efficacy of intervention measures
This archive contains materials required to replicate the analyses and reproduce the results presented in our study.
It is divided into four directories: the [`data`](#data) directory contains the genomic, geographic and sampling data used in our phylodynamic analyses, and also information regarding human-travel volume and public-health measures; the [`analyses`](#analyses) directory contains the XML scripts used to perform phylodynamic analyses using the `BEAST` software package; the [`scripts`](#scripts) directory contains the `R` scripts used to curate data and to postprocess output from analyses, and; the [`program`](#program) directory contain a modified version of the `BEAST` software package that implements various extensions required to perform some of the phylodynamic analyses in our study. The details of each subdirectory are described below.
Note that we also provide this archive as a [Dryad repository](link here).

## <a name="data"></a>Data
The `data` directory is divided into three subdirectories, one for each primary type of data used in our study, including:
* SARS-CoV-2 genomic sequence data ([`data/sequence_data`](#sequence_data)),
* human-travel data ([`data/sequence_data`](#travel_data)), and
* intervention-measure data ([`data/measure_data`](#measure_data)).

### <a name="sequence_data"></a>SARS-CoV-2 Sequence Data
Our study performed analyses of two SARS-CoV-2 genomic datasets: the *reduced dataset*, with 1271 sequences, and the *entire dataset*, with 2598 sequences.
The reduced dataset is based on all SARS-CoV-2 genomic sequences&mdash;collected during the early phase of the COVID-19 pandemic&mdash; that were deposited in [GISAID](https://www.gisaid.org) as of April 19 2020.
The entire dataset is based on all SARS-CoV-2 sequences collected during the early phase of the COVID-19 pandemic that were deposited in GISAID as of September 22.

GISAID policy prevents us from sharing the `FASTA` files for the SARS-CoV-2 genomic sequences comprising our reduced and entire datasets. We therefore provide the GISAID accession numbers for the SARS-CoV-2 sequences: those comprising the reduced dataset are listed in `data/sequence_data/gisaid_acknowledgement_table_041920_030820.tsv`, and sequences comprising the entire dataset arew listed in `data/sequence_data/gisaid_acknowledgement_table_092220_030820.tsv`.
The metadata for each sequence&mdash;including the sampling date and geographic location&mdash;are listed in `data/sequence_data/metainfo_041920.tsv` for the reduced dataset, and `data/sequence_data/metadata_092220.tsv` for the entire dataset.
`data/sequence_data/exclude_seqs.tsv`. The accession numbers of sequences that we excluded during the data curation stage (on the basis that they are duplicate sequences) are listed in (https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt)(modified from [Nextstrain ncov exclude list].

### <a name="travel_data"></a>Human-Travel Data
Our study included three types of human-travel data, including:
* data on the volume of global air travel ([`data/travel_data/global_flights`](#global_flights_data)),
* data on the volume of domestic travel in China ([`data/travel_data/china_domestic`](#china_domestic_data)), and
* data on the volume of domestic travel in the United States ([`data/travel_data/us_domestic`](#us_domestic_data)).

#### <a name="global_flights_data"></a>Global air-travel volume
Data on the daily volume of global air travel is located in `data/travel_data/global_flights/nflights_daily_byaircraft.csv`&mdash;sourced from FlightAware&mdash;, which contains the number of all commercial passenger flights (with the number of flights for each type of aircraft) on each day between December 30, 2019 to March 8, 2020.
We summarize data on the passenger occupancy (number of seats) in each type of aircraft in  `data/travel_data/global_flights/aircraft_nseats.csv`.
We then estimated the daily air-travel volume (assuming maximum capacity for all flights) as follows: first, we compute the product of the number of flights for each type of aircraft times the occupancy for that type of aircraft, and then summing over all aircraft types.

#### <a name="china_domestic_data"></a>Volume of domestic travel in China
We summarize data on daily mobility within China between January 1 to March 8, 2020 in `data/travel_data/china_domestic/baidu_mobility_index.csv`; we acquired these data on the daily domestic mobility index from [the Baidu Migration platform](https://qianxi.baidu.com/2020) sourced from the Baidu location-based services.

#### <a name="us_domestic_data"></a>Volume of domestic travel in the United States
We extracted the daily domestic mobility indices within the US from two publicly available datasets, [Apple COVID-19 mobility trends reports](https://covid19.apple.com/mobility) and [Google COVID-19 community mobility reports](https://www.google.com/covid19/mobility/).
`data/travel_data/us_domestic/apple_mobility_indices.csv` contains three types of mobility indices (walking, driving, and transit) provided by Apple, which are based on user direction requests in Apple Maps.
Each daily mobility index is the relative number of direction requests made on that day compared to the corresponding baseline volume on January 13, 2020 (these indices are only made available by Apple from this date onward to help study mobility trends during the COVID‑19 pandemic).
`data/travel_data/us_domestic/google_mobility_indices.csv` contains mobility index for each of six types of places (transit, retail and recreation, groceries and pharmacy, parks, workplaces, and residential), derived from Google location-based services.
These mobility indices are relative values, where the mobility index for a given day is computed as the number of visits on that day divided by the median number of visits on that day measured from January 3 to February 6, 2020 (these indices are only made available by Google from February 6 onward onward to help study mobility trends during the COVID‑19 pandemic).

### <a name="measure_data"></a>Measure Data
We focussed on two types of intervention measures enacted during the early phase of COVID-19 pandemic that involved China: targeted-containment measures with China (*i.e.*, international air-travel bans), and the domestic-mitigation measures within China.
`data/measure_data/international_airtravelban_withchina.csv` contains a collection of countries or territories that enacted international travel bans with China, as well as the associated initiation date and information source.
`data/measure_data/china_domestic.csv` contains a collection of provinces or cities of China that enacted domestic mitigation measures, as well as the associated implementation period, measure type, and information source.

## <a name="analyses"></a> Phylodynamic Analyses
The `analyses` subdirectory contains the `BEAST` XML scripts we used to preform the phylodyamic analyses in this study.
It is further divided into five subdirectories, corresponding to the five analysis steps we took in this study, including:
* the analyses that infer a dated phylogeny of the sampled viruses from their sequence alignment and sampling dates using the reduced dataset ([`analyses/step1_dated_phylogeny_inference`](#step1_dated_phylogeny_inference_analyses));
* the analyses that evaluate the fit of candidate biogeographic models to the reduced dataset ([`analyses/step2_model_evaluation`](#step2_model_evaluation_analyses));
* the analyses that jointly estimate the dispersal history and the dated phylogeny of SARS-CoV-2 using the entire dataset ([`analyses/step3_joint_analyses`](#step3_joint_analyses)), under the preferred model identified by the model-evaluation analyses;
* the ancillary analyses that estimate daily global dispersal rate of SARS-CoV-2 using the entire dataset ([`analyses/ancillary1_daily_rate`](#ancillary1_daily_rate)); and
* the ancillary analyses that assess the sensitivity of our results to incomplete and non-random sampling sampling of SARS-CoV-2 sequences, by replicating the series of analyses using the entire dataset with the reduced dataset ([`analyses/ancillary2_reduced_dataset`](#ancillary2_reduced_dataset)).

All these analyses were performed using a Bayesian phylogenetic program, `BEAST` (the specific version varies among types of analyses; see details below), with the `BEAGLE` library (compiled from [the `hmc-clock` branch, commit `dd36bf5`](https://github.com/beagle-dev/beagle-lib/tree/dd36bf5b8d88348c77a93eeeef0917d90df71a4f)) enabled to accelerate computation (both CPU and GPU were used).
The XML scripts provided here use `?` as placeholder for each sequence (as we are not allowed to publish the sequence data directly per the terms of use of GISAID); therefore to reproduce the analyses, one needs to download the sequences from GISAID using the GISAID accession numbers provided in the [`data/sequence_data`](#sequence_data) subdirectory to, then run the `R` scripts provided in the [`scripts/sequence_data_curation`](#sequence_data_curation_scripts) subdirectory to generate the curated sequence alignments we used in this study, and lastly replace the `?`s with the corresponding sequences to produce the exact XML scripts we used.
For better interpretation, the name of each sequence in the XML scripts comprises three parts, denoting the GISAID accession ID, sampling geographic area, and sampling date, respectively; for instance, the sequence listed under GISAID ID `EPI_ISL_402124` was sampled in China on December 30, 2019, so its name in the XML scripts is thus `402124_China_123019`.

### <a name="step1_dated_phylogeny_inference_analyses"></a>Step 1: Estimating the dated phylogeny of the reduced dataset
`data/step1_dated_phylogeny_inference/coalExp_ucln_posterior_run1.xml` is the XML script used to infer the dated phylogeny of the reduced SARS-CoV-2 dataset.
`data/step1_dated_phylogeny_inference/coalExp_ucln_posterior_MCC.tre` is the Maximum Clade Credibility (MCC) tree summarized from the inferred posterior distribution of dated phylogenies (negative branch lengths of this MCC tree, resulted from the summarizing process, were assigned a small positive value, 0.001 days).
These analyses were performed using [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1).

### <a name="step2_model_evaluation_analyses"></a>Step 2: Evaluating candidate biogeographic models using the reduced dataset
We explored a pool of candidate biogeographic models corresponding to possible combinations of:
* symmetric and asymmetric Q matrices;
* default and alternative priors on the total number of dispersal routes;
* default and alternative priors on the average dispersal rate, and;
* models where the Q matrices and average dispersal rates are piecewise constant over 1, 2 (before and after February 2 2020), or 4 (partitioned by January 12, February 2, and February 16 2020) pre-specified intervals.

We assessed the *relative fit* of these candidate biogeographic models by computing Bayes factors using the marginal likelihood estimated for each model.
The XML scripts for the power-posterior analyses that estimate marginal likelihood can be found in subdirectory `analyses/step2_model_evaluation/marginal_likelihood`.
We also assessed the *absolute fit* of each candidate biogeographic model using posterior-predictive simulation.
Subdirectory `analyses/step2_model_evaluation/constant/posterior_predictive` contains the XML scripts for the analyses inferring the posteriors, which were then used to perform posterior-predictive simulations.

We did these analyses using the reduced dataset.
To ensure numerical stability of the estimates, these geographic model-evaluation analyses were conditioned on the MCC tree inferred in the [previous step](#step1_dated_phylogeny_inference_analyses).
The analyses under the constant biogeographic models were performed using [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1), whereas the analyses under the piece-wise constant models were performed using a modified version of `BEAST` (see the [program section](#program) below for details).

### <a name="step3_joint_analyses"></a>Step 3: Jointly inferring SARS-CoV-2 phylogeny and dispersal history using the entire dataset
We estimated the interval-specific dispersal dynamics under the 4-interval piecewise-constant alternative-prior asymmetric geographic model (the best model identified in the [previous section](#step2_model_evaluation_analyses)) jointly with the dated phylogeny to accommodate phylogenetic uncertainty using the entire dataset.
The XML scripts for these joint analyses can be found in `analyses/step3_joint_analyses`.
The inferred posterior distribution under the 4-interval model using the entire dataset can be found in the Dryad repository (not available here due to file size limit).
These analyses were performed with the modified version of `BEAST`.

### <a name="ancillary1_daily_rate"></a>Ancillary analyses 1: Estimating daily global viral dispersal rates for the entire dataset
We estimated the daily global dispersal rate using the entire SARS-CoV-2 dataset under a more granular phylodynamic model that allows the average dispersal rate to vary from day to day (while keeping the other biogeographic model components identical to the joint inference) to assess the correlation between these estimates with the daily volume of global air travel during the early phase.
`analyses/ancillary1_daily_rate` contains the XML scripts for the analyses inferring the daily global dispersal rates.
In these analyses, we accommodated phylogenetic uncertainty by averaging over (a subsample of) the marginal posterior distribution of dated phylogenies inferred in the joint analyses; this distribution of trees can be found in the Dryad repository (not available here due to the file size limit).
These analyses were performed with the modified version of `BEAST`.

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

### <a name="sequence_data_curation_scripts"></a>Sequence Data Curation
The first step of sequence data curation is to filter the raw sequences downloaded from GISAID those that undermine the analyses due to various reasons (*e.g.*, a sequence that is incomplete or lacks the associated sampling time or location information).
The specific procedures differ slightly between the reduced dataset and the entire dataset, since we downloaded unaligned sequences from GISAID to generate the reduced dataset, but sequence alignment to generate the entire dataset.
We present both procedures here so that either way can be replicated.
For the reduced dataset (when the initial input are raw sequences), `scripts/sequence_data_curation/sequence_preprocess_beforealign.R` is run first as a preliminary filtration of the sequences, and then (after the alignment is inferred using the output of the preliminary filtration), `scripts/sequence_data_curation/sequence_process_aftealign.R` is run to filter the inferred alignment again (as some filtration conditions in this step require that site homology has been established).
For the entire dataset (when the initial input is sequence alignment), `scripts/sequence_data_curation/sequence_process_gisaidalignment.R` is run to filter it in one installment.
These scripts also trim the alignment after the filtration to generate nucleotide alignments that only keep the coding regions.

After the curation step is done, `scripts/sequence_data_curation/sequence_process_alignment_foranalysis.R` can be run to partition the entire alignment according to various partition schemes (*e.g.*, by codon positions and/or by gene regions), so that different substitution model among the partitions can be easily specified in the analyses.
`scripts/sequence_data_curation/discretetrait_writer.R` generates a data table containing the sampling age (in unit of days) and sampling geographic area (according to the area-grouping scheme defined in `scripts/sequence_data_curation/main_config.R`) information of each virus sample included in the processed sequence alignment.
Finally, `scripts/sequence_data_curation/generate_xml_data_head.R` reads in the partitioned sequence alignment (generated by `sequence_process_alignment_foranalysis.R`) and the sampling-information table (generated by `discretetrait_writer.R`) to generate the `taxa` part of a XML script that can replace the corresponding part of XML scripts (whose sequence data were represented by `?`s due to the terms of use of GISAID) provided in the `analyses` subdirectory to reproduce our analyses.

### <a name="history_simulation_scripts"></a>History Simulations
We simulated full dispersal history over the sampled dated phylogeny for two major purposes: inferring the posterior distribution of histories (through stochastic mapping) and assessing the adequacy of geographic models (through posterior-predictive simulation).
`scripts/history_simulation/history_simulator_functions.R` is the interface script that a user can execute the `history_simulator` routine within it to either perform stochastic mapping or posterior-predictive simulation (by setting the argument `conditional` to true or false).
Other scripts in this subdirectory contains subroutines that will be evoked by the interface routine.
One exception is `scripts/history_simulation/posterior_predictive_teststatistics_functions`, which provides functions to calculate posterior-predictive summary statistics (including the tipwise multinomial statistic and the parsimony statistic).

### <a name="parameter_summary_scripts"></a>Geographic Model Parameter Summary
`scripts/parameter_summary/get_BFs_counts.R` processes the `BEAST` parameter log file (and simulated history log file if available) to summarize the geographic model parameter estimates for both constant and piecewise constant model, including (interval-specific) pairwise dispersal routes, (interval-specific) global dispersal rate, and (interval-specific) number of pairwise dispersal events.
`scripts/parameter_summary/get_daily_events_nlieanges.R` processes the `BEAST` simulated history log file (or the tree log file) to summarize the number of dispersal events occurred on each dispersal route (or the number of active viral lineages) in each arbitrarily specified time interval.

## <a name="program"></a>Program
We included an executable `BEAST` program compiled from the `BEAST` source code we modified, which extends the machinery (that was designed for the constant model) to work for the piecewise constant model (including both phylogenetic likelihood computation and stochastic mapping).
Our modification of the source code was based on [the `master` branch, commit `c0e5f8d3a`](https://github.com/beast-dev/beast-mcmc/commit/c0e5f8d3af1f943bed5754d8e6466aab4eb776f4); the modified source code (from which the executable that we used in our analyses in this study was compiled) can be found as [the `epochal_model_extension` branch, commit `49e476e41`](https://github.com/jsigao/beast-mcmc/commit/49e476e4198d6328933d3d8fcaf3923374eaf5e2) of a forked `BEAST` source code repository.
As a quick example, this executable can be run by invoking:
```
java -jar beast.jar -working ./analyses/parameter_estimation/interval_specific/092220_030820/4interval_exphyper_asymmetric_poissonintermediate_orf1restcp_posterior_run1.xml
```
See [`BEAST` official tutorial website](http://beast.community/index.html) for further details about running `BEAST` analyses using a `Java` executable.
