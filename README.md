# Supplementary Archive for Phylodynamic insights on the early spread of the COVID-19 pandemic and the efficacy of mitigation measures
This archive contains the materials that are necessary and sufficient to replicate this study; it is divided into four subdirectories: [`data`](#data), [`analyses`](#analyses), [`scripts`](#scripts), and [`program`](#program), each described below.
This archive is also available as an [Dryad repository](link here).

## <a name="data"></a>Data
The `data` subdirectory is divided into three subdirectories, corresponding to the three types of data we used in this study, including:
* SARS-CoV-2 genomic sequence data (`data/sequence_data`),
* human travel data (`data/sequence_data`), and
* the information about measures enacted to contain COVID-19 (`data/measure_data`).

### Sequence Data
Two genomic sequence datasets were used in this study, one with 1271 sequences (*reduced dataset*) and the other with 2598 sequences (*entire dataset*).
The reduced dataset was curated on April 19, based on all available SARS-CoV-2 genomic sequences&mdash;that were collected during the early phase&mdash;from [GISAID](https://www.gisaid.org) as of that date; all the model-exploration analyses in this study used this dataset.
Sequences submitted to GISAID between April 20 and September 22 were later added to this dataset, resulting into a more comprehensive dataset with 2598 sequences; the final-stage analyses (the joint phylodynamic analysis and the daily dispersal rate analysis), where we inferred focal parameters and summaries (which we presented in the main text and drew our conclusions upon) under the preferred geographic model (identified with the reduced dataset), used this dataset.

Following the terms of use of GISAID, we can not share the `FASTA` files containing the exact sequences we used in our study here; instead, we provide the GISAID accession numbers of the used sequences: `data/sequence_data/gisaid_acknowledgement_table_041920_030820.tsv` for the reduced dataset, and `data/sequence_data/gisaid_acknowledgement_table_092220_030820.tsv` for the entire dataset.
The associated sampling time and location information can be found in `data/sequence_data/metainfo_041920.tsv` (for the reduced dataset) and `data/sequence_data/metadata_092220.tsv`.
`data/sequence_data/exclude_seqs.tsv` (modified from [Nextstrain ncov exclude list](https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt) contains name of the sequences that were considered to be duplicates of other sequences in the same dataset, which is used in the data curation step to remove those sequences.

### Travel Data
Three types of travel data were used in this study, including:
* global air-travel-volume data (`data/travel_data/global_flights`),
* China domestic travel-volume data (`data/travel_data/china_domestic`), and
* US domestic travel-volume data (`data/travel_data/us_domestic`).

#### Global air-travel-volume data
`data/travel_data/global_flights/nflights_daily_byaircraft.csv` is the global air-travel-volume data report acquired from FlightAware, containing the number of all commercial passenger flights (broken down into the number of each aircraft type) per day ranging from December 30 2019 to March 8 2020.
`data/travel_data/global_flights/aircraft_nseats.csv` contains the information we manually collected about the the number of seats availale on each type of aircraft.
Estimate of the daily number of air passengers (assuming full capacity on each aircraft) can be obtained by summing up the product of the number of flights under each aircraft type per day and the correponding number of seats.

#### China domestic travel-volume data
`data/travel_data/china_domestic/baidu_mobility_index.csv` contains the daily domestic mobility index of China (ranging from January 1 2020 to March 8 2020), acquired from [the Baidu Migration platform](https://qianxi.baidu.com/2020); this dataset is derived from Baidu location-based services.

#### US domestic travel-volume data
We extracted the daily domestic mobility indices of the US from two publicly available datasets, [Apple COVID-19 mobility trends reports](https://covid19.apple.com/mobility) and [Google COVID-19 community mobility reports](https://www.google.com/covid19/mobility/).
`data/travel_data/us_domestic/apple_mobility_indices.csv` contains three types of mobility indices (walking, driving, and transit) provided by Apple, which were produced based on the direction requests in Apple Maps.
Each mobility index of each day is a relative volume of directions requested made on that day compared to the corresponding baseline volume on January 13 2020 (so these indices are only available after that date).
`data/travel_data/us_domestic/google_mobility_indices.csv` contains mobility index for each of six types of places (transit, retail and recreation, groceries and pharmacy, parks, workplaces, and residential) derived from Google location-based services.
These mobility indices are also relative values where the mobility index of each day is computed as ratio of the number of visits on that day over the baseline value, which is the median value (of the correponding day of week) during January 3--February 6, 2020 (so these indices are only available after February 6).

### Measure Data
We focussed on two types of containment measures involving China enacted during the early phase of COVID-19 pandemic in this study: international travel bans with China and the domestic containment measures within China.
`data/measure_data/international_airtravelban_withchina.csv` contains a collection of countries or territories that enacted international travel bans with China, and the associated initiation date and information source.
`data/measure_data/china_domestic.csv` contains a collection of provinces or cities of China that enacted domestic containment measures, and the associated implementation period, measure type, and information source.

## <a name="analyses"></a> Analyses
The `analyses` subdirectory contains the `BEAST` XML scripts we used to preform the analyses in this study.
It is further divided into three subdirectories, corresponding to the three types of analyses we did in this study, including:
* the analyses that inferred a dated phylogeny of the sampled viruses from their sequence alignment and sampling dates (`analyses/dated_phylogeny_inference`),
* the analyses that evaluated the fit of candidate biogeographic models to the data (`analyses/model_exploration`), and
* the analyses that estimated the geographic model parameters and dispersal history of SARS-CoV-2 (`analyses/parameter_estimation`) under the preferred model identified in the previous step.

All these analyses were performed using a Bayesian phylogenetic program, `BEAST` (versions vary among types of analyses; see details below), with the `BEAGLE` library (compiled from [the `hmc-clock` branch, commit `dd36bf5`](https://github.com/beagle-dev/beagle-lib/tree/dd36bf5b8d88348c77a93eeeef0917d90df71a4f)) enabled to accerlate computation (both CPU and GPU were used).
The XML scripts provided here use `?` as placeholder for each sequence (as we are not allowed to publish the sequence data directly per the terms of use of GISAID); therefore to replicate the analyses, one needs to download the sequences from GISAID using the GISAID accession numbers provided in the `data/sequence_data` subdirectory to, then run the `R` scripts provided in the `scripts/sequence_data_curation` subdirectory to generate the curated sequence alignments we used in this study, and lastly replace the `?`s with them to produce the exact XML scripts we used.

### Estimating Dated Phylogeny from SARS-CoV-2 Genomes
We inferred the dated phylogeny of the SARS-CoV-2 samples we curated by running `data/dated_phylogeny_inference/coalExp_ucln_posterior_run1.xml` with [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1).
These analyses used the reduced dataset.
`data/dated_phylogeny_inference/coalExp_ucln_posterior_MCC.tre` is the Maximum-Clade Credibility (MCC) tree summarized from the inferred posterior distribution of dated phylogenies (negative branch lengths of this MCC tree, resulted from the summarizing process, were assigned a small positive value, 0.001 days).

### Evaluating Candidate Biogeographic Models
We performed rigorous Bayesian model-selection and model-checking analyses to identify the geographic model with the best fit to the data, which we used later to estimate geographic model parameters and dispersal history.
These model-evaluation analyses were carried out in two steps: we first compared various geographic (prior)models under a constant dispersal process (`analyses/model_exploration/constant`); conditioning on the preferred (prior)model, we then compared the preferred constant models with piecewise-constant biogeographic models (`analyses/model_exploration/piecewise_constant`) where the dispersal dynamics may vary across time intervals.
We did these analyses using the reduced dataset.
To ensure numerical stability of the estimates, these geographic model-evaluation analyses were conditioned on the MCC tree inferred in the previous step.

#### Evaluating geographic (prior)models under a constant dispersal process
We evaluated a suite of biogeographic (prior)models, including all combinations of:
* symmetric and asymmetric Q matrices;
* default and alternative priors on the number of dispersal routes; and
* default and alternative priors on the global dispersal rate.

We estimated the marginal likelihood for each of the eight candidate models to evaluate its relative fit to the data; the XML scripts can be found in subdirectory `analyses/model_exploration/constant/marginal_likelihood`
Subdirectory `analyses/model_exploration/constant/posterior_predictive` contains the XML scripts for the analyses inferring the posteriors, which were then used to perform posterior-predictive simulations to evaluate the absolute fit of each candidate model.
These analyses were performed using [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1).

#### Comparing constant with piecewise-constant biogeographic models
We explored two piecewise constant models, one with two time intervals (before and after February 2 2020) and the other with four time intervals (partitioned by January 12, February 2, and February 16 2020), where each interval has its own global viral dispersal rate and dispersal dynamics.
Specifically, we evaluated a suite of biogeographic (prior)models, including all combinations of:
* constant, 2- and 4-interval piecewise constant models; and
* default and alternative priors with asymmetric Q matrices (*i.e.*, the preferred prior models identified in the previous section).

We estimated the marginal likelihood for each of the four candidate models to evaluate its relative fit to the data; the XML scripts can be found in subdirectory `analyses/model_exploration/piecewise_constant/marginal_likelihood`
Subdirectory `analyses/model_exploration/piecewise_constant/posterior_predictive` contains the XML scripts for the analyses inferring the posteriors, which were then used to perform posterior-predictive simulations to evaluate the absolute fit of each candidate model.
These analyses were performed using a modified version of `BEAST` (see the program section below for details).

### Estimating Geographic Model Parameters and Dispersal History
After identifying the favored geographic model that adequately describes the dispersal process of SARS-CoV-2, we used it together with other model components to infer the dispersal dynamics and biogeographic history of SARS-CoV-2.
We first estimated the interval-specific dispersal dynamics under the 4-interval piecewise-constant alternative-prior asymmetric geographic model (best model identified the previous section) jointly with the dated phylogeny to accommodate phylogenetic uncertainty (`analyses/model_exploration/interval_specific`).
Then we further performed another sets of analyses to estimate the daily global dispersal rate (while keeping the four time intervals identical for the Q matrices), marginalizing over the posterior distribution of dated phylogenies inferred in the joint analysis (`analyses/model_exploration/daily_rate`).
We did these analyses using both the entire dataset as well as the reduced dataset to ensure that our focal results are robust sampling artifacts.
These analyses were all performed with the modified version of `BEAST`.


#### Jointly inferring SARS-CoV-2 phylogeny and dispersal history
`analyses/model_exploration/interval_specific/041920_030820` contains the XML scripts (including constant, 2-interval, and 4-interval models) for the reduced dataset, and the XML scripts (including constant and 4-interval models) for the entire dataset can be found in `analyses/model_exploration/interval_specific/092220_030820`.
Here we inferred the interval-specific dispersal dynamics and the dispersal history jointly with the dated phylogeny to allow the uncertainties in our estimates to be appropriately revealed; the inferred posterior distribution under the 4-interval model using the entire dataset can be found in the Dryad repository (not available here due to file size limit).

#### Estimating daily global viral dispersal rates
`analyses/model_exploration/daily_rate` contains the XML scripts for the analyses inferring the daily global dispersal rates.
Here we marginalized over (a subsample of) the posterior distribution of dated phylogenies inferred in the joint analysis; this distribution of trees can be found in the Dryad repository (not available here due to the file size limit).

## <a name="scripts"></a>Scripts
The `scripts` subdirectory contains the `R` scripts we used in this study; it is further divided into three subdirectories, corresponding to the three major types processing we did (either before or after running the `BEAST` analyses), including:
* the scripts that curate SARS-CoV-2 genomic sequences acquired from GISAID to generate the sequence alignment, sampling time, and sampling geographic area data to be used in the analyses (`scripts/sequence_data_curation`),
* the scripts that perform both stochastic mappings and posterior-predictive simulations using the inferred geographic model parameters and dated phylogenies (`scripts/history_simulation`), and
* the scripts that process the `BEAST` output files to produce summaries of geographic model parameters (`scripts/parameter_summary`).

### Sequence Data Curation
The first step of sequence data curation is to fliter the raw sequences downloaded from GISAID those that undermine the analyses due to various reasons (*e.g.*, a sequence that is incomplete or lacks the associated sampling time or location information).
The specific procedures differ between the reduced dataset and the entire dataset slightly, since we downloaded unaligned sequences from GISAID to generate the reduced dataset, but sequence alignment to generate the entire dataset.
We present both procedures here so that either way can be replicated.
For the reduced dataset (when the initial input are raw sequences), `scripts/sequence_data_curation/sequence_preprocess_beforealign.R` is run first as a preliminary filteration of the sequences, and then (after the alignment is inferred using the output of the preliminary filteration), `scripts/sequence_data_curation/sequence_process_aftealign.R` is run to filter the inferred alignment again (as some filteration conditions assume that site homology has been established).
For the entire dataset (when the initial input is sequence alignment), `scripts/sequence_data_curation/sequence_process_gisaidalignment.R` is run to filter it in one installment.
These scripts also trim the alignment after the filteration to generate nucleotide alignments that only keep the coding regions (as well as the translated amino-acid alignments).

After the curation step is done, `scripts/sequence_data_curation/sequence_process_alignment_foranalysis.R` can be run to partition the entire alignment according to various partition schemes (*e.g.*, by codon positions and/or by gene regions), so that different substitution model among the partitions can be easily specified in the analyses.
Lastly, `scripts/sequence_data_curation/discretetrait_writer.R` generates a data table containing the sampling age (in unit of days) and sampling geographic area (according to the area-grouping scheme defined in `scripts/sequence_data_curation/main_config.R`) information of each virus sample included in the processed sequence alignment.

### History Simulations
We simulated full dispersal history over the sampled dated phylogeny for two major purposes: inferring the posterior distribution of histories (through stochastic mapping) and assessing the adequacy of geographic models (through posterior-predictive simulation).
`scripts/sequence_data_curation/history_simulator_functions.R` is the interface script that a user can execute the `history_simulator` rountine within it to either perform stochastic mapping or posterior-predictive simulation (by setting the argument `conditional` true or false).
Other scripts in this subdirectory contains subroutines that will be evoked by the interface rountine.
One exception is `scripts/sequence_data_curation/posterior_predictive_teststatistics_functions`, which provides functions to calculate posterior-predictive summary statistics (including the tipwise multinomial statistic and the parsimony statistic).

### Geographic Model Parameter Summary
`scripts/sequence_data_curation/get_BFs_counts.R` processes the `BEAST` parameter log file and simulated history log file to summarize the geographic model parameter estimates for both constant and piecewise constant model, including (interval-specific) pairwise dispersal routes, (interval-specific) global dispersal rate, and (interval-specific) number of dispersal events.
`scripts/sequence_data_curation/get_daily_events_nlieanges.R` processes the `BEAST` simulated history log file (or the tree log file) to summarize the number of dispersal events occurred on each dispersal route (or the number of active viral lineages) in each arbitrarily specified time interval.

## <a name="program"></a>Program
We included an executable `BEAST` program compiled from the `BEAST` source code we modified to extend its machinery (that was designed for the constant model) to work for the piecewise constant model (including both phylogenetic likelihood computation and stochastic mapping) as well.
Our modification of the source code was based on [the `master` branch, commit `c0e5f8d3a`](https://github.com/beast-dev/beast-mcmc/commit/c0e5f8d3af1f943bed5754d8e6466aab4eb776f4); the modified source code (which the executable that we used in our analyses in this study was compiled from) can be found as [the `epochal_model_extension` branch, commit `49e476e41` of this forked `BEAST` source code repository](https://github.com/jsigao/beast-mcmc/commit/49e476e4198d6328933d3d8fcaf3923374eaf5e2).
