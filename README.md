# Supplementary Archive for: New Interval-Specific Phylodynamic Models Improve Inference of the Geographic History of Disease Outbreaks
This archive contains the necessary materials to replicate the analyses and reproduce the results of our study.
These materials are located within four directories: 
(1) the [`simulation_study`](#simulation_study) directory contains the `BEAST` XML scripts to specify the phylodynamic analyses that we performed to obtain parameter values that we subsequently used to simulate geographic datasets and also the scripts to analyze those simulated datasets; 
(2) the [`empirical_application`](#empirical_application) directory contains both the SARS-CoV-2 datasets that we use to illustrate our new methods and the XML scripts that specify the phylodynamic analyses of these empirical datasets that we performed using the `BEAST` software package; 
(3) the [`scripts`](#scripts) directory contains the `R` scripts that we used to curate empirical data and to process the output from our analyses, and; 
(4) the [`program`](#program) directory contain a modified version of the `BEAST` software package that implements various extensions required to perform some of the phylodynamic analyses in our study.
Each of these directories is detailed below.
Note that we also provide this archive as a [Dryad repository](https://datadryad.org/stash/share/qM_lKSgUqZ9jk3yKU5Fo5JurDb4iichKAnkt4PPh6Wg).

## <a name="simulation_study"></a>Simulation Study
The `simulation_study` directory contains the `BEAST` XML scripts to specify phylodynamic analyses of the empirical dataset that we performed to obtain realistic parameter values for simulating geographic datasets (within the [`generating_simulated_datasets`](#generating_simulated_datasets) subdirectory) and also the scripts for subsequently analyzing those simulated datasets (within the [`analyzing_simulated_datasets`](#analyzing_simulated_datasets) subdirectory).
These analyses condition on the MCC summary phylogeny inferred in [empirical application step 1](#step1_dated_phylogeny_inference_analyses).

### <a name="generating_simulated_datasets"></a> Generating Simulated Datasets
Our simulation study explored a region of parameter space that is centered on the joint posterior probability distribution of phylodynamic model parameters estimated from our empirical analyses.
To that end, we first used the `BEAST` XML scripts in the `generating_simulated_datasets` subdirectory to analyze our empirical dataset (the reduced SARS-CoV-2 dataset, see [the empirical-data section](#data)) under each of two models, <!-- $1\mu1\mathbf{Q}$ and $2\mu2\mathbf{Q}$,  -->
time-constant and interval-specific (2 intervals), and then centered the parameter values of our simulation on the resulting posterior median estimates of the corresponding parameters (using the `R` scripts provided in [the history-simulation-script section](#history_simulation_scripts)).

### <a name="analyzing_simulated_datasets"></a> Analyzing Simulated Datasets
<!-- We performed separate analyses of each simulated dataset under each of the two models, resulting in four true:inference model combinations (${1\mu1\mathbf{Q}{:}1\mu1\mathbf{Q}}$, ${2\mu2\mathbf{Q}{:}2\mu2\mathbf{Q}}$, ${1\mu1\mathbf{Q}{:}2\mu2\mathbf{Q}}$, and ${2\mu2\mathbf{Q}{:}1\mu1\mathbf{Q}}$). -->
We estimated the posterior distribution for each simulated dataset separately under each inference model by performing MCMC simulation using `BEAST`; we provide example XML scripts for this analysis in the subdirectory `analyzing_simulated_datasets/posterior_predictive`.
We used the resulting posterior estimates in two ways: 
(1) to evaluate the statistical behavior of each inference model (by computing the coverage probability and absolute error for the model parameters and the pairwise and total number of dispersal events), and; 
(2) to perform posterior-predictive simulation (using the `R` scripts provided in [the history-simulation-script section](#history_simulation_scripts)) to assess the *absolute fit* (adequacy) of each model.
We also assessed the *relative fit* of competing inference models by computing Bayes factors based on their marginal-likelihood estimates.
We estimated the marginal likelihood for each inference model by performing power-posterior analyses for each simulated dataset using `BEAST`; we provide example XML scripts for this analysis in the subdirectory `analyzing_simulated_datasets/marginal_likelihood`.

## <a name="empirical_application"></a>Empirical Application
The `empirical_application` directory contains the genomic, geographic and sampling data used in our empirical phylodynamic analyses (within the [`data`](#data) subdirectory) and the `BEAST` XML scripts to specify these phylodynamic analyses (within the [`analyses`](#analyses) subdirectory).

### <a name="data"></a>Data
The `data` subdirectory is divided into three subdirectories, one for each primary type of data used in our study:
* SARS-CoV-2 genomic sequence data ([`data/sequence_data`](#sequence_data)),
* air-travel data ([`data/travel_data`](#travel_data)), and
* intervention-measure data ([`data/measure_data`](#measure_data)).

#### <a name="sequence_data"></a>SARS-CoV-2 Sequence Data
Our study involves analyses of two SARS-CoV-2 genomic datasets: the *reduced dataset*, with 1271 sequences, and the *entire dataset*, with 2598 sequences.
The reduced dataset is based on all SARS-CoV-2 genomic sequences&mdash;collected during the early phase of the COVID-19 pandemic&mdash;that were deposited in [GISAID](https://www.gisaid.org) as of April 19, 2020.
The entire dataset is based on all SARS-CoV-2 sequences collected during the early phase of the COVID-19 pandemic that were deposited in GISAID as of September 22, 2020.

GISAID policy prohibits us from directly sharing the `FASTA` files for the SARS-CoV-2 genomic sequences used in our study.
Instead, we provide the corresponding GISAID accession numbers of those SARS-CoV-2 sequences: sequences comprising the reduced dataset are listed in `data/sequence_data/gisaid_acknowledgement_table_041920_030820.tsv`, and sequences comprising the entire dataset are listed in `data/sequence_data/gisaid_acknowledgement_table_092220_030820.tsv`.
The metadata for each sequence&mdash;including the sampling date and geographic location&mdash;are listed in `data/sequence_data/metainfo_041920.csv` for the reduced dataset, and `data/sequence_data/metadata_092220.tsv` for the entire dataset.
The accession numbers of sequences that we excluded during the data-curation stage (based on our filtering criteria) are listed in `data/sequence_data/exclude_seqs.tsv` (modified from [Nextstrain ncov exclude list](https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt)).

#### <a name="travel_data"></a>Air-Travel Data
We used the daily number of commercial passenger flights obtained from FlightAware as a proxy for the global passenger travel volume; this dataset is available in `data/travel_data/global_flights/nflights_daily_byaircraft.csv`.
This spreadsheet includes the number of all commercial passenger flights (with the number of flights for each type of aircraft) on each day between December 30, 2019 and March 8, 2020.
We summarize data on the passenger occupancy (number of seats) in each type of aircraft in `data/travel_data/global_flights/aircraft_nseats.csv`.
We then estimated the daily air-travel volume (assuming maximum capacity for all flights) as follows: first, we compute the product of the number of flights for each type of aircraft times the passenger occupancy for that type of aircraft, then we sum over all aircraft types.

#### <a name="measure_data"></a>Intervention Measures
We focused on two types of intervention measures involving China during the early phase of COVID-19 pandemic: targeted-containment measures (*i.e.*, international bans on air-travel with China), and the domestic-mitigation measures within China.
The table `data/measure_data/international_airtravelban_withchina.csv` contains a collection of countries or territories that imposed international travel bans with China, as well as the associated initiation date and information source.
The table `data/measure_data/china_domestic.csv` contains a collection of provinces or cities in China that imposed domestic-mitigation measures, as well as the associated implementation period, type of intervention measure, and information source.

### <a name="analyses"></a> Phylodynamic Analyses
The `analyses` subdirectory contains the XML scripts that we used to perform the empirical phylodyamic analyses using the `BEAST` software package.
We include a subdirectory for each stage of analyses that we performed in our study:
* analyses to infer a dated phylogeny for the reduced SARS-CoV-2 dataset ([`analyses/step1_dated_phylogeny_inference`](#step1_dated_phylogeny_inference_analyses));
* analyses to evaluate the fit of candidate biogeographic models to the reduced SARS-CoV-2 dataset ([`analyses/step2_model_evaluation`](#step2_model_evaluation_analyses));
* analyses to infer the phylogeny, divergence times, and biogeographic history for the entire SARS-CoV-2 dataset  ([`analyses/step3_joint_analyses`](#step3_joint_analyses)); and
* ancillary analyses to estimate daily global viral dispersal rates for the entire SARS-CoV-2 dataset ([`analyses/ancillary_daily_rate`](#ancillary_daily_rate)).

We performed all analyses using the `BEAST` software package (although the specific version varied among analyses, see details below), with the `BEAGLE` library (compiled from [the `hmc-clock` branch, commit `dd36bf5`](https://github.com/beagle-dev/beagle-lib/tree/dd36bf5b8d88348c77a93eeeef0917d90df71a4f)) enabled to accelerate both CPU and GPU computation.
Note that, because GISAID policy prohibits us from directly sharing the sequence data, the XML scripts that we include in this repository include `?` as placeholders for each sequence. 
In order to reproduce our analyses, it is therefore necessary to: 
(1) download the sequences from GISAID using the GISAID accession numbers that we provide in the [`data/sequence_data`](#sequence_data) subdirectory; 
(2) run the `R` scripts that we provide in the [`scripts/sequence_data_curation`](#sequence_data_curation_scripts) subdirectory to recreate our curated sequence alignments, and; 
(3) replace the `?`s with the corresponding sequences to replicate our XML scripts.
To improve clarity, the name of each sequence in the XML scripts comprises three parts, denoting the GISAID accession ID, sampling geographic area, and sampling date.
For example, the sequence listed under GISAID ID `EPI_ISL_402124` was sampled in China on December 30, 2019, so its name in the XML scripts is thus `402124_China_123019`.

#### <a name="step1_dated_phylogeny_inference_analyses"></a>Step 1: Estimating the dated phylogeny for the reduced SARS-CoV-2 dataset
We inferred the dated phylogeny for the reduced SARS-CoV-2 dataset using the XML script `data/step1_dated_phylogeny_inference/coalExp_ucln_posterior_run1.xml`.
We summarized the tree inferred from these analyses as a Maximum Clade Credibility (MCC) summary tree, which is located in `data/step1_dated_phylogeny_inference/coalExp_ucln_posterior_MCC.tre`.
We performed these analyses using [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1).

#### <a name="step2_model_evaluation_analyses"></a>Step 2: Evaluating candidate geographic models using the reduced SARS-CoV-2 dataset
We assessed both the *relative* and *absolute* fit of each of nine candidate geographic models to our reduced SARS-CoV-2 dataset.
These models assign interval-specific parameters&mdash;for the average rate of viral dispersal and/or relative rates of viral dispersal&mdash;to one, two, four, or five pre-specified time intervals.
<!-- These models assign interval-specific parameters&mdash;for the average rate of viral dispersal, $\mu$, and/or relative rates of viral dispersal, $\mathbf{Q}$&mdash;to one, two, four, or five pre-specified time intervals; *i.e.*, ${1\mu1\mathbf{Q}}$, ${1\mu2\mathbf{Q}}$, ${2\mu1\mathbf{Q}}$, ${2\mu2\mathbf{Q}}$, ${1\mu4\mathbf{Q}}$, ${4\mu1\mathbf{Q}}$, ${4\mu4\mathbf{Q}}$, ${5\mu5\mathbf{Q}}$, and ${5\mu5\mathbf{Q}^*}$.
We specified interval boundaries based on external information regarding events within the study period that might plausibly impact viral dispersal dynamics (*e.g.*, onset of international air-travel restrictions against China on February 2).
The final candidate model, ${5\mu5\mathbf{Q}^*}$, includes five arbitrary and uniform (bi-weekly) intervals. -->

We assessed the *relative fit* of competing geographic models by computing Bayes factors based on their marginal-likelihood estimates.
We estimated the marginal likelihood for each candidate geographic model by performing power-posterior analyses; these XML scripts are located in the subdirectory `analyses/step2_model_evaluation/marginal_likelihood`.
We also assessed the *absolute fit* of each candidate geographic model using posterior-predictive simulation.
We estimated the posterior distribution for each candidate geographic model by performing MCMC simulation; these XML scripts are located in the subdirectory `analyses/step2_model_evaluation/posterior_predictive`.

Our evaluation of candidate geographic models condition on the MCC summary phylogeny inferred in [step 1](#step1_dated_phylogeny_inference_analyses).
Our analyses under the time-constant geographic models were performed using [`BEAST` version 1.10.5](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.5pre1), whereas those under the interval-specific geographic models were performed using our modified version of `BEAST` (see the [program section](#program) for details).

#### <a name="step3_joint_analyses"></a>Step 3: Jointly inferring the phylogeny, divergence times and geographic history for the entire SARS-CoV-2 dataset
We inferred the phylogeny, divergence times and geographic history for the entire SARS-CoV-2 dataset.
The phylodynamic model that we used for these analyses included the previously selected geographic (sub)model (see [step 2](#step2_model_evaluation_analyses)).
We performed these joint phylodynamic analyses using our modified version of `BEAST`; the corresponding XML files are located in `analyses/step3_joint_analyses`.
We provide the inferred posterior distribution of phylodynamic model parameters in our Dryad repository (file-size constraints prevent us from including it in the GitHub repository).

#### <a name="ancillary_daily_rate"></a>Ancillary analyses: Estimating daily global viral dispersal rates for the entire SARS-CoV-2 dataset
We estimated the daily global viral dispersal rate for the entire SARS-CoV-2 dataset under a more granular phylodynamic model that allows the average dispersal rate to vary from day to day to assess the correlation between daily global viral dispersal rates and daily global air-travel volume during the early phase of the COVID-19 pandemic.
We performed these analyses using our modified version of `BEAST`; the corresponding XML scripts are located in `analyses/ancillary_daily_rate`.
To accommodate phylogenetic uncertainty, we performed these analyses while averaging over (a subsample of) the marginal posterior distribution of dated phylogenies inferred previously in step 3 ([`joint analyses`](#step3_joint_analyses)); we include this tree file in our Dryad repository (file-size constraints prevent us from including it in the GitHub repository).

## <a name="scripts"></a>Scripts
The `scripts` directory contains the `R` scripts that we used in this study to:
* curate SARS-CoV-2 genomic sequences sourced from GISAID to generate the sequence alignments, with the sampling dates and geographic locations of each sequence (see [`scripts/sequence_data_curation`](#sequence_data_curation_scripts)),
* (1) perform stochastic mapping, (2) perform posterior-predictive simulation, and (3) to simulate geographic datasets, using the inferred geographic-model parameters and dated phylogeny (see [`scripts/history_simulation`](#history_simulation_scripts)), and;
* process `BEAST` output files to generate summaries of various geographic-model parameters (see [`scripts/parameter_summary`](#parameter_summary_scripts)).
* compute marginal-likelihood estimates and posterior-predictive summary statistics to assess geographic model fit (see [`scripts/model_assessment`](#model_assessment_scripts)).

### <a name="sequence_data_curation_scripts"></a>Curating SARS-CoV-2 Sequence Data
Curation of the SARS-CoV-2 sequence data involved filtering the raw sequences downloaded from GISAID for various reasons (*e.g.*, because a sequence was incomplete and/or lacked the associated sampling dates or locality metadata).
Our curation procedure differed slightly between the reduced and entire SARS-CoV-2 datasets, as we downloaded unaligned sequences from GISAID to generate the reduced dataset, but downloaded aligned sequences to generate the entire dataset.
Specifically, processing of the reduced dataset (where we started from unaligned sequences), we first applied a preliminary filter using the script `scripts/sequence_data_curation/sequence_preprocess_beforealign.R`, and&mdash;after aligning the sequences&mdash;applied a secondary filter using the script `scripts/sequence_data_curation/sequence_process_aftealign.R`.
By contrast, we processed the entire dataset (where we started from aligned sequences) by applying a single filter using the script `scripts/sequence_data_curation/sequence_process_gisaidalignment.R`.
These scripts also trim the alignment after filtration to generate nucleotide alignments that only retain coding regions.

After the initial curation, we then defined data subsets within the curated alignments using the script `scripts/sequence_data_curation/sequence_process_alignment_foranalysis.R`. 
Defining subsets of the alignments&mdash;*e.g.*, subdividing sites in the alignment by codon position and/or by gene regions&mdash;allows us to specify independent substitution models for each subset of the alignment to accommodate variation in the substitution process across the SARS-CoV-2 genomes.
The script `scripts/sequence_data_curation/discretetrait_writer.R` generates a data table that contains the sampling date (in units of days) and sampling locality (according to the discrete geographic areas defined in `scripts/sequence_data_curation/main_config.R`) for each sequence in the alignment.
Finally, the script `scripts/sequence_data_curation/generate_xml_data_head.R` reads in the partitioned sequence alignment and the sampling-information table (*e.g.*, `data/sequence_data/metainfo_041920.csv`) to generate the `taxa` part of the XML scripts that are provided in the `analyses` directory.

### <a name="history_simulation_scripts"></a>Simulating Geographic Histories
We simulated dispersal histories over sampled phylogenies to: 
(1) infer the posterior distribution of geographic histories (using stochastic mapping); 
(2) assess the adequacy of geographic models (using posterior-predictive simulation), and; 
(3) generate geographic datasets for our simulation study (using forward-in-time simulation).
The script `scripts/history_simulation/history_simulator_functions.R` contains the main routine `history_simulator` to achieve these goals.
We performed stochastic mapping by setting the argument `conditional` to `true` (setting it to the default value, `false`, for the other two types of simulations).
We generated geographic datasets for our simulation study by setting the argument `true_value_type` to `median` as we used the posterior median estimates as the true parameter values (setting it to the default value, `sample`, for the other two types of inferences because they simulate histories by sampling from the joint posterior distribution).
Other scripts in this subdirectory include subroutines that are executed by `history_simulator_functions.R` and are thus not intended to be used directly.

### <a name="parameter_summary_scripts"></a>Summarizing Geographic Model Parameters of Phylodynamic Analyses
We used the script `scripts/parameter_summary/get_BFs_counts.R` to process the parameter log file output by `BEAST` (and simulated-history log file, if available) to summarize the geographic-model parameters for both time-constant and interval-specific geographic models.
These summaries include (interval-specific) support (as Bayes factor) for pairwise dispersal routes, (interval-specific) global average dispersal rate, (interval-specific) relative dispersal rates, and (interval-specific) number of pairwise dispersal events.
We used the script `scripts/parameter_summary/get_nevents_daily.R` to process a simulated-history log file output by `BEAST` to summarize the number of dispersal events that occurred over each dispersal route during each arbitrarily specified time interval (*e.g.*, each day).
The script `scripts/parameter_summary/simparams_summary_helpr.R` processes the output produced by `get_BFs_counts.R` for the simulated datasets to summarize the statistical behavior (*e.g.*, the coverage probability of the average dispersal rate or the total number of dispersal events) of a geographic model.

### <a name="model_assessment_scripts"></a>Assessing Geographic Model Fit
We used the `scripts/model_assessment/marginal_likelihood_calculator.R` script to compute a combined log marginal likelihood from replicate power-posterior log files produced by `BEAST`.
The script `scripts/model_assessment/posterior_predictive_summarystatistics.R` calculates (time-slice) posterior-predictive summary statistics (including the tip-wise multinomial statistic and the parsimony statistic) from the posterior-predictive distribution simulated by the `history_simulator` function provided in [`scripts/history_simulation`](#history_simulation_scripts).

## <a name="program"></a>Modified BEAST Program
We provide an executable for a modified version of the `BEAST` program; our modifications implement various extensions to enable analyses under interval-specific phylodynamic models, including likelihood computation and stochastic mapping.
Our modifications were made to the source code of [the `master` branch, commit `58c4e71`](https://github.com/beast-dev/beast-mcmc/commit/58c4e7154a8028d9c5095831a15aebdd7f74df3a); the modified source code (from which we compiled the executable used in several of our analyses) is located in [the `master` branch, commit `d1ee6a4`](https://github.com/jsigao/beast-mcmc/commit/d1ee6a495d2f3d24cdb52a85f5d622d449924b67) of a forked `BEAST` source-code repository.
Briefly, this executable can be run by invoking:
```
java -jar beast.jar -working ./empirical_application/analyses/step2_model_evaluation/marginal_likelihood/4mu4q_exphyper_asymmetric_poissonintermediate_powerposterior_run1.xml
```
See [`BEAST` official tutorial website](http://beast.community/index.html) for further details about running `BEAST` analyses using a `Java` executable.
