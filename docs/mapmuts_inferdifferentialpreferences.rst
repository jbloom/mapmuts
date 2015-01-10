.. _mapmuts_inferdifferentialpreferences.py:

========================================
mapmuts_inferdifferentialpreferences.py
========================================

**This program is no longer the state-of-the-art:** You are suggested to instead
 use the `dms_tools`_ program ``dms_inferdiffprefs``.

This script is designed to infer the differential preferences for each amino acid at each site when comparing two or more selection pressures. It does this using a Bayesian approach implemented via MCMC (Markov chain Monte Carlo).

Because this script uses MCMC, it takes a fairly long time to run -- perhaps as long as several days depending on your computer and the specific data set.

This script takes as input the *codoncounts.txt* files created by :ref:`mapmuts_parsecounts.py`, and so is designed to be run after you have completed the analysis with that script.

Dependencies
--------------
This script requires `pymc`_ for the Bayesian inference, which in turn requires `numpy`_. If you want the script to make plots, you also need `matplotlib`_. The use of `pymc`_ version 2.3 is recommended. This script may **not** work with versions 3.*, and may work suboptimally with earlier versions.

When would you use this script?
----------------------------------
If you just want to infer the site-specific amino-acid preferences by comparing a pre-selection mutant library to a post-selected sample, then you should probably instead use :ref:`mapmuts_inferpreferences.py`, which implements the algorithm described in `An experimentally determined evolutionary model dramatically improves phylogenetic fit`_ (in the Methods section entitled *Inference of the Amino Acid Preferences*).

This script is tailored to the case where you have several samples that have been selected in slightly different conditions, and you want to identify sites and amino acids that are selected differently in the different conditions. Most commonly, you would have already selected for basic functionality, and are now looking for differential functionality under various selections. This script identifies the differential preferences in these different conditions. The algorithm here is still similar in many respects to that described in `An experimentally determined evolutionary model dramatically improves phylogenetic fit`_, but now is elaborated to infer the differential preferences.

The samples that compose the data
------------------------------------
We assume that we are sequencing the following samples:

* *error_control*: This sample is used to measure the error rate associated with sequencing and sample processing. It will generally be the unmutated gene or a wildtype virus grown from unmutated gene, whichever better captures the errors that afflict the samples. For influenza, the errors (or unaccounted for mutations) generally come from mutations during viral replication, errors during reverse transcription, errors during PCR when preparing the sample, and sequencing errors. So here, the appropriate sample is generally wildtype virus grown from an unmutated gene.

* *starting_sample*: This is the sample that was the starting point for the selections. It could be your mutant library, but more typically it will be a mutant virus pool that you have already grown from your mutant library and are now subjecting to further selection. An implicit assumption of the way that the error correction with *error_control* is done is that at all sites this sample are still mostly wildtype (the mutations are relatively rare).

* *control_selection*: This is *starting_sample* passaged in your control condition. For instance, this might be virus that is passaged a second time in a control cell line or animal. An implicit assumption of the way that the error correction with *error_control* is done is that all sites in this sample are still mostly wildtype (the mutations are relatively rare). 

* *selection_1*: This is *starting_sample* passaged in the first selection condition. For instance, this might be virus that is passaged a second time under immune selection. An implicit assumption of the way that the error correction with *error_control* is done is that all sites in this sample are still mostly wildtype (the mutations are relatively rare).

* *selection_2*, ..., *selection_i*: These are an arbitrary number of additional samples like *selection_1* that are passaged in other selection conditions.

Let :math:`n^{\rm{error}}_{r,x}`, :math:`n^{\rm{start}}_{r,x}`, :math:`n^{\rm{control}}_{r,x}`, and :math:`n^{s_i}_{r,x}` be the observed number of counts of codon :math:`x` at site :math:`r` in the sequencing of *error_control*, *starting_sample*, *control_selection*, and *selection_i* samples, respectively. Define :math:`\overrightarrow{n}^{\rm{error}}_{r}` as the vector of :math:`n^{\rm{error}}_{r,x}` values for all codons :math:`x`, and make similar definitions for :math:`\overrightarrow{n}^{\rm{start}}_{r}`, :math:`\overrightarrow{n}^{\rm{control}}_{r}`, and :math:`\overrightarrow{n}^{s_i}_{r}`. Define :math:`N^{\rm{error}}_{r} = \sum_x n^{\rm{error}}_{r,x}` be the total number of counts of any codon at site :math:`r` in the *error_control* sample, and make similar definitions for :math:`N^{\rm{start}}_{r}`, :math:`N^{\rm{control}}_{r}` and :math:`N^{s_i}_{r}`.

The likelihoods in terms of unknown parameters
-------------------------------------------------
The *error_control* which is just sequencing of the wildtype gene. Let :math:`\rm{wt}_r` denote the wildtype codon at site :math:`r`. Let :math:`\epsilon_{r,x}` be the rate that site :math:`r` is read to be some mutant codon :math:`x \ne \rm{wt}_r`, and let :math:`\epsilon_{r,\rm{wt}_r} = 1 - \sum_{y \ne \rm{wt}_r} \epsilon_{r,x}`. Define :math:`\overrightarrow{\epsilon}_r` as the vector of :math:`\epsilon_{r,x}` values for all :math:`x`. The likelihood of observing :math:`\overrightarrow{n}_{\rm{error},r}` is 

.. math::   

   \Pr\left(\overrightarrow{n}^{\rm{error}}_{r} \mid N^{\rm{error}}_{r}, \overrightarrow{\epsilon}_r\right) = \operatorname{Multi}\left(\overrightarrow{n}^{\rm{error}}_{r}; N^{\rm{error}}_{r}, \overrightarrow{\epsilon}_r\right)

where :math:`\operatorname{Multi}` denotes the multinomial distribution.

The *starting_sample* gives the starting frequencies for the various variants before the selections. We assume that all sites are mostly the wildtype codon, and so account for errors from the wildtype codon to mutant codons (as quantified by the *error_control* sample), but do not account for errors from mutant codons to other codons. Let :math:`f_{r,x}` be the frequency of codon :math:`x` at site :math:`r` with :math:`\sum_x f_{r,x} = 1`, and let :math:`\overrightarrow{f}_r` be the vector of :math:`f_{r,x}` values for all codons :math:`x`. The likelihood of observing :math:`\overrightarrow{n}^{\rm{start}}_{r}` is 

.. math::   \Pr\left(\overrightarrow{n}^{\rm{start}}_{r} \mid N^{\rm{start}}_{r}, \overrightarrow{\epsilon}_r, \overrightarrow{f}_r\right) = \operatorname{Multi}\left(\overrightarrow{n}^{\rm{start}}_{r}; N^{\rm{start}}_{r}, \overrightarrow{\epsilon}_r + \overrightarrow{f}_r - \overrightarrow{\delta}_r\right)

where :math:`\overrightarrow{\delta}_r` is the vector with elements :math:`\delta_{x, \rm{wt}_r}` that has one in the element corresponding to the wildtype codon and zero for all othere elements.

The *control_selection* gives the frequencies after the *starting_sample* has been passaged in the control condition. These frequencies may change substantially depending on the selection pressures that generated *starting_sample* versus those acting on the *control_selection*. We assume that the selection acts entirely on the amino-acid sequence, such that all synonymous codons that encode the same amino-acid are equivalent. Let :math:`\mathcal{A}\left(x\right)` denote the amino acid encoded by codon :math:`x`. Let :math:`\pi_{r,a}` be the preference for amino-acid :math:`a` at site `r` in the *control_selection*, with the constraint that :math:`\sum_a \pi_{r,a} = 1`. Let :math:`\overrightarrow{\pi}_r` be the vector of :math:`\pi_{r,a}` values for all amino acids :math:`a`.
Define the vector-valued function :math:`\overrightarrow{\mathcal{C}}` as a mapping for the value of each codon :math:`x` to the corresponding value for its amino acid :math:`\mathcal{A}\left(x\right)`, such that :math:`\overrightarrow{\mathcal{C}}\left(\overrightarrow{\pi}_{r}\right) = \left(\ldots, \pi_{r,\mathcal{A}\left(x\right)}, \ldots\right)` for all codons :math:`x`. The expected frequency of codon :math:`x` at site :math:`r` in the control sample is defined in terms of the preferences as 
:math:`\frac{f_{r,x} \times \pi_{r,\mathcal{A}\left(x\right)}}{\sum_x f_{r,x} \times \pi_{r,\mathcal{A}\left(x\right)}}`. 
Therefore, if we again assume that all sites in the *control_selection* sample are mostly wildtype so that we only account for errors from the wildtype codon to some other codon, then the likelihood of observing :math:`\overrightarrow{n}^{\rm{control}}_{r}` is

.. math:: 

   \Pr\left(\overrightarrow{n}^{\rm{control}}_{r} \mid N^{\rm{control}}_{r}, \overrightarrow{\epsilon}_r, \overrightarrow{f}_r, \overrightarrow{\pi}_r\right) 
   = \operatorname{Multi}\left(\overrightarrow{n}^{\rm{control}}_{r}; N^{\rm{control}}_{r}, \overrightarrow{\epsilon}_r
   + \frac{\overrightarrow{f}_r \circ \overrightarrow{\mathcal{C}}\left(\overrightarrow{\pi}_r\right)}{\overrightarrow{f}_r \cdot \overrightarrow{\mathcal{C}}\left(\overrightarrow{\pi}_r\right)} - \overrightarrow{\delta}_r \right)

where :math:`\cdot` indicates the vector dot product and :math:`\circ` indicates the Hadamard (entry-wise) vector product.

We are most interested in how the other selections differ from the control selection -- in particular whether there is increased or decreased selection for some amino acids at certain sites. Let :math:`\Delta\pi^{s_i}_{r,a}` denote the difference in the preference for amino-acid :math:`a` at site :math:`r` for *selection_i* (denoted by :math:`s_i`). If *selection_i* is identical to the control selection, then we expect these differences in preferences to all be zero. Even when they are not zero, we have the constraint :math:`\sum_a \Delta\pi_{s_i}^{r,a} = 0`, since any increase in preference for one amino acid must be compensated by a decrease in preference for other amino acids. Let :math:`\overrightarrow{\Delta\pi}^{s_i}_{r}` be the vector of :math:`\Delta\pi^{s_i}_{r,a}` values for all amino acids :math:`a`. If we again assume that all sites in the *control_selection* sample are mostly wildtype so that we only account for errors from the wildtype codon to some other codon, then the likelihood of observing :math:`\overrightarrow{n}^{s_i}_{r}` is

.. math:: 

   \Pr\left(\overrightarrow{n}^{s_i}_{r} \mid N^{s_i}_{r}, \overrightarrow{\epsilon}_r, \overrightarrow{f}_r, \overrightarrow{\pi}_r, \overrightarrow{\Delta\pi}^{s_i}_{r}\right) 
   = \operatorname{Multi}\left(\overrightarrow{n}^{s_i}_{r}; N^{s_i}_{r}, \overrightarrow{\epsilon}_r 
   + \frac{\overrightarrow{f}_r \circ \overrightarrow{\mathcal{C}}\left(\overrightarrow{\pi}_r + \overrightarrow{\Delta\pi}^{s_i}_{r}\right)}
   {\overrightarrow{f}_r \cdot 
   \overrightarrow{\mathcal{C}}\left(\overrightarrow{\pi}_r + \overrightarrow{\Delta\pi}^{s_i}_{r}\right)} - \overrightarrow{\delta}_r \right)

Priors over the unknown parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We need to place priors over our unknown parameters, which are: :math:`\overrightarrow{\epsilon}_r`, :math:`\overrightarrow{f}_r`, :math:`\overrightarrow{\pi}_r`, and :math:`\overrightarrow{\Delta\pi}_{s_i,r}` for all selections :math:`s_i`.

The prior over :math:`\overrightarrow{\epsilon}_r` is defined exactly as in the 13th equation of `An experimentally determined evolutionary model dramatically improves phylogenetic fit`_. Briefly, the prior is

.. math::

   \Pr\left(\overrightarrow{\epsilon}_r\right) =
   \operatorname{Dir}\left(\overrightarrow{\epsilon}_r; n_{\rm{codons}} \cdot \sigma_{\epsilon} \cdot \overrightarrow{\alpha}_{\epsilon,r}\right)

where :math:`\operatorname{Dir}` is the `Dirichlet distribution`_, :math:`n_{\rm{codons}}` is the number of different codon identities (typically either 64 or 61, depending on whether stop codons are considered), :math:`\sigma_{\epsilon}` is the scalar `concentration parameter`_, and :math:`\overrightarrow{\alpha}_{\epsilon,r}` is the vector giving the average mutation rate for all codon mutations with that number of mutations over all sites in the *error_control* sample, exactly as defined in the description of the 13th equation of `An experimentally determined evolutionary model dramatically improves phylogenetic fit`_.

The prior over :math:`\overrightarrow{f}_r` is defined similarly to the prior over :math:`\overrightarrow{\mu}_r` defined in the 12th equation of `An experimentally determined evolutionary model dramatically improves phylogenetic fit`_. Specifically, the prior is

.. math::

   \Pr\left(\overrightarrow{f}_r\right) = 
   \operatorname{Dir}\left(\overrightarrow{f}_r; n_{\rm{codons}} \cdot \sigma_{f} \cdot \overrightarrow{\alpha}_{f,r}\right)

where :math:`\sigma_{f}` is the scalar `concentration parameter`_ and :math:`\overrightarrow{\alpha}_{f,r}` is the vector with elements :math:`\alpha_{f,r,x} = \frac{\overline{f}}{n_{\rm{codons}} - 1} + \delta_{x,\rm{wt}_r}\left(1 - \overline{f}\times\frac{n_{\rm{codons}}}{n_{\rm{codons}} -1}\right)` where :math:`\overline{f}` is the fraction of all reads in *start_sample* that do **not** give the wildtype codon :math:`\rm{wt}_r` averaged over all sites, and :math:`\delta_{x,\rm{wt_r}}` is the Kronecker-delta function which equal to one if codon :math:`x` is the wildtype codon :math:`\rm{wt}_r` and zero otherwise. In other words, the prior estimate is that all non-wildtype amino acids (mutations) are present at the average rate that any mutation is found over the entire gene. This is a reasonable uninformative prior.

The prior over :math:`\overrightarrow{\pi}_{r}` is defined exactly as in the 15th equation of `An experimentally determined evolutionary model dramatically improves phylogenetic fit`_. Specifically, the prior is

.. math::

   \Pr\left(\overrightarrow{\pi}_r\right) =
   \operatorname{Dir}\left(\overrightarrow{\pi}_r; \sigma_{\pi} \cdot \overrightarrow{1}\right)

where :math:`\sigma_{\pi}` is the scalar `concentration parameter`_ and :math:`\overrightarrow{1}` is a vector that is all ones. In other words, the prior is a symmetric `Dirichlet distribution`_, so the a prior estimate is that all amino acids are equally preferred at each site.

We also need to define priors over the differential preferences :math:`\overrightarrow{\Delta\pi}^{s_i}_r` under the alternative selections. We would like these priors to be such that the mean prior estimate for :math:`\Delta\pi^{s_i}_{r,a}` for all sites :math:`r` and amino acids :math:`a`; in other words, we would like the priors to not favor any tendency to estimate non-zero differential preferences. We also need the priors to enforce the conditions that :math:`\sum_a \Delta\pi^{s_i}_{r,a} = 0` and that :math:`0 \le \pi^{s_i}_{r,a} + \Delta\pi^{s_i}_{r,a} \le 1` for all :math:`a`. These conditions therefore require that the priors over :math:`\overrightarrow{\Delta\pi}^{s_i}_{r}` be defined in terms of :math:`\overrightarrow{\pi}_r`. Specifically, we define

.. math::

   \Pr\left(\overrightarrow{\Delta\pi}^{s_i}_r \mid \overrightarrow{\pi}_r\right) = \operatorname{Dir}\left(\overrightarrow{\Delta\pi}^{s_i}_r + \overrightarrow{\pi}_r; \sigma_{\Delta\pi} \cdot n_{\rm{aa}} \cdot \overrightarrow{\pi}_r \right) 

where :math:`\sigma_{\Delta\pi}` is the `concentration parameter`_ (assumed to be the same for all selections :math:`s_i`) and :math:`n_{\rm{aa}}` is the number of amino acids.

Inference with MCMC
-----------------------
The goal of the inference is to determine the differential preferences :math:`\overrightarrow{\Delta\pi}^{s_i}_r` for all selection pressures :math:`s_i`. In some cases we might also be interested in the preferences :math:`\overrightarrow{\pi}_r` for the *control_selection*.

The other parameters :math:`\overrightarrow{\epsilon_r}` and :math:`\overrightarrow{f}_r` are nuisance parameters, with values that must be calculated but are not of direct interest for the analyses performed by this script.

Finally, the `concentration parameter`_  used for each priors (:math:`\sigma_{\epsilon}`, :math:`\sigma_{f}`, :math:`\sigma_{\pi}`, and :math:`\sigma_{\Delta\pi}` are variables that must be specified a priori by the user of the script.

The posterior distributions over the parameters of interest are inferred by MCMC (Markov Chain Monte Carlo). Typically the values will be summarized by their posterior means. The sampling is done over the products of the likelihoods and priors specified in the sections above. Specifically, we sample from the following posterior:

.. math::

   \rm{posterior} \propto& 
   \Pr\left(\overrightarrow{n}^{\rm{error}}_r \mid N_r^{\rm{error}}, \overrightarrow{\epsilon}_r\right) \times \\
   & \Pr\left(\overrightarrow{n}^{\rm{start}}_{r} \mid N^{\rm{start}}_{r}, \overrightarrow{\epsilon}_r, \overrightarrow{f}_r\right) \times \\
   & \Pr\left(\overrightarrow{n}^{\rm{control}}_{r} \mid N^{\rm{control}}_{r}, \overrightarrow{\epsilon}_r, \overrightarrow{f}_r, \overrightarrow{\pi}_r\right) \times \\
   & \left(\prod_{s_i} \Pr\left(\overrightarrow{n}^{s_i}_{r} \mid N^{s_i}_{r}, \overrightarrow{\epsilon}_r, \overrightarrow{f}_r, \overrightarrow{\pi}_r, \overrightarrow{\Delta\pi}^{s_i}_{r}\right)\right) \times \\
   & \Pr\left(\overrightarrow{\epsilon}_r\right) \times \\
   & \Pr\left(\overrightarrow{f}_r\right) \times \\
   & \Pr\left(\overrightarrow{\pi}_r\right) \times \\
   & \left(\prod_{s_i} \Pr\left(\overrightarrow{\Delta\pi}^{s_i}_r \mid \overrightarrow{\pi}_r\right) \right)


Implementation of the MCMC
----------------------------
The posterior distributions for the differential preferences and preferences are estimated by MCMC using the `pymc`_ package. The script automatically tests for convergence using the *convergence* and *nruns* options. It is suggested that you use *nruns > 1* and *stepincrease > 1* so that you can try to ensure convergence. Substantial effort has been invested to make sure that by default the run will burn-in appropriately and set reasonable step sizes for the variables. The details of this implementation are documented in the source code. 

The log file created by this script will include information about whether the inference for each site converged according the criterium set by *convergence*. If it consistently fails to converge, you might consider increasing *nsteps*.

When examining convergence, the difference between the differential preferences and preferences inferred from different runs is quantified as the `Shannon-Jensen divergence`_. The logarithms are taken to the base 2, so this divergence can range from zero (identical inferred preferences) to one (completely independent inferred preferences). For the *control_selection*, the `Shannon-Jensen divergence`_ is computed using the preferences :math:`\pi_{r,a}`. For the other selections, the `Shannon-Jensen divergence`_ is computed using the sum of the preferences and differential preferences, :math:`\pi_{r,a} + \Delta\pi^{s_i}_{r,a}`.

Running the script
--------------------
To run this script from the prompt, first create a text `Input file`_ of the
format described below. Then simply run the script with the single text file as the argument, as in::

    mapmuts_inferdifferentialpreferences.py infile.txt


Input file
--------------
The input file is a text file with a series of *key* / *value* pairs. The required keys are indicated below. The values should not include spaces.

Lines beginning with # and empty lines are ignored.

Keys for the input file:

* *error_control* : This key should specify a file giving the counts of each codon at each position for the *error_control* sample; these are the :math:`\overrightarrow{n}^{\rm{error}}_r` values. This file should be in the format of a ``*_codoncounts.txt`` created by :ref:`mapmuts_parsecounts.py`. Typically, this sample would be from sequencing either an unmutated gene or a wildtype virus to estimate the baseline error rate.

* *starting_sample* : This key specifies a file giving the counts that are the :math:`\overrightarrow{n}^{\rm{start}}_r` values. This file is in the same format as for *error_control*. Typically, this sample would be the pool of functional mutant variants or viruses that were created from the mutant library and are now being subjected to the subsequent selections.

* *control_selection* : This key specifies a file giving the counts that are the :math:`\overrightarrow{n}^{\rm{control}}_r` values. This file is in the same format as for *error_control*. Typically, this sample would be the pool in *starting_sample* that has been subjected to the control growth condition, such as passage of a virus in control cells without any additional selection.

* *selection_???* : These keys specify the *selection_i* conditions used to look for differential selection relative to *control_selection*. There must be at least one of these keys, but there can be arbitrarily more. The actual key should take the form *selection_heat* or *selection_antibody* or *selection_* followed by any non-whitespace-containing string giving the name assigned to the additional selection. The name given to the selection condition in the results is whatever is the suffix that follows the underscore in the *selection_???* string. Typically, these samples would be *starting_sample* passaged in a condition that is being compared to *control_selection*. The values should be codon counts files in the same format as for *error_control*.

* *pi_concentration* : the `concentration parameter`_ of the symmetric `Dirichlet distribution`_  used as the prior over the preferences :math:`\pi_{r,a}` values. This is the parameter denoted as :math:`\sigma_{\pi}` above. Note that the concentration parameter here is taken to represent the value of each of the individual elements in the parameter vector. Choosing a value of *pi_concentration = 1* corresponds to a uniform prior over all of the possible vectors of *pi_concentration* values, and is probably what you want to choose unless you have a good reason to do otherwise. A value of *pi_concentration > 1* favors vectors where all *pi_concentration* values are similar, and a value of *pi_concentration < 1* favors vectors where some :math:`\pi_{r,a}` values are large and others are small.

* *epsilon_concentration* : the `concentration parameter`_ of the (non-symmetric) `Dirichlet distribution`_ used as a prior over the error rate observed in the *error_control* library. This is the parameter denoted as :math:`\sigma_{\epsilon}` above. The suggested value is one.

* *f_concentration* : the `concentration parameter`_ of the (non-symmetric) `Dirichlet distribution`_ used as a prior over the frequencies in the *starting_sample*. This is the parameter denoted as :math:`\sigma_{f}` above. The suggested value is one.

* *deltapi_concentration* : the `concentration parameter`_ of the (non-symmetric) `Dirichlet distribution`_ used as a prior over differential preferences. This is the parameter denoted as :math:`\sigma_{\Delta\pi}` above. The suggested value is one.

* *minvalue* : a number specifying the minimum value for the prior estimates for the error rates (:math:`\overrightarrow{\epsilon}_r`) determined from the overall library means. Having such a minimum value keeps any of the prior estimates from being set to zero or negative values. A reasonable value is *1e-7*.

* *seed* : integer seed for the random number generator. Runs with the same seed should generate exactly the same output. Otherwise the output may differ for different seeds as the MCMC uses random numbers.

* *nruns* : the number of independent MCMC runs. Must be at least one, but it is recommended that you use *nruns* equal to at least 2 (and probably 3), as this allows checks for MCMC convergence. When multiple runs are performed, the reported equilibrium preferences come from all runs.

* *nsteps* : the number of steps for each MCMC run. A reasonable starting value is probably 200000, although this may need to be increased if convergence is poor. Note that the number of burn-in steps are not explicitly specified in this script, but the program implements internal methods that should ensure adequate burn-in that is not included in the output. You will have to look at the source code if you want to understand that in detail. 

* *thin* : the posterior is only recorded every *thin* steps. If you want to record every value, then set *thin* to one. However, consecutive steps are usually correlated, so it is often preferable to set *thin* to a number greater than one. A suggested value is 200. *nsteps* must be a multiple of *thin*.

* *convergence* : the test applied to the MCMC runs to see if they converged. This option is only meaningful if *nruns* > 1 (as is recommended). For each pair of MCMC runs, we calculate the `Shannon-Jensen divergence`_ (logarithm base 2) between the preferences and differential preferences for that value of *r*. If the runs converge to identical values, this difference will be zero, whereas its maximum value for highly diverged runs is one. We test whether the actual divergence is less than the value specified by *convergence*. If the divergence is greater than *convergence* for any of the preferences or differential preferences, then what is done next depends on the values of *stepincrease*. A recommended value for *convergence* is 0.02. The log file contains information about whether the inference converged for each site.

* *stepincrease* : If the chain fails to converge for a given inference (according to the criterion of *convergence*), then the default action is to try to increase the number of steps to be equal to *stepincrease * nsteps*. We then test again for convergence for these longer MCMC runs. Typically increasing the step number might improve convergence. So *stepincrease* must be a number >= 1. If it is one, that is equivalent to not doing any *stepincrease*. If *stepincrease* is greater than one, we then try running the MCMC for *stepincrease * nsteps* steps for each run. If the chain still does not converge, nothing further is done, but we do use the inferences from the longer chain. A suggested value of *stepincrease* is 4.

* *MCMC_traces* is the directory to which we write plots showing the traces of the preferences and differential preferences as a function of the number of MCMC steps after burn-in and thinning. If you set this to the string *None*, then no such plots are created. Using this option requires `matplotlib`_. If this directory does not already exist, it is created.

* *preference_plots* is the directory to which we write plots showing the preferences and differential preferences for all amino acids at a site. Shown are the posterior mean and the median-centered 95% credible interval of the posterior. If you set this to the string *None*, then no such plots are created. Using this option requires `matplotlib`_. If this directory does not already exist, it is created.

* *outfileprefix* is a string giving the prefix for the `Output files`_ described below. If you do not want a prefix on the `Output files`_, then set this option to *None*.

* *ncpus* is an option that allows you to run the analyses using multiple CPUs. This will usually be helpful, because these runs can take a long time. If you set *ncpus* to 1, then we just use one CPU. But if you are working on a multi-threaded processor, then setting to a larger number will take advantage of the additional CPUs. Note however that if you set *ncpus* to a number larger than the number of spare CPUs, you will overload the processor -- so if you only have one CPU, set this to just 1.


Example input file
---------------------
Here is an example input file::

    # Input file for mapmuts_inferdifferentialpreferences.py
    error_control virus_codoncounts.txt
    starting_sample mutvirus_codoncounts.txt
    control_selection no_serum_codoncounts.txt
    selection_vaccinated_serum vaccinated_serum_codoncounts.txt
    selection_infected_serum infected_serum_codoncounts.txt
    pi_concentration 1.0
    epsilon_concentration 1.0
    f_concentration 1.0
    deltapi_concentration 1.0
    minvalue 1e-7
    seed 1
    nruns 3
    nsteps 200000
    thin 200
    convergence 0.02
    stepincrease 4
    MCMC_traces MCMC_traces/
    preference_plots preference_plots/
    outfileprefix None
    ncpus 4


Output files
----------------
Output files are created with the prefix specified by *outfileprefix* followed by suffixes shown below. These files have the following suffixes:

* ``inferdifferentialpreferences_log.txt`` : a text file logging the progress of this script. 

* ``preferences_control_selection.txt`` : a text file giving the MCMC-inferred preferences :math:`\pi_{r,a}` for each amino acid :math:`a` at each site :math:`r` for the *control_selection*. Stop codons (denoted by the * character) are included as possible identities. The file also gives the site entropy :math:`h_r = \sum_a \pi_{r,a} \log_2 \pi_{r,a}` in bits for each site. 

  The first line is a title line that begins with a # character and then lists the column headers separated by tabs. The remaining lines contain the data.

  The columns contain tab-separated values as follows:

    * The first column gives the residue number in 1, 2, ... numbering. The header for this column is *SITE*.
    
    * The second column gives the wildtype amino acid at the site (such as *N* or *R*). The header for this column is *WT_AA*.
    
    * The third column lists the site entropy :math:`h_r`. The header for this column is *SITE_ENTROPY*.
    
    * The remaining 21 columns list the preference :math:`\pi_{r,a}` values for each amino acid :math:`a` at this site :math:`r`. The amino acids are listed in alphabetical order according to their one letter codes, and then the final column gives the equilibrium preference for a stop codon (denoted by a * character). The headers for these columns are *PI_A*, *PI_C*, *PI_D*, ..., *PI_W*, *PI_Y*, *PI_**. Here is an example of a few lines::

        #SITE   WT_AA   SITE_ENTROPY    PI_A    PI_C    PI_D    PI_E    PI_F    PI_G    PI_H    PI_I    PI_K    PI_L    PI_M    PI_N    PI_P    PI_Q    PI_R    PI_S    PI_T    PI_V    PI_W    PI_Y    PI_*
        1   M   0.885694    0.000680047 0.00185968  8.48048e-06 0.00128217  0.0289946   0.000654589 0.0212144   0.0205306   0.0021597   0.00056726  0.876551    0.00130602  0.000789781 0.00144568  0.000414457 0.00071164  0.00126055  0.00167084  0.0352283   0.00176486  0.000905383
        2   A   1.1671  0.709575    0.00177747  3.79503e-05 0.00298854  0.00390524  0.00131246  0.00379599  0.000913006 0.00121497  0.000535209 0.00276335  0.00140061  0.000694694 0.0032134   0.00464972  0.000463924 0.000908368 0.254241    0.00335749  0.00137437  0.000876985
        3   N   3.2593  0.129992    0.000540239 4.2514e-06  0.00171887  0.0212062   0.0006086   0.102743    0.0431489   0.00524637  0.0771062   0.021714    0.0641921   0.00338082  0.0149898   0.0106522   0.13886 0.0247841   0.00607319  0.0500692   0.282389    0.000581555


* ``preferences_control_selection_credibleintervals_95.txt`` : a text file giving the median-centered 95% credible intervals on the :math:`\pi_{r,a}` values in ``preferences_control_selection.txt``. The first column is *SITE* and the remaining columns give the upper and lower bounds on the credible intervals for each of the :math:`\pi_{r,a}` values as comma-separated numbers (no spaces) in the same order that they are listed in ``preferences_control_selection.txt``. Here is an example::

        #SITE   PI_A_95cred PI_C_95cred PI_D_95cred PI_E_95cred PI_F_95cred PI_G_95cred PI_H_95cred PI_I_95cred PI_K_95cred PI_L_95cred PI_M_95cred PI_N_95cred PI_P_95cred PI_Q_95cred PI_R_95cred PI_S_95cred PI_T_95cred PI_V_95cred PI_W_95cred PI_Y_95cred PI_*_95cred
        1   1.49582e-05,0.00269598  7.2227e-05,0.00604361   4.37378e-06,1.17951e-05 3.26837e-05,0.00433571  0.0114465,0.0491834 1.82621e-05,0.00208322  0.00989126,0.0318659    0.00135281,0.0554989    4.87099e-05,0.00754094  1.41237e-05,0.00202465  0.821886,0.925528   3.69704e-05,0.00458556  1.94864e-05,0.0027967   3.89356e-05,0.00557027  8.44362e-06,0.00134572  1.76935e-05,0.0025792   3.51256e-05,0.00460279  3.68112e-05,0.00664131  0.0270458,0.0403275 5.19957e-05,0.00543017  2.3493e-05,0.00336828
        2   0.683813,0.729505   5.48957e-05,0.00581796  3.3971e-05,4.26404e-05  6.74445e-05,0.0109755   4.8097e-05,0.0145215    5.40342e-05,0.00379731  0.000377183,0.00853322  2.13206e-05,0.00337318  3.52239e-05,0.00407263  1.62767e-05,0.00193362  7.56385e-05,0.00950661  4.20842e-05,0.00481382  1.91134e-05,0.00253192  0.000126937,0.0090789   0.0036588,0.00587762    1.0427e-05,0.0016747    2.57371e-05,0.00350532  0.231977,0.267635   9.01372e-05,0.0116063   3.61905e-05,0.00479547  2.30435e-05,0.0032731
        3   0.122949,0.137744   1.51152e-05,0.00177613  1.65846e-07,1.00925e-05 0.000558957,0.0032513   0.015292,0.0269782  4.06388e-05,0.00157236  0.0925438,0.111993  0.0384322,0.0481046 0.00322823,0.00788021   0.0722691,0.0822557 0.0159843,0.0271257 0.0569958,0.0722337 0.000939467,0.00723918  0.0115211,0.0187761 0.00937067,0.0117235    0.135176,0.142903   0.0211448,0.0287984 0.00433432,0.00801638   0.0457743,0.0548462 0.2689,0.29917  2.91378e-05,0.00163938

* Files with suffixes of the form ``differentialpreferences_selection_???.txt`` where the ``???`` indicates the different selections. For instance, if the `Input file`_ specifies::

    selection_vaccinated_serum vaccinated_serum_codoncounts.txt
    selection_infected_serum infected_serum_codoncounts.txt

  as the alternative selections, then the following two files are created: ``differentialpreferences_selection_vaccinated_serum.txt`` and ``differentialpreferences_selection_infected_serum.txt``. These files give the differential preferences :math:`\Delta\pi^{s_i}_{r,a}` for each amino acid :math:`a` at each site :math:`r` for each selection :math:`s_i`. The files also give the root mean square (RMS) differential preference at each site for each selection. These are defined as :math:`\Delta\pi^{s_i}_{\rm{RMS},r} = \sqrt{\sum_a \left(\Delta\pi^{s_i}_{r,a}\right)^2}`. The minimal value of the RMS differential preference is zero, and the maximal value is 2.

  The first line is a title line that begins with a # character and then lists the column headers separated by tabs. The remaining lines contain the data.

  The columns contain tab-separated values as follows:

    * The first column gives the residue number in 1, 2, ... numbering. The header for this column is *SITE*.
    
    * The second column gives the wildtype amino acid at the site (such as *N* or *R*). The header for this column is *WT_AA*.
    
    * The third column lists the RMS differential preference for that site, :math:`\Delta\pi^{s_i}_{\rm{RMS},r}`. The header for this column is *RMS_dPI*.
    
    * The remaining 21 columns list the differentila preference :math:`\Delta\pi^{s_i}_{r,a}` values for each amino acid :math:`a` at this site :math:`r` for the selection pressure :math:`s_i` described in this file. The amino acids are listed in alphabetical order according to their one letter codes, and then the final column gives the equilibrium preference for a stop codon (denoted by a * character). The headers for these columns are *dPI_A*, *dPI_C*, *dPI_D*, ..., *dPI_W*, *dPI_Y*, *dPI_**. Here is the header::

        #SITE   WT_AA   RMS_dPI    dPI_A    dPI_C    dPI_D    dPI_E    dPI_F    dPI_G    dPI_H    dPI_I    dPI_K    dPI_L    dPI_M    dPI_N    dPI_P    dPI_Q    dPI_R    dPI_S    dPI_T    dPI_V    dPI_W    dPI_Y    dPI_*

* For each selection pressure, there is also a ``differentialpreferences_selection_???_credibleintervals_95.txt``. These are text files giving the median-centered 95% credible intervals on the :math:`\Delta\pi^{s_i}_{r,a}` values in the ``differentialpreferences_selection_???.txt`` files. The first column is *SITE* and the remaining columns give the upper and lower bounds on the credible intervals for each of the :math:`\Delta\pi^{s_i}_{r,a}` values as comma-separated numbers (no spaces). Here is the header::

        #SITE   dPI_A_95cred dPI_C_95cred dPI_D_95cred dPI_E_95cred dPI_F_95cred dPI_G_95cred dPI_H_95cred dPI_I_95cred dPI_K_95cred dPI_L_95cred dPI_M_95cred dPI_N_95cred dPI_P_95cred dPI_Q_95cred dPI_R_95cred dPI_S_95cred dPI_T_95cred dPI_V_95cred dPI_W_95cred dPI_Y_95cred dPI_*_95cred

* If the *preference_plots* option is being used, then plots are created showing the preferences :math:`\pi_{r,a}` for the *control_sample* and the differential preferences :math:`\Delta\pi_{r,a}^{s_i}` for each selection :math:`s_i`. These plots are in the subdirectory specified by *preference_plots*. Within this subdirectory, they have names beginning with *outfileprefix* and then having the suffix ``preferences_control_selection_1.pdf``, ``preferences_control_selection_2.pdf``, ... for residues 1, 2, ... for the *control_sample* preferences. For the differentail preferences, the names are ``differentialpreferences_selection_???_1.pdf``, ``differentialpreferences_selection_???_2.pdf``, ... where ``???`` gives the name of the selection pressure. These plots show the posterior mean and the median-centered 95% credible intervals.

* If the *MCMC_traces* option is being used, then plots are created showing the traces of the *control_sample* preferences and the differential preferences during the MCMC. These plots are in the subdirectory specified by *MCMC_traces*, and have names beginning with the prefix given by *outfileprefix* followed by the suffixes ``preferences_control_selection_1.pdf``, ``preferences_control_selection_2.pdf```, ... for the *control_sample* preferences, and followed by the suffix ``differentialpreferences_selection_???_1.pdf``, ``differentialpreferences_selection_???_2.pdf`` for the other selection pressures where ``???`` indicates the name of the selection pressure. These plots show traces of the values during the MCMC.


.. include:: weblinks.txt
