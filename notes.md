# The main idea

Given a VCF file, predict the population of individuals based on how their genome compares to other genomes.


## Input

The mandatory input is the **VCF file** (which can optionally be gzipped).

We need to allele frequency for the populations that we are trying to analyse. This can either be based on the `INFO` column of the VCF file (for superpopulations, we get the properties `AFR_AF`, `AMR_AF` and so on), or calculated from scratch. For this calculation, we need to know the population of each individual; this must be given with a **population** file.

Optionally, we may provide a file for the **individual identifiers** that we are interested in (the ones to classify); if not provided, all individuals are classified.

Optionally, we may provide a file for the **polymorphisms** of interest (the ones to base the calculations on); if not provided, all polymorphisms are considered.


## The procedure

For each relevant line in the VCF file (i.e, for each SNP):
1. Compute allele frequency for all populations. The reference allele is considered as `0` and all alternative alleles as `1`. Let these values be stored in a map $F$, which associates *populations* with their *allelic frequency* (we call $F_i$ the allele frequency for the population $i$).
    * *This may be a point for future development, since we are at the moment conflating all non-reference alleles (i.e. all alternative alleles) as equivalent...*
2. Convert individual genotypes into a numeric value $y$:
    * `0|0` becomes `0.0`;
    * `0|1` and `1|0` become `0.5`;
    * `1|1` becomes `1.0`.
3. For each individual, *do the math*. There are several algorithms that can be implemented here (see the section below).


## The *math* algorithms

There are several ways to compute a classification for each individual. They must all be based on a monoid structure (e.g. a list of population labels) assigned to each individual, which is appended (possibly multiple times) for each polymorphism; after the last polymorphism, this structure is used to compute a final classification label for the individual it is assigned to.

Here are a few implemented algorithms, each one used to compute the values to append to the monoid structure of an individual


### The "Mariana" way

1. The monoid structure is a list of population labels
2. For all distinct pairs of population $(p_1, p_2)$,:
    1. Calculate $D = |y-F_{p_1}| - |y-F_{p_2}|$
    2. If $D<0$, append $p_1$ to the list of labels; otherwise append $p_2$.
3. In the end, count the number of times each population $i$ appears in the list (call this value $c_i$) and assign to this individual the population with highest $c_i$.

There are two **extensions** to this method:
1. In step 2.2, and given a fixed threshold $\theta_D$, only append a label to the list if $|D| \ge \theta_D$.
2. In step 3, and given a fixed threshold $\theta_c$, only assign a classification to the individual if the difference between the maximum of all the $c_i$ values is at least $\theta_c$ units greater that the second maximum value.
    * Note that $\theta_c$ can be given in absolute values or relative to the highest $c_i$. The relative approach ensures that the same threshold can be applied to different amounts of polymorphisms


## The "João" way

1. The monoid here is also a list of population labels
2. For each population $i$, calculate $d_i = |y - F_i|$.
    1. Append to the list of labels the population that minimizes $d_i$.
3. As in the previous method, count the number of times each population $i$ appears in the list ($c_i$) and assign to this individual the population with highest $c_i$.

There are two **extensions** to this method, equivalent to the extensions in the previous one:
1. Given a fixed threshold $\theta_d$, step 2.1 only appends the label to the list if the difference between the minimum $d_i$ and the second lowest one is at least $\theta_d$.
2. In step 3, and given a fixed threshold $\theta_c$, only assign a classification to the individual if the difference between the maximum of all the $c_i$ values is at least $\theta_c$ units greater that the second maximum value.
    * Note that $\theta_c$ can be given in absolute values or relative to the highest $c_i$. The relative approach ensures that the same threshold can be applied to different amounts of polymorphisms
