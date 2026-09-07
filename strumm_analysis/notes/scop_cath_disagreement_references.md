**Yes, this is a very well-known and pervasive problem in structural bioinformatics.**

The phenomenon you are observing—where a PDB structure contains clear domains that align structurally to SCOP or CATH folds, yet lacks an entry in their official release tables or SIFTS mapping—stems from a combination of **curation lag**, **disagreement on domain boundaries**, and **unclassified regions (the "dark matter" of protein structures)**.

---

### Key Reasons for Unannotated Domains

* **Curation Lag and Release Delays:** Historical SCOP relied heavily on manual curation, leading to multi-year lag times between new PDB depositions and database updates. While CATH and SCOP/SCOPe use semi-automated pipelines today, manual checks are still required for novel structural variants or complex multi-domain proteins.
* **Domain Boundary Disagreements:** Automatic domain segmentation tools (such as PUU, PDP, or CATHEDRAL) frequently disagree on domain boundaries, especially for continuous vs. discontinuous domains. In many cases, one database will fail to split a multidomain chain, leaving whole domain regions unassigned or merged into a larger fragment (Csaba et al., 2009; Redfern et al., 2007).
* **Incomplete Domains and Expression Constructs:** Many PDB entries consist of truncated constructs, flexible linkers, or partial domain fragments designed for crystallization. Automated filters often mark these as low-confidence or non-standard, excluding them from automated classification pipelines.
* **Continuous vs. Discrete Fold Space:** Some protein structures fall into ambiguous, continuous regions of structural space (e.g., repeating $\alpha/\beta$ motifs) that do not fit neatly into discrete structural families (Redfern et al., 2007).

---

### Key Literature & References

* **Differences and Inconsistencies in Classification & Boundaries:**
Csaba et al. (2009) systematically analyzed domain assignments between SCOP and CATH, showing that large portions of the PDB suffer from inconsistent boundary assignments, resulting in large sets of domains present in one hierarchy but missing or split differently in the other.
* **Automated Boundary Prediction Limitations:**
Redfern et al. (2007) highlighted the challenge of domain assignment in multidomain proteins, showing that traditional automated algorithms frequently produce inconsistent boundary assignments when identifying recurrent folds in newly deposited structures.

---

**References**

Csaba, G., Birzele, F., & Zimmer, R. (2009). Systematic comparison of SCOP and CATH: a new gold standard for protein structure analysis. *BMC Structural Biology*, *9*, Article 23. [https://doi.org/10.1186/1472-6807-9-23](https://doi.org/10.1186/1472-6807-9-23)
Cited by: 122

Redfern, O. C., Harrison, A., Dallman, T., Pearl, F. M. G., & Orengo, C. A. (2007). CATHEDRAL: A fast and effective algorithm to predict folds and domain boundaries from multidomain protein structures. *PLoS Computational Biology*, *3*(11), e232. [https://doi.org/10.1371/journal.pcbi.0030232](https://doi.org/10.1371/journal.pcbi.0030232)
Cited by: 130