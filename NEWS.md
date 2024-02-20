### Package changes from previous lphom version 0.3.1-1

Compared to version 0.3.1-1 of lphom, the new version:

* includes two new arguments (`apriori` and `lambda`) in functions `lphom()`, `tslphom()`, `nslphom()` and `lclphom()`. The `apriori` argument allows to introduce in the inference process information from polls or some other a priori knowledge/intuition about the row standardized transfer probabilities. The `lambda` argument controls the weight the user assigns to the a priori information.
* offers an implementation of the bottom-up approach for ecological inference through the `rslphom()` function.
* expands the number of scenarios in terms of how entries of new voters and exits of previous voters are handled by the different functions. The new options include: `"ordinary"`, `"enriched"`, `"semifull"`, `"fullreverse"`, `"adjust1"` and `"adjust2"`.
* modifies the `lclphom()` algorithm by including the two possibilities of calculating the distance with the global matrix. The distance with the global matrix that is used to generate the local solution and the distance with the global matrix where the local solution integrates into. This is controlled with the `type.errors` argument. 
* incorporates to the outputs of the functions `lphom()`, `tslphom()`, `nslphom()`, `lclphom()`, `rslphom()`, `lphom_joint()`,
`tslphom_joint()` and `nslphom_joint()` the list `deterministic.bounds`. The object `deterministic.bounds` contains both the deterministic global and unit bounds that derive from the observed margins.
* reviews documentation and updates references.