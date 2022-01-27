## Cryptic cleavage sites in membrane proteins

![Workflow](ppt/spc_diag.png)

### Part A
1. Consider the entire human proteome's first 70aa and run [SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-4.1) with **TM** and **noTM** networks
2. Compare the **Y-score** in both the cases and highlight the proteins (scatter plot) that show a cleavage site in the **noTM** mode but not in the **TM** mode
3. The hunch is that these **cyptic cleavage sites** are present in proteins but otherwise inaccessible
4. However, upon a mutation at another site in the same protein, these **cyptic cleavage sites** become accessible (perhaphs due to the protein misfold caused by the mutation)
5. The revealation of these sites can lead to its post-translational cleavage by SPC (this can be a way of quality control)
7. For example, Andrea showed that a mutation C201R in Cx32 can lead to revealation of a **cyptic cleavage site** in first 70aa, which was otherwise inaccessible in the WT

### Part B
1. Same as Part A but rather than taking the first 70 aa, consider all Type II (N-term: Cytoplasmic) TMDs such that 5aa before the TMD + TMD itself + x aa is equal to 70

### Part C
1. Same as Part A but rather than taking the first 70 aa, consider all Type II (N-term: Cytoplasmic) TMDs such that the TMD itself + x aa is equal to 70

### webApp
The web application can be accessed [here](http://shiny.russelllab.org/spc/webApp/)

### To-do (post meeting on Oct 5, 2021) -- check Issues
1. Look for correlation b/w loop length vs. num of mutations
2. Look for known cleavage sites in [Topfind4.1](https://topfind.clip.msl.ubc.ca/)
3. Add more information in the tables (gene and position-level information) for PART B
4. Color the True Positive cases in the figure (should mostly be in the top-right quadrant)
5. Web-site development
6. Workflow figure
7. Manuscript writing
