### Scripts used in manual transposable element curation of NS clone Myzus persicae assembly

This pipeline was used for automated part of the manual transposable element curation to:

1. Cluster TEs from the RepeatModeler library to reduce redundancy
2. Search the genome for each sequence and select the best hits to use in manual curation
3. For each consensus sequence:
    1. Extend the coordinates of all hits to it
    2. Local align all extended hits to each other
    3. Trim off extended regions not seen in other hits
    4. Create multiple sequence alignment for manual curation
