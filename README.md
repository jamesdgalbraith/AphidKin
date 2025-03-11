### Scripts used in manual transposable element curation of NS clone Myzus persicae assembly

This pipeline was used for automated part of the manual transposable element curation to

1) Cluster TEs from the RepeatModeler library to reduce redundancy
2) Search the genome for each sequence and select the best hits to use in manual curation
3) For each consensus sequence
    a) Extend the coordinates of all hits to it
    b) Local align all extended hits to each other
    c) Trim off extended regions not seen in other hits
    d) Create multiple sequence alignment for manual curation
