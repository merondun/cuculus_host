# Evolutionary analysis of NDUFAF4 and W chromosome paralog

Assay these species and sexes for the presence of NDUFAF4, and then create a phylogeny based on the copies:

| Cuculus  canorus        | Female | Illumina | SRR11394165                          |
| ----------------------- | ------ | -------- | ------------------------------------ |
| Cuculus canorus         | Female | Illumina | SRR11531702                          |
| Cuculus canorus         | Female | Illumina | SRR11531718                          |
| Cuculus canorus         | Male   | Illumina | SRR11531726                          |
| Cuculus canorus         | Male   | Illumina | SRR11531741                          |
| Cuculus micropterus     | Male   | Illumina | SRR14117632                          |
| Cuculus micropterus     | Female | Illumina | SRR14117631                          |
| Cuculus poliocephalus   | Male   | Illumina | SRR14117570                          |
| Cuculus poliocephalus   | Female | Illumina | SRR14117568                          |
| Cuculus canorus         | Female | Illumina | SRR99999129                          |
| Clamator glandarius     | Female | HiFi     | SRR26807982                          |
| Coccyzus lansbergi      | Female | HiFi     | SRR26902471                          |
| Cuculus canorus         | Female | HiFi     | All bCucCan1.pri Reads; PRJNA1008121 |
| Dromococcyx pavoninus   | Female | HiFi     | SRR26905917                          |
| Geococcyx californianus | Female | Illumina | SRR9994302                           |
| Piaya cayana            | Female | Illumina | SRR9947006                           |

1... : Subsampling reads to 10GB and aligning them to the whole genome cuckoo reference
2: Extract consensus sequence over the chr3 and W chromosome copy
3: Align the nucleotide sequences, filter for gaps, and create a tree
4: Plot
