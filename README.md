# KBS_LTER_Crop_microbiome
Shared Space for Those working on the KBS LTER Microbiome project.

In this project we work on:

  -Microbiomes from rotational crop soil, roots, stems and leaves. This is the origin
  
  -Micriobiomes from conventional, no-till, and organic management systems. This is called manangement systems.
  
  -Microbiomes from different growth stages (soybean V2, R3, R6). This is growth stage
  
  -If you see a label in a mappingfile labelled "bar_label" that is a mix of the management system and growth stage in some form.
  

I would love to ask if anyone associated with this repository/ project would organize their uploads.
 
 We can follow the following methods for naming:
  
  1) R files (save where the name applies): <crop>_<year>_<origin>_<your initials>_<DAY.MONTH.YEAR.AM OR PM>.R
    
    2) Try and group files together when possible
      ex) 2022/Stems/Fungi/Roots/amplicons/demultiplex
      ex) 2022/Stems/Fungi/mapping_files
   
   3) ANNOTATE YOUR CODE WELL WITH WHAT YOU DO, WHAT PACKAGE IT IS UTILIZING, AND WHY THAT PACKAGE. IF YOU HAVE ISSUES, ANNOTATE WITH AT LEAST 3 ###!
   
   4) If you are saving figures you created, Please save the figure like: 
   
      <what kind of figure (ie. stacked_barplot)>_<origin>_<year>_<FIRST LINE OF CODE CORRESPONDING WITH MAKING GRAPH_LAST LINE OF CODE CORRESPONDING           WITH MAKING GRAPH>
    
    5) Pipeline folders: <step in pipeline (ie 1)>_<what it do>_program

      ex) Pipelines/clustered_otus_pipeline/1_demultiplexing_QIIME_myco_class_2022.sb
    
    6) Try and save one file with all the codes in one place if we can.
    
      ex) 2022/stems/R_codes -or- 2022/R_codes
   
   7) I have found for R, that having all the codes we need in one place is super helpful. I made a file called "Working_data" and it contained OTU              tables, mapping file, taxonomy file, otu_fasta and other files you need to work.
   
      ex) 2022/stems/working_codes -OR we can do- 2022/working_codes
     
     
As time goes on, we can find what works for us to share things. Between us collaborators, the MSU HPCC (a great way to combine files and download to upload here), and pushing things from Rstudio to here without having to download, we should be good!!!
