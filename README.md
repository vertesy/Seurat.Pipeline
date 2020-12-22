# Seurat Pipeline
Seurat Pipeline is an example single-cell analysis pipeline using Seurat and a set of custom packages developed for the [Knoblich Lab](https://www.imba.oeaw.ac.at/research/juergen-knoblich/team/).



### Dependencies

#### Published Libraries Used

- Seurat
- tidyverse
- cowplot
- future
- doMC
- tictoc


#### Custom Function Libraries Used
- [Seurat.utils](https://github.com/vertesy/Seurat.utils/)
- [Seurat.multicore](https://github.com/vertesy/Seurat.multicore) (optional)
- [CodeAndRoll](https://github.com/vertesy/CodeAndRoll)
- [ggExpressDev](https://github.com/vertesy/ggExpressDev)
- [MarkdownReports](https://github.com/vertesy/MarkdownReportsDev)


### Structure of the pipeline
There is one **wrapper Frame script** (“00.Frame.ProjectCode.etc.R”) in each projects. It loads **Daughter scripts**:
   1. Custom scripts (All in the GitHub repo of the project, e.g. [**SEO**](https://github.com/vertesy/SEO/))
      2. [**Parameters**](https://github.com/vertesy/SEO/blob/master/SEO2/Parameters.SEO2.R) script
      2. [**Genelists**](https://github.com/vertesy/SEO/blob/master/SEO2/Gene.Lists.SEO.R) scriptCustom **Elements** script
   
   2. Reusable scripts ([**Seurat.pipeline**](https://github.com/vertesy/Seurat.Pipeline/) repo)
   
      1. **Elements** scripts that are used in many projects
   
   3. Custom Package
   
      1. [MarkdownReports](https://github.com/vertesy/MarkdownReportsDev) 		→ See web for installation from Github.
   
   4. Custom Function Libraries
   
      - [Seurat.utils](https://github.com/vertesy/Seurat.utils/)
      - [Seurat.multicore](https://github.com/vertesy/Seurat.multicore) (optional)
      - [CodeAndRoll](https://github.com/vertesy/CodeAndRoll)
      - [ggExpressDev](https://github.com/vertesy/ggExpressDev)
      
        