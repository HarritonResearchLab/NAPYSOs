# NOTES

## TO-DO

### order for 1/23/2022
* figure 12 period recovery (and change text), fix figs 7 & 8, do supplimentary plots, choose examples, upload tar and examples 

- [x] go back and fix a lot of things...some of the new ids are formatted weird

- [x] fix references ?!? big O denotes to do
- [x] check all text references and all figures for things to change
    * text has been fixed, figs 5,7,8 need to be fixed 
- [x] change tabular information about them 

- [x] Produce final compiled results file (with changes from LAH and Miles)   
- [x] Make latex version of results for table 2 in paper  
- [x] add final results csv to overleaf    
- [x] update the period / timescales plots in figure 4 with updated periods / classes now    
- [x] update table 4 with new dist %s calculated on january final results    
- [x] update figure 6 with new data     
- [x] update figure 12??  
    * old file is located at `.\recovered\2.0\efforts\august21\aug23-29\froebrich\froebrich_periods.py`
- [x] check that figures with references qs haven't had them changed    
    * if they have fix them :)    
- [x] update infrared plots with new qm data 
- [x] O figure 13   
    * remake all the appendix plots with new data      
    * replace the old plots in fig 13 with new ones ... or choose new ones    
    * add zip file of all the plots into overleaf
    * file path to old script: `.\recovered\2.0\efforts\xi\Misc\supplimentary_figures\make_sup_plots.py`     
    * something got messed up with my code when I tried to do this on saturday? 
 
- [x] Email LAH   
    * should we still have table 3?  
    * should we add information about changing periods based on how they accounted for aliases // beat periods into the  paper?
    * ask about what we needed to change for appendix kepler plots


## OTHER

* **Schemes for QM Plots**:  
    * Colors and markers should be used as follows:  
        ```python  
        color_marker_dict={'p':['c0', 's'],  
                        'qps':['C1', 'o'], 
                        's':['C2','*'], 
                        'apd':['C3','D'], 
                        'qpd':['C4','h'], 
                        'b':['C5','d'], 
                        'l':['black','p'], 
                        'mp':['C6','P']} 
        ```

        where seaborn has been preivously importated as ```sns``` 

    * Edgecolors should be black with linewidths set to ```0.5```
