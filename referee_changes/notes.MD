# NOTES

## TO-DO
- [x] Produce final compiled results file (with changes from LAH and Miles)   
- [x] Make latex version of results for table 2 in paper  
- [x] add final results mrt to overleaf    
- [x] update the period / timescales plots in figure 4 with updated periods / classes now    
- [x] update table 4 with new dist %s calculated on january final results    
- [x] update figure 6 with new data     
- [ ] update figure 12??     
- [x] check that figures with references qs haven't had them changed    
    * if they have fix them :)    
- [x] update infrared plots with new qm data 
- [ ] figure 13   
    * remake all the appendix plots with new data      
    * replace the old plots in fig 13 with new ones    
    * add zip file of all the plots into overleaf     
 
- [ ] Email LAH   
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
