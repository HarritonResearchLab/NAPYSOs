# NOTES

## TO-DO
- [x] Produce final compiled results file (with changes from LAH and Miles)  
- [ ] Make latex version of results for table 2 in paper 
- [ ] add final results mrt to overleaf
- [ ] update the period / timescales plots in figure 4 with updated periods / classes now
- [ ] update table 4 with new dist %s calculated on january final results
- [ ] update figure 6 with new data
- [ ] update figure 12??
- [ ] figure 13
    * remake all the appendix plots with new data
    * replace the old plots in fig 13 with new ones
    * add zip file of all the plots into overleaf 
 
- [ ] Email LAH 
    * should we still have table 3? 
    * shoudl we add information about changing periods based on how they accounted for aliases // beat periods into the paper?  

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
    