#!/bin/bash
# INSTRUCTIONS + POINTERS
#sudo apt-get install fim # if fim is not already installed
#chmod u+x PATH_TO_THIS FILE 
# Execute that ^^^ prior to running this for the first time
# Run this file by pasting it's path into terminal session
# (when you're in the directory with the plots, that is)
# Note that this code was specially written for multiple class
# capabilities (QM classifications)

# ACTUAL CODE
echo Enter report file path: 
read report_file
touch $report_file
exec > $report_file
echo 'ID,first_class,second_class'
for i in $(ls); do
    fim $i
    read varname
    IFS=':' read -r -a array <<< "$i"
    filename="${array[2]}"
    ID=${filename/.png/""}
    read second # second class
    echo $ID,$varname,$second
    done

# ^^^ Adds column labels so you can read file into pandas df later