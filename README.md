# Primer explorer2

A small library to help seleting primers for K-Seq users .

Authors: P.Ziarsolo, V. Garcia-Carpintero, J.Blanca

This small library provides various scripts to help the user to select the best kmers for your experiment.


## Installation
    
    a) pip install primer_explorer2

    b) python setup.py install


### primer_xp.py
    
    It tries to select the best 10 primers for a K-Seq experiment and predicts the outcome of all 10 primer combinations.
    
    The outcome consist in two file:
        * An excel file with some stats for each primer combinations.
        * A database with the predicted pcr products for each pair combination. This can be use as input for other provided scripts.
               
    You can use the script to look for primers in different ways. 

    You can choose the length of the primer. The bigger is the length of the primer, less sites in the genome will find, thus generating less pcr products
    You can define regions in witch you prefer to have less primers. P.ej: heterocromating, repetitive...

### generate_prediction_beds.py        

    It generates a bed file for each combination of primer pair in the pcr_products pickle database. 

### draw_pcr_product_size_distribution.py

    It draws the insert size distribution of the predicted pcr products.
    
### count_gff_features.py

    Given the prediction in bed files. It counts how many predicted pcr_products are in the gff features. By default in only counts pcr_produts in exons and genes. It can be changed.
