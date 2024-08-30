# Bin Structure Classification Model

NC350 chromosomes 2, 7, and 10 were used to create the classification models because they all contained Expansions. 
The top priority for these models is identifying Expansions correctly.

For the purposes of model making, non-overlapping 20kb bins were used to reduce the sample size. However, all summary measures are made to 
be proportionate to the total # of monomers in a bin, so the size of the bin should not affect accuracy of the model.

## Starting data
Each bin was then converted to a network and revelant summary information was extracted: 
Model_Making.R (up to line 121)
For three chromosomes, the structure of the networks and monomer patterns were hand checked (lines 124-126), and logged in the model-building data in the Class 
column: network_summary_or_chr2_10_7_CLASS3.csv


## Model building
Two types of classification models were then built -- a decision tree using rpart and an LDA model using MASS. The model-building data was split into testing and
training data to assess the models. Both models were saved: 
lda_model.rds
tree_model.rds



LDA Model information
```
fit2_test<- readRDS("lda_model.rds")
fit2_test
#Call:
#lda(Class ~ dup_monomers + uncon_clust + prop_clust1 + prop_clust2 + 
#    modularity + mean_dist + avg_jacc + max_contract, data = train, 
#    na.action = "na.omit")

#Prior probabilities of groups:
#  Disorder        Exp        HOR      Order 
#0.39890710 0.04371585 0.32240437 0.23497268 

#Group means:
#         dup_monomers uncon_clust prop_clust1 prop_clust2 modularity mean_dist  avg_jacc max_contract
#Disorder    0.2843868   0.8626432   0.1774810  0.11423842  0.2270600  4.622747 0.8257814   0.11236594
#Exp         0.9391121   0.1584511   0.9783296  0.01379588  0.0000000  2.000000 0.9947765   0.87617808
#HOR         0.3300667   0.3842357   0.2106107  0.13248385  0.7349731  1.530723 0.8540105   0.06220439
#Order       0.2956020   0.2970516   0.5642663  0.08882801  0.2305042  2.208686 0.8898959   0.07653974

#Coefficients of linear discriminants:
#                      LD1         LD2         LD3
#dup_monomers  -1.01763235  0.16125755 -0.10658159
#uncon_clust   -4.47445894 -3.64050454  4.17873647
#prop_clust1   -2.17936862  5.24591431  8.08530854
#prop_clust2  -15.62437926  0.82750494 -3.26076198
#modularity    -0.30184930  0.94008521 -1.21792072
#mean_dist     -0.02667681  0.08569143  0.02510456
#avg_jacc       0.28535090 -1.80825421 -1.58049920
#max_contract  20.75385791 -7.05276978 -5.81356492

Proportion of trace:
   LD1    LD2    LD3 
0.7648 0.1485 0.0867 
