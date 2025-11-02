import os
import rdata
import pandas as pd

converted_data = rdata.read_rds("/Users/eltonugbogu/Documents/dsmz/data/VST_tcga_brca.rds")
brca_metadata = pd.read_csv("/Users/eltonugbogu/classification_of_dsmz_celllines/dsmz/scripts/R/TCGA_BRCA_project_subtype.csv")
print(dim(converted_data))
print(brca_metadata.head())

# subset converted_data
brca_subtype_expression_matrix = converted_data[brca_metadata['sample_barcode']] # x_data
print(dim(brca_subtype_expression_matrix))
print(brca_subtype_expression_matrix.head())

# save expression data without subtype labels
brca_wo_subtype_exp_matrix = converted_data[brca_subtype_expression_matrix.index != brca_metadata['sample_barcode']]
print(dim(brca_wo_subtype_exp_matrix))

# extract subtype labels in the same order as subtype_expression_data
brca_subtype_labels = brca_metadata['PAM50_Subtype'][brca_metadata['sample_barcode']]  # y_labels
print(brca_subtype_labels.head())

# add the labels as a new column to subtype_expression_data 
brca_subtype_expression_matrix['PAM50_Subtype'] = brca_subtype_labels
print(brca_subtype_expression_matrix.head())
brca_subtype_expression_matrix.to_csv("/Users/eltonugbogu/Documents/dsmz/data/VST_tcga_brca_subtype_expression_matrix.csv", index=False)

# boxplot of each subtype
for subtype in brca_subtype_labels.unique():
    plt.boxplot(brca_subtype_expression_matrix[brca_subtype_labels == subtype])
    plt.title(subtype)
    plt.show()

# scatter plot of each subtype
for subtype in brca_subtype_labels.unique():
    plt.scatter(brca_subtype_expression_matrix.index, brca_subtype_expression_matrix[brca_subtype_labels == subtype])
    plt.title(subtype)
    plt.show()  

# outliers detection per subtype using PCA
for subtype in brca_subtype_labels.unique():
    pca = PCA(brca_subtype_expression_matrix[brca_subtype_labels == subtype])
    plt.scatter(pca.explained_variance_ratio_)
    plt.title(subtype)
    plt.show()  
