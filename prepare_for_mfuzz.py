
# Load the provided pdui_result.txt file
import pandas as pd
pdui_df = pd.read_csv('~/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/01-pdui_and_motif_res/pdui_result.txt', sep="\t")

# Extract the relevant columns (PDUI values across time points) and set the Gene column as the index
mfuzz_matrix = pdui_df.set_index('Gene').iloc[:, 3:]

# Updated list of time points based on the provided column names
updated_time_points = ["proerythroblast", "early_basophilic", "late_basophilic", "polychromatic", "orthochromatic", "BFU", "CFU"]

# Calculate median for each time point across its samples
for time_point in updated_time_points:
    relevant_cols = [col for col in mfuzz_matrix.columns if time_point in col]
    mfuzz_matrix[time_point] = mfuzz_matrix[relevant_cols].median(axis=1, skipna=True)

# Keep only the median columns for each provided time point
updated_median_matrix = mfuzz_matrix[updated_time_points]

# Identify genes where the absolute difference between the median PDUI values of any two time points exceeds 0.15
updated_genes_to_keep = []

for index, row in updated_median_matrix.iterrows():
    keep_gene = False
    for i in range(len(updated_time_points)):
        for j in range(i+1, len(updated_time_points)):
            if abs(row[updated_time_points[i]] - row[updated_time_points[j]]) >= 0.15:
                keep_gene = True
                break
        if keep_gene:
            break
            
    if keep_gene:
        updated_genes_to_keep.append(index)

# Extract rows corresponding to these genes from the median matrix
updated_filtered_median_df = updated_median_matrix.loc[updated_genes_to_keep]

# Order of time points based on erythroid development
ordered_time_points = ["BFU", "CFU", "proerythroblast", "early_basophilic", "late_basophilic", "polychromatic", "orthochromatic"]

# Sort the dataframe based on the provided order
sorted_filtered_median_df = updated_filtered_median_df[ordered_time_points]

# Save the sorted filtered median matrix to a file
output_file_sorted = "02-diff_0.15_for_mfuzz_sorted.txt"
sorted_filtered_median_df.to_csv(output_file_sorted, sep="\t")
