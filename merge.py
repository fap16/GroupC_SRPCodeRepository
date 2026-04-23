import pandas as pd
import glob
import os

# use glob to find all 730 Cufflinks output files across subdirectories.
# ensures the script is scalable for any number of samples.
# 1. Define paths
path = "/scratch/alice/j/jvk3/camp_alice/results/*_cufflinks_V2/genes.fpkm_tracking"
output = "/scratch/alice/j/jvk3/camp_alice/camp_organoid_matrix_V2.csv"

# 2. Find all 730 files
files = glob.glob(path)
total_files = len(files)
print(f"Found {total_files} valid samples. Starting merge...")

# 3. Process files one by one
df_list = []
count = 0

for f in files:
    # extract the SRR ID from the folder name
    srr_id = f.split('/')[-2].replace('_cufflinks_V2', '')
    
    # read only gene_id and FPKM columns
    temp_df = pd.read_csv(f, sep='\t', usecols=['gene_id', 'FPKM'])
    
    # remove duplicate gene IDs 
    # to get past 'InvalidIndexError' 
    temp_df = temp_df.drop_duplicates(subset='gene_id')
    
    # rename FPKM column to the Sample ID for identification in the master matrix
    temp_df.columns = ['gene_id', srr_id]
    temp_df.set_index('gene_id', inplace=True)
    
    df_list.append(temp_df)
    
    # progress tracker
    count += 1
    if count % 50 == 0:
        print(f"Processed {count}/{total_files} samples...")

# 4. Stitch everything together
#'sort=False' maintains original gene ordering to reduce computation time.
print("Stitching all cells into one master matrix. This may take a minute...")
final_matrix = pd.concat(df_list, axis=1, sort=False)

# 5. Fill any missing values (NAs) with 0
# if gene found in one cell but not another, we assume 0 expression
final_matrix = final_matrix.fillna(0)

# 6. Save the final matrix
final_matrix.to_csv(output)
print(f"SUCCESS! Matrix saved to: {output}")
