import os

import pandas as pd
from sqlalchemy import create_engine

# Used for sqlalchemy engine services
sql_engine = create_engine("mysql+pymysql://root:Liuju83198431@localhost:3306/hemu_database")

# Final database sheet
final_TE_datasheet_df = pd.DataFrame()
# Counter for samples
samples_completed = 1

# Constructing dict, sample_id -> sample_tissue
sampleinfo_data = pd.read_csv("./sampleinfo.csv", low_memory=False, sep=",")
sampleinfo_df = pd.DataFrame(sampleinfo_data)
sample_tissue_dict = {}  # sample_id -> sample_tissue
for sample_entry in sampleinfo_df.index:
    sample_tissue_dict[sampleinfo_df.loc[sample_entry].sample_id] = str(sampleinfo_df.loc[sample_entry].sample_tissue).lower()
print("Dict construction completed. Length:", len(sample_tissue_dict.keys()), sep=" ")

for indv_file in os.listdir("./Normalized_TE"):  # ex. SRR0001_TE.csv

    if indv_file.endswith("_TE.csv"):
        TE_data = pd.read_csv("./Normalized_TE/" + indv_file, low_memory=False, sep=",")
        TE_df = pd.DataFrame(TE_data)

        # Delete single-copy TEs
        TE_df.drop(TE_df[TE_df['TE_id'].str.contains(pat='Seq', regex=False)].index, inplace=True)
        # Delete columns (axis=1)
        TE_df_final = TE_df.drop(['fam_length'], axis=1)
        # Insert sample_id at the second column
        TE_df_final.insert(loc=1, column='sample_id', value=indv_file.rstrip("_TE.csv"))

        # Insert tissue information at the last column
        try:
            tissue_type = str(sample_tissue_dict[indv_file.rstrip("_TE.csv")])
        except KeyError:
            print("Index not found for %s, considering tissue unknown." % indv_file.rstrip("_TE.csv"))
            tissue_type = "unknown"
        TE_df_final.insert(loc=6, column='tissue_type', value=tissue_type)

        final_TE_datasheet_df = pd.concat([final_TE_datasheet_df, TE_df_final])
        print("Completed %d samples" % samples_completed)
        samples_completed += 1
        #break  # Diag purposes

print("\n\nCreating database sheet..")
final_TE_datasheet_df.to_sql(name='zea_te_test', con=sql_engine, if_exists="replace")
print("All done")