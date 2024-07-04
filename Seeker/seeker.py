# Desired sample names should be in column "ID", with indices in "P7" and "P5"
# NOTE the "skiprows" argument of pd.read_excel() may need to be changed!

import pandas as pd
from numpy import nanmean
import sys
import re


input_file = sys.argv[1]
output_file = input_file.replace("xlsx", "csv")
fragment_file = input_file.replace(".xlsx", "_frag_size_summary.txt")

excel_file = pd.ExcelFile(input_file)

cdna_sheets = []
for s in excel_file.sheet_names:
  if re.search("cdna", s, flags = re.I):
    cdna_sheets.append(s)

cdna_tables = pd.read_excel(excel_file, sheet_name = cdna_sheets, skiprows = 3)


# Indices
indices_1 = {}
for key, df in cdna_tables.items(): 
  indices_1[key] = df[["ID", "P7", "P5"]]

indices_2 = {}
for key, df in indices_1.items():
  temp = df.apply(lambda row: row.P7 + "-" + row.P5, axis = 1)
  indices_2[key] = df.assign(i7_i5 = temp.values)[["ID", "i7_i5"]]

combined_index_tables = pd.concat(indices_2.values(), keys = indices_2.keys())

tidy_index_tables = combined_index_tables.rename_axis(["sheet", "row"]).reset_index().drop(["row"], axis = 1)


# Calculate fragment min, max, and mean
bp_tables = {}
for key, df in cdna_tables.items():
  bp_tables[key] = df["Aver. Bp"]
 
bp_values = pd.concat(bp_tables.values()).to_list()

summary_bp = pd.DataFrame(data = {
  "measure": ["min", "max", "mean"],
  "value": [min(bp_values), max(bp_values), nanmean(bp_values)]
})
print("\nFragment size stats, also saved to '" + fragment_file + "':")
print(summary_bp)


# Write outputs
tidy_index_tables.to_csv(output_file, index = False)
summary_bp.to_csv(fragment_file, sep = "\t", index = False)

