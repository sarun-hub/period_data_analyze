import pandas as pd
import os
from utils import getCompoundName 

files = os.listdir("data")
os.makedirs("data_0_00",exist_ok=True)

for file in files:
# Load the Excel file
    compound_name = getCompoundName(file) 
    excel_file = os.path.join("data",file)

    # If the Excel file has multiple sheets, specify the sheet name or index
    df = pd.read_excel(excel_file, sheet_name=0)

    # Save as CSV
    df.to_csv(f"data_0_00/{compound_name}.csv", index=False)  # index=False avoids writing row numbers