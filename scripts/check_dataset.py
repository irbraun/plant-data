import sys
import pandas as pd
import numpy as np
import glob
import os
import warnings
warnings.simplefilter('ignore')
sys.path.append("../../oats")
from oats.biology.dataset import Dataset




path = "../genes_texts_annots.csv"




# Recreating the dataset object from the saved csv file.
# Double checking for things that we know should be true about the descriptions if preprocessing works as intended.
df = Dataset(path).to_pandas()
descriptions = df["descriptions"].values
assert len([s for s in descriptions if "|" in s]) == 0
assert len([s for s in descriptions if "  " in s]) == 0







# Check again this time reading in the dataframe directly from the csv file.
# Double checking for things that we know should be true about the descriptions if preprocessing works as intended.
df = pd.read_csv(path)
descriptions = df["descriptions"].values
assert len([s for s in descriptions if "|" in s]) == 0
assert len([s for s in descriptions if "  " in s]) == 0




print("completed check of genes_texts_annots file")


