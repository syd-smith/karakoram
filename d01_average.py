# %%
"""
Author: Sydney Smith
Created: June 3, 2026
"""

from pathlib import Path
import sys


# ==================================
# - Establish Relative File Path - 
# ==================================

try:
    current_file_directory = Path(__file__).resolve().parent
except NameError:
    current_file_directory = Path().resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))

import from_savanna.nclcmaps as cmap
from preprocessor.data_access import get_data


noise = get_data('RAINNC', 7, 'd01', 'noise')
exp = get_data('RAINNC', 7, 'd01', 'MODISImproved')

# %%
for i in range(0, 5):
    domain_anom = exp[i].mean() - noise[i].mean()
    print(domain_anom)

# %%
