import pandas as pd
import numpy as np

# Read the file containing the CSS values calculated into a pandas dataframe
data = pd.read_csv('CSS.txt', sep='\t', quotechar='"')

# Column names are 'chromosome', 'position' and 'css'
# Sort by 'chromosome' and 'position' just to be sure
data = data.sort_values(['chromosome', 'position'])

# Get unique chromosomes
chromosomes = data['chromosome'].unique()

# Process each chromosome
for chromo in chromosomes:
    chromo_data = data[data['chromosome'] == chromo]

    # Start from the minimum position and end at the maximum
    min_pos = chromo_data['position'].min()
    max_pos = chromo_data['position'].max()

    # Start windowing
    start = min_pos
    output = pd.DataFrame()  # Initialize an empty dataframe for each chromosome
    while start <= max_pos:
        end = start + 20000 # 20 kb window
        
        # Get subset of data within window
        subset = chromo_data[(chromo_data['position'] >= start) & (chromo_data['position'] < end)]

        if not subset.empty:
            # Compute mean
            mean_css = subset['css'].mean()

            # Append to the output dataframe
            output = output.append({'chromosome': chromo, 'start': start, 'end': end, 'mean_css': mean_css}, ignore_index=True)
        
        # Move the window.
        start = end

    # Save the output for each chromosome to a separate file
    output.to_csv(f'output_chromosome_{chromo}.csv', index=False)
