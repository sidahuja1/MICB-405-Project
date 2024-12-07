# Here we will load the .pombase file from PombBase's website and convert it to a tsv. It is inititially 
# There is a lot of excess information here so we will remove the other columns other than the gene ID's and GO term accessions

# Locate the file downloaded from PomBases's Website
spombe_topGOIDs = 'C:/Users/brenn/Downloads/gene_association.pombase'

# Read and preview the contents
with open(spombe_topGOIDs, 'r') as file:
    data_preview = file.readlines()

# Check the first 10 rows to make sure it's what we are expecting
data_preview[:10]

# Here we are converting the GAF file to a topGO file
from collections import defaultdict

# Create a dictionary to preserve the gene to GO term mapping relationship
gene_to_go = defaultdict(list)

# Process each line in the file
for line in data_preview:
    # Make sure we skip metadata lines starting with '!'
    if line.startswith('!'):
        continue
    
    # Split the line into columns
    columns = line.strip().split('\t')
    
    # Extract the gene ID (column 2) and GO term (column 5)
    gene_id = columns[1]
    go_term = columns[4]
    
    # Append the GO term to the gene's list
    gene_to_go[gene_id].append(go_term)

# Create a list of lines for the topGO-compatible TSV format
topgo_lines = ["GeneID\tGO Terms"]
for gene, go_terms in gene_to_go.items():
    # Join GO terms with commas
    go_terms_str = ",".join(set(go_terms))  # Use set to remove duplicates
    topgo_lines.append(f"{gene}\t{go_terms_str}")

# Save the processed data to the MICB405 project folder
output_path = 'C:/Users/brenn/MICB405/finalproject/gene_association_topGO.tsv'
with open(output_path, 'w') as output_file:
    output_file.write("\n".join(topgo_lines))

output_path
