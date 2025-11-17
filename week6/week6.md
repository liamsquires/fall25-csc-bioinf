# Get Sample Data


```python
import os
os.makedirs("tempfiles", exist_ok=True)
os.makedirs("output", exist_ok=True)

# Set ALEVIN_FRY_HOME environment variable for simpleaf (use absolute path)
alevin_fry_home = os.path.abspath('.')
%env ALEVIN_FRY_HOME=$alevin_fry_home
```

    env: ALEVIN_FRY_HOME=/home/liam/Bioinformatics/fall25-csc-bioinf/week6



```python
!curl -L "https://github.com/f0t1h/3M-february-2018/raw/refs/heads/master/3M-february-2018.txt.gz" -o "tempfiles/3M-february-2018.txt.gz"
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0

      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0


     71 17.5M   71 12.4M    0     0  10.0M      0  0:00:01  0:00:01 --:--:-- 10.0M

    100 17.5M  100 17.5M    0     0  13.0M      0  0:00:01  0:00:01 --:--:-- 46.5M



```python
!gunzip -f tempfiles/3M-february-2018.txt.gz
```


```python
!tar -xf tempfiles/toy_read_ref_set.tar.gz -C tempfiles/
```

# Install Required Packages


```python
import sys

# Download and install simpleaf
!curl --proto '=https' --tlsv1.2 -LsSf https://github.com/COMBINE-lab/simpleaf/releases/download/v0.19.4/simpleaf-installer.sh | sh
!mv simpleaf-x86_64-unknown-linux-gnu/simpleaf simpleaf 2>/dev/null || true
!chmod +x simpleaf
!./simpleaf --version
```

    downloading simpleaf 0.19.4 x86_64-unknown-linux-gnu


    installing to /home/liam/.cargo/bin


      simpleaf
    everything's installed!


    simpleaf 0.19.4



```python
# Download and install salmon
!curl -L "https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz" -o salmon.tar.gz
!tar -xzf salmon.tar.gz
!chmod +x salmon-latest_linux_x86_64/bin/salmon

# Download and install alevin-fry
!curl -L "https://github.com/COMBINE-lab/alevin-fry/releases/download/v0.11.2/alevin-fry-x86_64-unknown-linux-gnu.tar.xz" -o alevin-fry.tar.xz
!tar -xf alevin-fry.tar.xz
!chmod +x alevin-fry-x86_64-unknown-linux-gnu/alevin-fry

# Add both to PATH for this session
import os
salmon_path = os.path.abspath('salmon-latest_linux_x86_64/bin')
alevin_fry_path = os.path.abspath('alevin-fry-x86_64-unknown-linux-gnu')
os.environ['PATH'] = f"{salmon_path}:{alevin_fry_path}:{os.environ['PATH']}"

!salmon --version
!alevin-fry --version
```

      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0

      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0

      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0


     43 88.9M   43 39.1M    0     0  36.1M      0  0:00:02  0:00:01  0:00:01 36.1M

    100 88.9M  100 88.9M    0     0  44.3M      0  0:00:02  0:00:02 --:--:-- 54.1M


      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0

      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0

      0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0


    100 1409k  100 1409k    0     0  4155k      0 --:--:-- --:--:-- --:--:-- 4155k


    salmon 1.10.0


    alevin-fry 0.11.2



```python
import scanpy as sc
import pyroe
import celltypist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor='white')
```

    /home/liam/Bioinformatics/fall25-csc-bioinf/.venv/lib/python3.12/site-packages/pyranges/__init__.py:18: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
      import pkg_resources
    /home/liam/Bioinformatics/fall25-csc-bioinf/.venv/lib/python3.12/site-packages/celltypist/classifier.py:11: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
      from scanpy import __version__ as scv


# Build Reference Index with Simpleaf


```python
!export ALEVIN_FRY_HOME=$(pwd) && export PATH="$(pwd)/salmon-latest_linux_x86_64/bin:$(pwd)/alevin-fry-x86_64-unknown-linux-gnu:$PATH" && ./simpleaf set-paths
```

    [2m2025-11-17T07:35:08.471049Z[0m [32m INFO[0m [2msimpleaf::utils::prog_utils[0m[2m:[0m could not find piscem executable, so salmon will be required.
    found `salmon` in the PATH at /home/liam/Bioinformatics/fall25-csc-bioinf/week6/salmon-latest_linux_x86_64/bin/salmon
    found `alevin-fry` in the PATH at /home/liam/Bioinformatics/fall25-csc-bioinf/week6/alevin-fry-x86_64-unknown-linux-gnu/alevin-fry


    [2m2025-11-17T07:35:08.491831Z[0m [32m INFO[0m [2msimpleaf::utils::prog_utils[0m[2m:[0m Could not find macs3 executable, peak calling cannot be peformed by simpleaf



```python
!export ALEVIN_FRY_HOME=$(pwd) && export PATH="$(pwd)/salmon-latest_linux_x86_64/bin:$(pwd)/alevin-fry-x86_64-unknown-linux-gnu:$PATH" && ./simpleaf index \
    --output simpleaf_index \
    --fasta tempfiles/toy_ref_read/toy_human_ref/fasta/genome.fa \
    --gtf tempfiles/toy_ref_read/toy_human_ref/genes/genes.gtf \
    --rlen 91 \
    --threads 2 \
    --no-piscem
```

    [2m2025-11-17T07:35:08.612451Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::indexing[0m[2m:[0m preparing to make reference with roers
    [2m2025-11-17T07:35:08.615028Z[0m [32m INFO[0m [2mgrangers::reader::gtf[0m[2m:[0m Finished parsing the input file. Found 3 comments and 2439 records.
    [2m2025-11-17T07:35:08.615494Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Built the Grangers object for 2439 records
    [2m2025-11-17T07:35:08.623963Z[0m [33m WARN[0m [2mgrangers::grangers_info[0m[2m:[0m The exon_number column contains null values. Will compute the exon number from exon start position .
    [2m2025-11-17T07:35:08.635366Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Found 2148 exon records from 271 transcripts.


    [2m2025-11-17T07:35:08.860523Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Wrote transcript sequences to output file.
    [2m2025-11-17T07:35:08.860546Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Processing intronic records.


    [2m2025-11-17T07:35:08.880050Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Found 1877 intronic records.
    [2m2025-11-17T07:35:08.880083Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Added flanking length to intronic records.
    [2m2025-11-17T07:35:08.884976Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Merged overlapping intronic records.


    [2m2025-11-17T07:35:09.016199Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Wrote intronic sequences to output file.


    [2m2025-11-17T07:35:09.021456Z[0m [32m INFO[0m [2mroers[0m[2m:[0m Done!
    [2m2025-11-17T07:35:09.022018Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::indexing[0m[2m:[0m salmon index cmd : /home/liam/Bioinformatics/fall25-csc-bioinf/week6/salmon-latest_linux_x86_64/bin/salmon index -k 31 -i simpleaf_index/index -t simpleaf_index/ref/roers_ref.fa --threads 2


    [2m2025-11-17T07:35:09.616503Z[0m [32m INFO[0m [2msimpleaf::utils::prog_utils[0m[2m:[0m command returned successfully (exit status: 0)
    [2m2025-11-17T07:35:09.616605Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::indexing[0m[2m:[0m Some("simpleaf_index/ref/gene_id_to_name.tsv")


# Quantify with Alevin-fry


```python
!export ALEVIN_FRY_HOME=$(pwd) && export PATH="$(pwd)/salmon-latest_linux_x86_64/bin:$(pwd)/alevin-fry-x86_64-unknown-linux-gnu:$PATH" && ./simpleaf quant \
    --reads1 tempfiles/toy_ref_read/toy_read_fastq/selected_R1_reads.fastq \
    --reads2 tempfiles/toy_ref_read/toy_read_fastq/selected_R2_reads.fastq \
    --threads 2 \
    --index simpleaf_index/index \
    --chemistry 10xv3 \
    --resolution cr-like \
    --expected-ori fw \
    --unfiltered-pl tempfiles/3M-february-2018.txt \
    --t2g-map simpleaf_index/index/t2g_3col.tsv \
    --output alevin_output
```

    [2m2025-11-17T07:35:09.736127Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::quant[0m[2m:[0m found uncompressed file
    [2m2025-11-17T07:35:09.736199Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::quant[0m[2m:[0m salmon alevin cmd : /home/liam/Bioinformatics/fall25-csc-bioinf/week6/salmon-latest_linux_x86_64/bin/salmon alevin --index simpleaf_index/index -l A -1 tempfiles/toy_ref_read/toy_read_fastq/selected_R1_reads.fastq -2 tempfiles/toy_ref_read/toy_read_fastq/selected_R2_reads.fastq --chromiumV3 --threads 2 -o alevin_output/af_map --sketch


    [2m2025-11-17T07:35:09.957041Z[0m [32m INFO[0m [2msimpleaf::utils::prog_utils[0m[2m:[0m command returned successfully (exit status: 0)
    [2m2025-11-17T07:35:09.957201Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::quant[0m[2m:[0m alevin-fry generate-permit-list cmd : /home/liam/Bioinformatics/fall25-csc-bioinf/week6/alevin-fry-x86_64-unknown-linux-gnu/alevin-fry generate-permit-list -i alevin_output/af_map -d fw -t 2 --unfiltered-pl tempfiles/3M-february-2018.txt --min-reads 10 -o alevin_output/af_quant


    [2m2025-11-17T07:35:11.254876Z[0m [32m INFO[0m [2msimpleaf::utils::prog_utils[0m[2m:[0m command returned successfully (exit status: 0)
    [2m2025-11-17T07:35:11.254908Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::quant[0m[2m:[0m alevin-fry collate cmd : /home/liam/Bioinformatics/fall25-csc-bioinf/week6/alevin-fry-x86_64-unknown-linux-gnu/alevin-fry collate -i alevin_output/af_quant -r alevin_output/af_map -t 2
    [2m2025-11-17T07:35:11.260238Z[0m [32m INFO[0m [2msimpleaf::utils::prog_utils[0m[2m:[0m command returned successfully (exit status: 0)
    [2m2025-11-17T07:35:11.260260Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::quant[0m[2m:[0m cmd : "/home/liam/Bioinformatics/fall25-csc-bioinf/week6/alevin-fry-x86_64-unknown-linux-gnu/alevin-fry" "quant" "-i" "alevin_output/af_quant" "-o" "alevin_output/af_quant" "-t" "2" "-m" "simpleaf_index/index/t2g_3col.tsv" "-r" "cr-like"
    [2m2025-11-17T07:35:11.264785Z[0m [32m INFO[0m [2msimpleaf::utils::prog_utils[0m[2m:[0m command returned successfully (exit status: 0)
    [2m2025-11-17T07:35:11.264891Z[0m [32m INFO[0m [2msimpleaf::simpleaf_commands::quant[0m[2m:[0m successfully copied the gene_name_to_id.tsv file into the quantification directory.


# Load Data into AnnData Object


```python
adata = pyroe.load_fry("alevin_output/af_quant", output_format='scRNA')
print(f"Loaded AnnData with shape: {adata.shape}")

# Convert Ensembl IDs to gene symbols for CellTypist compatibility
# Load the gene ID to gene name mapping
gene_mapping = pd.read_csv("simpleaf_index/ref/gene_id_to_name.tsv", sep="\t", header=None, names=["gene_id", "gene_name"])
gene_mapping = gene_mapping.set_index("gene_id")

# Map the gene IDs to gene names
adata.var["gene_ids"] = adata.var_names  # Keep original IDs
adata.var["gene_names"] = adata.var_names.map(lambda x: gene_mapping.loc[x, "gene_name"] if x in gene_mapping.index else x)
# Use gene symbols as the primary index
adata.var_names = adata.var["gene_names"]
adata.var_names.name = "gene_names"

print(f"Converted gene IDs to symbols. First 5 genes: {list(adata.var_names[:5])}")
adata
```

    USA mode: True
    Using pre-defined output format: scrna
    Will populate output field X with sum of counts frorm ['S', 'A'].
    Will combine ['U'] into output layer unspliced.


    Loaded AnnData with shape: (139, 20)
    Converted gene IDs to symbols. First 5 genes: ['NDFIP1', 'UBE2D2', 'APC', 'VCAN', 'UBE2B']





    AnnData object with n_obs × n_vars = 139 × 20
        obs: 'barcodes'
        var: 'gene_ids', 'gene_names'
        layers: 'unspliced'



# Quality Control and Preprocessing


```python
# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Filter cells and genes (adjusted for small toy dataset)
# For a toy dataset with only ~20 genes, we use much lower thresholds
sc.pp.filter_cells(adata, min_genes=3)  # Reduced from 200 to 3 for toy data
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.shape}")

# Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

# For this tiny toy dataset, mark all genes as highly variable
# This ensures that all genes are used in downstream analysis
adata.var['highly_variable'] = True

# Store the log-normalized data for CellTypist (before scaling)
# CellTypist requires log1p normalized expression to 10000 counts per cell
# IMPORTANT: Make a deep copy to preserve the state before scaling
adata.raw = adata.copy()

# Scale the data (will now scale all genes since they're all marked as highly variable)
# Scanpy's scale function handles zero-variance genes by centering them to 0
sc.pp.scale(adata, max_value=10)
```

    After filtering: (132, 19)



```python
# Perform PCA (adjusted for small toy dataset)
# With only ~19 genes, we can only compute a maximum of 18 principal components
n_comps = min(adata.n_vars - 1, 10)  # Use at most 10 PCs or n_genes-1, whichever is smaller
sc.tl.pca(adata, n_comps=n_comps, svd_solver='arpack')

# Compute the neighborhood graph
# Adjust n_pcs to match the number of components we actually computed
n_pcs_to_use = min(n_comps, 10)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs_to_use)

# Compute UMAP
sc.tl.umap(adata)
```

# Leiden Clustering


```python
# Perform Leiden clustering
sc.tl.leiden(adata, resolution=0.5)

print(f"Number of clusters: {len(adata.obs['leiden'].unique())}")
```

    Number of clusters: 3


# Visualization: Clustering Plot


```python
sc.pl.umap(adata, color=['leiden'], legend_loc='on data', title='Leiden Clustering', save='_leiden_clustering.png')
plt.show()
```

    WARNING: saving figure to file figures/umap_leiden_clustering.png



    
![png](week6_files/week6_22_1.png)
    


# Cell Type Annotation with CellTypist


```python
# Download and use a pre-trained model
celltypist.models.download_models(force_update=False)

# Show available models
print("Available models:")
print(celltypist.models.models_description())
```

    📂 Storing models in /home/liam/.celltypist/data/models


    ⏩ Skipping [1/58]: Immune_All_Low.pkl (file exists)


    ⏩ Skipping [2/58]: Immune_All_High.pkl (file exists)


    ⏩ Skipping [3/58]: Adult_COVID19_PBMC.pkl (file exists)


    ⏩ Skipping [4/58]: Adult_CynomolgusMacaque_Hippocampus.pkl (file exists)


    ⏩ Skipping [5/58]: Adult_Human_MTG.pkl (file exists)


    ⏩ Skipping [6/58]: Adult_Human_PancreaticIslet.pkl (file exists)


    ⏩ Skipping [7/58]: Adult_Human_PrefrontalCortex.pkl (file exists)


    ⏩ Skipping [8/58]: Adult_Human_Skin.pkl (file exists)


    ⏩ Skipping [9/58]: Adult_Human_Vascular.pkl (file exists)


    ⏩ Skipping [10/58]: Adult_Mouse_Gut.pkl (file exists)


    ⏩ Skipping [11/58]: Adult_Mouse_OlfactoryBulb.pkl (file exists)


    ⏩ Skipping [12/58]: Adult_Pig_Hippocampus.pkl (file exists)


    ⏩ Skipping [13/58]: Adult_RhesusMacaque_Hippocampus.pkl (file exists)


    ⏩ Skipping [14/58]: Adult_cHSPCs_Illumina.pkl (file exists)


    ⏩ Skipping [15/58]: Adult_cHSPCs_Ultima.pkl (file exists)


    ⏩ Skipping [16/58]: Autopsy_COVID19_Lung.pkl (file exists)


    ⏩ Skipping [17/58]: COVID19_HumanChallenge_Blood.pkl (file exists)


    ⏩ Skipping [18/58]: COVID19_Immune_Landscape.pkl (file exists)


    ⏩ Skipping [19/58]: Cells_Adult_Breast.pkl (file exists)


    ⏩ Skipping [20/58]: Cells_Fetal_Lung.pkl (file exists)


    ⏩ Skipping [21/58]: Cells_Human_Tonsil.pkl (file exists)


    ⏩ Skipping [22/58]: Cells_Intestinal_Tract.pkl (file exists)


    ⏩ Skipping [23/58]: Cells_Lung_Airway.pkl (file exists)


    ⏩ Skipping [24/58]: Developing_Human_Brain.pkl (file exists)


    ⏩ Skipping [25/58]: Developing_Human_Gonads.pkl (file exists)


    ⏩ Skipping [26/58]: Developing_Human_Hippocampus.pkl (file exists)


    ⏩ Skipping [27/58]: Developing_Human_Organs.pkl (file exists)


    ⏩ Skipping [28/58]: Developing_Human_Thymus.pkl (file exists)


    ⏩ Skipping [29/58]: Developing_Mouse_Brain.pkl (file exists)


    ⏩ Skipping [30/58]: Developing_Mouse_Hippocampus.pkl (file exists)


    ⏩ Skipping [31/58]: Fetal_Human_AdrenalGlands.pkl (file exists)


    ⏩ Skipping [32/58]: Fetal_Human_Pancreas.pkl (file exists)


    ⏩ Skipping [33/58]: Fetal_Human_Pituitary.pkl (file exists)


    ⏩ Skipping [34/58]: Fetal_Human_Retina.pkl (file exists)


    ⏩ Skipping [35/58]: Fetal_Human_Skin.pkl (file exists)


    ⏩ Skipping [36/58]: Healthy_Adult_Heart.pkl (file exists)


    ⏩ Skipping [37/58]: Healthy_COVID19_PBMC.pkl (file exists)


    ⏩ Skipping [38/58]: Healthy_Human_Liver.pkl (file exists)


    ⏩ Skipping [39/58]: Healthy_Mouse_Liver.pkl (file exists)


    ⏩ Skipping [40/58]: Human_AdultAged_Hippocampus.pkl (file exists)


    ⏩ Skipping [41/58]: Human_Colorectal_Cancer.pkl (file exists)


    ⏩ Skipping [42/58]: Human_Developmental_Retina.pkl (file exists)


    ⏩ Skipping [43/58]: Human_Embryonic_YolkSac.pkl (file exists)


    ⏩ Skipping [44/58]: Human_Endometrium_Atlas.pkl (file exists)


    ⏩ Skipping [45/58]: Human_IPF_Lung.pkl (file exists)


    ⏩ Skipping [46/58]: Human_Longitudinal_Hippocampus.pkl (file exists)


    ⏩ Skipping [47/58]: Human_Lung_Atlas.pkl (file exists)


    ⏩ Skipping [48/58]: Human_PF_Lung.pkl (file exists)


    ⏩ Skipping [49/58]: Human_Placenta_Decidua.pkl (file exists)


    ⏩ Skipping [50/58]: Lethal_COVID19_Lung.pkl (file exists)


    ⏩ Skipping [51/58]: Mouse_Dentate_Gyrus.pkl (file exists)


    ⏩ Skipping [52/58]: Mouse_Isocortex_Hippocampus.pkl (file exists)


    ⏩ Skipping [53/58]: Mouse_Postnatal_DentateGyrus.pkl (file exists)


    ⏩ Skipping [54/58]: Mouse_Whole_Brain.pkl (file exists)


    ⏩ Skipping [55/58]: Nuclei_Lung_Airway.pkl (file exists)


    ⏩ Skipping [56/58]: PaediatricAdult_COVID19_Airway.pkl (file exists)


    ⏩ Skipping [57/58]: PaediatricAdult_COVID19_PBMC.pkl (file exists)


    ⏩ Skipping [58/58]: Pan_Fetal_Human.pkl (file exists)


    👉 Detailed model information can be found at `https://www.celltypist.org/models`


    Available models:
                                          model  \
    0                        Immune_All_Low.pkl   
    1                       Immune_All_High.pkl   
    2                    Adult_COVID19_PBMC.pkl   
    3   Adult_CynomolgusMacaque_Hippocampus.pkl   
    4                       Adult_Human_MTG.pkl   
    5           Adult_Human_PancreaticIslet.pkl   
    6          Adult_Human_PrefrontalCortex.pkl   
    7                      Adult_Human_Skin.pkl   
    8                  Adult_Human_Vascular.pkl   
    9                       Adult_Mouse_Gut.pkl   
    10            Adult_Mouse_OlfactoryBulb.pkl   
    11                Adult_Pig_Hippocampus.pkl   
    12      Adult_RhesusMacaque_Hippocampus.pkl   
    13                Adult_cHSPCs_Illumina.pkl   
    14                  Adult_cHSPCs_Ultima.pkl   
    15                 Autopsy_COVID19_Lung.pkl   
    16         COVID19_HumanChallenge_Blood.pkl   
    17             COVID19_Immune_Landscape.pkl   
    18                   Cells_Adult_Breast.pkl   
    19                     Cells_Fetal_Lung.pkl   
    20                   Cells_Human_Tonsil.pkl   
    21               Cells_Intestinal_Tract.pkl   
    22                    Cells_Lung_Airway.pkl   
    23               Developing_Human_Brain.pkl   
    24              Developing_Human_Gonads.pkl   
    25         Developing_Human_Hippocampus.pkl   
    26              Developing_Human_Organs.pkl   
    27              Developing_Human_Thymus.pkl   
    28               Developing_Mouse_Brain.pkl   
    29         Developing_Mouse_Hippocampus.pkl   
    30            Fetal_Human_AdrenalGlands.pkl   
    31                 Fetal_Human_Pancreas.pkl   
    32                Fetal_Human_Pituitary.pkl   
    33                   Fetal_Human_Retina.pkl   
    34                     Fetal_Human_Skin.pkl   
    35                  Healthy_Adult_Heart.pkl   
    36                 Healthy_COVID19_PBMC.pkl   
    37                  Healthy_Human_Liver.pkl   
    38                  Healthy_Mouse_Liver.pkl   
    39          Human_AdultAged_Hippocampus.pkl   
    40              Human_Colorectal_Cancer.pkl   
    41           Human_Developmental_Retina.pkl   
    42              Human_Embryonic_YolkSac.pkl   
    43              Human_Endometrium_Atlas.pkl   
    44                       Human_IPF_Lung.pkl   
    45       Human_Longitudinal_Hippocampus.pkl   
    46                     Human_Lung_Atlas.pkl   
    47                        Human_PF_Lung.pkl   
    48               Human_Placenta_Decidua.pkl   
    49                  Lethal_COVID19_Lung.pkl   
    50                  Mouse_Dentate_Gyrus.pkl   
    51          Mouse_Isocortex_Hippocampus.pkl   
    52         Mouse_Postnatal_DentateGyrus.pkl   
    53                    Mouse_Whole_Brain.pkl   
    54                   Nuclei_Lung_Airway.pkl   
    55       PaediatricAdult_COVID19_Airway.pkl   
    56         PaediatricAdult_COVID19_PBMC.pkl   
    57                      Pan_Fetal_Human.pkl   
    
                                              description  
    0   immune sub-populations combined from 20 tissue...  
    1   immune populations combined from 20 tissues of...  
    2   peripheral blood mononuclear cell types from C...  
    3   cell types from the hippocampus of adult cynom...  
    4   cell types and subtypes (10x-based) from the a...  
    5   cell types from pancreatic islets of healthy a...  
    6   cell types and subtypes from the adult human d...  
    7            cell types from human healthy adult skin  
    8   vascular populations combined from multiple ad...  
    9   cell types in the adult mouse gut combined fro...  
    10   cell types from the olfactory bulb of adult mice  
    11          cell types from the adult pig hippocampus  
    12  cell types from the hippocampus of adult rhesu...  
    13  human circulating hematopoietic stem and proge...  
    14  human circulating hematopoietic stem and proge...  
    15  cell types from the lungs of 16 SARS-CoV-2 inf...  
    16  detailed blood cell states from 16 individuals...  
    17  immune subtypes from lung and blood of COVID-1...  
    18             cell types from the adult human breast  
    19    cell types from human embryonic and fetal lungs  
    20      tonsillar cell types from humans (3-65 years)  
    21  intestinal cells from fetal, pediatric (health...  
    22  cell populations from scRNA-seq of five locati...  
    23  cell types from the first-trimester developing...  
    24  cell types of human gonadal and adjacent extra...  
    25   cell types from the developing human hippocampus  
    26  cell types of five endoderm-derived organs in ...  
    27  cell populations in embryonic, fetal, pediatri...  
    28  cell types from the embryonic mouse brain betw...  
    29  cell types from the mouse hippocampus at postn...  
    30  cell types of human fetal adrenal glands from ...  
    31  pancreatic cell types from human embryos at 9-...  
    32  cell types of human fetal pituitaries from 7 t...  
    33  cell types from human fetal neural retina and ...  
    34        cell types from developing human fetal skin  
    35  cell types from eight anatomical regions of th...  
    36  peripheral blood mononuclear cell types from h...  
    37  cell types from scRNA-seq and snRNA-seq of the...  
    38  cell types from scRNA-seq and snRNA-seq of the...  
    39  cell types from the hippocampus of adult and a...  
    40  cell types of colon tissues from patients with...  
    41                 cell types from human fetal retina  
    42  cell types of the human yolk sac from 4-8 post...  
    43  endometrial cell types integrated from seven d...  
    44  cell types from idiopathic pulmonary fibrosis,...  
    45  cell types from the adult human anterior and p...  
    46  integrated Human Lung Cell Atlas (HLCA) combin...  
    47  cell types from different forms of pulmonary f...  
    48  cell types from first-trimester human placenta...  
    49  cell types from the lungs of individuals who d...  
    50  cell types from the dentate gyrus in perinatal...  
    51  cell types from the adult mouse isocortex (neo...  
    52  cell types from the mouse dentate gyrus in pos...  
    53        cell types from the whole adult mouse brain  
    54  cell populations from snRNA-seq of five locati...  
    55  cell types in the airway of paediatric and adu...  
    56  peripheral blood mononuclear cell types of pae...  
    57  stromal and immune populations from the human ...  



```python
# Use the Immune_All_Low model for annotation
# NOTE: This toy dataset has only 19 genes, so results may be limited
# For real datasets with thousands of genes, this would work much better
model = celltypist.models.Model.load(model='Immune_All_Low.pkl')

# Predict cell types
predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl', majority_voting=True)

# Transfer predictions to the original AnnData object
adata.obs['cell_type'] = predictions.predicted_labels.majority_voting
adata.obs['cell_type_predicted'] = predictions.predicted_labels.predicted_labels

print(f"Unique cell types: {adata.obs['cell_type'].unique()}")
```

    👀 Invalid expression matrix in `.X`, expect log1p normalized expression to 10000 counts per cell; will use `.raw.X` instead


    🔬 Input data has 132 cells and 19 genes


    🔗 Matching reference genes in the model


    🧬 10 features used for prediction


    ⚖️ Scaling input data


    🖋️ Predicting labels


    ✅ Prediction done!


    👀 Detected a neighborhood graph in the input object, will run over-clustering on the basis of it


    ⛓️ Over-clustering input data with resolution set to 5


    🗳️ Majority voting the predictions


    ✅ Majority voting done!


    Unique cell types: ['Epithelial cells', 'Fibroblasts', 'Double-positive thymocytes', 'Tcm/Naive helper T cells']
    Categories (4, object): ['Double-positive thymocytes', 'Epithelial cells', 'Fibroblasts', 'Tcm/Naive helper T cells']


# Visualization: Annotated Plot with Cell Types


```python
sc.pl.umap(adata, color=['cell_type'], title='Cell Type Annotation', save='_cell_type_annotation.png')
plt.show()
```

    WARNING: saving figure to file figures/umap_cell_type_annotation.png



    
![png](week6_files/week6_27_1.png)
    



```python
# Also plot both clustering and cell type side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

sc.pl.umap(adata, color=['leiden'], legend_loc='on data', title='Leiden Clustering', ax=ax1, show=False)
sc.pl.umap(adata, color=['cell_type'], title='Cell Type Annotation', ax=ax2, show=False)

plt.tight_layout()
plt.savefig('output/combined_clustering_annotation.png', dpi=150, bbox_inches='tight')
plt.show()
```


    
![png](week6_files/week6_28_0.png)
    


# Summary


```python
print("="*60)
print("SINGLE-CELL RNA-SEQ ANALYSIS COMPLETE")
print("="*60)
print(f"\nDataset shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"Number of Leiden clusters: {len(adata.obs['leiden'].unique())}")
print(f"Number of unique cell types identified: {len(adata.obs['cell_type'].unique())}")
print(f"\nCell type distribution:")
print(adata.obs['cell_type'].value_counts())
print("\n" + "="*60)
```

    ============================================================
    SINGLE-CELL RNA-SEQ ANALYSIS COMPLETE
    ============================================================
    
    Dataset shape: 132 cells × 19 genes
    Number of Leiden clusters: 3
    Number of unique cell types identified: 4
    
    Cell type distribution:
    cell_type
    Epithelial cells              103
    Tcm/Naive helper T cells       17
    Double-positive thymocytes      7
    Fibroblasts                     5
    Name: count, dtype: int64
    
    ============================================================

