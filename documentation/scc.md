#### Typical Commands for Shared Computing Cluster (SCC)

```bash
ssh <username>@scc4.bu.edu     # BU Kerberos ID and Password
qrsh -P <projectname>          # Must specify project for qsub
module load anaconda2          # Enable conda
source activate <environment>  # Activate conda environment
cd path/to/nextflow            # Go to nextflow project
./nextflow <workflow.nf>       # Run nextflow workflow
```