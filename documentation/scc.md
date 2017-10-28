#### Typical Commands for Shared Computing Cluster (SCC)

```bash
ssh <username>@scc4.bu.edu                   # BU Kerberos username and password
qrsh -P <projectname>                        # Must specify project for qsub
module load anaconda2                        # Enable conda
source activate <environment>                # Activate conda environment
cd path/to/nextflow                          # Go to nextflow project
./nextflow <workflow.nf> -c <params.config>  # Run nextflow workflow
```

*Resources for New Users*  
[Quick Start](http://www.bu.edu/tech/support/research/system-usage/scc-quickstart/)  
[Detailed Introduction](http://www.bu.edu/tech/files/2016/09/2016_fall-Tutorial-Intro-to-SCC.pdf)  
[VPN Setup for Linux](http://www.bu.edu/tech/services/cccs/remote/vpn/use/linux/)  
[VPN Setup for Mac](http://www.bu.edu/tech/services/cccs/remote/vpn/use/mac/)
