# Noise-Analysis, R code for gene expression noise analysis

By Tao Hu.

The code for the official implementation of paper **[Single cell transcriptomes reveal characteristics of miRNA in gene expression noise reduction](https://www.biorxiv.org/content/10.1101/465518v1)**. 

The whole code will be made publicly available shortly.

## Data download
The data can be download in [Google Driver](https://drive.google.com/open?id=1I2mDz1z-LQZl9FawbB8Q4Swf1C2DFsPG). If you want to use it, please put the data folder under the project path.

**Env Requirements:** 
- R 3.4.2 or higher.

## Run scripts for a demo

0. Download data in [Google Driver](https://drive.google.com/open?id=1I2mDz1z-LQZl9FawbB8Q4Swf1C2DFsPG), extract files and put the data folder in the project directory you defined. 

1. After the download, link your datasets to the current directory, like,
    ```
    cd ./project_path
    Rscript targetscan_mapping.R
    Rscript miRNA_exp_processing.R
    Rscript gene_noise_deconvlution.R
    Rscript noise_bias_test.R
    ```

## License
For academic use, this project is licensed under the 2-clause BSD License - see the LICENSE file for details. For commercial use, please contact the authors. 