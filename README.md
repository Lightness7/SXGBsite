# SXGBsite  
SXGBsite is a method to predict protein-ligand binding site residues, used PSSM-DCT and PSA to extract features containing sequence information, and a prediction model was constructed using SMOTE and XGBoost.  
  
## Install  
- Libsvm needs to be installed in MATLAB path to achieve feature normalization. 

    Please add 'X:\XXX\libsvm-3.21\matlab\' to the MATLAB path, where 'X' stands for path in the computer. Feature extraction mainly uses the 'line_map.m' function.  

- The python code in SXGBsite was written in python3.7. 

    The python packages of math, numpy, pandas, xgboost, scipy, matplotlib.pyplot, sklearn, imblearn and collections are required, and the version of packages has less impact on the running code.  
- Python environment build recommend Anaconda3 (download at: https://www.anaconda.com/distribution/), examples of how to install python package in Anaconda Prompt, as follows:  

    ```
    pip install numpy
    ```
    
    ```
    pip install sklearn
    ```
    etc.  
  
## Benchmark Datasets  
The benchmarks were constructed by Yu et al.[1], the PSSM files were obtained by PSI-BLAST, and the PSA files were obtained by Sann.  
  
[1] Yu, D.J.; Hu, J.; Yang, J.; Shen, H.B.; Tang, J.; Yang, J.Y. Designing template-free predictor for targeting protein-ligand binding sites with classifier ensemble and spatial clustering. IEEE/ACM Trans. Comput. Biol.Bioinform. 2013, 10, 994-1008.  
  
## Feature Extraction  
- The PSSM-DCT, PSSM-DWT, and PSA features used in this method and the corresponding MATLAB codes refer to the following two papers:  
  
  [2] Ding, Y.; Tang, J.; Guo, F. Identification of protein–ligand binding sites by sequence information and ensemble classifier. J. Chem. Inf. Model. 2017, 57, 3149-3161.  

  [3] Shen, C.; Ding, Y.; Tang, J.; Song, J.; Guo, F. Identification of DNA–protein binding sites through multi-scale local average blocks on sequence information. Molecules 2017, 22, 2079.  

  Thanks to their unselfish sharing of open source codes to help more researchers continue their research, including me.  

- The datasets after feature extraction are obtained by running 'Extract_Feature_Demo.m'. The input are PSSM files and PSA files of protein sequence. Different dataset feature extraction only needs to modify the input path, as follows:  

  Modify the path in line 8, 11, 17, 67, 102, 105, 111, 161, and 202 to the local path in 'Extract_Feature_Demo.m'.

- The datasets after feature extraction are downloaded at: 
  
  Feature data sets of metal ions  
  https://figshare.com/s/e51520c754d931649bab  

  Feature data sets of nucleotide, dna and heme  
  https://figshare.com/s/2f72d2f28614d472e95e  
  
## Prediction Model Construction  
The predictive model construction of this method only needs to run 'SXGBsite_demo.py' in a fully configured python environment. If you want to get the prediction results of different datasets, you only need to modify the path as follows:  

The results of different datasets are obtained by modifying the path in line 17 of 'SXGBsite_demo.py'.
```
data_path = "D:\Data\GTP\GTP_DCT_PSA.mat"
```
  
## Example  
Here, GTP and Mg datasets are taken as examples, and the following operations are performed:  
  
- Create a document 'Data' on the (D:) drive (or other, please choose by yourself) and put the downloaded feature datasets into this document. Run ‘SXGBsite_demo.py’ to get the results of SXGBsite on the GTP independent test set.  
  
- Modify line 17 of 'SXGBsite_demo.py'  

  ```
  data_path = "D:\Data\GTP\GTP_DCT_PSA.mat"
  ```  
  as follows:  

  ```
  data_path = "D:\Data\MG\MG_DCT_PSA.mat"
  ```  
  The results of SXGBsite on the Mg independent test set are available.
