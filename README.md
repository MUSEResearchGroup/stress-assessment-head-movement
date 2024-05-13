# Stress Assessment for Augmented Reality Applications based on Head Movement Features


In the repository, you can find the codes pertinent to the paper: 


A. Ferrarotti, S. Baldoni, M. Carli and F. Battisti, "Stress Assessment for Augmented Reality Applications Based on Head Movement Features," in IEEE Transactions on Visualization and Computer Graphics, doi: 10.1109/TVCG.2024.3385637.

## Repository structure
The repository contains the following items:

- **Codes**:
    - stats_ANOVA is a MATLAB code containing all the ANOVA tests performed. In addition, it allows to generate Figure 6(a) of the paper.
    - stats_WelchANOVA is a Python code containing the Welch's ANOVA tests performed of the STFT coefficients of the analyzed features.
    - nasa_results is a MATLAB code that allows to replicate the t-test performed to evaluate whether the NASA-TLX results of the Stroop Color Word Test and the Mental        Arithmetic test are statistically different.
    - phaseStress_classifier is a MATLAB code used to train and validate the classifiers of the proposed architecture.
    - testStressClassifier is a MATLAB code that allows to test the proposed architecture.


- **Functions**
    - getFeatures_fixed is a MATLAB function used to extract the recorded head movement data from the .txt files related to the users performing the Stroop Color Word         test.
    - getFeatures_fixed is a MATLAB function used to extract the recorded head movement data from the .txt files related to the users performing the Mental Arithmetic         test.
    - hypothesis_check is a MATLAB function used to verify the hypothesis required by the ANOVA tests for the selected head movement features.

      
- **Data**
  - **trained_models_unw** contains the Support Vector Machines models trained on the Stroop Color Word Test and that are used in the testStressClassifier code.

    
- **Folders:**
  - **file_paths** contains all the file paths pertinent to the recorded head movement data.
  - **Head Movements Dataset.zip** contains the complete dataset acquired for the study. Further details on its structure are given in the next section. **Plese make         sure to unzip the folder and place it in the same folder as the other codes to run them.**
  - **Stats-Stroop STFT** contains the STFT coefficients of the analyzed features that are used for statistical analysis in the code stats_WelchANOVA.py.
 
  ## Dataset

  
  
  ## Strumenti usati

  ## Come citare
