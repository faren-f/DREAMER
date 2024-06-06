# DREAMER: Exploring Common Mechanisms of Adverse Drug Reactions and Disease Phenotypes through Network-Based Analysis

Adverse drug reactions (ADRs) are a major concern in clinical healthcare, significantly affecting patient safety and drug development. The need for a deeper understanding of ADR mechanisms is crucial for improving drug safety profiles in drug design and drug repurposing. This study introduces DREAMER (Drug adverse REAction Mechanism ExplaineR), a novel network-based method for exploring the mechanisms underlying adverse drug reactions and disease phenotypes at a molecular level by leveraging a comprehensive knowledge graph obtained from various datasets. By considering drugs and diseases that cause similar phenotypes and investigating their commonalities regarding their impact on specific modules of the protein-protein interaction network, DREAMER can robustly identify protein sets associated with the biological mechanisms underlying ADRs and unravel the causal relationships that contribute to the observed clinical outcomes. Applying DREAMER to 649 ADRs, we identified the mechanism of action for 67 ADRs across multiple organ systems, e.g., ventricular arrhythmia, metabolic acidosis, and interstitial pneumonitis. In particular, DREAMER highlights the importance of GABAergic signaling and proteins of the coagulation pathways for personality disorders and intracranial hemorrhage, respectively. We further demonstrate the application of DREAMER in drug repurposing and propose sotalol (targeting KCNH2), ranolazine (targeting SCN5A, currently under clinical trial), and diltiazem (indicated drug targeting CACNA1C and SCN3A) as candidate drugs to be repurposed for cardiac arrest. In summary, DREAMER effectively detects molecular mechanisms underlying phenotypes emphasizing the importance of network-based analyses with integrative data for enhancing drug safety and accelerating the discovery of novel therapeutic strategies.

## Instructions

All the data are available in the `Data/` folder and all the scripts are available in `Scripts/` to reproduce the results.

### Requirements

- **R version 4.3.2**
- **RStudio version 2023.09.1+494**

#### R Packages:
- igraph version 2.0.3
- MASS version 7.3.60.0.1
- parallel version 4.3.2
- dplyr version 1.1.4
- clusterProfiler version 4.10.0
- org.Hs.eg.db version 3.18.0
- ggplot2 version 3.5.0
- openxlsx version 4.2.5.2
- biomaRt version 2.58.2
- ReactomePA version 1.46.0
- writexl version 1.5.0

The scripts that reproduce the ADR-related proteins, DP-related proteins, and ADR-DP-related proteins as well as controlling for confounding effects of indications of drugs are provided in the `scripts/main` folder.
The scripts that reproduce the figures and tables are provided in the `scripts/visualization` folder.
The scripts that reproduce the validations including comparing ADR-related proteins and DP-related proteins with the baseline method and holdout validation method are provided in the `scripts/validation` folder.


### Folder Structure

```
DREAMER
├── data/
│   ├── Knowledge_graph/
│   ├── preprocessed_graph/
│   └── other_resources/
├── scripts/
│   ├── main/
│   ├── validation/
│   └── visualization/
└── README.md
```

## How to Cite Us in the Future (Manuscript is under review)

If you use DREAMER in your research, please cite our work as follows:

Farzaneh Firoozbakht. (202.). DREAMER: Exploring Common Mechanisms of Adverse Drug Reactions and Disease Phenotypes through Network-Based Analysis. Journal Name, Volume(Issue), Page Numbers. DOI

@article{faren_dreamer,
  title={DREAMER: Exploring Common Mechanisms of Adverse Drug Reactions and Disease Phenotypes through Network-Based Analysis},
  author={Farzaneh Firoozbakht},
  journal={-},
  volume={-},
  number={-},
  pages={-},
  year={-},
  publisher={-},
  doi={-}
}

