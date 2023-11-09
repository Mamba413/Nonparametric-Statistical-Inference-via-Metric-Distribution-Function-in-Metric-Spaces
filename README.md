# Reproducible materials
This repository contains scripts to run the simulation described in [Nonparametric Statistical Inference via Metric Distribution Function in Metric Spaces](https://arxiv.org/abs/2107.07317). 

![](real_theory_bridge.png)

## Software

Please install `Ball` 1.3.13 to use the tests mentioned in the paper. 
> Note: download this package from [here](https://github.com/Mamba413/Nonparametric-Statistical-Inference-via-Metric-Distribution-Function-in-Metric-Spaces/blob/main/Ball_1.3.13.tar.gz). Some advanced features haven't been integrated into the `Ball` package (version 1.3.13) in R CRAN. We expect to make the `Ball` package on R CRAN support these new features in version 1.3.14.

## Reproducible code

We have organized our scripts to improve the reproducibility of experiments. Specifically, each R script file corresponds to a certain part of the results in the paper, which are listed below:

- `large_homogeneity_n.R` <-> Results in Figure 3A
- `large_joint_independence_n.R` <-> Results in Figure 3B
- `compare_homo_methods.R` <-> Results in Figure 4A
- `compare_indep_methods.R` <-> Results in Figure 4B
- `real_data_adni_analysis.R` <-> Results in Table 1
- `real_data_adhd200_analysis.R` <-> Results in Table 2

Notice that, before conducting these scripts, please modify the your_path in each script accordingly. 

## Citations

Please cite the following publications if you make use of the material here.

- Xueqin Wang, Jin Zhu, Wenliang Pan, Junzhu Zhu, and Heping Zhang (JASA, 2023+). Nonparametric statistical inference via metric distribution function in metric spaces. arXiv preprint arXiv:2107.07317.

The corresponding BibteX entries:

```
@article{wang2023nonparametric,
  title={Nonparametric statistical inference via metric distribution function in metric spaces},
  author={Wang, Xueqin and Zhu, Jin and Pan, Wenliang and Zhu, Junhao and Zhang, Heping},
  journal={Journal of American Statistical Association},
  year={2023+},
  doi={10.1080/01621459.2023.2277417},
}
```


## Contact
Please direct questions and comments to the [issues page](https://github.com/Mamba413/Nonparametric-Statistical-Inference-via-Metric-Distribution-Function-in-Metric-Spaces/issues).
