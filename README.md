### 
FLAT-FIELDING OF SOLAR Hα OBSERVATIONS BASED ON THE MAXIMUM CORRENTROPY CRITERION

#### Abstract
The flat-field CCD calibration method of Kuhn et al. (KLL) is an efficient method for flat-fielding. However, since
it depends on the minimum of the sum of squares error (SSE), its solution is sensitive to noise, especially non-
Gaussian noise. In this paper, a new algorithm is proposed to determine the flat field. The idea is to change the
criterion of gain estimate from SSE to the maximum correntropy. The result of a test on simulated data
demonstrates that our method has a higher accuracy and a faster convergence than KLL’s and Chae’s. It has been
found that the method effectively suppresses noise, especially in the case of typical non-Gaussian noise. And the
computing time of our algorithm is the shortest.


#### Usage
You can run the test_**.m file to understand how to use the code. After testing, it has been successfully executed on MATLAB 2018b version.
```angular2html
'test_adative.m' is the code for testing the convergence criterion of the algorithm.
'test_adative_sigma.m' is the code for the algorithm's adaptive sigma criterion.
'test_std_no_noise.m' is the code to assess the impact of non-Gaussian noise (salt-and-pepper noise) on the algorithm.
'test_std_noise.m' is the code to assess the impact of Gaussian noise on the algorithm.
'test_std_poisson.m' is the code to assess the impact of Poisson noise on the algorithm.
```
*Note* the code uploaded in this repository is the final version from the [paper](https://iopscience.iop.org/article/10.3847/0004-637X/827/2/137) and is considered stable (no longer maintained).


#### Citation
If you use this code in a scientific publication, I would appreciate citation/reference to [this paper](https://iopscience.iop.org/article/10.3847/0004-637X/827/2/137).
```
@article{Xu_2016,
doi = {10.3847/0004-637X/827/2/137},
url = {https://dx.doi.org/10.3847/0004-637X/827/2/137},
year = {2016},
month = {aug},
publisher = {The American Astronomical Society},
volume = {827},
number = {2},
pages = {137},
author = {Gao-Gui Xu and Sheng Zheng and Gang-Hua Lin and Xiao-Fan Wang},
title = {FLAT-FIELDING OF SOLAR Hα OBSERVATIONS BASED ON THE MAXIMUM CORRENTROPY CRITERION},
journal = {The Astrophysical Journal},
abstract = {The flat-field CCD calibration method of Kuhn et al. (KLL) is an efficient method for flat-fielding. However, since it depends on the minimum of the sum of squares error (SSE), its solution is sensitive to noise, especially non-Gaussian noise. In this paper, a new algorithm is proposed to determine the flat field. The idea is to change the criterion of gain estimate from SSE to the maximum correntropy. The result of a test on simulated data demonstrates that our method has a higher accuracy and a faster convergence than KLL’s and Chae’s. It has been found that the method effectively suppresses noise, especially in the case of typical non-Gaussian noise. And the computing time of our algorithm is the shortest.}
}

```

