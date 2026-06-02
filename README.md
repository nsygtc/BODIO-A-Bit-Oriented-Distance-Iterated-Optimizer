# BODIO: A Bit-Oriented Distance-Iterated Optimizer for Key Search with Bitwise Linear Leakage

[![Paper](https://img.shields.io/badge/Paper-Applied%20Soft%20Computing-blue)](https://doi.org/10.1016/j.asoc.2026.115531)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.asoc.2026.115531-blue)](https://doi.org/10.1016/j.asoc.2026.115531)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

This repository provides the official implementation of **BODIO**, a **Bit-Oriented Distance-Iterated Optimizer** for key search under bitwise linear leakage.

The associated paper has been published in *Applied Soft Computing*:

> **BODIO: A bit-oriented distance-iterated optimizer for key search with bitwise linear leakage**
> Ting Gao, Xuyang Chen, Yuanxing Liu, An Wang
> *Applied Soft Computing*, Volume 201, 2026, Article 115531
> DOI: [10.1016/j.asoc.2026.115531](https://doi.org/10.1016/j.asoc.2026.115531)

Official Elsevier Share Link:
https://authors.elsevier.com/a/1nB4F5aecT0-B6

---

## Overview

Key search with bitwise linear leakage aims to recover the correct cryptographic key from side-channel power traces. This problem is challenging because the search space is discrete, high-dimensional, and strongly affected by physical noise.

BODIO reformulates the key search problem as a high-dimensional binary optimization problem and introduces a bit-oriented optimization paradigm. Instead of treating a candidate key only as a whole block, BODIO explicitly optimizes individual key bits while evaluating the candidate key at the population level.

The optimizer is inspired by the binary Grey Wolf Optimizer and is redesigned for key search under bitwise linear leakage. It integrates distance-iterated population updates, hierarchical leader guidance, individual diversification, and one-bit-flip to improve robustness under noisy fitness evaluations.

---

## Main Features

* **Bit-oriented optimization** for key search under bitwise linear leakage.
* **Distance-iterated population update** inspired by binary Grey Wolf Optimization.
* **Hierarchical guidance** using multiple leader candidates.
* **Individual diversifier** to reduce premature convergence and maintain population diversity.
* **One-bit-flip** for local search after dynamic convergence.
* **Simulation and FPGA experiment support**, depending on available data.

---

## Repository Structure

```text
BODIO/
├── README.md
├── LICENSE
├── bodio/
│   ├── __init__.py
│   ├── optimizer.py
```

You may adjust the structure according to your own implementation.

---

## Data


---

## Citation

If you use this repository or find this work helpful, please cite our paper:

```bibtex
@article{gao2026bodio,
  title={BODIO: A bit-oriented distance-iterated optimizer for key search with bitwise linear leakage},
  author={Gao, Ting and Chen, Xuyang and Liu, Yuanxing and Wang, An},
  journal={Applied Soft Computing},
  pages={115531},
  year={2026},
  publisher={Elsevier}
}
```

---

## Paper Access

The official published article is available through Elsevier:

* DOI: https://doi.org/10.1016/j.asoc.2026.115531
* Elsevier Share Link: https://authors.elsevier.com/a/1nB4F5aecT0-B6

This repository does not host the published journal article PDF. Please use the official DOI or Elsevier Share Link to access the paper.

---

## License

The source code in this repository is released under the license specified in the `LICENSE` file.

Please note that the copyright of the published article belongs to Elsevier B.V. The article should be accessed through the official DOI or Elsevier Share Link.

---

## Acknowledgements

We thank the reviewers and editors of *Applied Soft Computing* for their constructive comments.

---

## Contact

For questions, issues, or possible collaborations, please open an issue in this repository or contact the authors via email listed in the paper.
