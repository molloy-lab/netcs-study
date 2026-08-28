# Evaluation, Selection, and Plotting Scripts

This directory contains the scripts used to evaluate simulated and
estimated phylogenetic trees and networks, select the CAMUS network used
for evaluation, and generate the draft figures for the NetCS experiments.

These scripts are **post-processing, evaluation, network-selection, and
plotting utilities**. They are not part of the NetCS inference algorithm
itself.

They are used to:

- identify networks containing degree-4 blobs,
- compute hybrid-ancestry errors between the true network and the
  estimated network,
- calculate hardwired cluster distance,
- compare hybrid placement within blobs between the true and estimated
  networks,
- compare circular orderings using Kendall-tau distance,
- compare the true tree of blobs (TOB) with an estimated TOB,
- select the CAMUS output network used for evaluation, and
- generate the draft figures from the processed experimental results.

The directory contains:

```text
scripts/
├── check_four_blob.jl
├── compare_hybrid_vec.jl
├── compare_network.jl
├── compare_two_networks_by_hybrid_and_circle_ordering.jl
├── compare_two_tree.py
├── select_network_by_hybrid_num.jl
├── figure2-plot.py
├── figure3-figure4-plot.py
├── figure5-plot.py
├── figure-6-plot.py
└── README.md