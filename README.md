**SS-PPCA:** Stochastic Subspace via Probabilistic Principal Component Analysis 
> A methodology for efficient probabilistic characterization of model-form error via stochastic reduced-order modeling.

Model-form error is ubiquitous in computational science and engineering. Its probabilistic characterization is a critical and challenging problem.

SS-PPCA introduces probabilistic modeling of subspaces. The stochastic ROMs constructed using the stochastic subspaces are then used to characterize model error.

The code demonstrates the performance of SS-PPCA for characterization of model error.

## Repository Contents

- `Ex1/` — Parametric linear static problem.
- `Ex2/` — Linear static problem for characterizing high-dimensional model (HDM) error.

## Requirements

- MATLAB R2023b or later (older versions might also work, but are not tested).
- Required Toolboxes:
  - Statistics and Machine Learning Toolbox
  - Optimization Toolbox

## Running the Experiments

1. Clone the repository:

    ```bash
    git clone https://github.com/UQUH/SS_PPCA.git
    cd SS_PPCA
    ```

2. Open MATLAB and set the repository folder as the working directory.

3. Run the scripts for different examples.

4. Results and figures will be saved automatically in the `results/` folder.

Each example (`Ex1`, `Ex2`) contains a separate `model` folder and method folder `SS-PPCA`:
- **SS-PPCA**: Stochastic Subspace via Probabilistic Principal Component Analysis.

## Citation

If you found SS-PPCA helpful, please cite the [following
paper](https://arxiv.org/abs/2504.19963):
```
@Article{yadav2025ss,
  author  = {Akash Yadav and Ruda Zhang},
  journal = {Computational Mechanics},
  title   = {Stochastic Subspace via Probabilistic Principal Component Analysis for Characterizing Model Error},
  year    = {2025},
  doi     = {10.1007/s00466-025-02701-6},
}
```

## The Team

The SS-PPCA method is developed by the [Uncertainty Quantification Lab](https://uq.uh.edu/group-members) at the University of Houston.

The primary contributors are:
- Akash Yadav
- Ruda Zhang
