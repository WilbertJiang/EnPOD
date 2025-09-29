# EnPOD
EnPOD (Ensemble Proper Orthogonal Decomposition with Galerkin Projection, EnPOD-GP) is a thermal simulation methodology developed for multi-core CPUs that balances accuracy, efficiency, and scalability. Unlike the earlier GPOD-GP approach, which generates a single global set of POD modes for the entire chip and thus requires extensive training to capture variations in all power sources, EnPOD constructs multiple independent POD models—one for each functional unit (FU). Each individual POD model is trained separately using temperature data induced by the FU’s own power source, and the full-chip thermal profile is then obtained through superposition. This ensemble strategy eliminates the need for global training power maps, significantly reduces training effort, and improves robustness and accuracy, particularly for processors with many cores. Compared to finite element simulations, EnPOD delivers orders-of-magnitude speedup while maintaining fine spatial and temporal resolution, making it well-suited for large-scale dynamic thermal analysis and practical runtime applications such as thermal-aware task scheduling and reliability management.
<p align="center">
  <img src="Image/enPOD.jpg" alt="Workflow of EnPOD" width="600">
</pr>
```
[1] Jiang L, Dowling A, Liu Y, Cheng M-C. Ensemble learning model for effective thermal simulation of multi-core CPUs. Integration. 2024;97:102201.
```
