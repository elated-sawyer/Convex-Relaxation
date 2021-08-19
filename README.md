# Convex-Relaxation

**Main.jl**: Call the relaxation models in **Relaxation Models.jl** to do experiments and save the outputs.

**Relaxation Models.jl**: Box-constrained Quadratic Programming Model, RLT Relaxatoin Model, DNN Relaxatoin Model, DNN-SOCP Relaxation Model, DNN Relaxatoin with Time Limitation Model are included.

**Data Files**: "Records for Iterations.csv", "Results for Solvers" obtained by the output of **Main.jl**.

**Plots.jl**: Plot figures based on data files, "Records for Iterations.csv", "Results for Solvers".

**Experimental Instances** are from Prof. Burer's repository--*sburer/BoxQP_instances*(link: https://github.com/sburer/BoxQP_instances)
