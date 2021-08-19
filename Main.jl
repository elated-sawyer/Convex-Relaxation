using JuMP
using Ipopt, MosekTools, Gurobi, CPLEX
using LinearAlgebra
using Random 
using TickTock
using CSV, DataFrames
using DelimitedFiles
include("Relaxation Models.jl")

#experiments on original problem
#solving the problem by IPOPT directly
for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/basic/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended2/")
    for file in files
        dat = readdlm(joinpath.(root, file))
        n = dat[1][1]
        c = Vector{Float64}(dat[2,:])
        Q = Matrix{Float64}(dat[3: size(dat)[1], :])
        #other parameters
        l = zeros(n)
        u = ones(n)
        R_QP = @timed QP(n, Q, c, l, u)
        CSV.write("EXP_Original_problem.csv", [(File_name = file, Time_QP = R_QP[2], objective_value_QP = R_QP[1][1],  termination_status_QP = R_QP[1][2])],  append = true)
    end
end

#experiments on RLT constraints problem
#solving the problem by Mosek directly
for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/basic/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended2/")
    for file in files
        dat = readdlm(joinpath.(root, file))
        n = dat[1][1]
        c = Vector{Float64}(dat[2,:])
        Q = Matrix{Float64}(dat[3: size(dat)[1], :])
        #other parameters
        l = zeros(n)
        u = ones(n)
        R_RLT = @timed RLT(n, Q, c, l, u)
        CSV.write("EXP_RLT_problem.csv", [(File_name = file, Time_RLT = R_RLT[2], objective_value_RLT = R_RLT[1][1],  termination_status_RLT = R_RLT[1][2])],  append = true)
    end
end

#experiments on SDP problem
#solving the problem by Mosek directly
for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/basic/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended2/")
    for file in files
        dat = readdlm(joinpath.(root, file))
        n = dat[1][1]
        c = Vector{Float64}(dat[2,:])
        Q = Matrix{Float64}(dat[3: size(dat)[1], :])
        #other parameters
        l = zeros(n)
        u = ones(n)
        R_DNP = @timed DNP(n, Q, c, l, u)
        CSV.write("EXP_DNP_problem.csv", [(File_name = file, Time_DNP = R_DNP[2], objective_value_DNP = R_DNP[1][1],  termination_status_DNP = R_DNP[1][2])],  append = true)
    end
end

#experiments on DNN-SOCP problem
#solving the problem by Mosek
for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/basic/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended2/")
    for file in files
        dat = readdlm(joinpath.(root, file))
        n = dat[1][1]
        c = Vector{Float64}(dat[2,:])
        Q = Matrix{Float64}(dat[3: size(dat)[1], :])
        #other parameters
        l = zeros(n)
        u = ones(n)
        #hyperparameter in SOCP
        limit = 500
        theta = -1e-8 #original -1e-8
        time_limit = 360 #time limitation for SOCP
        R_SOCP = @timed SOCP(n, Q, c, l, u, limit, theta, time_limit, file)
        CSV.write("EXP_SOCP_problem.csv", [(File_name = file, Time_SOCP = R_SOCP[2], objective_value_SOCP = R_SOCP[1][1], termination_status_SOCP = R_SOCP[1][2], iteration_times = R_SOCP[1][4],  min_eigval = R_SOCP[1][5])],  append = true)
    end
end

#experiments on SDP problem with time limitation
#solving the problem by Mosek directly
for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/basic/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended/")
#for (root, dirs, files) in walkdir("./Experiments/BoxQP_instances-master/extended2/")

    for file in files
        dat = readdlm(joinpath.(root, file))
        n = dat[1][1]
        c = Vector{Float64}(dat[2,:])
        Q = Matrix{Float64}(dat[3: size(dat)[1], :])
        #other parameters
        l = zeros(n)
        u = ones(n)
        R_DNP = @timed DNP_TimLim(n, Q, c, l, u)
        CSV.write("EXP_DNPTimLim_problem.csv", [(File_name = file, Time_DNP = R_DNP[2], objective_value_DNP = R_DNP[1][1],  termination_status_DNP = R_DNP[1][2], primal_status_DNP = R_DNP[1][3])],  append = true)
    end
end
