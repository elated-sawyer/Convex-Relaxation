using JuMP
using Ipopt, MosekTools, Gurobi, CPLEX
using LinearAlgebra
using Random 
using TickTock
using CSV, DataFrames
using DelimitedFiles

#Box-constrained Quadratic Programs
#Calculate optimial value or upper bound (feasible solution)
function QP(n, Q, c, l, u)
    QP = Model(Ipopt.Optimizer)
    set_silent(QP) #unset_silent(QP)
    @variable(QP, x[1:n])
    @objective(QP, Max, 0.5*x'*Q*x + c'*x)#Example-2
    @constraint(QP, x .>= l) #dot comparision
    @constraint(QP, x .<= u) #dot comparision
    optimize!(QP)
    return objective_value(QP), termination_status(QP), primal_status(QP)
end

#RLT constraints
function RLT(n, Q, c, l, u)
    RLT_constraints = Model(Mosek.Optimizer)
    unset_silent(RLT_constraints)
    @variable(RLT_constraints, X[1:n,1:n], Symmetric)
    @variable(RLT_constraints, x[1:n])
    @objective(RLT_constraints, Max, 0.5*dot(Q,X)+c'*x)
    @constraint(RLT_constraints, -x*u'-(x*u')'+X .>= -1)
    @constraint(RLT_constraints, x*u'-X .>= 0)
    @constraint(RLT_constraints, X.>=0)
    optimize!(RLT_constraints)
    return objective_value(RLT_constraints), termination_status(RLT_constraints), primal_status(RLT_constraints)
end

#Semidefinite Programming
#Doubly nonnegative relaxation
function DNP(n, Q, c, l, u)
    SemiDP = Model(Mosek.Optimizer)
    set_silent(SemiDP) #unset_silent(SemiDP)
    @variable(SemiDP, X[1:n,1:n], PSD)
    @variable(SemiDP, x[1:n])
    @objective(SemiDP, Max, 0.5*dot(Q,X)+c'*x)
    @constraint(SemiDP, -x*u'-(x*u')'+X .>= -1)
    @constraint(SemiDP, x*u'-X .>= 0)
    @constraint(SemiDP, X.>=0)
    @SDconstraint(SemiDP, [1 x'; x X]>=0)
    optimize!(SemiDP)
    return objective_value(SemiDP), termination_status(SemiDP), primal_status(SemiDP)
end

#Doubly nonnegative-Second-Order Cone Programming
function SOCP(n, Q, c, l, u, limit, theta, time_limit, file)
    tick()  #timer start
    model = Model(CPLEX.Optimizer)
    unset_silent(model)
    @variable(model, X[1:n,1:n],Symmetric)
    @variable(model, x[1:n])
    @objective(model, Max, 0.5*dot(Q,X)+c'*x)
    @constraint(model, -x*u'-(x*u')'+X .>= -1)
    @constraint(model, x*u'-X .>= 0)
    @constraint(model, X.>=0)
    @constraint(model, [0; 0] in SecondOrderCone())
    optimize!(model)
    x_val = value.(x)
    X_val = value.(X)
    X_val = (X_val+X_val')/2  #Symmetric
    Add_Con = 0
    CSV.write("EXP_SOCP_detail.csv", [(File_name = file, Time_SOCP = solve_time(model), objective_value_SOCP = objective_value(model), Adding_Constraints = Add_Con)],  append = true)
    v = eigvals(X_val-x_val*x_val')
    min_eigval = minimum(v)
    i_count = 1
    while min_eigval <= theta && i_count <= limit && peektimer() <= time_limit#check the time
        #Approximate the semidefinite constraints
        V = eigvecs(X_val-x_val*x_val')
        V_index = findall(v .<= theta) #get the index of negative eigenvalues
        for i in V_index
            d = V[:,i]  #eigenvectors
            @constraint(model, [(1+d'*X*d)/2; (1-d'*X*d)/2; d'*x] in SecondOrderCone())
        end
        optimize!(model)
        value.(x)
        value.(X)
        x_val = value.(x)
        X_val = value.(X)
        X_val = (X_val+X_val')/2  #Symmetric
        Add_Con = length(V_index)
        CSV.write("EXP_SOCP_detail.csv", [(File_name = file, Time_SOCP = solve_time(model), objective_value_SOCP = objective_value(model), Adding_Constraints = Add_Con)],  append = true)
        v = eigvals(X_val-x_val*x_val')
        min_eigval = minimum(v)
        i_count = i_count+1
    end
    tok()#timer stop
    return objective_value(model), termination_status(model), primal_status(model), i_count-1, min_eigval
end 

#Semidefinite Programming with time limitation
#Doubly nonnegative relaxation
function DNP_TimLim(n, Q, c, l, u)
    SemiDP = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_time_limit_sec(SemiDP, 360.0)
    @variable(SemiDP, X[1:n,1:n], PSD)
    @variable(SemiDP, x[1:n])
    @objective(SemiDP, Max, 0.5*dot(Q,X)+c'*x)
    @constraint(SemiDP, -x*u'-(x*u')'+X .>= -1)
    @constraint(SemiDP, x*u'-X .>= 0)
    @constraint(SemiDP, X.>=0)
    @SDconstraint(SemiDP, [1 x'; x X]>=0)
    optimize!(SemiDP)
    return objective_value(SemiDP), termination_status(SemiDP), primal_status(SemiDP)
end