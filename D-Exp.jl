
using JuMP
using Ipopt, MosekTools, Gurobi
using LinearAlgebra
using Random 
using TickTock

#Coefficients
Random.seed!(20);
n = 30 #dimension of variable x
Q_p = rand(-50:50,n,n)
Q = 0.5*(Q_p+Q_p')
c = rand(-50:50,n)
l = zeros(n)
u = ones(n)
#hyperparameter in SOCP
limit = 30 
theta = -1e-8 #original -1e-8
#time_limit = 5 #time limitation for SOCP
time_rate = 1 #set time limitation be the time_rate% of the consumed time of DNP 

#Box-constrained Quadratic Programs
#Calculate optimial value or upper bound (feasible solution)
function QP(Q, c, l, u)
    QP = Model(Ipopt.Optimizer)
    set_silent(QP) #unset_silent(QP)
    @variable(QP, x[1:n])
    @objective(QP, Min, x'*Q*x + c'*x)
    @constraint(QP, x .>= l) #dot comparision
    @constraint(QP, x .<= u) #dot comparision
    #print(QP)
    optimize!(QP)
    #println("objective value: ", objective_value(QP))
    #println("x: ", value.(x))
    #println("termination_status:", termination_status(QP))
    return objective_value(QP), termination_status(QP), primal_status(QP)#####################
end

#Semidefinite Programming
#Doubly nonnegative relaxation
function DNP(Q, c, l, u)
    #SemiDP = Model(with_optimizer(Mosek.Optimizer, MSK_DPAR_OPTIMIZER_MAX_TIME = 100.0))
    SemiDP = Model(Mosek.Optimizer)
    set_silent(SemiDP) #unset_silent(SemiDP)
    @variable(SemiDP, X[1:n,1:n], PSD)##########
    @variable(SemiDP, x[1:n])
    @objective(SemiDP, Min, dot(Q,X)+c'*x)
    @constraint(SemiDP, -x*u'-(x*u')'+X .>= -1) #matrix
    @constraint(SemiDP, x*u'-X .>= 0) #matrix
    @constraint(SemiDP, X.>=0) #matrix
    @SDconstraint(SemiDP, [1 x'; x X]>=0)
    #print(SemiDP)
    optimize!(SemiDP)
    #println("objective value:", objective_value(SemiDP))
    #println("x: ", value.(x))
    #println('X', value.(X))
    #println("termination_status:", termination_status(SemiDP))
    #println(primal_status(SemiDP))
    #solution_summary(SemiDP)
    return objective_value(SemiDP), termination_status(SemiDP), primal_status(SemiDP)#####################
end
#@time DNP(Q, c, l, u)

#Second-Order Cone Programming
#use SOCP constraints to approximate SDP constraints
#function SOCP(Q, c, l, u, limit, theta, time_limit)
function SOCP(Q, c, l, u, limit, theta, time_rate)#set time limitation be the time_rate% of the consumed time of DNP 
    time_limit = R_DNP[2]*time_rate
    tick()#timer start
    model = Model(Gurobi.Optimizer)
    set_silent(model) #unset_silent(model)
    set_time_limit_sec(model, 10.0) #unset_time_limit_sec(model)
    @variable(model, X[1:n,1:n],Symmetric) #symmetric########################
    @variable(model, x[1:n])
    @objective(model, Min, dot(Q,X)+c'*x)
    @constraint(model, -x*u'-(x*u')'+X .>= -1) #matrix
    @constraint(model, x*u'-X .>= 0) #matrix
    @constraint(model, X.>=0) #matrix
    optimize!(model)
    x_val = value.(x)
    X_val = value.(X)
    X_val = (X_val+X_val')/2###Symmetric#####
    min_eigval = eigvals(X_val-x_val*x_val')[1]
    i_count = 1
    while min_eigval <= theta && i_count <= limit && peektimer() <= time_limit#check the time
        #Approximate the semidefinite constraints
        v = eigvals(X_val-x_val*x_val')
        V = eigvecs(X_val-x_val*x_val')
        V_index = findall(v .<= theta)#get the index of negative eigenvalues
        for i in V_index
            d = V[:,i]#eigenvectors
            @constraint(model, [(1+d'*X*d)/2; (1-d'*X*d)/2; d'*x] in SecondOrderCone())
        end
        optimize!(model)
        value.(x)
        value.(X)
        x_val = value.(x)
        X_val = value.(X)
        X_val = (X_val+X_val')/2###Symmetric#####
        min_eigval = eigvals(X_val-x_val*x_val')[1]  
        i_count = i_count+1
    end
    tok()#timer stop
    #println("termination_status:", termination_status(model))
    #println("iteration times (i_count-1):", i_count-1)
    #print(model)
    #println("objective value:", objective_value(model))
    #println("x: ", value.(x))
    #println('X', value.(X))
    #println("minimum(eigvals(X_val-x_val*x_val'))=", eigvals(X_val-x_val*x_val')[1])
    return objective_value(model), termination_status(model), primal_status(model), i_count-1, eigvals(X_val-x_val*x_val')[1]#####################
end 

R_QP = @timed QP(Q, c, l, u)
R_DNP = @timed DNP(Q, c, l, u)
R_SOCP = @timed SOCP(Q, c, l, u, limit, theta, time_rate)
println("Time_QP:", R_QP[2])
println("Time_DNP:", R_DNP[2])
println("Time_SOCP:", R_SOCP[2])
println("objective_value(QP):", R_QP[1][1], "  termination_status(QP):", R_QP[1][2], "  primal_status(QP):", R_QP[1][3])
println("objective_value(DNP):", R_DNP[1][1], "  termination_status(DNP):", R_DNP[1][2], "  primal_status(DNP):", R_DNP[1][3])
println("objective_value(SOCP):", R_SOCP[1][1], "  termination_status(SOCP):", R_SOCP[1][2], "  primal_status(SOCP):", R_SOCP[1][3],
    "  iteration times (i_count-1):", R_SOCP[1][4], "  minimum(eigvals(X_val-x_val*x_val'))=",R_SOCP[1][5])
