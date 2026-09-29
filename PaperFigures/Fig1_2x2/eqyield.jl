# Compute equilibrium yield for 2x2 square, with given 
# volume fraction dens and given strong&weak bonds ds,dw
# weak bond is 34

using NonlinearSolve

function f4(u,p)
    c1 = u[1]   # concentration of 1,2
    c3 = u[2]   # concentration of 3,4
    ds = p[1]   # strong bond strength
    dw = p[2]   # weak bond strength
    rho = p[3]  # monomer concentration (of each individual monomer)
    es = exp(ds)
    ew = exp(dw)

    #println("es=",es,",ew=",ew,",rho=",rho,",c1=",c1,",c3=",c3)

    mass1 = c1 + es*c1^2 + es*c1*c3 +          # monomer+dimers
             2*es^2*c1^2*c3 + es*ew*c3^2*c1 +  # trimers
             es^3*ew*c1^2*c3^2 - rho;          # target  
    mass3 = c3 + ew*c3^2 + es*c1*c3 +          # monomer+dimers
            es^2*c1^2*c3 + 2*es*ew*c3^2*c1 +   # trimers
            es^3*ew*c1^2*c3^2 - rho;           # monomer+dimers
    #println("mass1=",mass1,",mass3=",mass3)
    return [mass1;mass3]
end

function f1(u,p)
    ed = exp(p[1])
    c1 = u[1]
    return c1 + 2*c1^2*ed + 3*c1^3*ed^2 + c1^4*ed^4 - p[2];
end

function ceq(ds,dw,rho=0.05)
    p = [ds,dw,rho];
    u0 = rho*[1,1];
    prob = NonlinearProblem(f4, u0, p)
    sol = solve(prob)
    yield = (exp(3*ds+dw)*sol.u[1]^2*sol.u[2]^2)/rho
    return yield;
end



#=
######  TESTING  ######
# test 1-component
rho = 0.05/4
del = 8
u1 = rho/10;
p1 = [del,rho]
prob1 = NonlinearProblem(f1, u1, p1)
sol1 = solve(prob1)
println(sol1.u)
y1=sol1.u^4*exp(4*del)/rho
println("yield (1-component)=",y1)

# multi-component
u0 = sol1.u*[1,1]
ds = 8;
dw = 4;
p = [ds,dw,rho];
prob = NonlinearProblem(f4, u0, p)
sol = solve(prob)
println(sol)
yield = (exp(3*ds+dw)*sol.u[1]^2*sol.u[2]^2)/rho
println("yield (4 components) = ",yield)
=#
