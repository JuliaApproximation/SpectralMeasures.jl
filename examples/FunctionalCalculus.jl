using ApproxFun, SpectralMeasures, Plots
gr(linewidth=3,legend=true,xlims=(0,100))

D = freejacobioperator() - I

## Discrete Diffusion equations
u0 = [zeros(20);ones(5);zeros(20);1.2*ones(3)]

# J1 = D
J1 = D
x1,U1 = eig(J1)
plot(pad(u0,100),ylims=(0,1.5),marker=2,label="t=0",xlabel="k",ylabel="u_k",title="Pure diffusion on the half-line")
  for t in [1,10,100,1000]
    ut1 = U1\(exp(x1*t)*(U1*u0))
    plot!(pad(ut1.coefficients,100),ylims=(0,1.5),marker=2,label="t=$(t)")
  end
  plot!()
Plots.savefig("purediffusion1.pdf")

J1n = full(J1[1:70,1:70])
plot(pad(u0,100),ylims=(0,1.5),marker=2,label="t=0",xlabel="k",ylabel="u_k",title="Pure diffusion on the half-line (n=70 finite section method)")
  for t in [1,10,100,1000]
    ut1 = expm(J1n*t)*pad(u0,70)
    push!(areasfs,sum(ut1))
    push!(energiesfs,sum(ut1.^2))
    plot!(pad(ut1,100),ylims=(0,1.5),marker=2,label="t=$(t)")
  end
  plot!()
Plots.savefig("purediffusion1finitesection.pdf")


# J2 = D + diag(v) where v is a reaction term
react = [zeros(35);-1.;-1.;-1]
J2 = D + SymTriOperator(react,[0.])
x2,U2 = eig(J2)
plot(pad(u0,100),marker=2,label="t=0",ylims=(0,1.5),xlabel="k",ylabel="u_k",title="Diffusion plus reaction on the half-line")
  for t in [1,10,100,1000]
    ut2 = U2\(exp(x2*t)*(U2*u0))
    plot!(pad(ut2.coefficients,100),ylims=(0,1.5),marker=2,label="t=$(t)")
  end
  plot!()
plot!(pad(-react,100),label="Reaction potential",color=:black,ylims=(0,1.5),linewidth=2)
Plots.savefig("reactiondiffusion.pdf")
eigs = discrete_eigs(J2)
plot(U2[1,1:100],marker=2,title="Stable eigenmodes of reaction-diffusion Jacobi operator",xlabel="k",ylabel="u_k",label="v_1")
  plot!(U2[2,1:100],marker=2,label="v_2")
  plot!(pad(react,100),color=:black,label="Reaction potential",yticks=-1.5:0.5:1.5)
Plots.savefig("stableeigenmodes.pdf")

# Pure fractional diffusion α = 0.65
J4 = D
x4,U4 = eig(J4)
α = 0.65
plot(pad(u0,100),ylims=(0,1.5),marker=2,label="t=0",xlabel="k",ylabel="u_k",title="Pure fractional diffusion on the half-line, \\alpha = 0.65")
  for t in [1,10,100,1000]
    h4 = Fun(x->exp(-(abs(x)^α)*t),space(x4))
    ut4 = U4\(h4*(U4*u0))
    plot!(pad(ut4.coefficients,100),ylims=(0,1.5),marker=2,label="t=$(t)")
  end
  plot!()
Plots.savefig("purefracdiffusion1.pdf")

# Pure fractional diffusion α = 0.85
J5 = D
x5,U5 = eig(J5)
α = 0.85
plot(pad(u0,100),ylims=(0,1.5),marker=2,label="t=0",xlabel="k",ylabel="u_k",title="Pure fractional diffusion on the half-line, \\alpha = 0.85")
  for t in [1,10,100,1000]
    h5 = Fun(x->exp(-(abs(x)^α)*t),space(x5))
    ut5 = U5\(h5*(U5*u0))
    plot!(pad(ut5.coefficients,100),ylims=(0,1.5),marker=2,label="t=$(t)")
  end
  plot!()
Plots.savefig("purefracdiffusion085.pdf")

# Pure fractional diffusion α = .25
J6 = D
x6,U6 = eig(J6)
α = .25
plot(pad(u0,100),ylims=(0,1.5),marker=2,label="t=0",xlabel="k",ylabel="u_k",title="Pure fractional diffusion on the half-line, \\alpha = 0.25")
  for t in [1,10,100,1000]
    h6 = Fun(x->exp(-(abs(x)^α)*t),space(x6))
    ut6 = U6\(h6*(U6*u0))
    plot!(pad(ut6.coefficients,100),ylims=(0,1.5),marker=2,label="t=$(t)")
  end
  plot!()
Plots.savefig("purefracdiffusion025.pdf")

# Fractional diffusion with reaction α = 0.85
react = [zeros(20);1.;1.;1.]
J7 = D + SymTriOperator(react,[0.])
discrete_eigs(J7)
x7,U7 = eig(J7)
α = 0.85
plot(pad(u0,100),ylims=(0,1.5),marker=2,label="t=0",xlabel="k",ylabel="u_k",title="Fractional diffusion with reaction on the half-line, \\alpha = 0.85")
#for t in [1;10;100;1000]
t = 1
h7 = Fun(x->exp(-(abs(x)^α)*t),space(x7)[3])
pt1 = space(x7)[1].points[1]
pt2 = space(x7)[2].points[1]
h7 = Fun(space(x7),[exp((pt1^α)*t);exp((pt2^α)*t);coefficients(h7)])

h7 = Fun(x->exp(-(abs(x)^α)*t),space(x7))

ut7 = U7\(h7*(U7*u0))
plot!(pad(ut7.coefficients,100),ylims=(0,1.5),marker=2,label="t=$(t)")
#end
plot!()
plot!(pad(react,100),color=:black,label="Reaction potential")
Plots.savefig("fracreactdiffusion085withpot.pdf")




SymTriPertToeplitz(V8+1.0,[.5],1.0,-.5)


# Discrete Schrodinger equations
V8 = [zeros(40);2*ones(5)]
  J8 = -D + SymTriOperator(V8,[0.])
  x8,U8 = eig(J8)
  λ8 = discrete_eigs(J8)
  n8 = length(λ8)
  Eigv8 = full(U8[1:n8,1:100])

Eigv8[1,40:100]

n8

V9 = [zeros(40);1*ones(3)]
  J9 = D + SymTriOperator(V9,[0.])
  x9,U9 = eig(J9)
  λ9 = discrete_eigs(J9)
  n9 = length(λ9)
  Eigv9 = full(U9[1:n9,1:100])

if n8 > 0
  plot(Eigv8[1,:]/4+λ8[1]-2,marker=2,title="Tunneling barrier potential and eigenstates",xlabel="k",ylabel="u_k",label="v_1")
  if n8 > 1
    for k in 2:n8
      plot!(Eigv8[k,:]/4+λ8[k]-2,marker=2,label="v_$(k)")
    end
  end
  plot!(pad(V8,100),color=:black,label="Potential")
end
Plots.savefig("tunnelbarriereig.pdf")

if n9 > 0
  plot(Eigv9[1,:]/4+λ9[1],marker=2,title="Double well potential and eigenstates (far)",xlabel="k",ylabel="u_k",label="v_1",ylims=(-.075,.4))
  if n9 > 1
    for k in 2:n9
      plot!((-1)^(k+1)*Eigv9[k,:]/4+λ9[k],marker=2,label="v_$(k)")
    end
  end
  plot!(.2-pad(V9,100),color=:black,label="Potential")
end
Plots.savefig("doublewelleigfar.pdf")

V8 = [zeros(40);2*ones(5)]
  J8 = -D + SymTriOperator(V8,[0.])
  n = 400
  uinit = pad!(exp(-(-28:.75:28).^2),n)
  Jn = full(J8[1:n,1:n])
t = 20
  ut8 = expm(im*Jn*t)*uinit
  plot(abs(ut8[1:100]),ylims=(-.5,1.5),marker=2,label="t=$(t) (abs)")
  plot!(pad(V8,100),ylims=(-.5,1.5),color=:black,label="Potential",xlabel="k",ylabel="u_k",title="No quantum tunneling on the half-line")

d,Q = eig(Jn)





for t in [0;5;10;20;30;40]
  ut8 = expm(im*Jn*t)*uinit
  plot(abs(ut8[1:100]),ylims=(-.5,1.5),marker=2,label="t=$(t) (abs)")
  plot!(pad(V8,100),ylims=(-.5,1.5),color=:black,label="Potential",xlabel="k",ylabel="u_k",title="No quantum tunneling on the half-line")
  Plots.savefig("noquantumtunnelingt=$(t).pdf")
end
plot!()


V8 = [zeros(10);-.5*ones(6);zeros(10);-.5*ones(6)]
  J8 = D + SymTriOperator(V8,[0.])
  n = 400
  uinit = pad!([zeros(15);.45],n)
  Jn = full(J8[1:n,1:n])
  t = 50
  ut8 = expm(im*Jn*t)*uinit
  plot(abs(ut8[1:100]),ylims=(-.5,1),marker=2,label="t=$(t)")
  plot!(pad(V8,100),ylims=(-.5,1),color=:black,label="Potential",xlabel="k",ylabel="u_k",title="Quantum tunneling on the half-line")
Plots.savefig("doublewellsimclose.pdf")


n = 1000
uinit = pad!([zeros(11);.05;.09;.09;.05],n)
  Jn = full(J9[1:n,1:n])
  plot(pad(V9,100),ylims=(0,1.5),marker=2,label="t=0",xlabel="k",ylabel="u_k",title="Quantum tunneling on the half-line")
  t = 1000000
  ut9 = expm(im*Jn*t)*uinit
  plot!(abs(ut9[1:100]),ylims=(0,1.5),marker=2,label="t=$(t)")
Plots.savefig("doublewellsimfar.pdf")


u0 = [zeros(10);1.;1.;1.;1.]
plot(pad(u0,100),marker=2)
for t = [5;10;15;20;25]
  ut = SchroEvo(u0U\(exp(im*x*t)*(U*u0))
  plot!(pad(abs(ut.coefficients),100),marker=2)
end
plot!()
