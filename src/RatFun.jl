import ApproxFun:evaluate,dimension
import Plots:plot,plot!

export RatFun, reciprocal

immutable RatFun{S1,T1,S2,T2}
    p::Fun{S1,T1}
    q::Fun{S2,T2}
end

function evaluate(r::RatFun,x)
    (r.p)(x)./(r.q)(x)
end

for op = (:*,:.*)
    @eval $op(r1::RatFun,r2::RatFun)=RatFun($op(r1.p,r2.p),$op(r1.q,r2.q))
    @eval $op(r::RatFun,a) = RatFun($op(r.p,a),r.q)
    @eval $op(a,r::RatFun) = RatFun($op(a,r.p),r.q)
end

reciprocal(r::RatFun) = RatFun(r.q,r.p)

(./)(r1::RatFun,r2::RatFun)=r1.*reciprocal(r2)
(./)(a,r::RatFun)=a.*reciprocal(r)
(./)(r::RatFun,a)=reciprocal(r).*a

(/)(r1::RatFun,r2::RatFun)=r1*reciprocal(r2)
(/)(a,r::RatFun)=a*reciprocal(r)
(./)(r::RatFun,a)=reciprocal(r)*a

for op = (:+,:.+,:-,:.-)
  @eval $op(r1::RatFun,r2::RatFun) = RatFun($op((r1.p.*r2.q),(r2.p.*r1.q)),r1.q.*r2.q)
end

plot(r::RatFun;grid=true,kwds...)=plot!(plot(grid=grid),r;kwds...)
plot!(f::RatFun;kwds...)=plot!(current(),r;kwds...)

plot(x::AbstractVector,r::RatFun;grid=true,kwds...)=plot!(plot(grid=grid),x,r;kwds...)
plot!(x::AbstractVector,f::RatFun;kwds...)=plot!(current(),x,r;kwds...)


# The padding in this function can be improved
# No support for functions with poles within the domain
function plotptsvals(r::RatFun)
    p = r.p
    q = r.q
    plen = length(p)
    qlen = length(q)
    if dimension(space(p)) == Inf && dimension(space(q)) == Inf
        p=pad(p,3plen+10*qlen+50)
        q=pad(q,3plen+10*qlen+50)
        r = RatFun(p,q)
    elseif dimension(space(p)) == Inf && dimension(space(q)) < Inf
        p=pad(p,3plen+10*dimension(space(q))+50)
        q=pad(q,3plen+10*dimension(space(q))+50)
        r = RatFun(p,q)
    elseif dimension(space(q)) == Inf && dimension(space(p)) < Inf
        q=pad(q,dimension(space(p))+10qlen+50)
        p=pad(p,dimension(space(p))+10qlen+50)
        r = RatFun(p,q)
    else
        p=pad(p,dimension(space(p))+dimension(space(q)))
        r = RatFun(p,q)
    end
    return points(r.p),values(r.p)./values(r.q)
end

plot!{S,T<:Real}(plt::Plots.Plot,r::RatFun{S,T};kwds...)=
                plot!(plt,plotptsvals(r)...;kwds...)

plot!{S,T<:Real}(plt::Plots.Plot,x::AbstractVector{T},r::RatFun{S,T};kwds...)=
                plot!(plt,x,evaluate(r,x);kwds...)


function plot!{S1<:PiecewiseSpace,S2<:PiecewiseSpace,T1<:Real,T2<:Real}(plt::Plots.Plot,r::RatFun{S1,T1,S2,T2},kwds...)
    vp = vec(r.p)
    vq = vec(r.q)
    c=plt.plotargs[:color_palette][plt.n+1]
    for k=1:length(vp)
        plot!(plt,RatFun(vp[k],vq[k]),color=c,kwds...)
    end
    plt
end

# For dirac space, we draw a dotted line extending to infinity
function plot!{S1<:DiracSpace,T1<:Real,S2<:PointSpace,T2<:Real}(plt::Plots.Plot,r::RatFun{S1,T1,S2,T2};kwds...)
    p = r.p
    q = r.q
    pts=space(p).points
    n=length(pts)
    ws=pad(p.coefficients./q.coefficients,length(pts))
    plt=plot!(plt,ones(2)*pts[1],[0,1]*ws[1];kwds...)
    c=plt.plotargs[:color_palette][plt.n]
    plot!(plt,ones(2)*pts[2:end]',[0,1]*ws[2:end]';color=c,kwds...)
    plot!(plt,ones(2)*pts',[1,2]*ws';color=c,linestyle=:dot,kwds...)
end

# for PointSpace, we draw just a line
function plot!{S1<:PointSpace,T1<:Real,S2<:PointSpace,T2<:Real}(plt::Plots.Plot,r::RatFun{S1,T1,S2,T2};kwds...)
    p = r.p
    q = r.q
    pts=space(p).points
    n=length(pts)
    ws=pad(p.coefficients./q.coefficients,length(pts))
    plt=plot!(plt,ones(2)*pts[1],[0,1]*ws[1];kwds...)
    c=plt.plotargs[:color_palette][plt.n]
    plot!(plt,ones(2)*pts[2:end]',[0,1]*ws[2:end]';color=c,kwds...)
end
