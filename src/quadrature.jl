"""
```matlab
function [x, w] = GaussLegendre_2(n)
% As you can see, the code consists of 2 blocks:
% 1: construct a symmetrical companion matrix
% 2: determine the (real) eigenvalues (i.e. the roots of the polynomial).
% It can produce the correct abscissas and weights, for any value n>=2.
%
% input: n ---- number of intergrate points
%
% output: x --- GaussLegendre intergrate point
%         w --- GaussLegendre intergrate weight
%
i = 1:n-1;
a = i./sqrt(4*i.^2-1);
CM = diag(a,1) + diag(a,-1);
[V L] = eig(CM);
[x ind] = sort(diag(L));
V = V(:,ind)';
w = 2 * V(:,1).^2;
return
end
```
------from [pasuka](http://bbs.fcode.cn/thread-97-1-1.html)
"""
function gauss_legendre(x::Int)
    x > 0 ? Nothing : throw(ArgumentError(" $x 个积分点？"))
    i = collect(1:x-1)  # 如果x为1时，i为Int64[]
    a = i ./ sqrt.(4 * i.^2 .- 1)
    CM = zeros(Float64, x, x)
    for i in 1:x-1
        CM[i, i+1] = CM[i+1, i] = a[i]
    end
    L, V = eigen(CM)
    w = 2 * V[:, 1].^2
    [L w]
end
function gauss_legendre(x::Int, y::Int)
    g1 = gauss_legendre(x)
    g2 = gauss_legendre(y)
    gx1, w1 = g1[:, 1], g1[:, 2]
    gx2, w2 = g2[:, 1], g2[:, 2]

    x_points = [gx1[i] for i in 1:x for j in 1:y]
    y_points = [gx2[j] for i in 1:x for j in 1:y]
    weights = [w1[i] * w2[j] for i in 1:x for j in 1:y]

    coords = vcat(x_points', y_points')
    return(coords, weights)
end
function gauss_legendre(x::Int, y::Int, z::Int)
    g1 = gauss_legendre(x)
    g2 = gauss_legendre(y)
    g3 = gauss_legendre(z)
    gx1, w1 = g1[:, 1], g1[:, 2]
    gx2, w2 = g2[:, 1], g2[:, 2]
    gx3, w3 = g3[:, 1], g3[:, 2]

    x_points = [gx1[i] for i in 1:x for j in 1:y for k in 1:z]
    y_points = [gx2[j] for i in 1:x for j in 1:y for k in 1:z]
    z_points = [gx3[k] for i in 1:x for j in 1:y for k in 1:z]
    weights = [w1[i] * w2[j] * w3[k] for i in 1:x for j in 1:y for k in 1:z]

    coords = vcat(x_points', y_points', z_points')
    return(coords, weights)
end

# 三角形积分点和权重
function tri_quad(order::Int)
    @assert order in 1:3 "$order 不在1~3范围内"
    coords =  [1/3  0.5  0    0.5  1/3   0.5  0    0.5  1     0     0
               1/3  0.5  0.5  0    1/3   0.5  0.5  0    0     1     0
               1/3  0    0.5  0.5  1/3   0    0.5  0    0     0     1]
    weights = [0.5  1/6  1/6  1/6  0.225 1/15 1/15 1/15 0.025 0.025 0.025]
    # 以上数据来自朱伯芳《有限单元法原理与应用（第四版）》P168
    s = [1 1
         2 4
         5 11]
    return(coords[:, s[order, 1]:s[order, 2]], weights[s[order, 1]:s[order, 2]])
end

# 四面体积分点和权重
function tet_quad(order::Int)
    @assert order in 1:3 "$order 不在1~3范围内"
    coords = [0.25 0.5854102 0.1381966 0.1381966 0.1381966 0.25 1/3  1/6  1/6  1/6
              0.25 0.1381966 0.5854102 0.1381966 0.1381966 0.25 1/6  1/3  1/6  1/6
              0.25 0.1381966 0.1381966 0.5854102 0.1381966 0.25 1/6  1/6  1/3  1/6
              0.25 0.1381966 0.1381966 0.1381966 0.5854102 0.25 1/6  1/6  1/6  1/3]
    weights = [1   0.25      0.25      0.25      0.25      0.8  0.45 0.45 0.45 0.45]
    # 以上数据来自朱伯芳《有限单元法原理与应用（第四版）》P168
    s = [1 1
         2 5
         6 10]
    return(coords[:, s[order, 1]:s[order, 2]], weights[s[order, 1]:s[order, 2]])
end

function quad_form(mode::Char, p::Int...)

    # 定义模式到函数的映射
    mode_map = Dict(
        'q' => (params) -> gauss_legendre(params[1], params[2]),
        'H' => (params) -> gauss_legendre(params[1], params[2], params[3]),
        't' => (params) -> tri_quad(params[1]),
        'T' => (params) -> tet_quad(params[1])
    )

    # 查找模式对应的函数
    func = get(mode_map, mode, nothing)

    if func === nothing
        error("The mode '$mode' is not valid. Use one of 't', 'q', 'H', or 'T'.")
    end

    # 参数数量检查
    valid_params = Dict(
        'q' => 2,
        'H' => 3,
        't' => 1,
        'T' => 1
    )

    required_params = valid_params[mode]
    length(p) == required_params || error("For mode '$mode', exactly $required_params parameters are required.")

    # 调用相应的函数
    return func(p)
end
