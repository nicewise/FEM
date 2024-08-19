# MandelNotation.jl
# encoding: utf-8
"""
定义了：
- Mandel.δ : 二阶单位张量
- Mandel.I : 对称四阶单位张量
- Mandel.J : 体积四阶单位张量
- Mandel.K : 偏差四阶单位张量
"""

"""
二阶对称张量的降阶表示为：
A = [ A11
      A22
      A33
    √2A23
    √2A31
    √2A12 ]
δ = [1
     1
     1
     0
     0
     0]
四阶张量U的降阶表示为：
U = [ U1111   U1122   U1133 √2U1123 √2U1131 √2U1112
      U2211   U2222   U2233 √2U2223 √2U2231 √2U2212
      U3311   U3322   U3333 √2U3323 √2U3331 √2U3312
    √2U2311 √2U2322 √2U2333  2U2323  2U2331  2U2312
    √2U3111 √2U3122 √2U3133  2U3123  2U3131  2U3112
    √2U1211 √2U1222 √2U1233  2U1223  2U1231  2U1212 ]
I = 1/2 (δ_ik δ_jl + δ_il δ_jk) ei⊗ej⊗ek⊗el
  = [ 1  0  0  0  0  0
      0  1  0  0  0  0
      0  0  1  0  0  0
      0  0  0  1  0  0
      0  0  0  0  1  0
      0  0  0  0  0  1 ]
J = 1/3 δ_ij δ_kl
  = [ 1  1  1  0  0  0
      1  1  1  0  0  0
      1  1  1  0  0  0
      0  0  0  0  0  0
      0  0  0  0  0  0
      0  0  0  0  0  0 ] / 3
"""
struct MandelNotation
    δ::SVector{6, Float64}
    I::SMatrix{6, 6, Float64}
    J::SMatrix{6, 6, Float64}
    K::SMatrix{6, 6, Float64}

    function MandelNotation()
        δ = @SVector [1.0, 1, 1, 0, 0, 0]
        I = @SMatrix [ 1  0  0  0  0  0.0
                       0  1  0  0  0  0
                       0  0  1  0  0  0
                       0  0  0  1  0  0
                       0  0  0  0  1  0
                       0  0  0  0  0  1 ]
        #J = δ*δ'/3
        #K = I - J
        J = SMatrix{6, 6}([ 1  1  1  0  0  0.0
                            1  1  1  0  0  0
                            1  1  1  0  0  0
                            0  0  0  0  0  0
                            0  0  0  0  0  0
                            0  0  0  0  0  0 ] / 3)
        K = SMatrix{6, 6}([ 2 -1 -1  0  0  0.0
                           -1  2 -1  0  0  0
                           -1 -1  2  0  0  0
                            0  0  0  3  0  0
                            0  0  0  0  3  0
                            0  0  0  0  0  3] / 3)
        return new(δ, I, J, K)
    end
end

const Mandel = MandelNotation()
