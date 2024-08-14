"""
定义了：
                   FE
                   |
           |-------|--------|
         fe2d             fe3d
           |                |
      |----|----|     |-----|-----|
     fet       feq   feT   feH   feP
      |         |     |     |     |
 |----|----|    |     |     |     |
fet1 fet2 fet3 feqx  feTx  feHx  fePx

- 三角形： triangle
- 四边形： quadrilateral
- 四面体： tetrahedron
- 六面体： hexahedron
- 三棱柱： prism
"""
abstract type FE end
abstract type fe2d <: FE end
abstract type fe3d <: FE end

abstract type fet <: fe2d end
abstract type feq <: fe2d end
abstract type feT <: fe3d end
abstract type feH <: fe3d end
abstract type feP <: fe3d end

abstract type fet1 <: fet end
abstract type fet2 <: fet  end
abstract type fet3 <: fet  end

abstract type feq1 <: feq end
abstract type feq2 <: feq end
abstract type feq3 <: feq end

abstract type feT1 <: feT end
abstract type feT2 <: feT end
abstract type feT3<: feT end

abstract type feH1 <: feH end
abstract type feH2 <: feH end
abstract type feH3 <: feH end

abstract type feP1 <: feP end
abstract type feP2 <: feP end
abstract type feP3 <: feP end

abstract type AbstractConstitutiveLaw{T<:FE} end
abstract type AbstractQuadraturePoint{T<:FE} end
abstract type AbstractElement{T<:FE} end
abstract type AbstractFem{T<:FE} end
