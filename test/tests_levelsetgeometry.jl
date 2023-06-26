using DrWatson, Test
@quickactivate :LevelSetSublimation

const LSS = LevelSetSublimation

dom1 = Domain(25, 24, 1.3, 2.0)
dom2 = Domain(45, 35, 1.2, 2.1)

ϕ1a = make_ϕ0(:circ, dom1)
reinitialize_ϕ_HCR!(ϕ1a, dom1)

# @test 