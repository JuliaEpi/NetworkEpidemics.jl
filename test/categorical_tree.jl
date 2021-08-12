
@testset "CategoricalTree is a collection" begin
    xs = collect(1:8)
    ct = CategoricalTree(xs)

    @test isempty(CategoricalTree([]))
    @test !isempty(ct)

    @test length(ct) == length(xs)

end

@testset "CategoricalTree is an indexable collection" begin
    xs = collect(1:8)
    ct = CategoricalTree(xs)
    
    @test firstindex(ct) == firstindex(xs)
    @test lastindex(ct) == lastindex(xs)
end