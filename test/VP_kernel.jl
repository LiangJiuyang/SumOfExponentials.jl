@testset "VP kernel accuracy, exp(-x^2)" begin
    f = x -> exp(-x^2)

    s, w = VP_sw(f, 4.0, 30, GaussParameter(100))
    x_array = [0.0:0.001:10.0...]
    error = [soe(x, s, w) - exp(-x^2) - w[1] for x in x_array]

    @test maximum(abs.(error)) < 1e-8
end