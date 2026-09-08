# TODO: make.

@testitem "dipole_field" begin
    using AURORA
    using AURORA: RE, μ₀, M
    using StaticArrays
    using LinearAlgebra

    ## ==================== Error handling ==================== ##
    @testset "Throws at origin" begin
        @test_throws ArgumentError dipole_field(0.0, 0.0, 0.0)
    end

    @testset "Throws inside Earth" begin
        @test_throws ArgumentError dipole_field(RE * 0.5, 0.0, 0.0)
    end

    @testset "Throws on Earth surface" begin
        @test_throws ArgumentError dipole_field(RE, 0.0, 0.0)
    end

    @testset "Throws outside valid range" begin
        @test_throws ArgumentError dipole_field(11 * RE, 0.0, 0.0)
    end

    @testset "Throws for wrong vector length" begin
        @test_throws ArgumentError dipole_field([RE * 2, 0.0])
        @test_throws ArgumentError dipole_field([RE * 2, 0.0, 0.0, 0.0])
    end

    ## ==================== Return type ==================== ##
    @testset "Returns SVector" begin
        B = dipole_field(2RE, 0.0, 0.0)
        @test B isa SVector{3, Float64}
    end

    ## ==================== Analytical values ==================== ##
    @testset "Equatorial field (z=0, on x-axis)" begin
        # At the magnetic equator on the x-axis, the dipole field is purely z-directed:
        # Bz = +(μ₀/4π) * M / r³,  Bx = By = 0
        r = 2.0 * RE
        B = dipole_field(r, 0.0, 0.0)

        Bz_analytical = (μ₀ / 4π) * M / r^3

        @test B[1] ≈ 0.0 atol=1e-20   # Bx = 0
        @test B[2] ≈ 0.0 atol=1e-20   # By = 0
        @test B[3] ≈ Bz_analytical rtol=1e-10
    end

    @testset "Field on magnetic axis (x=y=0)" begin
        # On the magnetic axis (x=y=0), the field is purely z-directed:
        # Bz = -2 * (μ₀/4π) * M / r³
        r = 2.0 * RE
        B = dipole_field(0.0, 0.0, r)

        Bz_analytical = -2 * (μ₀ / 4π) * M / r^3

        @test B[1] ≈ 0.0 atol=1e-20
        @test B[2] ≈ 0.0 atol=1e-20
        @test B[3] ≈ Bz_analytical rtol=1e-10
    end

    ## ==================== Symmetry ==================== ##
    @testset "Azimuthal symmetry (rotation around z-axis)" begin
        # Rotating around z-axis should not change |B|
        r = 2.0 * RE
        B1 = dipole_field(r, 0.0, RE)
        B2 = dipole_field(0.0, r, RE)  # 90° rotation around z

        @test norm(B1) ≈ norm(B2) rtol=1e-10
    end

    @testset "North-south symmetry" begin
        # Bx and By should be antisymmetric in z, Bz symmetric
        r = 2.0 * RE
        B_north = dipole_field(r, 0.0, RE)
        B_south = dipole_field(r, 0.0, -RE)

        @test B_north[1] ≈ -B_south[1] rtol=1e-10  # Bx antisymmetric
        @test B_north[2] ≈ -B_south[2] rtol=1e-10  # By antisymmetric
        @test B_north[3] ≈  B_south[3] rtol=1e-10  # Bz symmetric
    end

    ## ==================== Vector method ==================== ##
    @testset "Vector and scalar methods agree" begin
        B_scalar = dipole_field(2RE, 0.0, RE)
        B_vector = dipole_field([2RE, 0.0, RE])
        @test B_scalar ≈ B_vector
    end
end

@testitem "magnetic_basis" begin
    using AURORA
    using LinearAlgebra

    ## ==================== Orthonormality ==================== ##
    @testset "Basis vectors are unit length" begin
        B = [1.0, 2.0, 3.0]
        b, e1, e2, _ = magnetic_basis(B)

        @test norm(b)  ≈ 1.0 rtol=1e-10
        @test norm(e1) ≈ 1.0 rtol=1e-10
        @test norm(e2) ≈ 1.0 rtol=1e-10
    end

    @testset "Basis vectors are mutually orthogonal" begin
        B = [1.0, 2.0, 3.0]
        b, e1, e2, _ = magnetic_basis(B)

        @test dot(b, e1)  ≈ 0.0 atol=1e-10
        @test dot(b, e2)  ≈ 0.0 atol=1e-10
        @test dot(e1, e2) ≈ 0.0 atol=1e-10
    end

    ## ==================== b̂ direction ==================== ##
    @testset "b̂ is parallel to B" begin
        B = [1.0, 2.0, 3.0]
        b, _, _, B_mag = magnetic_basis(B)

        @test b ≈ B ./ B_mag rtol=1e-10
    end

    @testset "B_mag is correct" begin
        B = [3.0, 4.0, 0.0]
        _, _, _, B_mag = magnetic_basis(B)

        @test B_mag ≈ 5.0 rtol=1e-10
    end

    ## ==================== flip ==================== ##
    @testset "flip=true reverses b̂" begin
        B = [1.0, 2.0, 3.0]
        b_normal, _, _, _ = magnetic_basis(B)
        b_flipped, _, _, _ = magnetic_basis(B; flip=true)

        @test b_flipped ≈ -b_normal rtol=1e-10
    end

    @testset "flip=false keeps b̂ parallel to B" begin
        B = [1.0, 0.0, 0.0]
        b, _, _, _ = magnetic_basis(B; flip=false)

        @test dot(b, B) > 0   # b̂ should point same way as B
    end

    ## ==================== Error handling ==================== ##
    @testset "Throws for zero field" begin
        @test_throws ArgumentError magnetic_basis([0.0, 0.0, 0.0])
    end
end
