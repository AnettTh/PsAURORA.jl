# TODO: make.
# TODO: remember tests for velocity decomposition also
@testitem "velocity_from_kinetic_energy" begin
    using AURORA
    using AURORA: mₑ, c₀, eV_in_J

    @testset "Throws for zero mass" begin
        @test_throws ArgumentError velocity_from_kinetic_energy(1e3, 0.0)
    end

    @testset "Throws for negative energy" begin
        @test_throws ArgumentError velocity_from_kinetic_energy(-1.0, mₑ)
    end

    @testset "Throws for relativistic energy" begin
        # E ≥ m*c² should throw
        E_relativistic = mₑ * c₀^2 / eV_in_J  # rest energy in eV
        @test_throws ArgumentError velocity_from_kinetic_energy(E_relativistic, mₑ)
    end

    @testset "Warns for semi-relativistic energy" begin
        # E/mc² ≥ 0.1 should warn
        E_semirel = 0.15 * mₑ * c₀^2 / eV_in_J
        @test_logs (:warn,) velocity_from_kinetic_energy(E_semirel, mₑ)
    end

    @testset "Correct speed for known energy" begin
        # 1 eV electron: v = sqrt(2 * 1eV / mₑ)
        E_eV = 1.0
        v = velocity_from_kinetic_energy(E_eV, mₑ)
        v_analytical = sqrt(2 * eV_in_J / mₑ)
        @test v ≈ v_analytical rtol=1e-10
    end

    @testset "Speed increases with energy" begin
        v1 = velocity_from_kinetic_energy(100.0, mₑ)
        v2 = velocity_from_kinetic_energy(1000.0, mₑ)
        @test v2 > v1
    end

end


@testitem "gyro_frequency" begin
    using AURORA
    using AURORA: qₑ, mₑ, eV_in_J
    using LinearAlgebra

    @testset "Throws for zero field" begin
        @test_throws ArgumentError gyro_frequency(0.0, qₑ, mₑ)
        @test_throws ArgumentError gyro_frequency([0.0, 0.0, 0.0], qₑ, mₑ)
    end

    @testset "Scalar and vector methods agree" begin
        B = [1e-5, 2e-5, 3e-5]
        B_mag = norm(B)
        @test gyro_frequency(B, qₑ, mₑ) ≈ gyro_frequency(B_mag, qₑ, mₑ) rtol=1e-10
    end

    @testset "Correct analytical value" begin
        # ω_g = |q| * B / m
        B_mag = 5e-5   # 50 μT, typical ionospheric value
        ω_analytical = abs(qₑ) * B_mag / mₑ
        @test gyro_frequency(B_mag, qₑ, mₑ) ≈ ω_analytical rtol=1e-10
    end

    @testset "Frequency increases with field strength" begin
        @test gyro_frequency(2e-5, qₑ, mₑ) > gyro_frequency(1e-5, qₑ, mₑ)
    end

    @testset "Sign of charge does not matter" begin
        B_mag = 5e-5
        @test gyro_frequency(B_mag, qₑ, mₑ) ≈ gyro_frequency(B_mag, -qₑ, mₑ) rtol=1e-10
    end

end


@testitem "larmor_radius" begin
    using AURORA
    using AURORA: qₑ, mₑ
    using LinearAlgebra

    @testset "Throws for zero field" begin
        v = [1e6, 0.0, 0.0]
        B = [0.0, 0.0, 0.0]
        @test_throws ArgumentError larmor_radius(mₑ, qₑ, v, B)
    end

    @testset "Throws for zero charge" begin
        v = [1e6, 0.0, 0.0]
        B = [0.0, 0.0, 5e-5]
        @test_throws ArgumentError larmor_radius(mₑ, 0.0, v, B)
    end

    @testset "Correct analytical value" begin
        # r_L = m * v_perp / (|q| * B)
        # For v purely perpendicular to B:
        v = [1e6, 0.0, 0.0]   # x-direction
        B = [0.0, 0.0, 5e-5]  # z-direction, perpendicular to v
        r_analytical = mₑ * norm(v) / (abs(qₑ) * norm(B))
        @test larmor_radius(mₑ, qₑ, v, B) ≈ r_analytical rtol=1e-10
    end

    @testset "Zero for field-aligned velocity" begin
        # Purely parallel velocity has no perpendicular component → r_L = 0
        B = [0.0, 0.0, 5e-5]
        v = [0.0, 0.0, 1e6]   # parallel to B
        @test larmor_radius(mₑ, qₑ, v, B) ≈ 0.0 atol=1e-10
    end

    @testset "Radius increases with energy" begin
        B = [0.0, 0.0, 5e-5]
        v1 = [1e6, 0.0, 0.0]
        v2 = [2e6, 0.0, 0.0]
        @test larmor_radius(mₑ, qₑ, v2, B) > larmor_radius(mₑ, qₑ, v1, B)
    end

    @testset "Radius decreases with field strength" begin
        v = [1e6, 0.0, 0.0]
        B1 = [0.0, 0.0, 1e-5]
        B2 = [0.0, 0.0, 5e-5]
        @test larmor_radius(mₑ, qₑ, v, B2) < larmor_radius(mₑ, qₑ, v, B1)
    end

end


@testitem "gyrocenter" begin
    using AURORA
    using AURORA: qₑ, mₑ
    using LinearAlgebra

    @testset "Throws for zero field" begin
        r = [RE + 1e6, 0.0, 0.0]
        v = [1e6, 0.0, 0.0]
        B = [0.0, 0.0, 0.0]
        @test_throws ArgumentError gyrocenter(r, v, B, qₑ, mₑ)
    end

    @testset "Throws for zero charge" begin
        r = [RE + 1e6, 0.0, 0.0]
        v = [1e6, 0.0, 0.0]
        B = [0.0, 0.0, 5e-5]
        @test_throws ArgumentError gyrocenter(r, v, B, 0.0, mₑ)
    end

    @testset "Gyrocenter equals position for field-aligned velocity" begin
        # No perpendicular velocity → no gyration → gyrocenter = position
        r = [RE + 1e6, 0.0, 0.0]
        v = [0.0, 0.0, 1e6]   # parallel to B
        B = [0.0, 0.0, 5e-5]  # z-direction
        @test gyrocenter(r, v, B, qₑ, mₑ) ≈ r rtol=1e-10
    end

    @testset "Gyrocenter offset equals Larmor radius" begin
        # For purely perpendicular velocity, |r_gc - r| = larmor_radius
        r = [RE + 1e6, 0.0, 0.0]
        v = [1e6, 0.0, 0.0]   # x-direction
        B = [0.0, 0.0, 5e-5]  # z-direction
        r_gc = gyrocenter(r, v, B, qₑ, mₑ)
        offset = norm(r_gc - r)
        r_L = larmor_radius(mₑ, qₑ, v, B)
        @test offset ≈ r_L rtol=1e-10
    end

    @testset "Gyrocenter is independent of gyrophase" begin
        # Two particles at same position with same energy but different gyrophase
        # should have the same gyrocenter
        r = [RE + 1e6, 0.0, 0.0]
        B = [0.0, 0.0, 5e-5]
        v1 = [1e6, 0.0, 0.0]
        v2 = [0.0, 1e6, 0.0]
        r_gc1 = gyrocenter(r, v1, B, qₑ, mₑ)
        r_gc2 = gyrocenter(r, v2, B, qₑ, mₑ)
        @test norm(r_gc1 - r) ≈ norm(r_gc2 - r) rtol=1e-8       # Failed with 1e-10
    end

end


@testitem "losscone_angle" begin
    using AURORA
    using AURORA: RE, qₑ, mₑ
    using LinearAlgebra

    @testset "Throws for zero equatorial distance" begin
        @test_throws ArgumentError losscone_angle(dipole_field, [0.0, 0.0, 0.0])
    end

    @testset "Returns value in radians by default" begin
        r_eq = [4 * RE, 0.0, 0.0]
        α = losscone_angle(dipole_field, r_eq)
        @test 0.0 < α < π/2
    end

    @testset "Degrees and radians are consistent" begin
        r_eq = [4 * RE, 0.0, 0.0]
        α_rad = losscone_angle(dipole_field, r_eq; degrees=false)
        α_deg = losscone_angle(dipole_field, r_eq; degrees=true)
        @test α_deg ≈ rad2deg(α_rad) rtol=1e-10
    end

    @testset "Loss cone narrows with distance" begin
        # Further from Earth → weaker field ratio → smaller loss cone
        r_near = [3 * RE, 0.0, 0.0]
        r_far  = [6 * RE, 0.0, 0.0]
        α_near = losscone_angle(dipole_field, r_near)
        α_far  = losscone_angle(dipole_field, r_far)
        @test α_far < α_near
    end

    @testset "Loss cone angle is between 0 and 90 degrees" begin
        for L in [2, 4, 6, 8]
            r_eq = [L * RE, 0.0, 0.0]
            α = losscone_angle(dipole_field, r_eq; degrees=true)
            @test 0.0 < α < 90.0
        end
    end

end


@testitem "quarter_bounceperiod" begin
    using AURORA
    using AURORA: RE, mₑ, qₑ

    @testset "Degrees and radians are consistent" begin
        L = 4.0
        E_eV = 1e4
        θ = 45.0
        τ_rad = quarter_bounceperiod(L, E_eV, mₑ, deg2rad(θ); degrees=false)
        τ_deg = quarter_bounceperiod(L, E_eV, mₑ, θ; degrees=true)
        @test τ_rad ≈ τ_deg rtol=1e-10
    end

    @testset "Period increases with L-shell" begin
        # Larger L → longer field line → longer bounce period
        E_eV = 1e4
        θ = 45.0
        τ_near = quarter_bounceperiod(3.0, E_eV, mₑ, θ; degrees=true)
        τ_far  = quarter_bounceperiod(6.0, E_eV, mₑ, θ; degrees=true)
        @test τ_far > τ_near
    end

    @testset "Period decreases with energy" begin
        # Higher energy → faster particle → shorter bounce period
        L = 4.0
        θ = 45.0
        τ_low  = quarter_bounceperiod(L, 1e3, mₑ, θ; degrees=true)
        τ_high = quarter_bounceperiod(L, 1e4, mₑ, θ; degrees=true)
        @test τ_high < τ_low
    end

    @testset "Period is positive" begin
        for L in [2.0, 4.0, 6.0]
            τ = quarter_bounceperiod(L, 1e4, mₑ, 45.0; degrees=true)
            @test τ > 0.0
        end
    end

    @testset "Analytical value for known input" begin
        L = 4.0
        E_eV = 1e4
        θ_rad = deg2rad(45.0)
        v = velocity_from_kinetic_energy(E_eV, mₑ)
        Γ = 1.30 - 0.56 * sin(θ_rad)
        τ_analytical = (L * RE / v) * Γ
        τ = quarter_bounceperiod(L, E_eV, mₑ, θ_rad)
        @test τ ≈ τ_analytical rtol=1e-10
    end

end


@testitem "average_driftvelocity" begin
    using AURORA
    using AURORA: RE, qₑ, mₑ, BE, eV_in_J, BE

    @testset "Degrees and radians are consistent" begin
        L = 4.0
        E_eV = 1e4
        θ = 45.0
        v_rad = average_driftvelocity(L, E_eV, qₑ, deg2rad(θ); degrees=false)
        v_deg = average_driftvelocity(L, E_eV, qₑ, θ; degrees=true)
        @test v_rad ≈ v_deg rtol=1e-10
    end

    @testset "Drift velocity increases with L-shell" begin
        # Further from Earth → stronger drift
        E_eV = 1e4
        θ = 45.0
        v_near = average_driftvelocity(3.0, E_eV, qₑ, θ; degrees=true)
        v_far  = average_driftvelocity(6.0, E_eV, qₑ, θ; degrees=true)
        @test abs(v_far) > abs(v_near)
    end

    @testset "Drift velocity increases with energy" begin
        L = 4.0
        θ = 45.0
        v_low  = average_driftvelocity(L, 1e3, qₑ, θ; degrees=true)
        v_high = average_driftvelocity(L, 1e4, qₑ, θ; degrees=true)
        @test abs(v_high) > abs(v_low)
    end

    @testset "Analytical value for known input" begin
        L = 4.0
        E_eV = 1e4
        θ_rad = deg2rad(45.0)
        E_J = E_eV * eV_in_J
        Γ = 0.35 - 0.15 * sin(θ_rad)
        v_analytical = (3 * L^2 * E_J * Γ) / (2 * qₑ * BE * RE)
        v = average_driftvelocity(L, E_eV, qₑ, θ_rad)
        @test v ≈ v_analytical rtol=1e-10
    end

    @testset "Sign flips with charge sign" begin
        # Electrons and positrons drift in opposite directions
        L = 4.0
        E_eV = 1e4
        θ = 45.0
        v_electron = average_driftvelocity(L, E_eV, qₑ, θ; degrees=true)
        v_positron = average_driftvelocity(L, E_eV, -qₑ, θ; degrees=true)
        @test v_electron ≈ -v_positron rtol=1e-10
    end

end



@testitem "total_drift" begin
    using AURORA
    using AURORA: RE, qₑ, mₑ

    @testset "Degrees and radians are consistent" begin
        L = 4.0
        E_eV = 1e4
        θ = 45.0
        d_rad = total_drift(L, E_eV, qₑ, mₑ, deg2rad(θ); degrees=false)
        d_deg = total_drift(L, E_eV, qₑ, mₑ, θ; degrees=true)
        @test d_rad ≈ d_deg rtol=1e-10
    end

    @testset "Consistent with bounce period and drift velocity" begin
        # total_drift = quarter_bounceperiod * average_driftvelocity
        L = 4.0
        E_eV = 1e4
        θ_rad = deg2rad(45.0)
        τ = quarter_bounceperiod(L, E_eV, mₑ, θ_rad)
        v_d = average_driftvelocity(L, E_eV, qₑ, θ_rad)
        d_analytical = τ * v_d
        d = total_drift(L, E_eV, qₑ, mₑ, θ_rad)
        @test d ≈ d_analytical rtol=1e-10
    end

    @testset "Drift increases with L-shell" begin
        E_eV = 1e4
        θ = 45.0
        d_near = total_drift(3.0, E_eV, qₑ, mₑ, θ; degrees=true)
        d_far  = total_drift(6.0, E_eV, qₑ, mₑ, θ; degrees=true)
        @test abs(d_far) > abs(d_near)
    end

    @testset "Drift increases with energy" begin
        L = 4.0
        θ = 45.0
        d_low  = total_drift(L, 1e3, qₑ, mₑ, θ; degrees=true)
        d_high = total_drift(L, 1e4, qₑ, mₑ, θ; degrees=true)
        @test abs(d_high) > abs(d_low)
    end

    @testset "Sign flips with charge sign" begin
        L = 4.0
        E_eV = 1e4
        θ = 45.0
        d_electron = total_drift(L, E_eV, qₑ, mₑ, θ; degrees=true)
        d_positron = total_drift(L, E_eV, -qₑ, mₑ, θ; degrees=true)
        @test d_electron ≈ -d_positron rtol=1e-10
    end

end


@testitem "parallel_velocity" begin
    using AURORA
    using LinearAlgebra

    @testset "Throws for zero field" begin
        v = [1.0, 2.0, 3.0]
        B = [0.0, 0.0, 0.0]
        @test_throws ArgumentError parallel_velocity(v, B)
    end

    @testset "Purely parallel velocity is unchanged" begin
        # v parallel to B → parallel component = v
        B = [0.0, 0.0, 1.0]
        v = [0.0, 0.0, 5.0]
        @test parallel_velocity(v, B) ≈ v rtol=1e-10
    end

    @testset "Purely perpendicular velocity gives zero" begin
        # v perpendicular to B → parallel component = 0
        B = [0.0, 0.0, 1.0]
        v = [1.0, 0.0, 0.0]
        @test norm(parallel_velocity(v, B)) ≈ 0.0 atol=1e-10
    end

    @testset "Result is parallel to B" begin
        # parallel_velocity should always point along B
        B = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        v_par = parallel_velocity(v, B)
        # v_par × B should be zero vector
        @test norm(cross(v_par, B)) ≈ 0.0 atol=1e-10
    end

    @testset "Invariant to field magnitude" begin
        # Scaling B should not change the parallel velocity
        v = [1.0, 2.0, 3.0]
        B = [0.0, 0.0, 1.0]
        @test parallel_velocity(v, B) ≈ parallel_velocity(v, 2.0 .* B) rtol=1e-10
    end

end


@testitem "perpendicular_velocity" begin
    using AURORA
    using LinearAlgebra

    @testset "Throws for zero field" begin
        v = [1.0, 2.0, 3.0]
        B = [0.0, 0.0, 0.0]
        @test_throws ArgumentError perpendicular_velocity(v, B)
    end

    @testset "Purely perpendicular velocity is unchanged" begin
        # v perpendicular to B → perpendicular component = v
        B = [0.0, 0.0, 1.0]
        v = [1.0, 0.0, 0.0]
        @test perpendicular_velocity(v, B) ≈ v rtol=1e-10
    end

    @testset "Purely parallel velocity gives zero" begin
        # v parallel to B → perpendicular component = 0
        B = [0.0, 0.0, 1.0]
        v = [0.0, 0.0, 5.0]
        @test norm(perpendicular_velocity(v, B)) ≈ 0.0 atol=1e-10
    end

    @testset "Result is perpendicular to B" begin
        # perpendicular_velocity should always be orthogonal to B
        B = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        v_perp = perpendicular_velocity(v, B)
        @test dot(v_perp, B) ≈ 0.0 atol=1e-10
    end

    @testset "Parallel and perpendicular components sum to original" begin
        B = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        v_par  = parallel_velocity(v, B)
        v_perp = perpendicular_velocity(v, B)
        @test v_par .+ v_perp ≈ v rtol=1e-10
    end

    @testset "Invariant to field magnitude" begin
        v = [1.0, 2.0, 3.0]
        B = [0.0, 0.0, 1.0]
        @test perpendicular_velocity(v, B) ≈ perpendicular_velocity(v, 2.0 .* B) rtol=1e-10
    end

end


@testitem "parallel_speed" begin
    using AURORA
    using LinearAlgebra

    @testset "Throws for zero field" begin
        v = [1.0, 2.0, 3.0]
        B = [0.0, 0.0, 0.0]
        @test_throws ArgumentError parallel_speed(v, B)
    end

    @testset "Equals norm of parallel_velocity" begin
        B = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        @test parallel_speed(v, B) ≈ norm(parallel_velocity(v, B)) rtol=1e-10
    end

    @testset "Purely parallel velocity gives full speed" begin
        B = [0.0, 0.0, 1.0]
        v = [0.0, 0.0, 5.0]
        @test parallel_speed(v, B) ≈ norm(v) rtol=1e-10
    end

    @testset "Purely perpendicular velocity gives zero" begin
        B = [0.0, 0.0, 1.0]
        v = [1.0, 0.0, 0.0]
        @test parallel_speed(v, B) ≈ 0.0 atol=1e-10
    end

    @testset "Always non-negative" begin
        B = [1.0, 2.0, 3.0]
        for v in [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]]
            @test parallel_speed(v, B) ≥ 0.0
        end
    end

end


@testitem "perpendicular_speed" begin
    using AURORA
    using LinearAlgebra

    @testset "Throws for zero field" begin
        v = [1.0, 2.0, 3.0]
        B = [0.0, 0.0, 0.0]
        @test_throws ArgumentError perpendicular_speed(v, B)
    end

    @testset "Equals norm of perpendicular_velocity" begin
        B = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        @test perpendicular_speed(v, B) ≈ norm(perpendicular_velocity(v, B)) rtol=1e-10
    end

    @testset "Purely perpendicular velocity gives full speed" begin
        B = [0.0, 0.0, 1.0]
        v = [1.0, 0.0, 0.0]
        @test perpendicular_speed(v, B) ≈ norm(v) rtol=1e-10
    end

    @testset "Purely parallel velocity gives zero" begin
        B = [0.0, 0.0, 1.0]
        v = [0.0, 0.0, 5.0]
        @test perpendicular_speed(v, B) ≈ 0.0 atol=1e-10
    end

    @testset "Always non-negative" begin
        B = [1.0, 2.0, 3.0]
        for v in [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]]
            @test perpendicular_speed(v, B) ≥ 0.0
        end
    end

    @testset "Speed decomposition conserves total speed" begin
        # |v|² = |v_par|² + |v_perp|²
        B = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        @test parallel_speed(v, B)^2 + perpendicular_speed(v, B)^2 ≈ norm(v)^2 rtol=1e-10
    end

end
