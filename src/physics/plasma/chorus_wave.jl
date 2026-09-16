struct ChorusWave
    B_w0    :: Float64    # peak wave amplitude [T]
    ω       :: Float64    # wave frequency [rad/s]
    θ       :: Float64    # wave normal angle [rad] (0 = ducted)
    #λ0      :: Float64    # source latitude [rad]
    #Δλ      :: Float64    # latitudinal extent [rad]
    #ducted  :: Bool       # ducted or nonducted
end


wave_amplitude(wave::ChorusWave, λ) = nothing


wave_normal_angle(wave::ChorusWave, λ, plasma) = nothing
