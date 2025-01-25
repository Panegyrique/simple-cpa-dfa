using Statistics
using NPZ
using Plots
using ProgressBars
using Printf


function plot_traces(data; num_traces::Int=5, chunk_size::Int=10000)
    p = plot(size=(1200, 600), legend=false)
    
    x = 1:size(data, 2)

    for i in 1:min(num_traces, size(data, 1))
        trace_data = zeros(size(data, 2))
        
        for start_idx in 1:chunk_size:size(data, 2)
            end_idx = min(start_idx + chunk_size - 1, size(data, 2))
            trace_data[start_idx:end_idx] = data[i, start_idx:end_idx]
        end
        
        plot!(p, x, trace_data)
    end
    
    plot!(p, 
        xlabel="t", 
        ylabel="V", 
        grid=true,
    )

    display(p)
end


function get_sbox()
    [
        # 0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
        0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76, # 0
        0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, # 1
        0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15, # 2
        0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75, # 3
        0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, # 4
        0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf, # 5
        0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8, # 6
        0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, # 7
        0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73, # 8
        0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb, # 9
        0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, # A
        0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08, # B
        0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a, # C
        0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, # D
        0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf, # E
        0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16  # F
    ]
end


function hamming_weight(x::UInt8)
    count_ones(x)
end


function theoretical_power_leakage(hypothesis, plaintext, sbox, hamming_weights)
    index = (hypothesis ⊻ plaintext) + 1
    sbox_value = sbox[index]
    return hamming_weights[sbox_value + 1]
end


function optimized_pearson_correlation(X, Y)
    X_centered = X .- mean(X)
    Y_centered = Y .- mean(Y)
    numerator = sum(X_centered .* Y_centered)
    denominator = sqrt(sum(X_centered.^2) * sum(Y_centered.^2))
    denominator == 0 ? 0.0 : numerator / denominator
end


function correlation_power_analysis(plaintexts, traces; 
                              number_step = 1)
    
    number_traces, number_points = size(traces)
    sbox = get_sbox()
    hamming_weights = [hamming_weight(UInt8(n)) for n in 0x00:0xff] # precompute hamming weights
    
    recover_key = zeros(UInt8, 16)
    step_size = number_traces ÷ number_step
    bytes_guesses_evolution = zeros(16, 256, number_step)
    correlation_over_time = zeros(16, 256, number_points)

    pbar = ProgressBar(1:16*256)

    for byte_index in 1:16
        for hypothesis in 0x00:0xff
            set_postfix(pbar, Byte=@sprintf("%d/16", byte_index), Hypothesis=@sprintf("%d / 256", hypothesis + 1))
            update(pbar)

            theorical_power = [
                theoretical_power_leakage(hypothesis, plaintexts[trace_index, byte_index], 
                                          sbox, hamming_weights) 
                for trace_index in 1:number_traces
            ]
            
            for i in 1:number_step
                num_traces_to_use = i * step_size
                
                point_correlation = [
                    optimized_pearson_correlation(
                        theorical_power[1:num_traces_to_use], 
                        traces[1:num_traces_to_use, point_index]
                    ) 
                    for point_index in 1:number_points
                ]
                
                bytes_guesses_evolution[byte_index, hypothesis+1, i] = maximum(abs.(point_correlation))
                correlation_over_time[byte_index, hypothesis+1, :] = point_correlation
            end
        end
        
        final_correlations = bytes_guesses_evolution[byte_index, :, end]
        recover_key[byte_index] = argmax(final_correlations) - 1
    end
    
    return recover_key, bytes_guesses_evolution, correlation_over_time
end


function plot_correlation(octet, correlation_over_time, recover_key)
    num_hypothesis, num_points = size(correlation_over_time[octet, :, :])
    
    x = 1:num_points
    best_hypothesis = recover_key[octet] + 1
    
    p1 = plot(size=(1200, 600), legend=false)
    p2 = plot(size=(1200, 600), legend=true)
    
    for hypothesis in 1:num_hypothesis
        y = correlation_over_time[octet, hypothesis, :]
        
        plot!(p1, x, y)
        
        if hypothesis != best_hypothesis
            plot!(p2, x, y, color=:gray, alpha=0.5, label="")
        end
    end

    y_best = correlation_over_time[octet, best_hypothesis, :]
    plot!(p2, x, y_best, color=:green, linewidth=2, label="Best Hypothesis: $best_hypothesis")
    
    plot!(p1, xlabel="Time", ylabel="Correlation", title="All Hypotheses for Byte $octet", grid=true)
    plot!(p2, xlabel="Time", ylabel="Correlation", title="Best Hypothesis for Byte $octet", grid=true, 
          left_margin=10Plots.mm, bottom_margin=10Plots.mm)

    display(plot(p1, p2, layout=@layout([a{0.5h}; b{0.5h}]), size=(1200, 1300)))
end


# Load data
plaintext = npzread("./pts.npy")
println("Inputs shape: ", size(plaintext))
traces = npzread("./traces.npy")
println("Traces shape: ", size(traces))

# Plot raw traces
# plot_traces(traces, num_traces=50, chunk_size=50_000)

# Launch CPA
recover_key, _, correlation_over_time = correlation_power_analysis(plaintext, traces)

# Print result
hex_values = ["0x" * string(Int(val), base=16, pad=2) for val in recover_key]
ascii_string = join([Char(val) for val in recover_key])
println("\nRecovered Key (Hex): ", hex_values)
println("Recovered Key (ASCII): ", ascii_string)

# Plot all correlation
for byte in 1:16
    plot_correlation(byte, correlation_over_time, recover_key)
end