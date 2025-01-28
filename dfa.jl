using Statistics
using NPZ
using Plots
using ProgressBars
using Printf


function get_inv_sbox()
    [
        #0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
        0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb, # 0
        0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb, # 1
        0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e, # 2
        0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25, # 3
        0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92, # 4
        0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84, # 5
        0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06, # 6
        0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b, # 7
        0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73, # 8
        0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e, # 9 
        0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b, # A
        0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4, # B
        0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f, # C
        0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef, # D
        0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61, # E
        0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d  # F
    ]
end


function compute_error_distribution(cts, fcts)
    num_cipher, num_bytes = size(cts)
    error_distributions = zeros(Int, num_bytes, 256, num_cipher)
    inv_sbox = get_inv_sbox()
    
    for cipher_index in 1:num_cipher
        for byte_index in 1:num_bytes
            for hypothesis in 0:255
                c = cts[cipher_index, byte_index]
                c_fault = fcts[cipher_index, byte_index]
                idx1 = Int((hypothesis ⊻ c)) + 1
                idx2 = Int((hypothesis ⊻ c_fault)) + 1
                e = inv_sbox[idx1] ⊻ inv_sbox[idx2]
                error_distributions[byte_index, hypothesis + 1, cipher_index] = e
            end
        end
    end
    
    return error_distributions
end


function compute_error_entropy(error_distributions)
    num_bytes, num_hypothesis, num_samples = size(error_distributions)
    entropies = zeros(Float64, num_bytes, num_hypothesis)
    
    for byte_index in 1:num_bytes
        for hypothesis in 1:num_hypothesis
            errors = error_distributions[byte_index, hypothesis, :]
            error_counts = Dict{Int,Int}()
            for error in errors
                error_counts[error] = get(error_counts, error, 0) + 1
            end
            probabilities = [count/num_samples for (_, count) in error_counts]
            entropy = 0.0
            for p in probabilities
                if p > 0
                    entropy -= p * log2(p)
                end
            end
            entropies[byte_index, hypothesis] = entropy
        end
    end
    
    return entropies
end


function get_subkey(entropies)
    num_bytes, _ = size(entropies)
    results = zeros(num_bytes)
    
    for byte_index in 1:num_bytes
        best_key = argmin(entropies[byte_index, :]) - 1
        results[byte_index] = best_key
    end
    
    return results
end


function plot_hypothesis(byte_index, entropies, subkey)
    x = 0:255
    entropy_values = entropies[byte_index, :]
    subkey_values = subkey[byte_index]
    hex_value = "0x" * string(Int(subkey_values), base=16, pad=2)
    
    min_entropy_index = argmin(entropy_values)
    
    p = plot(size=(1200, 600), legend=true)
    
    plot!(p,
        x[1:end .!= min_entropy_index],
        entropy_values[1:end .!= min_entropy_index],
        seriestype=:scatter,
        color=:grey,
        label="",
        xlabel="Hypothesis", 
        ylabel="Entropy", 
        title="Entropy for each byte hypothesis $byte_index",
        grid=true,
        left_margin=10Plots.mm, 
        bottom_margin=10Plots.mm
    )
    
    plot!(p,
        [x[min_entropy_index]],
        [entropy_values[min_entropy_index]],
        seriestype=:scatter,
        color=:green,
        label="Best hypothesis: $hex_value",
    )

    savefig(p, "entropy_byte_$byte_index.png")

    display(p)
end


# Load data
cts = npzread("./data/dfa/sucess/cts.npy")
println("Ciphertexts shape: ", size(cts))
fcts = npzread("./data/dfa/sucess/fcts.npy")
println("Faulty ciphertexts shape: ", size(fcts))

# # DFA
# error_distribution = compute_error_distribution(cts, fcts)
# entropies = compute_error_entropy(error_distribution)
# subkey = get_subkey(entropies)

# # Print results
# hex_values = ["0x" * string(Int(val), base=16, pad=2) for val in subkey]
# ascii_string = join([Char(val) for val in subkey])
# println("\nRecovered SubKey (Hex): ", hex_values)
# println("Recovered SubKey (ASCII): ", ascii_string)

# # Plot of all entropies
# for byte in 1:16
#     plot_hypothesis(byte, entropies, subkey)
# end