# In this code ,we plot the IPR of the AA model

using LinearAlgebra, Plots

# parameters
L = 987
omega = (sqrt(5) - 1) / 2
phi = pi / 3
lambda_all = 0:0.05:4
Data = []

for lambda in lambda_all
    H = zeros(L, L)

    # Construct Hamiltonian matrix
    for i = 1:L - 1
        H[i, i + 1] = 1
        H[i + 1, i] = 1      
        H[i, i] = lambda * cos(2 * π * omega * i + phi)
    end
    H[L, L] = lambda * cos(2 * π * omega * L + phi)

    # Eigenvalues and eigenvectors
    E, Ev = eigen(H)

    # Calculate IPR (Inverse Participation Ratio)
    Ipr_all = 0
    for m = 1:L
        psi = Ev[:, m]
        Ipr_all += sum((psi .* conj(psi)) .* (psi .* conj(psi)))
    end
    Ipr_all /= L
    push!(Data, Ipr_all)
end

# Plot the results
plot(lambda_all, Data, xlabel="λ ", ylabel="IPR ", 
     title="IPR vs λ", legend=false, lw=2, color=:blue)

# Save the plot to a file
savefig("IPR_vs_lambda.png")  # Save as PNG file

