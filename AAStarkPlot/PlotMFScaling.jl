# In this code ,we plot the scaling of the IPR in the fractal dimension
using LinearAlgebra, Plots
using LsqFit  # 引入拟合库

# Parameters
L_all = [144, 233, 377, 610, 987, 1597]
omega = (sqrt(5) - 1) / 2
phi = π / 3
lambda = 2
Data = []

# Calculate IPR and normalized values
for L in L_all
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
    push!(Data, -log(Ipr_all) / log(L))
end

# Prepare x and y for fitting
x = [1 / log(L) for L in L_all]
y = Data

# Define linear model
model(x, p) = p[1] .+ p[2] .* x  # p[1] 是截距，p[2] 是斜率

# Perform the fit
initial_params = [0.0, 1.0]  # 初始参数猜测
fit_result = curve_fit(model, x, y, initial_params)
fitted_y = model(x, fit_result.param)

# Plotting
plot(x, y, marker=:o, label="Data", xlabel="1/log(L)", ylabel="Normalized IPR", title="Scaling of IPR", grid=true,ylim=(0, 1))
plot!(x, fitted_y, label="Linear Fit", color=:red, linestyle=:dash)

# Extend x-axis to 0
x_extended = [0.0; x]
fitted_y_extended = model(x_extended, fit_result.param)
plot!(x_extended, fitted_y_extended, label="Extrapolated Fit", color=:blue)


# Save the plot
savefig("Scaling_Fit.png")

