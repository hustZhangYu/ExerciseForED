# we plot the eigenstates of stark ladder model

using LinearAlgebra, Plots

L=100

H=zeros(L,L)
F=0.01
for i=1:L-1
    H[i,i+1]=1
    H[i+1,i]=1
    H[i,i]=F*i
end
H[L,L]=F*L

E,Ev=eigen(H)

# Select a particular eigenvector (e.g., 100th)
psi = Ev[50, :]

# Plot the probability distribution |ψ|^2
plot(1:L, psi .* conj(psi), 
    xlabel="Index (i)", 
    ylabel="|ψ(i)|²", 
    # title="Probability Distribution of Eigenstate", 
    color=:viridis,  # Color map choice, you can change it
    linewidth=2,      # Adjust line width for better visibility
    legend=false,     # Remove legend if not needed
    grid=true,        # Add grid for better readability
    xlims=(1, L),     # Set x-axis limits to match the size of the system
    ylims=(0, maximum(psi .* conj(psi)) * 1.1),  # Set y-axis limit to show the maximum value slightly above the peak
    fontsize=18,      # Set the font size for the labels and title
    tickfontsize=12   # Set the font size for the tick labels
)

# Save the figure
savefig("plot1.png")