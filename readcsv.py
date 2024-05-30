import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file into a DataFrame
file_path = 's817.csv'  # Replace with your actual file path
data = pd.read_csv(file_path)

# Plot each column
plt.figure(figsize=(10, 6))
for column in data.columns[1:]:  # Skip the first column if it's an index or identifier
    plt.plot(data.iloc[:, 0], data[column], label=column)

# Add labels and title
plt.xlabel(data.columns[0])  # X-axis label using the first column's header
plt.ylabel('Values')
plt.title('CSV Data Plot')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
