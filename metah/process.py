import os
import pandas as pd

directory = 'C:/Users/ASUS/Desktop/SO/projeto-so-mei-2/metah' 

def process_files(directory):
    file_means = []

    for filename in os.listdir(directory):
        if filename.startswith("GA_results"):
            file_path = os.path.join(directory, filename)
            data = pd.read_csv(file_path, delim_whitespace=True, header=None)
            averages = data.mean(axis=0)
            overall_mean = averages.mean() 
            file_means.append((filename, overall_mean, averages))

            # print(f"Filename: {filename}")
            # for col in range(len(averages)):
                # print(f"Average of column {col + 1}: {averages[col]:.2f}")
            # print(f"Overall mean: {overall_mean:.2f}")
            #print('---------------------------------------')
    
    # Sort files by their overall mean values
    file_means.sort(key=lambda x: x[1])
    
    # Sort files by their overall mean values
    for i in range(len(file_means)):
        filename, overall_mean, averages = file_means[i]
        print(f"Filename: {filename}")
        print(f"Smallest mean value: {overall_mean:.2f}")
        for col in range(len(averages)):
            print(f"Average of column {col + 1}: {averages[col]:.2f}")
        print('---------------------------------------')

# Call the function
process_files(directory)






