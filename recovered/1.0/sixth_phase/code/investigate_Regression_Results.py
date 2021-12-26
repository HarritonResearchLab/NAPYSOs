#Import(s)
import pandas as pd
import matplotlib.pyplot as plt

#Action
results_file = '/home/thaddaeus/FMU/HRL/LAH/fifth_phase/work/cmd_regression_results.csv'

df = pd.read_csv(results_file)

slopes = list(df['ols_slope'])
r_sqs = list(df['ols_r_sq'])

good_slopes = []
good_r_sqs = []

for item in r_sqs:
    if item > 0.8:
        item_index = r_sqs.index(item)

        good_r_sqs.append(item)

        good_slopes.append(slopes[item_index])

#Plot
plt.hist(good_slopes,color='cornflowerblue')
plt.xlabel('Slope (r^2 > 0.8)')
plt.ylabel('Frequency')
plt.show()