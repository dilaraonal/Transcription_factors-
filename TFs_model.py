import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# file read using pandas
TF_binding_data = pd.read_csv("TFbindings.txt", delimiter="\t")

# The mean and standard deviation of expression levels were calculated
Expression_level = TF_binding_data['Expression Level']
Mean_expression = Expression_level.mean()
std_dev_expression = Expression_level.std()

print(f"Mean of Expression Levels: {Mean_expression}")
print(f"Standard Deviation of Expression Levels: {std_dev_expression}")

#Total binding score was calculated, its mean and standard deviation were calculated
TF_binding_data['Total_Binding_Score'] = TF_binding_data.iloc[:, 1:-1].sum(axis=1)
print(TF_binding_data[['Gene', 'Total_Binding_Score']])

#The mean and standard deviation of total binding scores were calculated
total_binding_scores = TF_binding_data['Total_Binding_Score'] # Binding scores were assigned to a variable for ease of use in future operations
mean_total_binding = total_binding_scores.mean()
std_total_binding = total_binding_scores.std()

print(f"Mean of Total Binding Score: {mean_total_binding}")
print(f"Standard Deviation of Total Binding Score: {std_total_binding}")

#define a function to obtain histograms
def plot_histograms():
    
    plt.figure(figsize=(12,6)) # created a figure for the plot
    #Two plots were drawn in the same figure
    plt.subplot(121) # 1 is the number of columns, 2 is the number of rows, 1 at the end represents the first part of the divided figure.
    plt.hist(Expression_level,15, facecolor='blue', edgecolor='black', alpha=0.7) 
    plt.title(('Histogram of Expression Levels'))
    plt.xlabel('Expression Level')
    plt.ylabel('Frequency')
    #Mean and standard deviation information was added to the histogram
    plt.text(0.95, 0.95, f'Mean: {Mean_expression:.2f}\nStd Dev: {std_dev_expression:.2f}', transform=plt.gca().transAxes,fontsize=10, va='top', ha='right') 
    plt.grid(True) #added grid view on chart
    
    plt.subplot(122) #2nd plot in the created figure
    plt.hist(total_binding_scores,15, facecolor='green', edgecolor='black', alpha=0.7)
    plt.title(('Histogram of Total Binding Scores'))
    plt.xlabel('Total Binding Score')
    plt.ylabel('Frequency')
    #Mean and standard deviation information was added to the histogram
    plt.text(0.95, 0.95, f'Mean: {mean_total_binding:.2f}\nStd Dev: {std_total_binding:.2f}',transform=plt.gca().transAxes, fontsize=10, va='top', ha='right')
    plt.grid(True)
    
    plt.show()



# a plot was created to show the relationship between binding scores and expression levels by creating a linear equation
def linear_equation(TF_data):
    
    #A linear regression model was applied between binding scores and expression levels using the polyfit function
    regression_model = np.polyfit(total_binding_scores, Expression_level, 1) 
    #katsayıları belirlendi
    b = regression_model[0] 
    a = regression_model[1]

    #Estimates are made for each point with coefficients
    TF_data['Predicted_Expression'] = a + b * total_binding_scores
    # column is added for errors between actual values and prediction
    TF_data['Error'] = Expression_level - TF_data['Predicted_Expression']
    print(TF_data[['Gene', 'Expression Level', 'Predicted_Expression', 'Error']])

    plt.figure(figsize=(10, 6))
    #Scatter plot was created using binding scores and expression levels
    plt.scatter(total_binding_scores, Expression_level, color='blue', label='True Expression', marker='o')
    plt.plot(total_binding_scores, TF_data['Predicted_Expression'], color='red', label='Predictions', linewidth=2)
    plt.scatter(total_binding_scores, TF_data['Error'], color='green', label='Errors', marker='*', s=50)
    
    #linear regression equation created and added to the plot
    equation = f'Expression Level = {a:.2f} + {b:.2f} * Total Binding Score'
    plt.text(0.5, 0.95, equation, transform=plt.gca().transAxes, fontsize=10, va='top', ha='center')
    #Plot was created with the necessary titles
    plt.xlabel('Total Binding Score')
    plt.ylabel('Expression Level / Prediction Error')
    plt.title('Expression Level vs. Total Binding Score')
    plt.grid(True)
    plt.legend()

    plt.show()

plot_histograms()
linear_equation(TF_binding_data)
