import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

def normalize(df, cols):
    return (df[cols] - df[cols].mean()) / df[cols].std()

# x = current estrogen, progesterone, all other inputs
# y = estrogen(t+1)
def preprocess(dataset, hormone_t, hormone_tplus1):
    # write new column estrogen(t+1)
    dataset[hormone_tplus1] = dataset[hormone_t].shift(-1)
    # find the last row of each user id and set the estrogen(t+1) and progesterone(t+1) to the first day's value
    for user_id in dataset['user id'].unique():
        first_row = dataset[dataset['user id'] == user_id].iloc[0]
        dataset.loc[(dataset['user id'] == user_id) & (dataset['day'] == dataset['cycle length']), hormone_tplus1] = first_row[hormone_t]

    # split by numeric user id
    user_num = dataset['user id'].str.extract(r'user_(\d+)', expand=False).astype(int)
    dataset = dataset.drop(columns=['user id']) # for lin reg
    train_x = dataset[user_num <= 60]
    test_x = dataset[(user_num > 60) & (user_num <= 80)]
    val_x = dataset[user_num > 80]

    # normalize
    hormone_cols = ["estradiol (E2)", "estrone (E1)", "progesterone", "testosterone", "HCG", hormone_tplus1]
    train_x[hormone_cols] = normalize(train_x, hormone_cols)
    test_x[hormone_cols] = normalize(test_x, hormone_cols)
    val_x[hormone_cols] = normalize(val_x, hormone_cols)

    # separate into x and y
    train_y = train_x.pop(hormone_tplus1)
    test_y = test_x.pop(hormone_tplus1)
    val_y = val_x.pop(hormone_tplus1)

    return train_x, train_y, test_x, test_y, val_x, val_y


def predict(train_x, train_y):
    # implement lin reg
    model = LinearRegression()
    model.fit(train_x, train_y)
    y_pred = model.predict(train_x)
    # return y_pred

    # print loss curve, mae, rmse, r2
    plt.figure(figsize=(8,6)) 
    plt.scatter(train_x['estradiol (E2)'], train_y, color='blue', label='Data Points') 
    plt.plot(train_x['estradiol (E2)'], y_pred, color='red', linewidth=2, label='Regression Line') 
    plt.title('Linear Regression')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.grid(True)
    plt.show()

# print accuracy, precision, recall, F1, ROC-AUC

# do the same for progesterone(t+1)

dataset = pd.read_csv("simulated_hormone_cycles.csv")
x_train, y_train, x_test, y_test, x_val, y_val = preprocess(dataset, 'estradiol (E2)', 'estrogen(t+1)')
predict(x_train, y_train)

x_train, y_train, x_test, y_test, x_val, y_val = preprocess(dataset, 'progesterone', 'progesterone(t+1)')
predict(x_train, y_train)

