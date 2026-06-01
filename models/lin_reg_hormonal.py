import pandas as pd
from sklearn.linear_model import LinearRegression
import sklearn.metrics
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
    train_x = dataset[user_num <= 80]
    test_x = dataset[user_num > 80]

    # normalize
    hormone_cols = ["estradiol (E2)", "estrone (E1)", "progesterone", "testosterone", "HCG", hormone_tplus1]
    train_x[hormone_cols] = normalize(train_x, hormone_cols)
    test_x[hormone_cols] = normalize(test_x, hormone_cols)

    # separate into x and y
    train_y = train_x.pop(hormone_tplus1)
    test_y = test_x.pop(hormone_tplus1)

    return train_x, train_y, test_x, test_y


def train(train_x, train_y):
    # implement lin reg
    model = LinearRegression()
    model.fit(train_x, train_y)
    return model

def predict(model, test_x, test_y):
    y_pred = model.predict(test_x)
    r2 = model.score(test_x, test_y)
    mae = sklearn.metrics.mean_absolute_error(test_y, y_pred)
    rmse = sklearn.metrics.root_mean_squared_error(test_y, y_pred)
    print(f"r2 = {r2}")
    # plot mae
    plt.figure(figsize=(8, 5))
    plt.bar(['MAE'], [mae], color='skyblue')
    plt.ylabel("Mean Absolute Error")
    plt.title("Mean Absolute Error (MAE) on Test Set")
    plt.show()
    # plot rmse
    plt.figure(figsize=(8, 5))
    plt.bar(['RMSE'], [rmse], color='lightpink')
    plt.ylabel("RMSE")
    plt.title("RMSE on Test Set")
    plt.show()


dataset = pd.read_csv("simulated_hormone_cycles.csv")
x_train, y_train, x_test, y_test, x_val, y_val = preprocess(dataset, 'estradiol (E2)', 'estrogen(t+1)')
model_estrogen = train(x_train, y_train)
predict(model_estrogen, x_test, y_test)

x_train, y_train, x_test, y_test, x_val, y_val = preprocess(dataset, 'progesterone', 'progesterone(t+1)')
model_progesterone = train(x_train, y_train)
predict(model_progesterone, x_test, y_test)

