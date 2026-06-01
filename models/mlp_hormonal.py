import pandas as pd
import sklearn.metrics
import matplotlib.pyplot as plt
import torch
from torch import nn

# mlp
# compare metrics, training time, overfitting behavior with lin reg

def normalize(df, cols):
    return (df[cols] - df[cols].mean()) / df[cols].std()

# x = current estrogen, progesterone, all other inputs
# y = diagnosis (normal or abnormal)
def preprocess(dataset):
    # no class imbalance -- 48 / 52 normal:abnormal
    # normal = dataset[dataset['diagnosis_normal'] == True]
    # abnormal = dataset[(dataset['diagnosis_PCOS'] == True) | (dataset['diagnosis_fibroids'] == True)]
 
    # divide into train, test, validation by user id
    user_num = dataset['user id'].str.extract(r'user_(\d+)', expand=False).astype(int)
    train_x = dataset[user_num <= 60].copy()
    val_x = dataset[(user_num > 60) & (user_num <= 80)].copy()
    test_x = dataset[user_num > 80].copy()
    datasets = [train_x, val_x, test_x]

    for ds in datasets:
        normal_dim = ds[ds['diagnosis_normal'] == True].shape[0]
        abnormal_dim = ds[(ds['diagnosis_PCOS'] == True) | (ds['diagnosis_fibroids'] == True)].shape[0]
        ratio = normal_dim / (normal_dim + abnormal_dim)
        if ratio <= 0.35 or ratio >= 0.65:
            print(f"ratio = {ratio}")

    # normalize
    hormone_cols = ["estradiol (E2)", "estrone (E1)", "progesterone", "testosterone", "HCG"]
    train_x.loc[:, hormone_cols] = normalize(train_x, hormone_cols)
    val_x.loc[:, hormone_cols] = normalize(val_x, hormone_cols)
    test_x.loc[:, hormone_cols] = normalize(test_x, hormone_cols)

    # separate into x and y: 1 = abnormal (PCOS or fibroids), 0 = normal
    def diagnosis_label(df):
        return ((df['diagnosis_PCOS']) | (df['diagnosis_fibroids'])).astype(int)

    train_y = diagnosis_label(train_x)
    val_y = diagnosis_label(val_x)
    test_y = diagnosis_label(test_x)

    drop_cols = [
        'diagnosis_PCOS', 'diagnosis_fibroids', 'diagnosis_normal', 'diagnosis_nan',
        'user id',
    ]
    train_x = train_x.drop(columns=drop_cols)
    val_x = val_x.drop(columns=drop_cols)
    test_x = test_x.drop(columns=drop_cols)

    return train_x, train_y, test_x, test_y, val_x, val_y

def build_model(input_dim: int) -> nn.Module:
    # Simple tabular MLP.
    # Using BCEWithLogitsLoss means the output layer does NOT need a sigmoid.
    return nn.Sequential(
        nn.Linear(input_dim, 64),
        nn.ReLU(),
        nn.Linear(64, 32),
        nn.ReLU(),
        nn.Linear(32, 1),
    )


def to_tensor(df, dtype=torch.float32):
    # DataFrames can contain a mix of bool/ints/floats; force a single numeric dtype
    # so torch doesn't see an object array from mixed columns.
    return torch.as_tensor(df.to_numpy(dtype="float32"), dtype=dtype)


def train(X_train, y_train, X_val, y_val, epochs: int = 50):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    X_train_t = to_tensor(X_train, dtype=torch.float32)
    y_train_t = torch.as_tensor(y_train.to_numpy(), dtype=torch.float32).view(-1, 1)
    X_val_t = to_tensor(X_val, dtype=torch.float32)
    y_val_t = torch.as_tensor(y_val.to_numpy(), dtype=torch.float32).view(-1, 1)

    model = build_model(input_dim=X_train_t.shape[1]).to(device)

    criterion = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    history = {"train_loss": [], "val_loss": []}
    best_state = None
    best_val = float("inf")

    batch_size = 32
    n = X_train_t.shape[0]
    for epoch in range(epochs):
        model.train()
        # Shuffle indices each epoch.
        perm = torch.randperm(n)

        train_loss_epoch = 0.0
        num_batches = 0

        for start in range(0, n, batch_size):
            idx = perm[start : start + batch_size]
            xb = X_train_t[idx].to(device)
            yb = y_train_t[idx].to(device)

            optimizer.zero_grad()
            logits = model(xb)
            loss = criterion(logits, yb)
            loss.backward()
            optimizer.step()

            train_loss_epoch += loss.item()
            num_batches += 1

        train_loss_epoch /= max(1, num_batches)

        # Validation
        model.eval()
        with torch.no_grad():
            val_logits = model(X_val_t.to(device))
            val_loss_epoch = criterion(val_logits, y_val_t.to(device)).item()

        history["train_loss"].append(train_loss_epoch)
        history["val_loss"].append(val_loss_epoch)

        if val_loss_epoch < best_val:
            best_val = val_loss_epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)

    return model.to("cpu"), history


def evaluate(model, X_split, y_split, split_name="test"):
    model.eval()
    X_split_t = to_tensor(X_split, dtype=torch.float32)
    y_split_np = y_split.to_numpy(dtype="int32").ravel()

    with torch.no_grad():
        logits = model(X_split_t)
        y_score = torch.sigmoid(logits).cpu().numpy().ravel()

    y_pred = (y_score >= 0.5).astype(int)

    accuracy = sklearn.metrics.accuracy_score(y_split_np, y_pred)
    precision = sklearn.metrics.precision_score(y_split_np, y_pred, zero_division=0)
    recall = sklearn.metrics.recall_score(y_split_np, y_pred, zero_division=0)
    f1 = sklearn.metrics.f1_score(y_split_np, y_pred, zero_division=0)
    roc_auc = sklearn.metrics.roc_auc_score(y_split_np, y_score)

    print(f"{split_name} set metrics:")
    print(
        f"  accuracy={accuracy:.4f}, precision={precision:.4f}, "
        f"recall={recall:.4f}, f1={f1:.4f}, roc_auc={roc_auc:.4f}"
    )
    print(
        f"{split_name} confusion matrix:\n{sklearn.metrics.confusion_matrix(y_split_np, y_pred)}"
    )


def plot_loss_curve(history):
    plt.figure(figsize=(8, 5))
    plt.plot(history["train_loss"], label="train loss")
    plt.plot(history["val_loss"], label="val loss")
    plt.title("MLP Loss Curve")
    plt.xlabel("Epoch")
    plt.ylabel("BCEWithLogitsLoss")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()


def main():
    dataset = pd.read_csv("simulated_hormone_cycles.csv")
    X_train, y_train, X_test, y_test, X_val, y_val = preprocess(dataset)
    model, history = train(X_train, y_train, X_val, y_val, epochs=50)
    torch.save(model, 'model.pth')
    evaluate(model, X_val, y_val, split_name="val")
    evaluate(model, X_test, y_test, split_name="test")
    plot_loss_curve(history)


if __name__ == "__main__":
    main()