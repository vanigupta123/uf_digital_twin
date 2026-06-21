import pandas as pd
import sklearn
import torch
import matplotlib.pyplot as plt
from torch import nn
import json

class Normalizer:
    def __init__(self):
        self.mean = None
        self.std = None

    def fit(self, mean, std):
        self.mean = mean
        self.std = std

def normalize(df, cols, normalizer):
    if normalizer.mean is None:
        normalizer.fit(df[cols].mean(), df[cols].std())
    return (df[cols] - normalizer.mean) / normalizer.std

def preprocess(dataset):
    dataset["pain_level"] = dataset["pain_level"].astype(float)
    dataset.drop(columns=['labels', 'downsampled_shape', 'age_group'], inplace=True)
    user_num = dataset['patient_id'].astype(int)
    train_x = dataset[user_num <= 180].copy()
    val_x = dataset[(user_num > 180) & (user_num <= 240)].copy()
    test_x = dataset[user_num > 240].copy()
    datasets = [train_x, val_x, test_x]

    for ds in datasets:
        no_fibroids = ds[ds['fibroid_present'] == True].shape[0]
        fibroids_present = ds[(ds['fibroid_present'] == True) | (ds['fibroid_present'] == True)].shape[0]
        ratio = no_fibroids / (no_fibroids + fibroids_present)
        print(f"ratio = {ratio}")

    # separate into x and y
    train_y = train_x.pop("fibroid_present")
    val_y = val_x.pop("fibroid_present")
    test_y = test_x.pop("fibroid_present")

    normalize_cols = ["cycle_length_days", "symptom_duration_months", "pain_level", "ferritin_proxy"]
    normalizer = Normalizer()
    train_x[normalize_cols] = train_x[normalize_cols].fillna(train_x[normalize_cols].median())
    val_x[normalize_cols] = val_x[normalize_cols].fillna(train_x[normalize_cols].median())
    test_x[normalize_cols] = test_x[normalize_cols].fillna(train_x[normalize_cols].median())    
    train_x.loc[:, normalize_cols] = normalize(train_x, normalize_cols, normalizer)
    val_x.loc[:, normalize_cols] = normalize(val_x, normalize_cols, normalizer)
    test_x.loc[:, normalize_cols] = normalize(test_x, normalize_cols, normalizer)

    normalize_cols.append("num_fibroids")
    normalize_cols.append("fibroid_volume_ratio")
    normalize_cols.append("ferritin_proxy")
    stats = {
        "medians": train_x[normalize_cols].median().to_dict(),
        "means": normalizer.mean.to_dict(),
        "stds": normalizer.std.to_dict(),
        "feature_columns": list(train_x.columns)
    }
    json.dump(stats, open("models/preprocessing_stats.json", "w"))
    return train_x, train_y, test_x, test_y, val_x, val_y


def train(X_train, y_train, X_val, y_val, epochs: int = 50):
    device = torch.device("mps" if torch.mps.is_available() else "cpu")

    X_train_t = to_tensor(X_train, dtype=torch.float32)
    y_train_t = torch.as_tensor(y_train.to_numpy(), dtype=torch.float32).view(-1, 1)
    X_val_t = to_tensor(X_val, dtype=torch.float32)
    y_val_t = torch.as_tensor(y_val.to_numpy(), dtype=torch.float32).view(-1, 1)
    print(X_train_t.shape[1])
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
    dataset = pd.read_csv("umd_data_categorical.csv")
    X_train, y_train, X_test, y_test, X_val, y_val = preprocess(dataset)
    print(X_train.isnull().sum())
    print(X_train.shape)
    model, history = train(X_train, y_train, X_val, y_val, epochs=50)
    torch.save(model.state_dict(), 'model.pth')
    evaluate(model, X_val, y_val, split_name="val")
    evaluate(model, X_test, y_test, split_name="test")
    plot_loss_curve(history)

main()