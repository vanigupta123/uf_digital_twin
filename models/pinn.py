import pandas as pd
import sklearn.metrics
import matplotlib.pyplot as plt
import torch
from torch import nn
import torch.autograd as autograd
# mlp
# compare metrics, training time, overfitting behavior with lin reg
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
 
    # divide into train, test, validation by user id
    cols = ["patient_id","t_days","fibroid_volume_ratio","num_fibroids","patient_weight","pain_level","age_group_encoded","treatment_type","effective_rate"]
    dataset["treatment_type"] = dataset["treatment_type"].map({"none": 0, "hormonal": 1, "surgery": 2})
    user_num = dataset["patient_id"].astype(int)
    train_x = dataset[user_num <= 100].copy()
    val_x = dataset[(user_num > 100) & (user_num <= 200)].copy()
    test_x = dataset[user_num > 200].copy()
    datasets = [train_x, val_x, test_x]

    # normalize
    normalize_cols = ["t_days", "num_fibroids", "patient_weight", "pain_level", "effective_rate"]
    normalizer = Normalizer()
    train_x[normalize_cols] = train_x[normalize_cols].fillna(train_x[normalize_cols].median())
    val_x[normalize_cols] = val_x[normalize_cols].fillna(train_x[normalize_cols].median())
    test_x[normalize_cols] = test_x[normalize_cols].fillna(train_x[normalize_cols].median())   
    train_x.loc[:, normalize_cols] = normalize(train_x, normalize_cols, normalizer)
    val_x.loc[:, normalize_cols] = normalize(val_x, normalize_cols, normalizer)
    test_x.loc[:, normalize_cols] = normalize(test_x, normalize_cols, normalizer)


    train_y = train_x.pop("fibroid_volume_ratio")
    val_y = val_x.pop("fibroid_volume_ratio")
    test_y = test_x.pop("fibroid_volume_ratio")
    return train_x, train_y, test_x, test_y, val_x, val_y

def build_model(input_dim: int) -> nn.Module:
    return nn.Sequential(
        nn.Linear(input_dim, 64),
        nn.Tanh(),
        nn.Linear(64, 64), # tapering like (64, 32) -> (32, 1) makes sense for classifiers but not for PINN
        nn.Tanh(),
        nn.Linear(64, 64), # in PINN, uniform layers are better because gradient signals stay more accurate
        nn.Tanh(),
        nn.Linear(64, 1), # because we're learning a continuous function that needs to stay smooth and expressive
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
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    history = {"train_loss": [], "val_loss": []}
    best_state = None
    best_val = float("inf")

    batch_size = 32 # samples per gradient update
    n = X_train_t.shape[0] # total number of training samples
    lambda_reg = 0.5

    for epoch in range(epochs):
        model.train()
        # shuffle indices of training data in each epoch so model sees different batch order each time
        perm = torch.randperm(n)

        train_loss_epoch = 0.0
        num_batches = 0
        for start in range(0, n, batch_size): # sets up batching of training data
            idx = perm[start : start + batch_size]
            xb = X_train_t[idx].to(device)
            yb = y_train_t[idx].to(device)

            optimizer.zero_grad()
            # residual = dV_dt - (r * V * (1 - (V/K)))
            t = xb[:, 1].requires_grad_(True)
            xb_with_grad_t = xb.clone() # rebuild xb with grad-enabled t
            xb_with_grad_t[:, 1] = t
            V = model(xb_with_grad_t)
            r = xb[:, -1].unsqueeze(1)
            dV_dt = autograd.grad(V, t, grad_outputs=torch.ones_like(V), create_graph=True)[0].unsqueeze(1) # (batch_size, 1)
            residual = dV_dt - (r * V * (1 - V)) # K = 1.0 based on preprocessing
            loss = nn.MSELoss()(V, yb) + lambda_reg * torch.mean(residual ** 2) # nn.MSELoss()(residual, 0) replaced with torch.mean(residual ** 2)
            loss.backward()
            optimizer.step()

            train_loss_epoch += loss.item()
            num_batches += 1

        train_loss_epoch /= max(1, num_batches)

        # Validation
        model.eval()
        with torch.no_grad():
            val_loss_epoch = nn.MSELoss()(model(X_val_t.to(device)), y_val_t.to(device)).item()

        history["train_loss"].append(train_loss_epoch)
        history["val_loss"].append(val_loss_epoch)

        if val_loss_epoch < best_val:
            best_val = val_loss_epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)

    return model.to("cpu"), history


def evaluate(model, X_split, y_split, split_name):
    X_split_t = to_tensor(X_split, dtype=torch.float32)
    y_split_t = to_tensor(y_split, dtype=torch.float32).view(-1, 1) # should be same shape as logits

    with torch.no_grad():
        logits = model(X_split_t)

    mse = nn.MSELoss()(logits, y_split_t)
    mae = nn.L1Loss()(logits, y_split_t)
    r2 = 1 - mse / torch.var(y_split_t)

    print(f"{split_name} set metrics:")
    print(
        f"mse={mse:.4f}, mae={mae:.4f}, r2={r2:.4f}"
    )


def plot_loss_curve(history):
    plt.figure(figsize=(8, 5))
    plt.plot(history["train_loss"], label="train loss")
    plt.plot(history["val_loss"], label="val loss")
    plt.title("PINN Loss Curve")
    plt.xlabel("Epoch")
    plt.ylabel("PINN Loss")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

def main():
    dataset = pd.read_csv("pinn_fibroid_growth.csv")
    X_train, y_train, X_test, y_test, X_val, y_val = preprocess(dataset)
    model, history = train(X_train, y_train, X_val, y_val, epochs=50)
    torch.save(model, 'model.pth')
    evaluate(model, X_val, y_val, split_name="val")
    evaluate(model, X_test, y_test, split_name="test")
    plot_loss_curve(history)


if __name__ == "__main__":
    main()