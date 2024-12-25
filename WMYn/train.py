import os
import glob
import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import random
import numpy as np
import torch.nn.functional as F

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

class LearnableActivation(nn.Module):
    def __init__(self):
        super(LearnableActivation, self).__init__()
        self.alpha = nn.Parameter(torch.tensor(1.0))  
        self.beta = nn.Parameter(torch.tensor(0.0))  

    def forward(self, x):
        return self.alpha * (x * torch.sigmoid(x)) + self.beta 

class KanLikeLayer(nn.Module):
    def __init__(self, in_features, out_features, bias=True):
        super(KanLikeLayer, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.learnable_activation = LearnableActivation()
        self.weight = nn.Parameter(torch.Tensor(out_features, in_features))
        if bias:
            self.bias = nn.Parameter(torch.Tensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        nn.init.kaiming_uniform_(self.weight, a=nn.init.calculate_gain('sigmoid'))
        if self.bias is not None:
            fan_in, _ = nn.init._calculate_fan_in_and_fan_out(self.weight)
            bound = 1 / np.sqrt(fan_in)
            nn.init.uniform_(self.bias, -bound, bound)

    def forward(self, x):
        linear_output = F.linear(x, self.weight, self.bias)
        return self.learnable_activation(linear_output)

class ScaledFeatureExtractor(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(ScaledFeatureExtractor, self).__init__()
        self.fc1 = KanLikeLayer(input_dim, hidden_dim)
        self.fc2 = KanLikeLayer(hidden_dim, output_dim)

    def forward(self, x):
        x = self.fc1(x)
        x = self.fc2(x)
        return F.relu(x)


class WMY(nn.Module):
    def __init__(self, input_dim, output_dim, num_layers, num_heads, hidden_dim, dropout, batch_size):
        super(WMY, self).__init__()
        self.mlp = ScaledFeatureExtractor(input_dim, hidden_dim, hidden_dim)
        self.encoder = Encoder(hidden_dim, num_layers, num_heads, dropout)
        self.output_layer = KanLikeLayer(hidden_dim, output_dim)
        self.mlp_output = ScaledFeatureExtractor(batch_size, hidden_dim, 1)

    def forward(self, x):
        x = self.mlp(x)
        x = self.encoder(x)
        x = self.output_layer(x)
        x = x.squeeze(1)
        x = x.transpose(0, 1)
        x = self.mlp_output(x)
        x = x.transpose(0, 1)
        return x


class Encoder(nn.Module):
    def __init__(self, input_dim, num_layers, num_heads, dropout):
        super(Encoder, self).__init__()
        self.layers = nn.ModuleList([EncoderLayer(input_dim, num_heads, dropout) for _ in range(num_layers)])
        self.norm = nn.LayerNorm(input_dim)

    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
        x = self.norm(x)
        return x


class EncoderLayer(nn.Module):
    def __init__(self, input_dim, num_heads, dropout):
        super(EncoderLayer, self).__init__()
        self.self_attn = nn.MultiheadAttention(input_dim, num_heads, dropout=dropout)
        self.dropout1 = nn.Dropout(dropout)
        self.norm1 = nn.LayerNorm(input_dim)
        self.fc = KanLikeLayer(input_dim, input_dim)
        self.relu = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(input_dim)

    def forward(self, x):
        x = x.permute(1, 0, 2)
        attn_output, _ = self.self_attn(x, x, x)
        x = x + self.dropout1(attn_output)
        x = self.norm1(x)

        x = self.fc(x)
        x = self.relu(x)
        x = self.dropout2(x)
        x = x + self.norm2(x)

        return x


def load_data(csv_file):
    df = pd.read_csv(csv_file)
    labels = df.columns[1:].tolist()
    data = df.iloc[:, 1:].values
    data = data.T
    return labels, data


def load_target_vector(csv_file):
    df = pd.read_csv(csv_file, header=None, skiprows=1)
    target_vector = df.iloc[:, 1:].values.squeeze()
    return target_vector


def preprocess_data(data):
    data_tensor = torch.tensor(data, dtype=torch.float32)
    data_tensor = data_tensor.unsqueeze(1)
    return data_tensor


def set_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    random.seed(seed)
    np.random.seed(seed)


class CustomLRScheduler:
    def __init__(self, optimizer, lr_steps, patience=10):
        self.optimizer = optimizer
        self.lr_steps = lr_steps
        self.step_index = 0
        self.patience = patience
        self.counter = 0
        self.best_loss = float('inf')

    def step(self, current_loss):
        if current_loss < self.best_loss:
            self.best_loss = current_loss
            self.counter = 0
        else:
            self.counter += 1

        if self.counter >= self.patience:
            if self.step_index < len(self.lr_steps) - 1:
                self.step_index += 1
                new_lr = self.lr_steps[self.step_index]
                for param_group in self.optimizer.param_groups:
                    param_group['lr'] = new_lr
                self.counter = 0

    def get_last_lr(self):
        return [group['lr'] for group in self.optimizer.param_groups]


def main(csv_add, target_vec, weight_name, final_output_name):
    set_seed(42)

    train_csv_file = csv_add
    labels, train_data = load_data(train_csv_file)
    target_vector_csv = target_vec
    target_vector = load_target_vector(target_vector_csv)

    input_dim = train_data.shape[1]
    output_dim = len(target_vector)
    hidden_dim = 512
    batch_size = len(labels)

    model = WMY(input_dim, output_dim, num_layers=6, num_heads=8, hidden_dim=hidden_dim, dropout=0.1,
                batch_size=batch_size)

    train_data_tensor = preprocess_data(train_data)
    target_vector_tensor = torch.tensor(target_vector, dtype=torch.float32).unsqueeze(0)

    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.01)

    lr_steps = [0.01, 0.01, 0.001, 0.001, 0.0001]
    scheduler = CustomLRScheduler(optimizer, lr_steps, patience=500)

    num_epochs = 3000
    loss_threshold = 40
    best_loss = float('inf')

    for epoch in range(num_epochs):
        model.train()
        optimizer.zero_grad()
        output = model(train_data_tensor)
        loss = criterion(output, target_vector_tensor)
        loss.backward()
        optimizer.step()

        scheduler.step(loss.item())

        if loss.item() < best_loss:
            best_loss = loss.item()
            torch.save(model.state_dict(), weight_name)

        if loss.item() < loss_threshold:
            break

    best_model = WMY(input_dim, output_dim, num_layers=6, num_heads=8, hidden_dim=hidden_dim, dropout=0.1,
                     batch_size=batch_size)
    best_model.load_state_dict(torch.load(weight_name))

    best_model.eval()
    with torch.no_grad():
        prediction = best_model(train_data_tensor)
        prediction = prediction.squeeze(0).numpy()

    gt_df = pd.read_csv(target_vector_csv, header=None)
    prediction_df = pd.DataFrame(prediction)

    gt_df.iloc[1:, 1] = prediction_df.iloc[:, 0].clip(lower=0)
    gt_df.to_csv(final_output_name, index=False, header=False)


def batch_process(data_folder, project_folder):
    os.makedirs(os.path.join(project_folder, "result"), exist_ok=True)
    os.makedirs(os.path.join(project_folder, "weights"), exist_ok=True)
    csv_files = glob.glob(os.path.join(data_folder, "*.csv"))

    for csv_file in csv_files:
        file_name = os.path.basename(csv_file)
        if "_GT" in file_name or "prediction" in file_name:
            continue

        target_file = os.path.join(data_folder, file_name.replace(".csv", "_GT.csv"))

        if not os.path.exists(target_file):
            continue

        weight_name = os.path.join(project_folder, "weights", f"{os.path.splitext(file_name)[0]}_weights.pth")
        final_output_name = os.path.join(project_folder, "result", f"{os.path.splitext(file_name)[0]}_prediction.csv")

        main(csv_file, target_file, weight_name, final_output_name)


if __name__ == "__main__":
    data_folder = r" "  
    output_folder = r" " 
    batch_process(data_folder, output_folder)
