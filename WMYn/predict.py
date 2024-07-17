import os
import shutil
import torch
import numpy as np
import pandas as pd
from model import WMY, preprocess_data, load_data
from pathlib import Path

def parse_parameters_from_filename(filename):
    parameters = filename.split('_')[0].split()
    try:
        input_dim = int(parameters[1])
        output_dim = int(parameters[2])
        return input_dim, output_dim
    except (IndexError, ValueError):
        print("Invalid parameters found in filename.")
        return None, None

def WYMn_predict(csv_add, weight_name):
    train_csv_file = csv_add
    labels, train_data = load_data(train_csv_file)

    input_dim = train_data.shape[1]
    output_dim = 2468
    hidden_dim = 512  
    batch_size = len(labels)
    model = WMY(input_dim, output_dim, num_layers=6, num_heads=8, hidden_dim=hidden_dim, dropout=0.1, batch_size=batch_size)
    model.load_state_dict(torch.load(weight_name))

    train_data_tensor = preprocess_data(train_data)

    with torch.no_grad():
        prediction = model(train_data_tensor)
        prediction = prediction.squeeze(0).numpy()  

    return prediction

def normalize_columns(column):
    col_max = np.max(column)
    if col_max == 0:
        return np.zeros_like(column)
    normalized_column = column / col_max
    normalized_column[normalized_column < 0] = 0
    return normalized_column

def process_files(data_folder, weight_folder, gt_folder, prediction_folder, mode):
    no_weight_folder = os.path.join(data_folder, "no_weights")

    os.makedirs(prediction_folder, exist_ok=True)
    os.makedirs(no_weight_folder, exist_ok=True)

    for data_filename in os.listdir(data_folder):
        if data_filename.endswith(".csv") and "_GT" not in data_filename:
            weight_filename = f"{os.path.splitext(data_filename)[0]}_weights.pth"
            gt_filename = f"{os.path.splitext(data_filename)[0]}_GT.csv"

            data_filepath = os.path.join(data_folder, data_filename)
            weight_filepath = os.path.join(weight_folder, weight_filename)
            gt_filepath = os.path.join(gt_folder, gt_filename)

            if os.path.exists(weight_filepath):
                prediction = WYMn_predict(data_filepath, weight_filepath)
                if prediction is not None:
                    normalized_prediction = normalize_columns(prediction)

                    if os.path.exists(gt_filepath):
                        try:
                            gt_df = pd.read_csv(gt_filepath)
                            if len(gt_df) != len(normalized_prediction):
                                print(f"Mismatch in number of rows between GT file and prediction for {data_filename}.")
                                continue
                            gt_df.iloc[:, 1] = normalized_prediction
                            combined_data = gt_df

                            output_filename = os.path.join(prediction_folder, f"{os.path.splitext(data_filename)[0]}_prediction_out.csv")
                            combined_data.to_csv(output_filename, index=False)
                            print(f"Prediction saved to {output_filename}")

                        except ValueError as e:
                            print(f"Error reading GT file {gt_filepath}: {e}")
                            continue
                    else:
                        print(f"GT file not found for {data_filename}. Prediction could not be saved.")
            else:
                shutil.copy(data_filepath, os.path.join(no_weight_folder, data_filename))
                print(f"Weight file not found for {data_filename}. Please use the training module.")

if __name__ == "__main__":
    data_path =  " "
    project_folder = Path(" ")
    mode = "neg"  # "neg" or "pos"

    if mode not in ["neg", "pos"]:
        print("Invalid mode. Please choose 'neg' or 'pos'.")
        exit(1)

    weight_folder = os.path.join(project_folder, mode)
    if mode == "neg":
        gt_folder = os.path.join(project_folder, "neg", "neg_gt")
    else:
        gt_folder = os.path.join(project_folder, "pos", "pos_gt")

    prediction_folder = os.path.join(project_folder, "predictions")

    if os.path.isfile(data_path):
        data_folder = os.path.dirname(data_path)
        process_files(data_folder, weight_folder, gt_folder, prediction_folder, mode)
    elif os.path.isdir(data_path):
        process_files(data_path, weight_folder, gt_folder, prediction_folder, mode)
    else:
        print("Invalid data path. Please provide a valid file or directory.")
