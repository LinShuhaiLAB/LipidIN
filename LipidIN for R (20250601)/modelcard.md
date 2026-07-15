---
language:
- en  # Assuming English is the primary language for any text processing or documentation
license: apache-2.0  # Example license; adjust according to your use
library_name: pytorch  # The model is implemented in PyTorch
tags:
- multitask-learning
- weak-supervision
- attention
- transformer
datasets:
- custom-dataset  # Since you mentioned custom CSV files, you can name your dataset here
metrics:
- mse  # Mean Squared Error, as used in the loss function

model-index:
- name: WMY_Model  # You can replace this with the actual name of your model
  results:
  - task:
      type: multitask-learning  # Specify the nature of the task, e.g., "multitask-learning"
      name: Prediction of Target Vectors  # Describe the task briefly
    dataset:
      type: custom-dataset  # Name your dataset accordingly
      name: Custom Contour Data  # Give your dataset a descriptive name
      config: default  # If there's a specific configuration, mention it; otherwise, default works
      split: train  # Specify which data split was used for evaluation
    metrics:
      - type: mse  # Mean Squared Error
        value: {metric_value}  # Replace with the actual metric value
        name: Training MSE  # Optional; describe the metric result
    source:  
      name: Custom Evaluation  # Name of the evaluation source or method
      url: {source_url}  # Provide a URL if applicable
---

## WMY Model Card

### Model Description
The WMY model is designed for multitask learning using weak supervision techniques. It employs a deep neural network architecture with multiple layers of multi-head self-attention mechanisms and feed-forward networks (FFN). The model incorporates custom layers such as `KanLikeLayer` and `ScaledFeatureExtractor`, and is optimized using a custom learning rate scheduler.

### Intended Use
This model is intended for tasks that require multitask learning and weak supervision, especially in contexts where input data is heterogeneous, involving complex features such as specimen sampling site information, cytological images, or other structured data.

### Training Data
The model was trained on a custom dataset that includes input features and corresponding target vectors, represented as CSV files. The dataset consists of coordinate data and labels extracted from these files.

### Training Procedure
The training process involves 3000 epochs with early stopping based on a loss threshold (MSE < 40). The model is trained using the Adam optimizer with a learning rate schedule defined by the `CustomLRScheduler` class, which adjusts the learning rate based on the loss function.

### Evaluation Results
The model’s performance was evaluated using Mean Squared Error (MSE) on a validation dataset. The best model was saved based on the lowest validation loss achieved during training.

### Limitations and Future Work
While the model shows promising results in multitask learning contexts, it has limitations in terms of scalability and generalization to other types of data. Future work could involve fine-tuning the model on a broader range of datasets and exploring additional regularization techniques.

### License
The model is released under the Apache 2.0 license.

### Contact Information
For any questions or issues, please contact the developers at [your_contact_info].
