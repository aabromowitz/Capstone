# Imports
from darts.models import TFTModel
from darts.metrics import rmse
# import torch

def train_tft(train, val, test, exog_train, exog_val, exog_test):
    """Trains a TFT model and returns validation/test loss and variable importance."""

    # Define model
    model = TFTModel(
        input_chunk_length=12,  
        output_chunk_length=6,  
        hidden_size=32,
        lstm_layers=1,
        num_attention_heads=4,
        dropout=0.1,
        batch_size=32,
        n_epochs=50,
        add_relative_index=True,
        add_encoders={"cyclic": {"past": ["month"]}},  # Example encoder
        likelihood=None,  # For probabilistic forecasting, set a likelihood
        optimizer_kwargs={"lr": 0.001},
        random_state=42,
    )

    # Train model
    model.fit(
        series=train,
        past_covariates=exog_train,
        val_series=val,
        val_past_covariates=exog_val,
        verbose=True,
    )

    # Evaluate model
    val_pred = model.predict(n=len(val), past_covariates=exog_val)
    test_pred = model.predict(n=len(test), past_covariates=exog_test)
    
    val_loss = rmse(val, val_pred)
    test_loss = rmse(test, test_pred)

    # Get variable importance
    importance = model.interpretation()["variable_importance"]

    return model, val_loss, test_loss, importance