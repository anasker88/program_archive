import gc
import os
from functools import partial

import pandas as pd
import plotly.express as px
import torch
from datasets import load_dataset
from sae_lens import SAE, HookedSAETransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from tqdm import tqdm
from transformer_lens.utils import tokenize_and_concatenate

# 1. モデルをロードする
# デフォルトでは、GPT-2のような事前トレーニング済みのTransformerモデルをロードできます。

device = (
    "cuda"
    if torch.cuda.is_available()
    else "mps" if torch.mps.is_available() else "cpu"
)
model = HookedSAETransformer.from_pretrained("gpt2-small")

sae, cfg_dict, sparsity = SAE.from_pretrained(
    release="gpt2-small-res-jb",  # <- Release name
    sae_id="blocks.7.hook_resid_pre",  # <- SAE id (not always a hook point!)
    device=device,
)

sae.W_dec


class modifier:
    def __init__(self, model, sae, cfg_dict, replacement_values):
        self.model = model
        self.sae = sae
        self.cfg_dict = cfg_dict
        self.replacement_values = replacement_values
        self.counter = 0

    def modify_specific_feature(self, activations, hook, target_feature_ids):
        """
        activations: 元のアクティベーション（shape: (batch_size, seq_len, hidden_size)）
        hook: Transformerlensのフック情報
        target_feature_ids: 変更したいfeatureのインデックスのリスト
        """
        # activationsの対象のfeatureだけを変更
        modified_activations = activations.clone()
        for target_feature_id in target_feature_ids:
            modified_activations[:, :, target_feature_id] = self.replacement_values[
                0, self.counter
            ]
        self.counter += 1
        return modified_activations


def generate_with_feature_modification(
    model, sae, prompt, target_feature_idx, replacement_values
):
    input_ids = model.to_tokens(prompt, prepend_bos=sae.cfg.prepend_bos)
    # initialize modifier
    mod = modifier(model, sae, sae.cfg, replacement_values)


# テスト例
prompt = "Once upon a time"

# 特定のfeatureインデックス（例: 20115）
target_feature_idx = 20

# 各トークンごとに指定する値を作成 (batch_size=1, seq_len=10)
replacement_values = (
    torch.linspace(0.1, 1.0, steps=10).unsqueeze(0).to(model.cfg.device)
)

# 生成
steered_text = generate_with_feature_modification(
    model, sae, prompt, target_feature_idx, replacement_values
)
print("Text with modified feature activations:")
print(steered_text)


# ramdom forestによる2クラス分類
def run_random_forest(st_act, anti_st_act):
    """
    st_act: stereotypical promptのアクティベーション
    anti_st_act: anti-stereotypical promptのアクティベーション
    """
    # st_actとanti_st_actを結合
    activations = torch.cat([st_act, anti_st_act], dim=0)
    # ラベル
    labels = torch.cat([torch.ones(st_act.shape[0]), torch.zeros(anti_st_act.shape[0])])
    # ランダムフォレスト
    clf = RandomForestClassifier()
    clf.fit(activations.cpu().numpy(), labels.cpu().numpy())
    # 正答率
    accuracy = clf.score(activations.cpu().numpy(), labels.cpu().numpy())
    # 特徴量の重要度
    feature_importance = clf.feature_importances_

    return accuracy, feature_importance


# SVMによる2クラス分類
def run_svm(st_act, anti_st_act):
    """
    st_act: stereotypical promptのアクティベーション
    anti_st_act: anti-stereotypical promptのアクティベーション
    return: 正答率,法線ベクトル
    """
    # st_actとanti_st_actを結合
    activations = torch.cat([st_act, anti_st_act], dim=0)
    # ラベル
    labels = torch.cat([torch.ones(st_act.shape[0]), torch.zeros(anti_st_act.shape[0])])
    # SVM
    clf = SVC(kernel="linear")
    clf.fit(activations.cpu().numpy(), labels.cpu().numpy())
    # 正答率
    accuracy = clf.score(activations.cpu().numpy(), labels.cpu().numpy())
    # 法線ベクトル
    normal_vector = clf.coef_

    return accuracy, normal_vector
