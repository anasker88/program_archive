import gc
import os
from typing import Tuple

import pandas as pd
import plotly.express as px
import torch
from datasets import load_dataset
from sae_lens import SAE, HookedSAETransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from tqdm import tqdm
from transformer_lens.utils import tokenize_and_concatenate


def find_token_combinations(tokens, target_word):
    """Given a list of tokens and a target word,
    find combinations of consecutive tokens that form the target word.

    Args:
        tokens (list): List of tokens.
        target_word (str): Target word to find in the list of tokens.

    Returns:
        list of tuples: List of tuples where each tuple contains the start and end index of the tokens forming the target word.
    """
    results = []
    for start_idx in range(len(tokens)):
        # skip empty tokens
        if not tokens[start_idx]:
            continue
        combined_word = ""
        for end_idx in range(start_idx, len(tokens)):
            combined_word += tokens[end_idx]
            # Check if the combined word matches the target word
            if combined_word == target_word:
                results.append((start_idx, end_idx))
            # Break if the combined word is longer than the target word
            if len(combined_word) > len(target_word):
                break
    return results


def normalize_and_find_target(tokens, target_word):
    """Normalize tokens and find the target word in the list of tokens.

    Args:
        tokens (list): List of tokens.
        target_word (str): Target word to find in the list of tokens.

    Returns:
        list of tuples: List of tuples where each tuple contains the start and end index of the tokens forming the target word.
    """
    # Normalize tokens
    tokens = [token.replace(" ", "").lower() for token in tokens]
    # Find target word in the list of tokens
    result = find_token_combinations(tokens, target_word)
    if len(result) > 0:
        return result
    # If target word is not found, try with 's' appended
    result = find_token_combinations(tokens, target_word + "s")
    if len(result) > 0:
        return result
    # If target word is not found, try with 't' appended
    result = find_token_combinations(tokens, target_word + "t")
    if len(result) > 0:
        return result
    print(f"Target word '{target_word}' not found in the tokens.")
    return [(0, 0)]


def get_activations(
    model: HookedSAETransformer,
    sae: SAE,
    prompts: list[str],
    target_index: list[tuple[int, int]],
    width: int = 5,
) -> torch.Tensor:
    """
    Args:
        model: Model
        sae: SAE
        prompts: Prompts
        target_index: Target index
    Returns:
        activations: Max activations of features at target index(shape: (batch_size, feature_size))
    """
    activations = torch.zeros(len(prompts), sae.cfg.d_sae, device=device)
    for i in range((len(prompts) + 4) // width):
        start_index = width * i
        end_index = min(width * (i + 1), len(prompts))
        _, cache = model.run_with_cache_with_saes(
            prompts[start_index:end_index],
            saes=[sae],
            # stop_at_layer=sae.cfg.hook_layer + 1,
            # names_filter=[sae.cfg.hook_name + ".hook_sae_acts_post"],
        )
        cache = cache[sae.cfg.hook_name + ".hook_sae_acts_post"]
        # get the max activations of features among target indices
        for i in range(start_index, end_index):
            for j in range(target_index[i][0], target_index[i][1] + 1):
                activations[i] = torch.max(activations[i], cache[i - start_index, j, :])
        del cache
        torch.cuda.empty_cache()
    return activations


def remove_same_pairs(
    st_act: torch.Tensor, anti_st_act: torch.Tensor
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Args:
        st_act: Activations of stereotypical prompt(shape: (batch_size, feature_size))
        anti_st_act: Activations of anti-stereotypical prompt(shape: (batch_size, feature_size))
    Returns:
        st_act: Activations of stereotypical prompt without same pairs(shape: (batch_size, feature_size))
        anti_st_act: Activations of anti-stereotypical prompt without same pairs(shape: (batch_size, feature_size))
    """
    # remove same pairs
    new_st_act = st_act[~torch.all(st_act == anti_st_act, dim=1)]
    new_anti_st_act = anti_st_act[~torch.all(st_act == anti_st_act, dim=1)]
    return new_st_act, new_anti_st_act


# evaluate by diff
def evaluate_by_diff(
    st_act: torch.Tensor, anti_st_act: torch.Tensor
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Args:
        st_act: Activations of stereotypical prompt(shape: (batch_size, feature_size))
        anti_st_act: Activations of anti-stereotypical prompt(shape: (batch_size, feature_size))
    Returns:
        diff_sum: Difference of activations
        vals: Top 5 features
        inds: Top 5 feature indices
        avg_diff: Average of the absolute difference
    """
    diff = st_act - anti_st_act
    # get sum
    diff_sum = diff.sum(dim=0)
    # top 5 features
    vals, inds = torch.topk(diff_sum, 5)
    # average of the absolute difference
    avg_diff = diff_sum.abs().mean()
    return diff_sum, vals, inds, avg_diff


# 2-class classification by random forest
def run_random_forest(
    st_act: torch.Tensor, anti_st_act: torch.Tensor
) -> Tuple[float, torch.Tensor]:
    """
    Args:
        st_act: Activations of stereotypical prompt(shape: (batch_size, feature_size))
        anti_st_act: Activations of anti-stereotypical prompt(shape: (batch_size, feature_size))
    Returns:
        accuracy: Accuracy
        feature_importance: Feature importance
    """
    # concatenate st_act and anti_st_act
    activations = torch.cat([st_act, anti_st_act], dim=0)
    # set labels
    labels = torch.cat([torch.ones(st_act.shape[0]), torch.zeros(anti_st_act.shape[0])])
    # random forest
    clf = RandomForestClassifier()
    clf.fit(activations.cpu().numpy(), labels.cpu().numpy())
    accuracy = clf.score(activations.cpu().numpy(), labels.cpu().numpy())
    feature_importance = torch.tensor(clf.feature_importances_)

    return accuracy, feature_importance


# 2-class classification by SVM
def run_svm(
    st_act: torch.Tensor, anti_st_act: torch.Tensor
) -> Tuple[float, torch.Tensor]:
    """
    Args:
        st_act: Activations of stereotypical prompt(shape: (batch_size, feature_size))
        anti_st_act: Activations of anti-stereotypical prompt(shape: (batch_size, feature_size))
    Returns:
        accuracy: Accuracy
        normal_vector: Normal vector
    """
    # concatenate st_act and anti_st_act
    activations = torch.cat([st_act, anti_st_act], dim=0)
    # set labels
    labels = torch.cat([torch.ones(st_act.shape[0]), torch.zeros(anti_st_act.shape[0])])
    # SVM
    clf = SVC(kernel="linear")
    clf.fit(activations.cpu().numpy(), labels.cpu().numpy())
    accuracy = clf.score(activations.cpu().numpy(), labels.cpu().numpy())
    normal_vector = torch.tensor(clf.coef_)

    return accuracy, normal_vector


# sae_release = "llama_scope_lxm_8x"
# sae_id = "l{}m_8x"
sae_release = "gpt2-small-res-jb"
sae_id = "blocks.{}.hook_resid_pre"
output_dir = f"out/{sae_release}"
os.makedirs(output_dir, exist_ok=True)

# initialize
torch.set_grad_enabled(False)
device = (
    "cuda"
    if torch.cuda.is_available()
    else "mps" if torch.mps.is_available() else "cpu"
)
print(f"Device: {device}")

# load LLM
# model = HookedSAETransformer.from_pretrained("meta-llama/Llama-3.1-8B", device=device)
model = HookedSAETransformer.from_pretrained("gpt2-small", device=device)


# load dataset
df = pd.read_parquet(
    "hf://datasets/McGill-NLP/stereoset/intrasentence/validation-00000-of-00001.parquet"
)
# extract gender bias data
df = df[df["bias_type"] == "gender"]
st_prompt = []
anti_st_prompt = []
st_target_index = []
anti_st_target_index = []
for i, row in df.iterrows():
    # check labels
    # 0: anti-stereotypical, 1: stereotypical,2:unrelated
    labels = row["sentences"]["gold_label"]
    for label, sentence in zip(labels, row["sentences"]["sentence"]):
        if label == 0:
            anti_st_prompt.append(sentence)
            index = normalize_and_find_target(
                model.to_str_tokens(sentence), row["target"]
            )
            anti_st_target_index.append(index[0])
        elif label == 1:
            st_prompt.append(sentence)
            index = normalize_and_find_target(
                model.to_str_tokens(sentence), row["target"]
            )
            st_target_index.append(index[0])
        else:
            continue

# show the first 5 prompts
for i in range(5):
    print(f"Stereotypical prompt {i}: {st_prompt[i]}")
    print(f"Anti-stereotypical prompt {i}: {anti_st_prompt[i]}")
    print(f"Target index {i}: {st_target_index[i]}")
    print(f"Target index {i}: {anti_st_target_index[i]}")


layer_num = 32
bar = tqdm(range(9, 10))
with open(f"{output_dir}/output.txt", "w") as f:
    for layer in bar:
        bar.set_description(f"Loading SAE for layer {layer}")
        # load SAE
        sae, cfg_dict, sparsity = SAE.from_pretrained(
            release=sae_release,  # <- Release name
            sae_id=sae_id.format(layer),  # <- SAE id (not always a hook point!)
            device=device,
        )
        _, cache = model.run_with_cache_with_saes(
            st_prompt,
            saes=[sae],
            stop_at_layer=sae.cfg.hook_layer + 1,
            names_filter=[sae.cfg.hook_name + ".hook_sae_acts_post"],
        )
        cache = cache[sae.cfg.hook_name + ".hook_sae_acts_post"]
        # get activations
        st_act = get_activations(model, sae, st_prompt, st_target_index)
        anti_st_act = get_activations(model, sae, anti_st_prompt, anti_st_target_index)
        # remove same pairs
        st_act, anti_st_act = remove_same_pairs(st_act, anti_st_act)
        print(f"sample size: {st_act.shape[0]}")
        # evaluate by diff
        diff_sum, vals, inds, avg_diff = evaluate_by_diff(st_act, anti_st_act)
        # 2-class classification by random forest
        accuracy_rf, feature_importance_rf = run_random_forest(st_act, anti_st_act)
        # 2-class classification by SVM
        accuracy_svm, normal_vector_svm = run_svm(st_act, anti_st_act)
        # save results
        f.write(f"Layer {layer}\n")
        f.write(f"Diff sum: {diff_sum}\n")
        f.write(f"Top 5 features: {vals}\n")
        f.write(f"Top 5 feature indices: {inds}\n")
        f.write(f"Average of the absolute difference: {avg_diff}\n")
        f.write(f"Accuracy by random forest: {accuracy_rf}\n")
        vals, inds = torch.topk(feature_importance_rf, 10)
        f.write(f"Top 5 features by random forest: {vals}\n")
        f.write(f"Top 5 feature indices by random forest: {inds}\n")
        f.write(f"Accuracy by SVM: {accuracy_svm}\n")
        vals, inds = torch.topk(normal_vector_svm, 10)
        f.write(f"Top 5 features by SVM: {vals}\n")
        f.write(f"Top 5 feature indices by SVM: {inds}\n")
        # cosine similarity
        diff_svm_rf = torch.cosine_similarity(normal_vector_svm, feature_importance_rf)
        f.write(f"Cosine similarity between SVM and random forest: {diff_svm_rf}\n")
        f.write("\n")
