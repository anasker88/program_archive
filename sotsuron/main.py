import gc
import os
from functools import partial

import pandas as pd
import plotly.express as px
import torch
from datasets import load_dataset
from sae_lens import SAE, ActivationsStore, HookedSAETransformer
from sklearn.svm import SVC
from tqdm import tqdm
from transformer_lens.utils import tokenize_and_concatenate

# open output file(make if not exists)
f = open("out/result.txt", "w")

# initialize
torch.set_grad_enabled(False)
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Device: {device}")

# load LLM
model = HookedSAETransformer.from_pretrained("meta-llama/Llama-3.1-8B", device=device)

# load dataset
df = pd.read_parquet(
    "hf://datasets/McGill-NLP/stereoset/intrasentence/validation-00000-of-00001.parquet"
)
# extract gender bias data
df = df[df["bias_type"] == "gender"]
st_prompt = []
anti_st_prompt = []
target_index = []
for i, row in df.iterrows():
    print(row["target"])
    # get the last token of the target
    i = model.get_token_position(row["target"], row["context"])
    # check labels
    # 0: anti-stereotypical, 1: stereotypical,2:unrelated
    labels = row["sentences"]["gold_label"]
    for label, sentence in zip(labels, row["sentences"]["sentence"]):
        if label == 0:
            anti_st_prompt.append(sentence)
        elif label == 1:
            st_prompt.append(sentence)
        else:
            continue

layer_num = 32
bar = tqdm(range(layer_num))
for layer in bar:
    bar.set_description(f"Loading SAE for layer {layer}")
    # load SAE
    sae, cfg_dict, sparsity = SAE.from_pretrained(
        release="llama_scope_lxr_8x",  # <- Release name
        sae_id=f"l{layer}r_8x",  # <- SAE id (not always a hook point!)
        device=device,
    )
    # load activation store
    activation_store = ActivationsStore.from_sae(
        model=model,
        sae=sae,
        streaming=True,
        # fairly conservative parameters here so can use same for larger
        # models without running out of memory.
        store_batch_size_prompts=8,
        train_batch_size_tokens=4096,
        n_batches_in_buffer=32,
        device=device,
    )
    # get max activations
    sae.cfg.normalize_activations
    diff = torch.zeros(cfg_dict["d_sae"], device=device)
    width = 5
    bar.set_description("Running SAE with prompts")
    for i in range(len(st_prompt) // width):
        start_index = width * i
        end_index = min(width * (i + 1), len(st_prompt))
        _, st_cache = model.run_with_cache_with_saes(
            st_prompt[start_index:end_index],
            saes=[sae],
        )
        _, anti_st_cache = model.run_with_cache_with_saes(
            anti_st_prompt[start_index:end_index], saes=[sae]
        )
        batch_target_index = target_index[start_index:end_index]
        for i, target_pos in enumerate(batch_target_index):
            diff += (
                st_cache[cfg_dict["hook_name"] + ".hook_sae_acts_post"][
                    i, target_pos + 1, :
                ]
                - anti_st_cache[cfg_dict["hook_name"] + ".hook_sae_acts_post"][
                    i, target_pos + 1, :
                ]
            )
        # キャッシュをクリア
        del st_cache
        del anti_st_cache
        torch.cuda.empty_cache()
    # let's print the top 5 features and how much they fired
    vals, inds = torch.topk(
        diff,
        5,
    )
    f.write(f"Layer {layer}:\n")
    for val, ind in zip(vals, inds):
        f.write(f"  Feature {ind}: {val:.2f}\n")
    # visualize the activations
    px.line(
        diff.cpu().numpy(),
        title=f"Difference in activations for layer {layer}",
        labels={"index": "Feature", "value": "Activation"},
    ).write_html(f"out/activation_diff_{layer}.html")
    # average of the absolute difference
    diff = torch.abs(diff)
    avg_diff = diff.mean()
    f.write(f"Average difference: {avg_diff:.2f}\n")

f.close()
