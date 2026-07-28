"""ESM masked-marginal scoring for MAPT missense variants."""

from __future__ import annotations

from .constants import AMINO_ACIDS, MAPT_441_SEQUENCE
from .variants import MissenseVariant


def load_esm_model(model_name: str):
    try:
        import esm  # type: ignore
        import torch  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "ESM scoring requires optional dependencies. Install with "
            'pip install -e ".[esm]" or use environment.yml.'
        ) from exc

    loader = getattr(esm.pretrained, model_name, None)
    if loader is None:
        raise ValueError(f"Unknown fair-esm model loader: {model_name}")
    model, alphabet = loader()
    model.eval()
    return model, alphabet, torch


def masked_marginal_scores(
    model_name: str,
    sequence: str = MAPT_441_SEQUENCE,
    device: str | None = None,
) -> list[dict[str, object]]:
    """Compute alt minus wild-type log-probability at each MAPT position."""
    model, alphabet, torch = load_esm_model(model_name)
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)

    batch_converter = alphabet.get_batch_converter()
    _, _, tokens = batch_converter([("MAPT_441", sequence)])
    tokens = tokens.to(device)

    aa_to_token = {aa: alphabet.get_idx(aa) for aa in AMINO_ACIDS}
    mask_idx = alphabet.mask_idx
    rows: list[dict[str, object]] = []

    with torch.no_grad():
        for zero_based, wt_aa in enumerate(sequence):
            token_position = zero_based + 1
            masked = tokens.clone()
            masked[0, token_position] = mask_idx
            out = model(masked, repr_layers=[], return_contacts=False)
            logits = out["logits"][0, token_position]
            log_probs = torch.log_softmax(logits, dim=-1)
            wt_logp = float(log_probs[aa_to_token[wt_aa]].item())

            for mut_aa in AMINO_ACIDS:
                if mut_aa == wt_aa:
                    continue
                mut_logp = float(log_probs[aa_to_token[mut_aa]].item())
                variant = MissenseVariant(position=zero_based + 1, wt_aa=wt_aa, mut_aa=mut_aa)
                rows.append(
                    {
                        "variant_id": variant.variant_id,
                        "protein_change": variant.protein_change,
                        "position": variant.position,
                        "wt_aa": variant.wt_aa,
                        "mut_aa": variant.mut_aa,
                        "model": model_name,
                        "esm_wt_logp": wt_logp,
                        "esm_mut_logp": mut_logp,
                        "esm_llr": mut_logp - wt_logp,
                    }
                )

    return rows

