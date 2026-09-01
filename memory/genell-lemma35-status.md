---
name: genell-lemma35-status
description: GenEll Lemma 3.5 は証明済みだが条つき。条なし化の障害と正しい道が確定している。
metadata:
  type: project
---

**2026-09-01（第 944-1091、148 ブロック）で `Lemma 3.5` の高さ不等式が通った。**
`Found/GaloisRep/Lemma35Unconditional.lean` の `lemma_3_5_height_ineq`、
外部入力なし、`#print axioms` は 3 公理。

★ただし**原典に無い仮説 2 つ**を置いている:

- `l ≠ 2` —— Vélu の `±` の対（`veluWFull` の `/2`）が不動点を持たないため
- `d + 1 < l` —— 第 1044 が悪い素点で `p ∤ l` を出すのに使う

☆`Lemma 3.7` の**条件 (a) の枝では吸収される**（`l ≥ 100·d ≥ 100`）。
第 1088 `htFalt_le_of_condA_lcyclic` で実際に吸収した。
★しかし**条件 (b) の枝は `2 ≤ l` だけ**で `(†)` を消費するので吸収されない。
このため `.src` は条つきのまま（指標 18/24）。第 1087 で条なしにしたのは誤りで第 1089 で撤回。

## 条なし化の障害（実測して確定、第 1090-1091）

☆逃げ道は 2 つとも潰した:

1. **「商体 `K` に移す」は通らない** —— Tate 級数は `q`-進収束を要するが、
   `K` では `q` が可逆なので `IsHausdorff (q) K` が破れる。`R[1/l]` も `p ∣ l` では `K`。
2. **`l` の可逆性は本質的** —— μ-等級付きの核 `tateXterm_zeta_eq_poly`
   （`MuGraded.lean:362`）の係数 `tateXtermC` が `Ring.inverse ((l : R))^2` を含む。

★**正しい道は分母払い**。第 277（`CollDenomFree.lean`）が共線性で使った技法
（`u ≡ 1` で `IsUnit` が破れるので行ごとに分母を払う）を μ-等級付きの計算に当て、
`l^k × (恒等式)` の形で `R` の中で証明し、最後に `K` で `l^k` を割る。
在庫の `prod_one_sub_pow_erase`（`∏(1−ζ^i) = l`）が `1−ζ^i` 側を `l` の冪に集約する。
☆規模 40-80 ブロック。

★**群の水準は既に無条件**——`lemma_3_2` (ii)（`Found/GenEll/Lemma32.lean`）が
`(Lˣ/q^ℤ)/μ_l ≃ Lˣ/(q^l)^ℤ` と `vAdd v (q^l) = l·vAdd v q` を `IsUnit (l)` なしで与える。
止まっているのは「群の商」と「Vélu の曲線」を繋ぐ段だけ。

詳細は `ResearchPaper/mathlib-gap.json` の `lemma35BareCost20260901`。
関連: [[genell-track-b]]、[[frdi-split-nonisotropic-not-derivable]]（在庫を先に引く教訓）。
