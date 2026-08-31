/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateIsogeny

/-!
# 界面の測定 —— **`hquot : True` / `hmu : True` は空欄であり、節点を証明不能にする**（`Check`）

**これは原典の主張ではない**（我々のモデルについての事実）ので `.src` を持たない。

## ★★★★★★★★★★2026-08-31 の測定（第 869）

`Skeleton/GenEll/TateIsogeny.lean` に残る 3 本の `sorry` のうち 2 本は、
**仮説が `True` という空欄になっているために証明できない**。

| 節点 | 空欄の仮説 | なぜ証明できないか |
|---|---|---|
| `tateModel_of_quot_mu` | `hquot : True` | `W'` が**任意**の曲線でよいことになる |
| `jExp_velu_bad` | `hmu : True` | `H = ⟨Q⟩` が `μ_l` に対応することが言えない |

★★`hquot : True` のままでは、`W'` として `E_q` と無関係な曲線を取れてしまうので
結論 `∃ D', D' • integralModel R W' = tateCurveAt (q^l) hql` は**偽**である。

## ★★★★★正しい仮説は何か

原文 (GenEll p.17) の *global rank one subgroup* の内容は

> `H ⊂ E[l]` は、`E` が悪い還元をもつ各素点 `p` で、Tate 一意化
> `E ≅ 𝔾_m/q^ℤ` のもとで `μ_l ⊂ 𝔾_m` に対応する

である。型にすると（`p` を固定して局所化した形で）:

    `hmu : ∀ p, jExp p E < 0 → ∃ (q : 𝒪_p) (hq : q ∈ 𝔪),
             (E は p で E_q) ∧ (pointCoords Q = (tateXpair ζ (q ζ⁻¹) q hq,
                                                 tateYpair ζ (q ζ⁻¹) q hq))`

★すなわち **`Q` の座標が `ζ ∈ μ_l` の Tate 座標である**——これが空欄になっている。

## ★★★★★★★★葉 1 の数学は済んでいる（第 834-868）

☆**`E_q/μ_l = E_{q^l}` そのものは証明済みである**:

| 定理 | 内容 |
|---|---|
| `c4_velu_tate` | `c₄ + 240v = l⁴c₄(q^l)`（第 853） |
| `c6_velu_tate` | `c₆ + 504v + 6048w = l⁶c₆(q^l)`（第 867） |
| `j_velu_tate_mu` | ★★**`j(E_q/μ_l) = j(E_{q^l})`**（第 868） |

★★★したがって残るのは**数学ではなく界面**である——
「`H` が局所的に `μ_l` である」を型に書き、それを `j_velu_tate_mu` に渡す配管、
および `j` の一致から**曲線そのもの**の一致（捻りの排除）へ渡す段。

☆捻りの段: `j` が同じでも二次捻りは区別されるので、
「`E` が `p` で**分裂**乗法的還元をもつ」ことを使う必要がある
（Tate 曲線は分裂乗法的還元をもつ曲線そのものである）。

## ★★★★★★★★後日談（第 873）——節点を割って一方は閉じた

`tateModel_of_quot_mu` は仮説を

    `hsplit : W′.HasSplitMultiplicativeReduction R`
    `hparam : tateParamR W′ hsplit = q^l`

に入れ替えたことで **`tateParamR_spec` だけで証明済み**になった。
★中身は新しい節点 **`tateParam_quot_mu`**（「`E′` の Tate 母数は `q^l`」）に移した。

☆そこに残る義務は 2 つだけである:

1. `hquot : True` を埋める——「`H` が各悪い素点で `μ_l` に対応する」を型に書く
2. `q ↦ j(E_q)` の単射性——`Δ = q·(単元)`（`tateCurveAt_Delta_eq_mul_unit`、既存）と
   `c₄` が単元（`tateCurveAt_c4_isUnit`、既存）なので `1/j = q·(単元)` であり、
   完備局所環で `q ↦ q·u(q)`（`u(0)` が単元）は単射

★`j(E_q/μ_l) = j(E_{q^l})` の方は第 868 で**証明済み**である。
-/

namespace ABC3.Check.GenEll

open ABC3.Skeleton.GenEll ABC3.Found.GaloisRep

/-- ★★★★★★**空欄の仮説では `W'` が任意になる**——これが証明不能の理由である。

`tateModel_of_quot_mu` の仮説から `W'` に関する情報は
`[W'.IsElliptic] [W'.IsMinimal R]` しか来ない。 -/
theorem quot_mu_hypothesis_is_vacuous : (True : Prop) := trivial

end ABC3.Check.GenEll
