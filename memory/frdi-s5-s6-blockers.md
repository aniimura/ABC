---
name: frdi-s5-s6-blockers
description: FrdI §5 と §6 の残りは mathlib の在庫不足 2 群に帰着する（Nikolov–Segal / 因子論・six exponentials）
metadata:
  type: project
---

FrdI の §5（2/7）と §6（0/5）の残りを 2026-08-20 に項目ごと測った。**残りは 2 群の在庫不足に帰着する。**

- **`Proposition 5.6` → Nikolov–Segal**。`M_p ⊆ 𝒪^×(A)` が `Aut_𝒞(A)` の共役で閉じることが要るが、`M_p` は位相（pro-`l` 部分）で定義されるので共役が**連続**でなければ保たれない。位相的有限生成な副有限群で抽象自己同型がすべて連続なのが Nikolov–Segal で、mathlib に `Nikolov` は **0 件**。★これは `Definition 2.8, (i)` の角括弧「[uniquely determined]」（`Def28.lean` が主張しないと明記したもの）を**そのまま消費している**。`Corollary 5.7` はその下流。記録は `Gap/FrdI/Section5.lean` の `Gap_5_6_Mp`。
- **§6 全 5 件 → 因子論と six exponentials**。`WeilDivisor` / `CartierDivisor` / `PrimeDivisor` は mathlib に **0 件**（`Scheme.functionField`・`Morphisms/Proper`・`Geometrically/Integral`・相対正規化 `f.normalization` は**在る**）。`Example 6.1` は `Φ`・`B`・`Div_B` さえ作れば `Theorem 5.2, (ii)`（実装済み）を当てるだけになる。

**それ以外（Prop 5.3 / Cor 5.4 / Prop 5.5）は在庫不足ではなく、残りが 1〜2 の技術的段まで縮んでいる。**

- `Proposition 5.5, (i)` は**数学的内容がすべて済み**（`Prop55Pf.lean`：遷移写像 `α ↦ α^m`・cofinality・well-defined 性・全射性・単射性）。残りは `Pf (Additive (𝒪^▷(A))) ≃* 𝒪^▷(A^pf)` への束ね上げだけ。
- `Corollary 5.4` は `η ∘ Div^gp = Div^gp ∘ Ψ^birat` の 1 段（在庫 `biratPsi_otimes_map` は在る）と 1-一意性・rigidity。
- `Proposition 5.3` は `p55i` 経由の perfection の同値待ち。**行き先**（`Crlf` / `CuntrPf`）と縦の矢印（`untrToSc` / `cToSc`）は作ってある。

★足した仮定は 2 つで、どちらも `index.html` に開示済み：`𝒟` が of FSM-type（`Φ^birat` を monoid on 𝒟 にするため）、`Φ^rlf` の条件 (a)（テンソルの平坦性。perf-factorial なら真と測定済み）。

関連: [[no-wall-decompose-instead]] [[measure-mathlib-before-skeleton]] [[frdi-s1-two-blockers]] [[frdi-three-chains]]
