---
name: frdi-prop21-quantifier-false
description: [FrdI] Prop 2.1 (iii) は原文どおり「d 固定」で読むと偽。意図された「すべての d」版は証明済み。
metadata:
  type: project
---

[FrdI] `Proposition 2.1, (iii)`「`𝒞` が perfect 型 ⟺ `Ψ` が圏同値」は、
原文が冒頭で `Let d ∈N≥1.` と `d` を固定しているとおりに読むと **偽**である。
機械検証済み(`lean/ABC3/Check/FrdI/Prop21QuantifierGap.lean`、2026-08-17):

* `d = 1` では次数 1 の Frobenius 型射がつねに同型なので `Ψ₁ ≅ 𝟭`。
  よって **perfect 型でない Frobenioid すべてが反例**(`prop_2_1_iii_fixed_degree_false_general`)。
* `d ≥ 2` でも同じ。`Φ = ℤ/2`(定数)では 3 倍が全単射・2 倍が非全射なので、
  `Ψ₃` は圏同値だが `n = 2` で perfect が破れる(`prop_2_1_iii_fixed_degree_false_deg_three`)。

**Why:** `perfect object`(Definition 1.2, (iv))は**すべての `n ∈ ℕ≥1`** を量化するのに、
`Ψ_d` が与えるのは `n = d` の場合だけだから。次数は乗法的なので 1 つの `d` から他の `n` は作れない。

**How to apply:** 意図された主張「すべての `d` で `Ψ_d` が圏同値 ⟺ perfect 型」は
`Found/FrdI/Prop21.lean` の `prop_2_1_iii` で**証明済み**なので、後続で使うときはそちらを引く。
`Proposition 2.1` には `.src` を付けていない(`Gap/FrdI/Section2.lean` の `Gap_2_1_iii`、③ sourceGap)。
★これは数学の誤りではなく**書き方の不備**と見るのが妥当。[[frdi-s1-two-blockers]] の
Prop 1.14(仮定の落ち)とは種類が違う。
