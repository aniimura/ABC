---
name: lean-fullsubcat-procedures
description: 充満部分圏(Istr P など)を Lean で扱うときの3手順と、その例外条件
metadata:
  type: feedback
---

`Istr P`(= `ObjectProperty.FullSubcategory`)の射は `InducedCategory.Hom` に包まれ、
型が簡約されないため `rw` が motive を組めない。移送補題を書くときは**最初に**次を課す:

1. **`Istr P` を含む等式に `rw` を使わない** —— `calc` + `congrArg` + `InducedCategory.hom_ext` で書く。
2. **`inv` を書かない** —— 同型は `InducedCategory.isoMk` で持ち上げ、逆射は `IsIso.out` から明示的に取る。
3. **`include F in` を最初に書く**(`F` が主張に現れない補題では自動包含されない)。

**Why:** 原因を知ることと再発を止めることは別。手順に落とした後、8条連続で構造的修正 0 になった。

**How to apply:** 書き始める前に自分に課す。書いた後に直すのでは遅い。

## ★手順2の例外条件(2026-08-15、`istr_preStepSpan` で確定)

**主張そのものが `@inv _ _ _ _ (Base φ) hφ.2` を含む場合は避けられない。**
そのときは `IsIso.eq_inv_comp` を**インスタンスを明示して**当てる
(`(@IsIso.eq_inv_comp _ _ _ _ _ f hf _ _).mp/.mpr`)。
証明の**中で**使う逆射は依然として `IsIso.out` から取り、`inv` は書かない。

★親の予想「`𝒟` の射については `inv` を使ってよい」は**外れ**だった ——
`𝒟` 側の逆射も `.out` から取れて `inv` は不要。例外は「主張に現れるとき」だけである。

手順は固定するものではなく、反例が出たら条件を精密化するもの。関連: [[abc3-plan-two-track]]
