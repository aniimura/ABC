---
name: frdi-1uniqueness-rigidity-patterns
description: [FrdI] の「1-一意性」と「rigidity」は毎回同じ 2 つの型で落ちる。手で共役を組まないこと。
metadata:
  type: project
---

**[FrdI] の §3 以降に頻出する「1-一意な関手が図式に入る」「合成関手は rigid」は、
毎回同じ 2 つの型で落ちる。**

## ★型 1 —— 1-一意性は「**恒等になる合成を 1 つ見つける**」だけで決まる

| 定理 | 恒等になる合成 |
|---|---|
| `Theorem 3.4, (i)` | `ι ⋙ isotropification ≅ 𝟭`(`isotropificationRestrictIso`、`Proposition 1.9, (v)`) |
| `Proposition 3.11, (iii)` | `toProdN ⋙ fst = 𝟭`(★**`rfl`**) |

★どちらも **`Equivalence.invFunIdAssoc`**(mathlib)で **4 行**。
★随伴の一意性(`Adjunction.leftAdjointUniq`)を経由すると、
**右随伴どうしの同型を作る所で whisker 計算が要って詰まる**(2026-08-17 に実測、差し戻し)。

## ★★型 2 —— rigidity は「**その操作が関手(圏同値)であることを言う**」

`Ψ` が圏同値で `F` が rigid なら `Ψ ⋙ F` も rigid
(`isRigidFunctor_comp_of_isEquivalence`、`Thm34.lean`)。

★★**手で共役を組んで自然性を証明しようとすると詰まる**(2026-08-17 に実測、差し戻し)。
`Functor.fun_inv_map` で余単位に開いたあと `F.map` を跨いだ結合の書き換えが
`rw`/`simp` のどちらでも収束しない。

★★★**`Equivalence.congrLeft` を使う** —— `F ↦ Ψ ⋙ F` **自体が関手圏の圏同値**なので、
充満忠実性(`preimageIso` ＋ `map_preimage`)から `Aut(Ψ ⋙ F)` の元は `Aut(F)` の元の像。
**自然性は `congrLeft` が既に持っている。実測 8 行。**

★併せて `isRigidFunctor_of_iso`(rigidity は関手の同型で移る)も要る ——
**1-可換性が「2 つの合成関手は同型」を与えるので、片方で示せば両方に届く**
(原文の「each of the composite functors」の中身)。

## ★どちらにも共通する教訓

**「手で書く」より「既にある構造(関手・圏同値・随伴)に載せる」方が短い。**
★詰まったら「この操作は何かの関手か? 圏同値か?」を先に問う。

関連: [[frdi-split-nonisotropic-not-derivable]](在庫検索の規則)。
