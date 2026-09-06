---
name: interface-forall-too-strong
description: Interface の構造体に場が足りないと、それを ∀ で量化した Skeleton の主張が原典より強くなるか偽になる。2026-09-06 に 2 件確定
metadata:
  type: project
---

**`Interface/` の構造体を `∀` で量化した `Skeleton/` の主張は、
構造体の場が足りないと「原典より強い」か「偽」になる。**

2026-09-06 に 2 件見つかった。どちらも**証明を書こうとして初めて分かった**。

## 1. [pGC] Proposition 1.2（D13、10 例目の退化）

`Skeleton/PGC/Section1.lean` の `residueCard_and_degree_recoverable (RD : ResidueCardinality p)`。
`ResidueCardinality` の場は `card` / `isPrimePow` / `card_congr` の 3 つで、
**`card K = 剰余体の実際の濃度` という場が無い**。

Lean で同値が証明された（`Check/PGC/Prop12ForallRD.lean`、sorry 0）:

```
(∀ RD, (residueCardAndDegreeObject RD).RecoverableFromAbsGal) ↔
  ∀ {K K'}, (K.absGal ≃ₜ* K'.absGal) → Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
```

★右辺は**原典 Introduction が明示的に偽だと述べている命題**
（"the Grothendieck Conjecture cannot hold in the naive sense"、
Historical Remark の "only to discover that this was, in fact, false"、引用[8] = Jarden–Ritter）。

**★2026-09-05 の修理（`card_congr` を足した）が原因**。安い反例を消した代わりに、
`card` が「同型類の任意の関数」でよいままになり、主定理より強い命題になった。

## 2. [CorrHyp] Theorem 6.1（D17。★こちらは「偽」）

`thm_6_1 (D : HyperbolicCurveData)`。`HyperbolicCurveData` は
**Prop 値の公理フィールドが `Gamma_isDiscrete` 1 本だけ**で、
`Aut` / `idAut` / `IsGenericallyScheme` は無制約のデータ。
⇒ `Aut _ := Bool`, `idAut _ := true` と取れば**偽**。

**Why:** どちらも「まだ構成できていないから `Interface` に仮説として置く」という設計から来ている。
だが**構成ができた後もそのまま**にすると、`∀`（その構造体）が別の主張に化ける。
Prop 1.2 の場合、本物の構成 `realResidueCardinality` は
2026-09-05（第 1012）に既にできていた——**仮説に取る理由がもう無かった**。

**How to apply:**

- ★`Interface` の構造体を `∀` で量化する `Skeleton` の主張を見たら、
  **その構造体の場だけで結論が出るか**を先に確かめる。出ないなら statement が強すぎる。
- ★**「本物」が既に `Found/` にできていないか**を確認する。できていれば
  `∀ RD` をやめて**その項に固定する**のが正しい（Prop 1.2 では
  `residueCard_and_degree_recoverable_real` が 2 行で通った）。
- ★**消費者がいない宣言ほど長く生き残る**。D17 は消費者 0 かつ
  **ビルドの import 閉包の外**（`lake build ABC3` が CorrHyp をコンパイルしていない可能性）
  だったので誰も文句を言わなかった。[[orchestration-graph-first]] の
  `unwired.mjs --dead` がこの型を機械で出す。
- 退化の一覧は `lean/ABC3/Check/**/*Degenerate.lean`（2026-09-06 時点で 10 例 + 候補 1）。
  ★種別が 3 つある: (a) 落とすと**自明**になる、(b) 落とすと**偽**になる、
  (c) ★**素朴な修復が偽／強すぎる主張を作る**（D10・D13・D16）。
