---
name: genell-bd-class-direction
description: [GenEll] Def 1.2 (ii) の記号 ≲ の向きが、同じ論文の Thm 2.1 と [IUTchIV] Cor 2.3 と食い違う(3 箇所を PDF 目視、子が独立確認)。
metadata:
  type: reference
---

**[GenEll] Definition 1.2, (ii)(物理 p.5)に印字された `≲` の向きは、
同じ論文の `Theorem 2.1`(p.11)および [IUTchIV] `Corollary 2.3`(p.54)での用法と逆である。**

| 出典 | 目視した内容 |
|---|---|
| [GenEll] p.5, Def 1.2, (ii) | `α ≲ β` ⟺ ∃C, **`β(x) − α(x) ≤ C`**(= `β − α` が**上**に有界) |
| [GenEll] p.11, Thm 2.1, (i) | 表題が **“Effective Mordell/ABC/Vojta Conjecture”**、式は `ht ≲ (1+ε)(log-diff + log-cond)` |
| [IUTchIV] p.54, Cor 2.3 | 同じ式に「i.e., … `(1+ε)(…) − ht` **is bounded below** by a constant」。`[cf. [GenEll], Definition 1.2, (ii)]` と明示的に引く |

**Why:** abc/Vojta の内容は `ht ≤ (1+ε)(…) + C`(= 差が**下**に有界)。
Thm 2.1 の表題と IUTchIV の言い換えは互いに整合し、abc とも整合する。
Def 1.2 (ii) の印字だけが逆になる。**最も無理のない読みは「α と β が入れ替わった誤植」だが、
これは推論であって観測ではない。**

★**読字は独立に裏が取れている**——文脈を持たない子に**結論を伏せて**同じ 3 箇所を読ませ
(260 dpi 全頁 ＋ 該当箇所 500 dpi 拡大)、判定は「**逆**」。
`GAP` 型(否定的)の判定なので `OK` より強い証拠である([[challenger-audit-without-context]])。

**How to apply:** `lean/ABC3/Found/GenEll/BDClass.lean` は**どちらも型に出している**——
`BDle`(印字どおり)/ `BDge`(abc の向き)/ `bdle_ne_bdge`(反例で両者が別物と示す)。
★**abc の主張を書くときは `BDge` を使う**。この選択は機械可読になっている。
`PLAN.md` §5 の類型(①モデル化 ②未構築 ③飛躍)に**第 4 の類型「原典の誤植」が要る**——
飛躍と違い falsifier ではなく、**別の箇所での同じ記号の用法**が証拠になる。

関連: [[genell-track-b]] / [[abc3-plan-two-track]]
