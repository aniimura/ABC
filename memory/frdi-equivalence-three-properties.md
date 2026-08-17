---
name: frdi-equivalence-three-properties
description: [FrdI] の「商圏/部分圏が Frobenioid」を示すとき、圏同値の 3 性質は毎回別々の道具で落ちる。行き先が前順序圏なら充満性は存在問題に潰れる。
metadata:
  type: project
---

**[FrdI] `Definition 1.3, (iii), (d)` の 2 本の圏同値
(`coaPreUnderEquiv` / `coaPreOverEquiv`)は、`𝒞^istr`・`𝒞^un-tr`・`𝒞^pf` の
どれでも同じ形で落ちる。**

## ★★行き先 `Order(Φ(A))` は**前順序圏**である

★hom は `MLe` という **Prop** なので `Subsingleton`。したがって:

| 性質 | 何に帰着するか |
|---|---|
| **忠実性** | 添字圏の hom が **subsingleton** であること |
| **充満性** | ★**「射が 1 本あればよい」**——`F.map f = g` は `Subsingleton.elim` で自動 |
| **本質的全射性** | 既にある側(`𝒞^istr` など)の証人を押し出すだけ |

★★**充満性を「`map` の等式まで作る」問題だと読むと無駄に重くなる。**

## ★添字圏の hom が subsingleton になる理由は 2 つ(別々)

- **コスライス `_A(𝒞^coa-pre)`** —— `Z.hom` が **epi**(圏が totally epimorphic)
- **スライス `(𝒞^coa-pre)_A`** —— `W.hom` が **mono**(pre-step は mono、`Definition 1.3, (v), (a)`)

★**同じ「2 本の脚を別々の理由で消す」形が `Theorem 3.4, (i)` の rigidity
段 C にも出る** —— pre-step の span `B₀ ←φ— X —ψ→ A` で、
**左脚 `φ` は mono、右脚 `ψ` は epi**。

## ★★商圏へ持ち上げるときの Lean の手

`𝒞^istr → 𝒞^un-tr` のような**全射的な関手**で持ち上げてから性質を移すとき、
★**対象を先に `obtain ⟨…⟩ := Z` で分解して `hom` を局所変数にする**。
そうすれば `hζ : map ζ = zh` を **`subst` できる**ので、
`inv (Base …)` のような**インスタンス引数を含む項**の書き換えで詰まらない。
★これをやらないと `rw` の motive が壊れる(2026-08-17 に実測)。

関連: [[frdi-1uniqueness-rigidity-patterns]]、[[lean-fullsubcat-procedures]]。
