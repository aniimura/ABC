---
name: frdi-split-nonisotropic-not-derivable
description: 「導けない」と書く前に自分の在庫を検索する。[FrdI] の分裂の非 isotropic 拡張で一度これを誤った。
metadata:
  type: feedback
---

**「導けない」と記録する前に、まず自分たちが既に作った補題を検索する。**

**Why:** 2026-08-17、[FrdI] `Proposition 2.5, (iii)` の「分裂は `A` が isotropic でなくても
使える」(物理 p.49)について、★**「`Definition 2.3, (a), (b)` からは導けない」と測定して
コミットした**。★★**それは誤りだった。**

見落としていたのは `Definition 1.3, (iii), (b)`:
> If A′ → A is a co-angular pre-step of C, then any morphism A′ → A is co-angular.

★**恒等射は co-angular な pre-step** なので `A′ = A` に当てると
**自己射はすべて co-angular**——isotropy は要らない。これで
`Proposition 2.2, (iii)`(`Div` が等しい ⟺ `𝒪^×` 倍だけ違う)が任意の `A` で成り立ち
(`otri_div_eq_iff'`、`lean/ABC3/Found/FrdI/Prop25.lean`)、詰まっていた
「`u'` が像に入るか」は**迂回できる**——`u'` を持ち上げず、`Div x = Div t` から
直接 `x = t ≫ u` を得ればよい。

★★★**しかもその補題は既に `Prop18.lean` に `isCoAngular_of_endo` として存在した。**
自分たちの在庫を使い損ねたのである。

**How to apply:**
- 「原文のこの段が導けない」と書きたくなったら、★**まず `grep` で自分の在庫を探す**。
  次に**原文の同じ節の他の条**(ここでは `Definition 1.3` の (iii)(b))を読み直す。
- 測定の記録は安いが、**誤った測定は後の判断を歪める**。撤回も同じ場所に書くこと。
- 残作業(2026-08-17 時点): 非 isotropic な `A` での `τ(A)` を `τ(A^istr)` の
  引き戻しとして確定させること。原文の `τ` は `(𝒞^istr)^lin` 上の部分関手なので、
  非 isotropic での `τ(A)` は**引き戻しが定義**である。我々は `τ` を全対象で
  与えているため、この対応を `IsCharacteristicSplitting` の条件として書く必要がある。
  これが済めば `Prop25iii.lean` の `IsOfIsotropicType` を外す道が開く。
- 関連: [[frdi-prop21-quantifier-false]] / [[lean-build-check-discipline]]。
