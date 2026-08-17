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

**★2 度目(2026-08-17)**: `Theorem 3.4, (i)` の 1-可換性で
`isotropification` の四角形が要ったので `isotropification_map_comp` を**作った**。
★★**`Prop19.lean` に `isotropification_square` として既にあった**——同じ主張である。
撤去して既存を使った。★1 度目は「導けない」と誤り、2 度目は「作り直した」。
**同じ根**(在庫を検索しない)から出た別の症状である。

**★2026-08-17、同じ日に 4 回起きた。**

| # | 何をした | 在庫にあったもの |
|---|---|---|
| 1 | 「`Definition 2.3` から導けない」と誤測定 | `isCoAngular_of_endo`(`Prop18.lean`) |
| 2 | `isotropification_map_comp` を書いた | `isotropification_square`(`Prop19.lean`) |
| 3 | 「`Φ.map` の単射性は無いので原文の隠れ段だ」と測定しかけた | `MonoidOn.charInj` / `map_injective`(構造フィールド) |
| 4 | `Fsmff.lean` を新規作成した | `not_isFSMI_endo_of_isOfFSMFFType`(`CategoryVocabulary.lean:681`) |

★★**4 では `grep` を実行していた**。にもかかわらず外した。
検索語が `IsFSMIChain` / `IsOfFSMFFType`(**Lean の識別子**)だったからである。
★★★**このリポジトリの在庫は日本語の docstring で索引されている** ——
`FSMI 自己射` で引けば 1 発で当たった。

**How to apply:**
- ★**補題を書き始める前に `grep` する**。「導けないと書く前」だけでは足りない。
- ★★★**検索語は 3 種類を必ず併用する**:
  1. **Lean の識別子**(`degFr`, `≫ hullMap`)
  2. **主張の形**(`IsIso (P.Base`, `≫ ?f = ?g ≫`)
  3. ★★**日本語の言い回し**(`自己射`, `単射`, `一意`, `交換`, `伝わる`)
     —— **3 が一番よく当たる**。docstring が日本語だからである。
- 「原文のこの段が導けない」と書きたくなったら、★**まず `grep` で自分の在庫を探す**。
  次に**原文の同じ節の他の条**(ここでは `Definition 1.3` の (iii)(b))を読み直す。
- 測定の記録は安いが、**誤った測定は後の判断を歪める**。撤回も同じ場所に書くこと。
- 残作業(2026-08-17 時点): 非 isotropic な `A` での `τ(A)` を `τ(A^istr)` の
  引き戻しとして確定させること。原文の `τ` は `(𝒞^istr)^lin` 上の部分関手なので、
  非 isotropic での `τ(A)` は**引き戻しが定義**である。我々は `τ` を全対象で
  与えているため、この対応を `IsCharacteristicSplitting` の条件として書く必要がある。
  これが済めば `Prop25iii.lean` の `IsOfIsotropicType` を外す道が開く。
- 関連: [[frdi-prop21-quantifier-false]] / [[lean-build-check-discipline]]。
