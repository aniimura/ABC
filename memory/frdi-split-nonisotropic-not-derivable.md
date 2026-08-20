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

**★★★2026-08-17、同じ誤りが mathlib に対しても 2 回起きた。**

これまでは「**自分たちの**在庫を検索しない」という症状だった。
★★同じ日に、**mathlib の在庫**に対しても 2 回外した:

| # | 何を書いたか | mathlib に実際にあったもの | なぜ外したか |
|---|---|---|---|
| 5 | 「`Cartier` は mathlib に 0 件だから底が抜けている」 | `IdealSheafData.comap` で足りた | ★**名前で測った**(`Cartier` の語で検索) |
| 6 | 「一般のイデアルの拡大の重複度公式は無い」 | `emultiplicity_map_eq_ramificationIdx_mul` | ★★`ramificationIdx` で検索し、**`emultiplicity` で検索しなかった** |

★★★**mathlib は「数学の概念名」ではなく「形式化された補題の名前」で索引されている。**
`Cartier` という語が 0 件でも、**その中身**(イデアル層の引き戻し)は在る。
`ramificationIdx` を含む補題を全部見ても、**同じ内容が `emultiplicity` の語で**述べられている。

**How to apply:**
- ★**補題を書き始める前に `grep` する**。「導けないと書く前」だけでは足りない。
- ★★★**mathlib を検索するときは「概念名」と「量の名前」の両方で引く**:
  - 概念名(`Cartier`, `Arakelov`, `moving lemma`)—— ★**当たらないことが多い**
  - 量の名前(`emultiplicity`, `ramificationIdx`, `comap`, `map`)—— ★★こちらが本体
  - ★**ディレクトリを直接見る**(`ls Mathlib/AlgebraicGeometry/IdealSheaf/`)——
    ファイル名は概念で切られているので、これが一番速い
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

## ★★★2026-08-18: 同じセッションで **7 回**外した

FrdI §3/§4 を進める 1 セッションで、「無い」と見積もったものが
在庫にあった例が 7 件出た:

1. `endo_isCoAngular`(`Prop110.lean`)—— 自己射の co-angular 性
2. `IsOfStandardType` / `IsOfRationallyStandardType`(`Def31`/`Def45`)
3. `Remark 3.1.1`(`Remark311.lean` に丸ごとあった)
4. `isDivIdentity_of_isBaseIdentity`(`MorphismTypes.lean:530`)
5. `isGroupLikeObj_of_baseIso`(`Prop110.lean:2733`)
6. `prop_1_10_i_exists`(任意の射版。3 成分に割ってから気づいた)
7. `Proposition 1.4, (v)` —— ★**`FrobenioidCore` のフィールド**
   (`preStepMono` / `preStepFactor` / `preStepFactor'`)。
   名前つき定理を grep していて見落とした

★★**共通する失敗の形**: grep のパターンを
「自分が付けそうな名前」で探している。
★★★**対策**: 実装する前に、目標の命題を
 の**複数の語彙**で grep する。
★★★**成果として定理名だけでなく、**
- **構造のフィールド**(`FrobenioidCore` の 21 条など)
- **任意の射版**(場合分けした後で探すと遅い)
も探すこと。
それでも見つからなければ**書いてみて Lean の重複エラーに聞く**のが確実である
(5 件中 2 件は重複エラーで気づいた)。

★2026-08-19、**8 件目と 9 件目**:
8. `Definition A.1` / `Proposition A.2`(FrdI 付録)—— 逆に**在庫に無かった**ことを
   grep で確認してから書いた。これは正しい手順。
9. `psiBase` / `psiBaseUniq` —— `Prop311.lean` に既にあった
   (`Proposition 3.11, (iii)` の `Φ = 0` の特別な場合)。
   ★★**単一ファイルの `lake build` では衝突が出ない**。
   `lake build`(全体)まで走らせて初めて
   「environment already contains」で落ちた。
   ★対策に追記: **新しい名前を作る前に `Grep(pattern: "def <name>\b", path: "lean/ABC3")`**、
   かつ**節点を閉じる前に全体ビルドを 1 回通す**こと。

★2026-08-19、**10 件目**: `nfMap` が linear / base-iso / pre-step を保つことを
手で書き下したが、`Prop21.lean` に `nfMap_preStep` / `nfMap_frobType` /
`nfMap_pullBack` / `prop_2_1_ii_degFr` があり、さらに `Prop110.lean` に
**`prop_1_10_i_coAngular_of`** まであった(全部で 5 本の重複)。
★★**検索語を `naiveFrob` / `prop_2_1` に限ったのが原因**。
`nfMap_*` という実際の名前で引いていれば 1 回で当たった。

**対策の追記**: 「性質 P が構成 X で保たれる」を書く前に、**2 通りで引く**:
1. `prop_<原典番号>_.*_of`(原典由来の移送補題の命名規則)
2. `<構成の接頭辞>_<性質>`(例: `nfMap_preStep`、`biratMap_comp`)
★構成の接頭辞は定義の名前(`nfMap`・`biratMap`・`istrMap` など)であって、
論文の名前(`naiveFrob`)ではない。

★2026-08-19、**11 件目**: `Remark 4.5.1`(standard 型は `𝒞^istr` で保たれる)を
手で書き下したが、**`Rmk451.lean:227` に `istr_standardType` として完成していた**
(6 フィールドすべて)。

**対策の追記(検索の第 3 の手)**: **原典の番号でファイルを引く**。
`Remark 4.5.1` なら `Rmk451.lean`、`Proposition 1.10` なら `Prop110.lean`。
★本プロジェクトはこの命名規則で徹底されているので、
**実装しようとする命題の番号でファイルを探すのが最も確実**である。
★これで検索は 3 手: (1) `prop_<番号>_.*_of`、(2) `<構成の接頭辞>_<性質>`、
(3) **原典番号のファイル名**。

★2026-08-19、**12 件目**: `isIsotropicHull_comp_iso`(isotropic hull に同型を後置しても hull)を
書き下したが、**`Prop19.lean:238` に同名で在った**。
★原因: `isotropification_*` という**接頭辞でしか検索しなかった**。
`Prop19.lean` は `Proposition 1.9`(等長化)そのものなので、
★★**hull まわりの補題はまず `Prop19.lean` を通しで読む**。

★★11 件目の対策(原典番号でファイルを引く)は正しかったが、
**ファイルを開いても中を通しで読まなければ意味が無い** —— grep の接頭辞が外れると落ちる。
★対策の追記: 原典番号のファイルを見つけたら、**`grep '^theorem \|^def '` で宣言一覧を出す**。
