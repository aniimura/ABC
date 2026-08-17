import ABC3.Found.GenEll.PullbackBase

/-!
# [GenEll] Definition 1.1 —— **複素共役との両立**とアルキメデス側の底変換(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

## ★★★原文が要求する `ι_X` との両立の**理由**

原文 `Definition 1.1` はエルミート計量が**複素共役 `ι_X` と両立する**ことを
要求するが、★理由を書いていない。

★★★**理由は「高さの底変換不変性」である。**

`w : InfinitePlace K` を `F` に制限すると `v` になるが、
★★**`v.embedding` は `w.embedding|_F` **または その複素共役** に等しい**
(mathlib `NumberField.InfinitePlace.embedding_mk_eq`:
`embedding (mk φ) = φ ∨ embedding (mk φ) = conjugate φ`)。

★★★したがって Green 関数が複素共役で不変でなければ、
`archPoint` の値が定義体 `F` の取り方に依ってしまう。

## ★★これで「原文が書いていない理由」の 2 例目

| 原文の要求 | 理由 | 場所 |
|---|---|---|
| 次数を `[F:ℚ]` で正規化 | `≲` の定数を `F` に依らなくするため | `ArchBound.lean` |
| 計量が `ι_X` と両立 | ★★高さの**底変換不変性**のため | ★本ファイル |

★★どちらも **Gap(原文の飛躍)ではない**——原文は正しく、理由を省いていただけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★複素共役 -/

/-- ★**`Spec ℂ` の複素共役**。 -/
noncomputable def conjSpec : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of ℂ) :=
  Spec.map (CommRingCat.ofHom (starRingEnd ℂ))

/-- ★★**ℂ-点の複素共役**(原文の `ι_X`)。 -/
noncomputable def conjPoint {X : Scheme.{0}} (p : complexPoints X) : complexPoints X :=
  conjSpec ≫ p

/-- ★★**Green 関数が複素共役で不変**(原文「compatible with `ι_X`」)。

★★★これが無いと高さは定義体の取り方に依る(冒頭の docstring を参照)。 -/
def IsConjInvariant {X : Scheme.{0}} (g : GreenFn X) : Prop :=
  ∀ p : complexPoints X, g (conjPoint p) = g p

/-- ★**定数関数は複素共役不変**(非空虚性の witness)。 -/
theorem isConjInvariant_const {X : Scheme.{0}} (c : ℝ) :
    IsConjInvariant (fun _ : complexPoints X => c) := fun _ => rfl

/-- ★**複素共役不変な Green 関数の和も不変**(テンソル積で保たれる)。 -/
theorem IsConjInvariant.add {X : Scheme.{0}} {g h : GreenFn X}
    (hg : IsConjInvariant g) (hh : IsConjInvariant h) :
    IsConjInvariant (fun p => g p + h p) := fun p => by
  simp only [hg p, hh p]

/-! ## ★★★ℂ-点の底変換 -/

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]

/-- ★★**`archSpecMap` は底変換と合成できる**。

`Spec` は反変なので、`archSpecMap K w ≫ Spec.map φ = Spec.map (φ ≫ …)` になる。 -/
theorem archSpecMap_comp (w : InfinitePlace K)
    (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K)) :
    archSpecMap K w ≫ Spec.map φ
      = Spec.map (φ ≫ CommRingCat.ofHom (archRingHom K w)) := by
  rw [archSpecMap, ← Spec.map_comp]

/-- ★★★**埋め込みが制限で一致すれば、ℂ-点も一致する**。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★★`w.embedding|_F = v.embedding` のとき、`x_K` の `w` での ℂ-点は
`x_F` の `v` での ℂ-点に等しい。
★★★**一致しない場合は複素共役だけずれる**——そこで `IsConjInvariant` が効く。 -/
theorem archPoint_comp {X : Scheme.{0}} (v : InfinitePlace F) (w : InfinitePlace K)
    (xF : specRingOfIntegers F ⟶ X)
    (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K))
    (hφ : φ ≫ CommRingCat.ofHom (archRingHom K w)
      = CommRingCat.ofHom (archRingHom F v)) :
    archPoint (Spec.map φ ≫ xF) w = archPoint xF v := by
  rw [archPoint, archPoint, ← Category.assoc, archSpecMap_comp, hφ, archSpecMap]

/-! ## ★出典の紐付け(`.src`)

★条つきである。底変換不変性そのものには
「`w.embedding|_F` が `v.embedding` か その共役か」の場合分けと、
`degNormalized_baseChange` との接続が残っている。 -/

def IsConjInvariant.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(計量と複素共役 ι_X の両立——条件の定式化のみ)",
    sectionId := "genell-def-1-1-i" }

def archPoint_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(ℂ-点の底変換——埋め込みが一致する場合のみ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
