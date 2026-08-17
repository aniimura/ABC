import ABC3.Found.GenEll.IdealStalk

/-!
# [GenEll] Proposition 1.4 の足場 —— 茎による**引き戻し**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> tensor product of arithmetic line bundles on X. The

## ★★★なぜこれが `Proposition 1.4, (i)` の鍵なのか

`ht` の加法性 `ht_{L̄⊗M̄} = ht_{L̄} + ht_{M̄}` を**構成から**示すには、
引き戻しが積を保つ必要がある:

  `x_F^*(D + E) = x_F^* D + x_F^* E`

★**mathlib 実測(2026-08-17)**: `Scheme.IdealSheafData.comap` に
**積の保存は無い**(`IdealSheaf/Functorial.lean` を `mul` で検索して 0 件)。
そもそも `comap I f := (pullback.fst f I.subschemeι).ker` と
**ファイバー積の核**で定義されているので、積との両立は自明でない。

★★**茎で引き戻すと出る。** `f.stalkMap x` は**環準同型**なので、
`Ideal.map_mul` 1 本で積を保つ。

## ★★★「x が D を通らない」という仮定の正体

`Proposition 1.4` は `x_F` が因子を通らないことを要求する。
★本ファイルの `isEffectiveCartierAt_stalkPullback_of_ne_bot` が示すのは:

  **`Spec 𝓞_F` の茎は整域なので、引き戻したイデアルが `≠ ⊥` でありさえすれば
  自動的に有効 Cartier になる。**

★★整域では**非零元はすべて非零因子**だからである。
つまり「通らない」という幾何的条件は、**`≠ ⊥` という 1 本の条件**に潰れる。

## ★mathlib 実測(2026-08-17)

`Scheme.Hom.stalkMap f x : Y.presheaf.stalk (f x) ⟶ X.presheaf.stalk x`
(`Scheme.lean:237`)は在る。`stalkMap_id` / `stalkMap_comp` も在る。
★**イデアルの茎ごとの引き戻しは無い**——本ファイルで作る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

/-! ## ★点ごとの有効 Cartier 性 -/

/-- ★**点ごとの有効 Cartier 性** —— イデアルが**非零因子 1 つ**で生成される。

★`IsEffectiveCartierStalk` を各点に分解したもの。
★★引き戻しは `IdealSheafData` を返さない(茎のイデアルを返す)ので、
この分解が無いと引き戻し先で有効 Cartier 性を述べられない。 -/
def IsEffectiveCartierAt {R : Type*} [CommRing R] (I : Ideal R) : Prop :=
  ∃ f : R, I = Ideal.span {f} ∧ f ∈ nonZeroDivisors R

/-- ★分解は定義的に等しい。 -/
theorem isEffectiveCartierStalk_iff {X : Scheme} (I : X.IdealSheafData) :
    IsEffectiveCartierStalk I ↔ ∀ x : X, IsEffectiveCartierAt (I.stalk x) :=
  Iff.rfl

/-- ★点ごとの有効 Cartier 性は積で閉じる。 -/
theorem IsEffectiveCartierAt.mul {R : Type*} [CommRing R] {I J : Ideal R}
    (hI : IsEffectiveCartierAt I) (hJ : IsEffectiveCartierAt J) :
    IsEffectiveCartierAt (I * J) := by
  obtain ⟨f, hf, hfn⟩ := hI
  obtain ⟨g, hg, hgn⟩ := hJ
  exact ⟨f * g, by rw [hf, hg, Ideal.span_singleton_mul_span_singleton],
    Submonoid.mul_mem _ hfn hgn⟩

/-- ★全体は点ごとの有効 Cartier。 -/
theorem isEffectiveCartierAt_top {R : Type*} [CommRing R] :
    IsEffectiveCartierAt (⊤ : Ideal R) :=
  ⟨1, (Ideal.span_singleton_one).symm, Submonoid.one_mem _⟩

/-! ## ★★茎による引き戻し -/

variable {X Y Z : Scheme}

/-- ★★**イデアル層の茎ごとの引き戻し**。

`f : X ⟶ Y` と `I : Y.IdealSheafData` に対し、`x : X` での茎を
`f.stalkMap x : Y.presheaf.stalk (f x) ⟶ X.presheaf.stalk x` で押し出す。

★mathlib に無い(2026-08-17 実測)。
★★`comap` と違って**環準同型による像**なので、代数的な性質がそのまま通る。 -/
noncomputable def stalkPullback (I : Y.IdealSheafData) (f : X ⟶ Y) (x : X) :
    Ideal (X.presheaf.stalk x) :=
  (I.stalk (f.base x)).map (f.stalkMap x).hom

/-- ★★★**引き戻しは積を保つ** —— `x^*(I·J) = x^*I · x^*J`。

原文 (GenEll p.3):
> tensor product of arithmetic line bundles on X. The

★★これが `Proposition 1.4, (i)`(`ht` の加法性)の代数的な芯である。
★`comap` では取れなかった(mathlib に積の保存が無く、
ファイバー積の核という定義から作るのも自明でない)。 -/
theorem stalkPullback_mul (I J : Y.IdealSheafData) (f : X ⟶ Y) (x : X) :
    stalkPullback (I * J) f x = stalkPullback I f x * stalkPullback J f x := by
  rw [stalkPullback, Scheme.IdealSheafData.stalk_mul, Ideal.map_mul, stalkPullback,
    stalkPullback]

/-- ★**空因子の引き戻しは空因子**。 -/
@[simp] theorem stalkPullback_top (f : X ⟶ Y) (x : X) :
    stalkPullback (⊤ : Y.IdealSheafData) f x = ⊤ := by
  rw [stalkPullback, Scheme.IdealSheafData.stalk_top, Ideal.map_top]

/-- ★**恒等射による引き戻しは何もしない**。 -/
@[simp] theorem stalkPullback_id (I : X.IdealSheafData) (x : X) :
    stalkPullback I (𝟙 X) x = I.stalk x := by
  rw [stalkPullback, Scheme.Hom.stalkMap_id]
  -- ★`CommRingCat.Hom.hom (𝟙 _)` と `RingHom.id` は定義的にのみ等しいので `erw`
  erw [Ideal.map_id]
  rfl

/-- ★**引き戻しは合成と両立する**(関手性)。 -/
theorem stalkPullback_comp (I : Z.IdealSheafData) (f : X ⟶ Y) (g : Y ⟶ Z) (x : X) :
    stalkPullback I (f ≫ g) x
      = (stalkPullback I g (f.base x)).map (f.stalkMap x).hom := by
  rw [stalkPullback, stalkPullback, Ideal.map_map, ← CommRingCat.hom_comp,
    ← Scheme.Hom.stalkMap_comp]
  rfl

/-! ## ★★★「x が D を通らない」の正体 -/

/-- ★★★**整域の茎では「引き戻しが `≠ ⊥`」だけで有効 Cartier になる**。

`Proposition 1.4` の「`x_F` は因子 `D` を通らない」という幾何的な仮定は、
★★**`Spec 𝓞_F` の茎が整域である**ことにより
**`x^* D ≠ ⊥` という 1 本の条件**に潰れる。
整域では**非零元はすべて非零因子**だからである。

★`Spec 𝓞_F` の茎は `𝓞_F` の局所化なので整域であり、この仮定は満たされる。 -/
theorem isEffectiveCartierAt_stalkPullback_of_ne_bot
    {x : X} [IsDomain (X.presheaf.stalk x)] {I : Y.IdealSheafData} {f : X ⟶ Y}
    (hI : IsEffectiveCartierStalk I) (hne : stalkPullback I f x ≠ ⊥) :
    IsEffectiveCartierAt (stalkPullback I f x) := by
  obtain ⟨a, ha, -⟩ := hI (f.base x)
  refine ⟨(f.stalkMap x).hom a, ?_, ?_⟩
  · rw [stalkPullback, ha, Ideal.map_span, Set.image_singleton]
  · refine mem_nonZeroDivisors_of_ne_zero ?_
    intro hc
    exact hne (by rw [stalkPullback, ha, Ideal.map_span, Set.image_singleton, hc,
      Ideal.span_singleton_eq_bot.2 rfl])

/-- ★**引き戻しは有効 Cartier 性を積ごと保つ** —— 2 つが `≠ ⊥` なら積も。

★★これが `ht` の加法性を「構成から」出すための最後の環である:
`x^*(D+E)` が有効 Cartier であり、かつ `x^*D · x^*E` に等しい。 -/
theorem isEffectiveCartierAt_stalkPullback_mul
    {x : X} [IsDomain (X.presheaf.stalk x)] {I J : Y.IdealSheafData} {f : X ⟶ Y}
    (hI : IsEffectiveCartierStalk I) (hJ : IsEffectiveCartierStalk J)
    (hIne : stalkPullback I f x ≠ ⊥) (hJne : stalkPullback J f x ≠ ⊥) :
    IsEffectiveCartierAt (stalkPullback (I * J) f x)
      ∧ stalkPullback (I * J) f x = stalkPullback I f x * stalkPullback J f x := by
  refine ⟨?_, stalkPullback_mul I J f x⟩
  rw [stalkPullback_mul]
  exact (isEffectiveCartierAt_stalkPullback_of_ne_bot hI hIne).mul
    (isEffectiveCartierAt_stalkPullback_of_ne_bot hJ hJne)

/-! ## ★出典の紐付け(`.src`) -/

def stalkPullback_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Proposition 1.4, (i)(茎の水準で引き戻しが積を保つことのみ)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
