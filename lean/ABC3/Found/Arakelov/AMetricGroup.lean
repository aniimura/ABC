/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricIso
import ABC3.Meta.Claim

/-!
# 算術直線束の**結合律と単位律**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★原文の `APic(X)` が群であることの欄

    `IsIsometry ((L * M) * N) (L * (M * N)) (α_ …)`   （`isIsometry_mul_assoc`）
    `IsIsometry (L * 1) L (ρ_ L.sheaf)`               （`isIsometry_mul_one`）

★どちらも `isIsometry_mul`（`AMetricIso.lean`）と同じ機構である
——3 つ（2 つ）が同時に自明になるチャートへ降り、`IsTensorOf` で `h` を分解し、
`transUnit_pullTriv` で遷移単元が変わらないことを使う。
★★最後に残るのは実数の等式（`(a*b)*c = a*(b*c)`、`a*1 = a`）だけである。

## ★★★★★★★★自明化の側の 2 つの恒等式

    `pullTriv (α_ A B C) V (tensorTriv eA (tensorTriv eB eC)) = tensorTriv (tensorTriv eA eB) eC`
    `pullTriv (ρ_ A) V eA = tensorTriv eA (baseTriv X V)`

★どちらも **`rfl` ではない**（2026-08-28 実測）。結合律の側は 4 段:

| 段 | 根拠 |
|---|---|
| 関手のモノイダル性 | `Functor.OplaxMonoidal.oplax_associativity` |
| 結合子の自然性 | `MonoidalCategory.associator_naturality` |
| 三角形恒等式 | `MonoidalCategory.triangle` |
| 単位子の一致 | `MonoidalCategory.unitors_equal` |

★★`simp` に任せると `λ_` を関手の lax 構造へ展開して**逆方向へ進む**ので、
テンソル積を分解する 4 つの書き換え（`tensorHom_comp_*`）を作って手で並べた。

## ★残っている段（明示）

★★**逆元（双対計量）** と **商そのもの**である。
★★★`Definition 1.1` の項目全体には (ii) の `deg_F`
（台帳 `arakelov-degF-finite-places`）も要る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-! ## ★テンソル積を分解する 4 つの書き換え -/

theorem tensorHom_comp_srcWhisker {C : Type*} [Category C] [MonoidalCategory C]
    {P Q R S T : C} (f : P ⟶ Q) (g : R ⟶ S) (h : S ⟶ T) :
    f ⊗ₘ (g ≫ h) = (P ◁ g) ≫ (f ⊗ₘ h) := by
  rw [← MonoidalCategory.id_tensorHom, MonoidalCategory.tensorHom_comp_tensorHom,
    Category.id_comp]

theorem tensorHom_comp_tgtWhisker {C : Type*} [Category C] [MonoidalCategory C]
    {P Q R S T : C} (f : P ⟶ Q) (g : R ⟶ S) (h : S ⟶ T) :
    f ⊗ₘ (g ≫ h) = (f ⊗ₘ g) ≫ (Q ◁ h) := by
  rw [← MonoidalCategory.id_tensorHom, MonoidalCategory.tensorHom_comp_tensorHom,
    Category.comp_id]

theorem tensorHom_comp_srcWhiskerRight {C : Type*} [Category C] [MonoidalCategory C]
    {P Q R S T : C} (f : P ⟶ Q) (g : Q ⟶ R) (h : S ⟶ T) :
    (f ≫ g) ⊗ₘ h = (f ▷ S) ≫ (g ⊗ₘ h) := by
  rw [← MonoidalCategory.tensorHom_id, MonoidalCategory.tensorHom_comp_tensorHom,
    Category.id_comp]

theorem tensorHom_comp_tgtWhiskerRight {C : Type*} [Category C] [MonoidalCategory C]
    {P Q R S T : C} (f : P ⟶ Q) (g : Q ⟶ R) (h : S ⟶ T) :
    (f ≫ g) ⊗ₘ h = (f ⊗ₘ h) ≫ (g ▷ T) := by
  rw [← MonoidalCategory.tensorHom_id, MonoidalCategory.tensorHom_comp_tensorHom,
    Category.comp_id]

/-- ★★★★★★★★**結合律の自明化版**。

★機構は 4 段:
`oplax_associativity`（関手のモノイダル性）→ `associator_naturality` →
`triangle` → `unitors_equal`。 -/
theorem pullTriv_associator {A B C : X.PresheafOfModules} (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (eC : (restrictPresheafFunctor X V).obj C ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullTriv (α_ A B C) V (tensorTriv eA (tensorTriv eB eC))
      = tensorTriv (tensorTriv eA eB) eC := by
  ext1
  show ((restrictPresheafFunctor X V).mapIso (α_ A B C)).hom
      ≫ (tensorTriv eA (tensorTriv eB eC)).hom
    = (tensorTriv (tensorTriv eA eB) eC).hom
  simp only [tensorTriv, restrictPresheafTensor, Iso.trans_hom, Iso.symm_hom, tensorIso_hom,
    Functor.Monoidal.μIso_inv, Functor.mapIso_hom, Category.assoc]
  rw [tensorHom_comp_srcWhisker eA.hom
      (Functor.OplaxMonoidal.δ (restrictPresheafFunctor X V) B C) _,
    tensorHom_comp_srcWhiskerRight
      (Functor.OplaxMonoidal.δ (restrictPresheafFunctor X V) A B) _ eC.hom]
  rw [Category.assoc, Category.assoc, ← Functor.OplaxMonoidal.associativity_assoc]
  congr 2
  rw [tensorHom_comp_tgtWhisker eA.hom (eB.hom ⊗ₘ eC.hom)
      (λ_ (𝟙_ (PresheafModulesOn X V))).hom,
    tensorHom_comp_tgtWhiskerRight (eA.hom ⊗ₘ eB.hom)
      (λ_ (𝟙_ (PresheafModulesOn X V))).hom eC.hom]
  rw [Category.assoc, Category.assoc, ← MonoidalCategory.associator_naturality_assoc]
  congr 1
  rw [MonoidalCategory.triangle_assoc, ← MonoidalCategory.unitors_equal]

/-- ★★★★**右単位子は「基準の自明化とのテンソル」である**。

★機構は mathlib の `right_unitality_hom` ＋ `unitors_equal` ＋ `rightUnitor_naturality`。 -/
theorem pullTriv_rightUnitor {A : X.PresheafOfModules} (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullTriv (ρ_ A) V eA = tensorTriv eA (baseTriv X V) := by
  ext1
  show ((restrictPresheafFunctor X V).mapIso (ρ_ A)).hom ≫ eA.hom
    = (tensorTriv eA (baseTriv X V)).hom
  simp only [tensorTriv, baseTriv, restrictPresheafTensor, restrictPresheafUnit,
    Iso.trans_hom, Iso.symm_hom, tensorIso_hom, Functor.Monoidal.μIso_inv,
    Category.assoc, MonoidalCategory.tensorHom_def']
  rw [MonoidalCategory.unitors_equal, MonoidalCategory.rightUnitor_naturality]
  have hru : Functor.OplaxMonoidal.δ (restrictPresheafFunctor X V) A (𝟙_ X.PresheafOfModules)
        ≫ ((restrictPresheafFunctor X V).obj A
            ◁ (Functor.Monoidal.εIso (restrictPresheafFunctor X V)).inv)
        ≫ (ρ_ ((restrictPresheafFunctor X V).obj A)).hom
        ≫ eA.hom
      = (restrictPresheafFunctor X V).map (ρ_ A).hom ≫ eA.hom :=
    Functor.OplaxMonoidal.right_unitality_hom_assoc
      (F := restrictPresheafFunctor X V) A eA.hom
  rw [hru]
  rfl

/-- ★★★★★★★★**単位律** `L̄ ⊗ Ō_X ≅ L̄`（等長）。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★`Ō_X` の基準ノルムは基準の自明化で `1` なので、テンソル積の `h` は変わらない。 -/
theorem isIsometry_mul_one (L : AMetric X) :
    IsIsometry (L * 1) L (ρ_ L.sheaf) := by
  intro V e p hp
  have hten : (L * 1).metric.h V (tensorTriv e (baseTriv X V)) p
      = L.metric.h V e p * (1 : AMetric X).metric.h V (baseTriv X V) p :=
    isTensorOf_tensor L.triv (1 : AMetric X).triv L.metric (1 : AMetric X).metric V e
      (baseTriv X V) p hp
  have hone : (1 : AMetric X).metric.h V (baseTriv X V) p = 1 :=
    structLocalMetric_h_baseTriv X V p hp
  have hpull : pullTriv (ρ_ L.sheaf) V e = tensorTriv e (baseTriv X V) :=
    pullTriv_rightUnitor V e
  have hstep : (L * 1).metric.h V (pullTriv (ρ_ L.sheaf) V e) p
      = (L * 1).metric.h V (tensorTriv e (baseTriv X V)) p :=
    congrArg (fun t => (L * 1).metric.h V t p) hpull
  exact hstep.trans (hten.trans (by rw [hone, mul_one]))

/-- ★★**したがって同値類の上で `[L̄] · 1 = [L̄]`**。 -/
theorem isometric_mul_one (L : AMetric X) : Isometric (L * 1) L :=
  ⟨ρ_ L.sheaf, isIsometry_mul_one L⟩

/-! ## ★★★★★★★★★結合律 -/

/-- ★**3 つが同時に自明になるチャート**。 -/
structure TripleChart (A B C : X.PresheafOfModules) (V : X.Opens)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) where
  /-- 小さい開集合。 -/
  W : X.Opens
  /-- `V` の中にある。 -/
  hWV : W ≤ V
  /-- `p` を含む。 -/
  hpW : p ⁻¹ᵁ W = ⊤
  /-- `A` の自明化。 -/
  eA : (restrictPresheafFunctor X W).obj A ≅ 𝟙_ (PresheafModulesOn X W)
  /-- `B` の自明化。 -/
  eB : (restrictPresheafFunctor X W).obj B ≅ 𝟙_ (PresheafModulesOn X W)
  /-- `C` の自明化。 -/
  eC : (restrictPresheafFunctor X W).obj C ≅ 𝟙_ (PresheafModulesOn X W)

/-- ★★**3 重チャートも取れる**——2 重を 2 回取って細める。 -/
theorem nonempty_tripleChart {A B C : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B) (hC : IsLocallyTrivial X C)
    (V : X.Opens) (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) :
    Nonempty (TripleChart A B C V p) := by
  obtain ⟨c₁⟩ := nonempty_tensorChart hA hB V p hp
  obtain ⟨c₂⟩ := nonempty_tensorChart (hA.tensor hB) hC c₁.W p c₁.hpW
  exact ⟨⟨c₂.W, c₂.hWV.trans c₁.hWV, c₂.hpW,
    trivialOfLe A c₂.hWV c₁.eA, trivialOfLe B c₂.hWV c₁.eB, c₂.eB⟩⟩

/-- ★★★★★★★★★**結合律** `(L̄ ⊗ M̄) ⊗ N̄ ≅ L̄ ⊗ (M̄ ⊗ N̄)`（等長）。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は `isIsometry_mul` と同じ——3 つが同時に自明になるチャートへ降り、
`IsTensorOf` で `h` を分解し、`transUnit_pullTriv` で遷移単元が変わらないことを使う。
★★最後に残るのは実数の結合律 `(a*b)*c = a*(b*c)` である。 -/
theorem isIsometry_mul_assoc (L M N : AMetric X) :
    IsIsometry ((L * M) * N) (L * (M * N)) (α_ L.sheaf M.sheaf N.sheaf) := by
  intro V f p hp
  obtain ⟨c⟩ := nonempty_tripleChart L.triv M.triv N.triv V p hp
  have hpW := c.hpW
  set g := trivialOfLe (L * (M * N)).sheaf c.hWV f with hgdef
  have hfac : ‖evalOn p c.W hpW (transUnit (L * (M * N)).sheaf c.W
      (tensorTriv c.eA (tensorTriv c.eB c.eC)) g)‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p c.W hpW (isUnit_transUnit _ c.W _ _))
  -- 右辺: `L ⊗ (M ⊗ N)`
  have hR : (L * (M * N)).metric.h c.W g p
      * ‖evalOn p c.W hpW (transUnit (L * (M * N)).sheaf c.W
          (tensorTriv c.eA (tensorTriv c.eB c.eC)) g)‖
      = L.metric.h c.W c.eA p * (M.metric.h c.W c.eB p * N.metric.h c.W c.eC p) := by
    rw [(L * (M * N)).metric.compat c.W (tensorTriv c.eA (tensorTriv c.eB c.eC)) g p hpW]
    have h1 : (L * (M * N)).metric.h c.W (tensorTriv c.eA (tensorTriv c.eB c.eC)) p
        = L.metric.h c.W c.eA p * (M * N).metric.h c.W (tensorTriv c.eB c.eC) p :=
      isTensorOf_tensor L.triv (M * N).triv L.metric (M * N).metric c.W c.eA
        (tensorTriv c.eB c.eC) p hpW
    have h2 : (M * N).metric.h c.W (tensorTriv c.eB c.eC) p
        = M.metric.h c.W c.eB p * N.metric.h c.W c.eC p :=
      isTensorOf_tensor M.triv N.triv M.metric N.metric c.W c.eB c.eC p hpW
    rw [h1, h2]
  -- 左辺: `(L ⊗ M) ⊗ N`
  have e2 := ((L * M) * N).metric.compat c.W
    (pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W (tensorTriv c.eA (tensorTriv c.eB c.eC)))
    (pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W g) p hpW
  have e1 : transUnit ((L * M) * N).sheaf c.W
        (pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W (tensorTriv c.eA (tensorTriv c.eB c.eC)))
        (pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W g)
      = transUnit (L * (M * N)).sheaf c.W (tensorTriv c.eA (tensorTriv c.eB c.eC)) g :=
    transUnit_pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W _ _
  rw [e1] at e2
  have eassoc : pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W
      (tensorTriv c.eA (tensorTriv c.eB c.eC))
      = tensorTriv (tensorTriv c.eA c.eB) c.eC :=
    pullTriv_associator c.W c.eA c.eB c.eC
  have e3 : ((L * M) * N).metric.h c.W (pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W
        (tensorTriv c.eA (tensorTriv c.eB c.eC))) p
      = L.metric.h c.W c.eA p * (M.metric.h c.W c.eB p * N.metric.h c.W c.eC p) := by
    have h1 : ((L * M) * N).metric.h c.W (tensorTriv (tensorTriv c.eA c.eB) c.eC) p
        = (L * M).metric.h c.W (tensorTriv c.eA c.eB) p * N.metric.h c.W c.eC p :=
      isTensorOf_tensor (L * M).triv N.triv (L * M).metric N.metric c.W
        (tensorTriv c.eA c.eB) c.eC p hpW
    have h2 : (L * M).metric.h c.W (tensorTriv c.eA c.eB) p
        = L.metric.h c.W c.eA p * M.metric.h c.W c.eB p :=
      isTensorOf_tensor L.triv M.triv L.metric M.metric c.W c.eA c.eB p hpW
    refine (congrArg (fun t => ((L * M) * N).metric.h c.W t p) eassoc).trans ?_
    refine h1.trans ?_
    rw [h2, mul_assoc]
  have hkey : ((L * M) * N).metric.h c.W (pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W g) p
      = (L * (M * N)).metric.h c.W g p :=
    mul_right_cancel₀ hfac ((e2.trans e3).trans hR.symm)
  have h1 : ((L * M) * N).metric.h V (pullTriv (α_ L.sheaf M.sheaf N.sheaf) V f) p
      = ((L * M) * N).metric.h c.W
        (trivialOfLe ((L * M) * N).sheaf c.hWV
          (pullTriv (α_ L.sheaf M.sheaf N.sheaf) V f)) p :=
    (((L * M) * N).metric.restrict c.hWV
      (pullTriv (α_ L.sheaf M.sheaf N.sheaf) V f) p hpW).symm
  have h2 : (L * (M * N)).metric.h V f p = (L * (M * N)).metric.h c.W g p :=
    ((L * (M * N)).metric.restrict c.hWV f p hpW).symm
  have h3 : ((L * M) * N).metric.h c.W
        (trivialOfLe ((L * M) * N).sheaf c.hWV
          (pullTriv (α_ L.sheaf M.sheaf N.sheaf) V f)) p
      = ((L * M) * N).metric.h c.W (pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.W g) p :=
    congrArg (fun t => ((L * M) * N).metric.h c.W t p)
      (trivialOfLe_pullTriv (α_ L.sheaf M.sheaf N.sheaf) c.hWV f)
  exact h1.trans (h3.trans (hkey.trans h2.symm))

/-- ★★**したがって同値類の上で結合律が成り立つ**。 -/
theorem isometric_mul_assoc (L M N : AMetric X) : Isometric ((L * M) * N) (L * (M * N)) :=
  ⟨α_ L.sheaf M.sheaf N.sheaf, isIsometry_mul_assoc L M N⟩

/-! ### ★出典の紐付け(`.src`) -/

def pullTriv_rightUnitor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(右単位子は基準の自明化とのテンソルであること)",
    sectionId := "genell-def-1-1-i" }

def pullTriv_associator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(結合子の自明化版)",
    sectionId := "genell-def-1-1-i" }

def isIsometry_mul_assoc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(結合律 (L̄ ⊗ M̄) ⊗ N̄ ≅ L̄ ⊗ (M̄ ⊗ N̄)——等長)",
    sectionId := "genell-def-1-1-i" }

def isIsometry_mul_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(単位律 L̄ ⊗ Ō_X ≅ L̄——等長)",
    sectionId := "genell-def-1-1-i" }

def isIsometry_mul_one.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isTensorOf_tensor(構成した計量がテンソル積の計量であること)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isTensorOf_tensor") 3,
    .citation "[ABC3]" "structLocalMetric_h_baseTriv(Ō_X は基準の自明化で h = 1)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.structLocalMetric_h_baseTriv") 3,
    .citation "[mathlib]" "Functor.OplaxMonoidal.right_unitality_hom"
      (.inMathlib "CategoryTheory.Functor.OplaxMonoidal.right_unitality_hom") 3,
    .implicitStep
      ("★★★★**結合律も 2026-08-28 に入った**——" ++
       "pullTriv (α_ A B C) V (tensorTriv eA (tensorTriv eB eC)) " ++
       "= tensorTriv (tensorTriv eA eB) eC は rfl ではなかった(2026-08-28 実測)。" ++
       "★機構は 4 段(oplax_associativity ＋ associator_naturality ＋ " ++
       "triangle ＋ unitors_equal)で、simp は λ_ を関手の lax 構造へ展開して" ++
       "逆方向へ進むので、テンソル積を分解する 4 つの書き換えを作って手で並べた") 3,
    .implicitStep
      ("★★★逆元(双対計量)と商そのものもまだである。" ++
       "★Definition 1.1 の項目全体には (ii) の deg_F も要る") 3 ]

end ABC3.Found.Arakelov
