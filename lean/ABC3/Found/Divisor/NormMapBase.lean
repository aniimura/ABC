/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormalSections
import Mathlib.AlgebraicGeometry.Normalization
import Mathlib.AlgebraicGeometry.Pullbacks

/-!
# 底が動くときの相対正規化の射(鎖 `normalize` の `normalization-universal-normal`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism

## ★★★★★★測定の訂正(2026-08-25)—— `normalizationDesc` は**使えた**

台帳には長らく

> mathlib の `normalizationDesc` は「すでに `T` を経由する分解が与えられていれば
> 正規化がそこへ降りる」という**向きが逆**の形なので、そのままでは循環する

と書いてあったが、★**底変換(pullback)を 1 つ挟めば出る**。

```
Ψ := normalizationDesc f₂ ⟨lift⟩ (pullback.fst ψ f₁.fromNormalization) ⋯
       ≫ pullback.snd ψ f₁.fromNormalization
```

すなわち

* `f₁.fromNormalization : Y₁[L] ⟶ V₁` は**整射**(mathlib の instance)
* 整射は**底変換で保たれる**ので `pullback.fst ψ f₁.fromNormalization : T ⟶ V₂` も整射
* `f₂ : Spec L₂ ⟶ V₂` は `T` を経由する(仮定 `f₂ ≫ ψ = ε ≫ f₁` がちょうど lift の条件)
* よって `normalizationDesc` が `V₂[L₂] ⟶ T` を与え、`pullback.snd` で `V₁[L]` へ降りる

★★**正規性も支配性も要らない**。`Theorem 6.2, (i)` が要求するのはこの射だけなので、
「正規整スキームからの支配射は相対正規化を経由する」という**一般形は不要**である。

## ★★★これで消えた見積り

台帳は本節点を **★大**(「mathlib の `relativeGluingData` をもう一度書く」)と
見積もっていたが、実際は **15 行**だった。★★同じ型の見立て違いは **7 回目**である。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory Limits ABC3.Meta

universe u

variable {V₁ V₂ SpecL SpecL₂ : Scheme.{u}}

/-! ## ★1. 射の構成 -/

/-- ★★★★★★**底が動くときの相対正規化の射** ——
`ψ : V₂ ⟶ V₁` と `ε : Spec L₂ ⟶ Spec L` が両立すれば `V₂[L₂] ⟶ V₁[L]` が定まる。

★中身は `normalizationDesc` ＋ **底変換**の 2 本だけ。
`fromNormalization` が整射で、整射は底変換で保たれるのが要点である。 -/
noncomputable def normMapOfBase (f₁ : SpecL ⟶ V₁) (f₂ : SpecL₂ ⟶ V₂)
    [QuasiCompact f₁] [QuasiSeparated f₁] [QuasiCompact f₂] [QuasiSeparated f₂]
    (ψ : V₂ ⟶ V₁) (ε : SpecL₂ ⟶ SpecL) (hcompat : f₂ ≫ ψ = ε ≫ f₁) :
    f₂.normalization ⟶ f₁.normalization :=
  haveI : IsIntegralHom (pullback.fst ψ f₁.fromNormalization) := inferInstance
  Scheme.Hom.normalizationDesc f₂
      (pullback.lift f₂ (ε ≫ f₁.toNormalization) (by
        rw [Category.assoc, Scheme.Hom.toNormalization_fromNormalization]; exact hcompat))
      (pullback.fst ψ f₁.fromNormalization) (by rw [pullback.lift_fst])
    ≫ pullback.snd ψ f₁.fromNormalization

/-- ★★**`V₁` の上での可換性** —— 原文の四角形。 -/
@[reassoc] theorem normMapOfBase_fromNormalization (f₁ : SpecL ⟶ V₁) (f₂ : SpecL₂ ⟶ V₂)
    [QuasiCompact f₁] [QuasiSeparated f₁] [QuasiCompact f₂] [QuasiSeparated f₂]
    (ψ : V₂ ⟶ V₁) (ε : SpecL₂ ⟶ SpecL) (hcompat : f₂ ≫ ψ = ε ≫ f₁) :
    normMapOfBase f₁ f₂ ψ ε hcompat ≫ f₁.fromNormalization
      = f₂.fromNormalization ≫ ψ := by
  haveI : IsIntegralHom (pullback.fst ψ f₁.fromNormalization) := inferInstance
  rw [normMapOfBase, Category.assoc, ← pullback.condition, ← Category.assoc,
    Scheme.Hom.normalizationDesc_comp]

/-- ★★**生成点の側での可換性**。 -/
@[reassoc] theorem toNormalization_normMapOfBase (f₁ : SpecL ⟶ V₁) (f₂ : SpecL₂ ⟶ V₂)
    [QuasiCompact f₁] [QuasiSeparated f₁] [QuasiCompact f₂] [QuasiSeparated f₂]
    (ψ : V₂ ⟶ V₁) (ε : SpecL₂ ⟶ SpecL) (hcompat : f₂ ≫ ψ = ε ≫ f₁) :
    f₂.toNormalization ≫ normMapOfBase f₁ f₂ ψ ε hcompat
      = ε ≫ f₁.toNormalization := by
  haveI : IsIntegralHom (pullback.fst ψ f₁.fromNormalization) := inferInstance
  rw [normMapOfBase, ← Category.assoc, Scheme.Hom.toNormalization_normalizationDesc,
    pullback.lift_snd]

/-! ## ★2. 一意性 -/

/-- ★★★**2 つの可換性を満たす射は一意** ——
底変換の普遍性と `normalization.hom_ext` を合わせる。 -/
theorem normMapOfBase_unique (f₁ : SpecL ⟶ V₁) (f₂ : SpecL₂ ⟶ V₂)
    [QuasiCompact f₁] [QuasiSeparated f₁] [QuasiCompact f₂] [QuasiSeparated f₂]
    (ψ : V₂ ⟶ V₁) (ε : SpecL₂ ⟶ SpecL)
    (Ψ Ψ' : f₂.normalization ⟶ f₁.normalization)
    (hΨ : Ψ ≫ f₁.fromNormalization = f₂.fromNormalization ≫ ψ)
    (hΨ' : Ψ' ≫ f₁.fromNormalization = f₂.fromNormalization ≫ ψ)
    (htop : f₂.toNormalization ≫ Ψ = ε ≫ f₁.toNormalization)
    (htop' : f₂.toNormalization ≫ Ψ' = ε ≫ f₁.toNormalization) :
    Ψ = Ψ' := by
  haveI : IsIntegralHom (pullback.fst ψ f₁.fromNormalization) := inferInstance
  haveI : IsAffineHom (pullback.fst ψ f₁.fromNormalization) := inferInstance
  set L1 := pullback.lift (f := ψ) (g := f₁.fromNormalization)
    f₂.fromNormalization Ψ hΨ.symm with hL1
  set L2 := pullback.lift (f := ψ) (g := f₁.fromNormalization)
    f₂.fromNormalization Ψ' hΨ'.symm with hL2
  have key : L1 = L2 := by
    refine Scheme.Hom.normalization.hom_ext f₂ L1 L2 (pullback.fst ψ f₁.fromNormalization)
      ?_ ?_ ?_
    · refine pullback.hom_ext ?_ ?_
      · rw [Category.assoc, hL1, pullback.lift_fst, Category.assoc, hL2, pullback.lift_fst]
      · rw [Category.assoc, hL1, pullback.lift_snd, Category.assoc, hL2, pullback.lift_snd,
          htop, htop']
    · rw [hL1, pullback.lift_fst]
    · rw [hL2, pullback.lift_fst]
  have h1 : L1 ≫ pullback.snd ψ f₁.fromNormalization = Ψ := by rw [hL1, pullback.lift_snd]
  have h2 : L2 ≫ pullback.snd ψ f₁.fromNormalization = Ψ' := by rw [hL2, pullback.lift_snd]
  rw [← h1, ← h2, key]

/-! ## ★3. 関手性 —— `thm62-i-pull` が要求する 3 本 -/

/-- ★★**恒等**。 -/
theorem normMapOfBase_id {V SpecL : Scheme.{u}} (f : SpecL ⟶ V)
    [QuasiCompact f] [QuasiSeparated f] (h : f ≫ 𝟙 V = 𝟙 SpecL ≫ f) :
    normMapOfBase f f (𝟙 V) (𝟙 SpecL) h = 𝟙 (f.normalization) := by
  refine normMapOfBase_unique f f (𝟙 V) (𝟙 SpecL) _ _ ?_ ?_ ?_ ?_
  · exact normMapOfBase_fromNormalization _ _ _ _ _
  · rw [Category.id_comp, Category.comp_id]
  · exact toNormalization_normMapOfBase _ _ _ _ _
  · rw [Category.comp_id, Category.id_comp]

/-- ★★**合成**。 -/
theorem normMapOfBase_comp {V₁ V₂ V₃ SpecL₁ SpecL₂ SpecL₃ : Scheme.{u}}
    (f₁ : SpecL₁ ⟶ V₁) (f₂ : SpecL₂ ⟶ V₂) (f₃ : SpecL₃ ⟶ V₃)
    [QuasiCompact f₁] [QuasiSeparated f₁] [QuasiCompact f₂] [QuasiSeparated f₂]
    [QuasiCompact f₃] [QuasiSeparated f₃]
    (ψ₁ : V₂ ⟶ V₁) (ε₁ : SpecL₂ ⟶ SpecL₁) (h₁ : f₂ ≫ ψ₁ = ε₁ ≫ f₁)
    (ψ₂ : V₃ ⟶ V₂) (ε₂ : SpecL₃ ⟶ SpecL₂) (h₂ : f₃ ≫ ψ₂ = ε₂ ≫ f₂)
    (h₃ : f₃ ≫ (ψ₂ ≫ ψ₁) = (ε₂ ≫ ε₁) ≫ f₁) :
    normMapOfBase f₂ f₃ ψ₂ ε₂ h₂ ≫ normMapOfBase f₁ f₂ ψ₁ ε₁ h₁
      = normMapOfBase f₁ f₃ (ψ₂ ≫ ψ₁) (ε₂ ≫ ε₁) h₃ := by
  refine normMapOfBase_unique f₁ f₃ (ψ₂ ≫ ψ₁) (ε₂ ≫ ε₁) _ _ ?_ ?_ ?_ ?_
  · rw [Category.assoc, normMapOfBase_fromNormalization, ← Category.assoc,
      normMapOfBase_fromNormalization, Category.assoc]
  · exact normMapOfBase_fromNormalization _ _ _ _ _
  · rw [← Category.assoc, toNormalization_normMapOfBase, Category.assoc,
      toNormalization_normMapOfBase, ← Category.assoc]
  · exact toNormalization_normMapOfBase _ _ _ _ _

/-- ★★★★★**自然性の四角形** —— `thm62-i-pull` の関手性の背骨。

```
V₂[M₂] --π_M--> V₁[M]
   |              |
  n₂              n₁
   v              v
V₂[L₂] --π_L--> V₁[L]
```

★4 本とも `normMapOfBase` の実例なので、一意性(`normMapOfBase_unique`)で 1 度に出る。 -/
theorem normMapOfBase_naturality {V₁ V₂ SpecL SpecL₂ SpecM SpecM₂ : Scheme.{u}}
    (f₁ : SpecL ⟶ V₁) (f₂ : SpecL₂ ⟶ V₂) (g₁ : SpecM ⟶ V₁) (g₂ : SpecM₂ ⟶ V₂)
    [QuasiCompact f₁] [QuasiSeparated f₁] [QuasiCompact f₂] [QuasiSeparated f₂]
    [QuasiCompact g₁] [QuasiSeparated g₁] [QuasiCompact g₂] [QuasiSeparated g₂]
    (ψ : V₂ ⟶ V₁)
    (εL : SpecL₂ ⟶ SpecL) (hL : f₂ ≫ ψ = εL ≫ f₁)
    (εM : SpecM₂ ⟶ SpecM) (hM : g₂ ≫ ψ = εM ≫ g₁)
    (aL : SpecM ⟶ SpecL) (h₁ : g₁ ≫ 𝟙 V₁ = aL ≫ f₁)
    (aL₂ : SpecM₂ ⟶ SpecL₂) (h₂ : g₂ ≫ 𝟙 V₂ = aL₂ ≫ f₂)
    (hcube : aL₂ ≫ εL = εM ≫ aL) :
    normMapOfBase f₂ g₂ (𝟙 V₂) aL₂ h₂ ≫ normMapOfBase f₁ f₂ ψ εL hL
      = normMapOfBase g₁ g₂ ψ εM hM ≫ normMapOfBase f₁ g₁ (𝟙 V₁) aL h₁ := by
  refine normMapOfBase_unique f₁ g₂ ψ (aL₂ ≫ εL) _ _ ?_ ?_ ?_ ?_
  · rw [Category.assoc, normMapOfBase_fromNormalization, ← Category.assoc,
      normMapOfBase_fromNormalization, Category.assoc, Category.id_comp]
  · rw [Category.assoc, normMapOfBase_fromNormalization, ← Category.assoc,
      normMapOfBase_fromNormalization, Category.assoc, Category.comp_id]
  · rw [← Category.assoc, toNormalization_normMapOfBase, Category.assoc,
      toNormalization_normMapOfBase, ← Category.assoc]
  · rw [← Category.assoc, toNormalization_normMapOfBase, Category.assoc,
      toNormalization_normMapOfBase, ← Category.assoc, hcube]

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def normMapOfBase.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 支配射に沿った V[L] の引き戻し",
    sectionId := "frdi-thm-6-2" }

def normMapOfBase.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.normalizationDesc(★底変換を挟めば使える)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.normalizationDesc") 110,
    .citation "[mathlib]" "IsIntegralHom fromNormalization ＋ 整射は底変換で保たれる"
      (.inMathlib "AlgebraicGeometry.IsIntegralHom") 110,
    .citation "[mathlib]" "スキームの fiber 積"
      (.inMathlib "AlgebraicGeometry.Scheme.pullback") 110,
    .derivation
      "V₂ ×_{V₁} V₁[L] は V₂ の上で整。仮定 f₂ ≫ ψ = ε ≫ f₁ がちょうど lift の条件になる" 110,
    .implicitStep
      "★原文は仮定 (a) から因子と有理函数の引き戻しを「(i)」の一言で置いている" 110 ]

def normMapOfBase_unique.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 正規化への射は V₁ 上で一意",
    sectionId := "frdi-thm-6-2" }

def normMapOfBase_unique.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.normalization.hom_ext"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.normalization.hom_ext") 110,
    .derivation "底変換の普遍性で pullback への lift の一意性に帰着させる" 110 ]

def normMapOfBase_naturality.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — V[L] の引き戻しの自然性",
    sectionId := "frdi-thm-6-2" }

def normMapOfBase_naturality.needs : List ProofObligation :=
  [ .citation "[ABC3]" "normMapOfBase_unique(一意性)"
      (.inProject "ABC3" "ABC3.Found.Divisor.normMapOfBase_unique") 110,
    .derivation "四角形の 4 辺がどれも normMapOfBase の実例なので、一意性 1 本で出る" 110 ]

def normMapOfBase_comp.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — V[L] の引き戻しの合成",
    sectionId := "frdi-thm-6-2" }

def normMapOfBase_comp.needs : List ProofObligation :=
  [ .citation "[ABC3]" "normMapOfBase_unique"
      (.inProject "ABC3" "ABC3.Found.Divisor.normMapOfBase_unique") 110 ]

def normMapOfBase_id.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — V[L] の引き戻しの恒等",
    sectionId := "frdi-thm-6-2" }

def normMapOfBase_id.needs : List ProofObligation :=
  [ .citation "[ABC3]" "normMapOfBase_unique"
      (.inProject "ABC3" "ABC3.Found.Divisor.normMapOfBase_unique") 110 ]

end ABC3.Found.Divisor
