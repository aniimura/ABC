/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Thm62Pull
import ABC3.Found.Divisor.Ex61Model
import ABC3.Found.FrdI.Thm52Change

/-!
# `Theorem 6.2, (i)` —— 幾何のデータの射から `Ψ : 𝒞₁ ⥤ 𝒞₂` を作る

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism

## ★このファイルで閉じること

`Thm62Pull.lean` で作った 5 つの部品(`phiPullEff` / `bPullBaseAdd` /
`pullCoeff_square` / `ffMap_square` / `pullCoeff_weilDivOfFn`)を
`ModelDataHomOver` の 5 フィールドへ**束ねる**。

| 節 | 内容 |
|---|---|
| ★1 | `Φ^gp ≃ Γ` は `restrictEff` と両立する(`effSubGpHom_gpMap_restrictEff`) |
| ★2 | `divCompat` の中身(`divBHom_bPullBase`) |
| ★3 | 自然性 3 本(`phiPullEff_square` / `bPullBaseAdd_square` / `divBGeom_bPullBaseAdd`) |
| ★4 | 入力データ `GeomPullDatum` |
| ★5 | ★★★★★★`ModelDataHomOver` と `Ψ : 𝒞₁ ⥤ 𝒞₂` |

★★**要は `hpull`(原文の仮定 (a))の対偶**である ——
「`D_{L₂}` の外の点は `D_L` の外へ行く」が
因子側の「引き戻しの台が `D_{L₂}` に収まる」と
有理函数側の「そこでは `u` が茎の単元だから `ord = 0` が保たれる」を**同時に**与える。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI ABC3.Meta

universe u

/-! ## ★1. `Φ^gp ≃ Γ` は `restrictEff` と両立する -/

/-- ★`restrictEff` は `effSubIncl` を通すと元の射そのもの。 -/
theorem effSubIncl_restrictEff {S T : Type u} {Γ₁ : AddSubgroup (S →₀ ℤ)}
    {Γ₂ : AddSubgroup (T →₀ ℤ)} (φ : Γ₁ →+ Γ₂)
    (hnn : ∀ a : Γ₁, 0 ≤ (a : S →₀ ℤ) → 0 ≤ ((φ a : Γ₂) : T →₀ ℤ))
    (a : effSub Γ₁) :
    effSubIncl Γ₂ (restrictEff φ hnn a) = φ (effSubIncl Γ₁ a) := Subtype.ext rfl

/-- ★★★**`Φ^gp ≃ Γ` は `restrictEff` と両立する** ——
`divCompat` が `Φ^gp` の側で書かれているので要る。 -/
theorem effSubGpHom_gpMap_restrictEff {S T : Type u} {Γ₁ : AddSubgroup (S →₀ ℤ)}
    {Γ₂ : AddSubgroup (T →₀ ℤ)} (φ : Γ₁ →+ Γ₂)
    (hnn : ∀ a : Γ₁, 0 ≤ (a : S →₀ ℤ) → 0 ≤ ((φ a : Γ₂) : T →₀ ℤ))
    (z : Gp (effSub Γ₁)) :
    effSubGpHom Γ₂ (gpMap _ (restrictEff φ hnn) z) = φ (effSubGpHom Γ₁ z) := by
  obtain ⟨a, b, rfl⟩ := gp_eq_sub z
  simp only [map_sub, gpMap_toGp, effSubGpHom, gpLift_toGp, effSubIncl_restrictEff]

/-! ## ★2. `divCompat` の中身 -/

section DivCompat

variable {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
  {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
  [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
  [IsGalois V₁.functionField Kbar₁] [IsGalois V₂.functionField Kbar₂]
  [∀ L : FinSub V₁.functionField Kbar₁, IsLocallyNoetherian (normObj V₁ L)]
  [∀ L : FinSub V₁.functionField Kbar₁, CompactSpace (normObj V₁ L)]
  [∀ L : FinSub V₂.functionField Kbar₂, IsLocallyNoetherian (normObj V₂ L)]
  [∀ L : FinSub V₂.functionField Kbar₂, CompactSpace (normObj V₂ L)]

variable (DK₁ : Set (PrimeDivisorPt V₁)) (DK₂ : Set (PrimeDivisorPt V₂))
  (L : FinSub V₁.functionField Kbar₁) (FL : FinSub V₂.functionField Kbar₂)
  (π : normObj V₂ FL ⟶ normObj V₁ L) [IsDominant π]
  (hdim : ∀ w : PrimeDivisorPt (normObj V₂ FL),
    ringKrullDim ((normObj V₁ L).presheaf.stalk (π.base w.1)) ≤ 1)
  (hpull : ∀ w : PrimeDivisorPt (normObj V₂ FL), w ∉ DLSet V₂ DK₂ FL →
    ∀ hc : IsCodimOnePt (normObj V₁ L) (π.base w.1),
      (⟨π.base w.1, hc⟩ : PrimeDivisorPt (normObj V₁ L)) ∉ DLSet V₁ DK₁ L)

omit [IsGalois V₁.functionField Kbar₁] [IsGalois V₂.functionField Kbar₂] in
/-- ★★★★★★**`div` は底が動く引き戻しと両立する** —— `divCompat` の中身。

★中身は 1 行 —— **主因子は `⊤` の上で `u` 自身を局所方程式に持つ**ので、
`pullCoeff_eq` をそのまま当てればよい。 -/
theorem divBHom_bPullBase
    (x : Additive (BSubgroup V₁ DK₁ L (normObj_isNormalScheme V₁ L))) :
    divBHom V₂ DK₂ FL (bPullBaseAdd DK₁ DK₂ L FL π hdim hpull x)
      = phiPullBase DK₁ DK₂ L FL π hdim hpull (divBHom V₁ DK₁ L x) := by
  refine Subtype.ext (Finsupp.ext fun s => ?_)
  have hq : (((Additive.toMul x : BSubgroup V₁ DK₁ L _)
      : ((normObj V₁ L).functionField)ˣ) : (normObj V₁ L).functionField) ≠ 0 :=
    (Additive.toMul x : BSubgroup V₁ DK₁ L _).1.ne_zero
  have hDU : ∀ v : PrimeDivisorPt (normObj V₁ L),
      (v : normObj V₁ L) ∈ (⊤ : (normObj V₁ L).Opens) →
      toWeilOnDL V₁ DK₁ L ((divBHom V₁ DK₁ L x : cartierOnDL V₁ DK₁ L _)
          : (DLSet V₁ DK₁ L) →₀ ℤ) v
        = ordPt (normObj V₁ L) (normObj_isNormalScheme V₁ L) v
          (((Additive.toMul x : BSubgroup V₁ DK₁ L _)
            : ((normObj V₁ L).functionField)ˣ) : (normObj V₁ L).functionField) := by
    intro v _
    rw [show ((divBHom V₁ DK₁ L x : cartierOnDL V₁ DK₁ L _) : (DLSet V₁ DK₁ L) →₀ ℤ)
        = Finsupp.subtypeDomain (· ∈ DLSet V₁ DK₁ L)
          (weilDivOfFn (normObj_isNormalScheme V₁ L) hq) from rfl,
      embDomain_subtypeDomain V₁ DK₁ L _ (fun w hw => by
        rw [weilDivOfFn_apply]
        exact (Additive.toMul x).2 w hw)]
    rfl
  show ordPt (normObj V₂ FL) (normObj_isNormalScheme V₂ FL) s.1
      (((Additive.toMul (bPullBaseAdd DK₁ DK₂ L FL π hdim hpull x) : BSubgroup V₂ DK₂ FL _)
        : ((normObj V₂ FL).functionField)ˣ) : (normObj V₂ FL).functionField)
    = pullCoeff π (normObj_isNormalScheme V₂ FL)
      (normObj_isNormalScheme V₁ L) (divBHom V₁ DK₁ L x).2 s.1
  rw [pullCoeff_eq π (normObj_isNormalScheme V₂ FL) (normObj_isNormalScheme V₁ L)
    (divBHom V₁ DK₁ L x).2 s.1 (hdim s.1) hq (Set.mem_univ _) hDU]
  rw [bPullBaseAdd_val]
  rfl

end DivCompat

/-! ## ★3. 自然性 3 本 -/

section Squares

variable {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
  {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
  [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
  [IsGalois V₁.functionField Kbar₁] [IsGalois V₂.functionField Kbar₂]
  [∀ L : FinSub V₁.functionField Kbar₁, IsLocallyNoetherian (normObj V₁ L)]
  [∀ L : FinSub V₁.functionField Kbar₁, CompactSpace (normObj V₁ L)]
  [∀ L : FinSub V₂.functionField Kbar₂, IsLocallyNoetherian (normObj V₂ L)]
  [∀ L : FinSub V₂.functionField Kbar₂, CompactSpace (normObj V₂ L)]
  (DK₁ : Set (PrimeDivisorPt V₁)) (DK₂ : Set (PrimeDivisorPt V₂))

variable {L M : FinSub V₁.functionField Kbar₁} {FL FM : FinSub V₂.functionField Kbar₂}
  (g : L ⟶ M) (Fg : FL ⟶ FM)
  (πL : normObj V₂ FL ⟶ normObj V₁ L) [IsDominant πL]
  (πM : normObj V₂ FM ⟶ normObj V₁ M) [IsDominant πM]
  (hdimL : ∀ w : PrimeDivisorPt (normObj V₂ FL),
    ringKrullDim ((normObj V₁ L).presheaf.stalk (πL.base w.1)) ≤ 1)
  (hpullL : ∀ w : PrimeDivisorPt (normObj V₂ FL), w ∉ DLSet V₂ DK₂ FL →
    ∀ hc : IsCodimOnePt (normObj V₁ L) (πL.base w.1),
      (⟨πL.base w.1, hc⟩ : PrimeDivisorPt (normObj V₁ L)) ∉ DLSet V₁ DK₁ L)
  (hdimM : ∀ w : PrimeDivisorPt (normObj V₂ FM),
    ringKrullDim ((normObj V₁ M).presheaf.stalk (πM.base w.1)) ≤ 1)
  (hpullM : ∀ w : PrimeDivisorPt (normObj V₂ FM), w ∉ DLSet V₂ DK₂ FM →
    ∀ hc : IsCodimOnePt (normObj V₁ M) (πM.base w.1),
      (⟨πM.base w.1, hc⟩ : PrimeDivisorPt (normObj V₁ M)) ∉ DLSet V₁ DK₁ M)

variable (hkq₁ : IsKQCartier V₁ DK₁
    (fun (L : FinSub V₁.functionField Kbar₁) _ => normObj_isNormalScheme V₁ L))
  (hkq₂ : IsKQCartier V₂ DK₂
    (fun (L : FinSub V₂.functionField Kbar₂) _ => normObj_isNormalScheme V₂ L))

omit [IsGalois V₁.functionField Kbar₁] [IsGalois V₂.functionField Kbar₂] in
/-- ★★★★★★**`phiNat`** —— `Φ` の側の自然性。

```
V₂[FM] --π_M--> V₁[M]
   |              |
normMap V₂ Fg   normMap V₁ g
   v              v
V₂[FL] --π_L--> V₁[L]
```

★中身は `pullCoeff_square`(引き戻しの関手性を 2 回)。
★★2 本の `pullCoeff_congr` は「`D_L` 添字へ切り詰めてから戻すと元に戻る」
(`toWeilOnDL_phiPull` / `toWeilOnDL_phiPullBase`)を通すためのもの。 -/
theorem phiPullEff_square
    (hdimLM : ∀ w : PrimeDivisorPt (normObj V₂ FM),
      ringKrullDim ((normObj V₁ L).presheaf.stalk
        ((normMap V₂ Fg ≫ πL).base w.1)) ≤ 1)
    (hsq : normMap V₂ Fg ≫ πL = πM ≫ normMap V₁ g)
    (x : effSub (cartierOnDL V₁ DK₁ L (normObj_isNormalScheme V₁ L))) :
    phiPullEff DK₁ DK₂ M FM πM hdimM hpullM
        ((cartierDatumGeom V₁ DK₁ hkq₁).mapHom g.op x)
      = (cartierDatumGeom V₂ DK₂ hkq₂).mapHom Fg.op
        (phiPullEff DK₁ DK₂ L FL πL hdimL hpullL x) := by
  haveI : IsDominant (πM ≫ normMap V₁ g) := hsq ▸ inferInstance
  refine Subtype.ext (Finsupp.ext fun s => ?_)
  have hA : toWeilOnDL V₁ DK₁ M
        ((phiPull V₁ DK₁ g (effSubIncl _ x) : cartierOnDL V₁ DK₁ M _)
          : (DLSet V₁ DK₁ M) →₀ ℤ)
      = cartierPullback (normMap V₁ g) (normObj_isNormalScheme V₁ M)
        (normObj_isNormalScheme V₁ L) (effSubIncl _ x).2 (normMap_hdim g) :=
    Finsupp.ext fun v => toWeilOnDL_phiPull V₁ DK₁ g (effSubIncl _ x) v
  have hB : toWeilOnDL V₂ DK₂ FL
        ((phiPullBase DK₁ DK₂ L FL πL hdimL hpullL (effSubIncl _ x)
          : cartierOnDL V₂ DK₂ FL _) : (DLSet V₂ DK₂ FL) →₀ ℤ)
      = cartierPullback πL (normObj_isNormalScheme V₂ FL)
        (normObj_isNormalScheme V₁ L) (effSubIncl _ x).2 hdimL :=
    Finsupp.ext fun v =>
      toWeilOnDL_phiPullBase DK₁ DK₂ L FL πL hdimL hpullL (effSubIncl _ x) v
  show pullCoeff πM (normObj_isNormalScheme V₂ FM) (normObj_isNormalScheme V₁ M)
      (phiPull V₁ DK₁ g (effSubIncl _ x)).2 s.1
    = pullCoeff (normMap V₂ Fg) (normObj_isNormalScheme V₂ FM)
        (normObj_isNormalScheme V₂ FL)
      (phiPullBase DK₁ DK₂ L FL πL hdimL hpullL (effSubIncl _ x)).2 s.1
  rw [pullCoeff_congr πM (normObj_isNormalScheme V₂ FM) (normObj_isNormalScheme V₁ M)
      (phiPull V₁ DK₁ g (effSubIncl _ x)).2
      (hA ▸ (phiPull V₁ DK₁ g (effSubIncl _ x)).2) hA s.1,
    pullCoeff_congr (normMap V₂ Fg) (normObj_isNormalScheme V₂ FM)
      (normObj_isNormalScheme V₂ FL)
      (phiPullBase DK₁ DK₂ L FL πL hdimL hpullL (effSubIncl _ x)).2
      (hB ▸ (phiPullBase DK₁ DK₂ L FL πL hdimL hpullL (effSubIncl _ x)).2) hB s.1]
  exact (pullCoeff_square (normMap V₂ Fg) πL πM (normMap V₁ g)
    (normObj_isNormalScheme V₂ FM) (normObj_isNormalScheme V₂ FL)
    (normObj_isNormalScheme V₁ M) (normObj_isNormalScheme V₁ L) (effSubIncl _ x).2
    hdimL (normMap_hdim g) s.1 (normMap_hdim Fg s.1) (hdimM s.1) (hdimLM s.1)
    (hsq ▸ hdimLM s.1) hsq).symm

omit [IsGalois V₁.functionField Kbar₁] [IsGalois V₂.functionField Kbar₂]
  [∀ L : FinSub V₁.functionField Kbar₁, CompactSpace (normObj V₁ L)]
  [∀ L : FinSub V₂.functionField Kbar₂, CompactSpace (normObj V₂ L)] in
/-- ★★★★★**`bmonNat`** —— `B` の側の自然性。★中身は `ffMap_square` 1 本。 -/
theorem bPullBaseAdd_square
    (hsq : normMap V₂ Fg ≫ πL = πM ≫ normMap V₁ g)
    (x : Additive (BSubgroup V₁ DK₁ L (normObj_isNormalScheme V₁ L))) :
    bPullBaseAdd DK₁ DK₂ M FM πM hdimM hpullM (bHom V₁ DK₁ g x)
      = bHom V₂ DK₂ Fg (bPullBaseAdd DK₁ DK₂ L FL πL hdimL hpullL x) := by
  haveI : IsDominant (πM ≫ normMap V₁ g) := hsq ▸ inferInstance
  refine Additive.toMul.injective (Subtype.ext (Units.ext ?_))
  show ffMap πM (ffMap (normMap V₁ g)
      (((Additive.toMul x : BSubgroup V₁ DK₁ L _)
        : ((normObj V₁ L).functionField)ˣ) : (normObj V₁ L).functionField))
    = ffMap (normMap V₂ Fg) (ffMap πL
      (((Additive.toMul x : BSubgroup V₁ DK₁ L _)
        : ((normObj V₁ L).functionField)ˣ) : (normObj V₁ L).functionField))
  exact (ffMap_square (normMap V₂ Fg) πL πM (normMap V₁ g) hsq _).symm

/-- ★★★★★★**`divCompat`** —— `Div_B` との両立(`Φ^gp` の側で書いたもの)。 -/
theorem divBGeom_bPullBaseAdd
    (x : (bmonGeom V₁ DK₁).val (Opposite.op L)) :
    divBGeom V₂ DK₂ hkq₂ (Opposite.op FL)
        (bPullBaseAdd DK₁ DK₂ L FL πL hdimL hpullL x)
      = gpMap _ (phiPullEff DK₁ DK₂ L FL πL hdimL hpullL)
        (divBGeom V₁ DK₁ hkq₁ (Opposite.op L) x) := by
  refine phiGpHomC_injective (cartierDatumGeom V₂ DK₂ hkq₂) finSubOp_isOfFSMType ?_
  rw [phiGpHomC_divBGeom]
  show _ = effSubGpHom (cartierOnDL V₂ DK₂ FL (normObj_isNormalScheme V₂ FL))
    (gpMap _ (restrictEff (phiPullBase DK₁ DK₂ L FL πL hdimL hpullL) _)
      (divBGeom V₁ DK₁ hkq₁ (Opposite.op L) x))
  rw [effSubGpHom_gpMap_restrictEff]
  show _ = phiPullBase DK₁ DK₂ L FL πL hdimL hpullL
    (phiGpHomC (cartierDatumGeom V₁ DK₁ hkq₁) finSubOp_isOfFSMType
      (divBGeom V₁ DK₁ hkq₁ (Opposite.op L) x))
  rw [phiGpHomC_divBGeom]
  exact divBHom_bPullBase DK₁ DK₂ L FL πL hdimL hpullL x

end Squares

/-! ## ★4. 入力データ

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism
-/

/-- ★★★★★★**[FrdI] `Theorem 6.2, (i)` の入力データ**(幾何版)——
原文の「次のデータからなると思ってよい」をそのまま構造体にしたもの。

★★底のスキームの射 `ψ : V₂ ⟶ V₁` と、底の圏の関手 `L ↦ K₂(L)`、
各 `L` に対する体の埋め込み `L ↪ K₂(L)` の族(と、その自然性)である。

★★★`hpull` が**原文の仮定 (a)** —— 「`D_{K₁}` の引き戻しは `D_{K₂}` に入る」。
これの対偶が、因子側と有理函数側の両方を同時に支える。

★逸脱の記録: 原文は「余次元が保たれる」を暗黙に使うので、
`hdimLe`(`π_L` は余次元を上げない)を明示の仮定として置いた。
同じ `V` の中の `normMap` については `ringKrullDim_stalk_normMap_le` で証明済みであり、
底が動く場合も going-down から同様に出るが、`ψ` に何も仮定しない一般形では成り立たない。 -/
structure GeomPullDatum {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
    (Kbar₁ Kbar₂ : Type u) [Field Kbar₁] [Field Kbar₂]
    [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
    [∀ L : FinSub V₁.functionField Kbar₁, IsLocallyNoetherian (normObj V₁ L)]
    [∀ L : FinSub V₂.functionField Kbar₂, IsLocallyNoetherian (normObj V₂ L)]
    [∀ L : FinSub V₂.functionField Kbar₂, CompactSpace (normObj V₂ L)]
    (DK₁ : Set (PrimeDivisorPt V₁)) (DK₂ : Set (PrimeDivisorPt V₂))
    (psi : V₂ ⟶ V₁) [IsDominant psi] where
  /-- 底の圏の関手 `L ↦ K₂(L)`。 -/
  base : FinSub V₁.functionField Kbar₁ ⥤ FinSub V₂.functionField Kbar₂
  /-- 各 `L` に対する体の埋め込み `L ↪ K₂(L)`。 -/
  emb : ∀ L : FinSub V₁.functionField Kbar₁,
    CommRingCat.of L.toIF ⟶ CommRingCat.of (base.obj L).toIF
  /-- 埋め込みは `ψ` の関数体への作用と両立。 -/
  hemb : ∀ L : FinSub V₁.functionField Kbar₁,
    CommRingCat.ofHom (algebraMap V₁.functionField L.toIF) ≫ emb L
      = ffMap psi ≫ CommRingCat.ofHom (algebraMap V₂.functionField (base.obj L).toIF)
  /-- 埋め込みの族は自然。 -/
  hnat : ∀ {L M : FinSub V₁.functionField Kbar₁} (g : L ⟶ M),
    emb L ≫ CommRingCat.ofHom (FinSub.hom (base.map g)).toRingHom
      = CommRingCat.ofHom (FinSub.hom g).toRingHom ≫ emb M
  /-- ★`π_L` は余次元を上げない。 -/
  hdimLe : ∀ (L : FinSub V₁.functionField Kbar₁) (y : normObj V₂ (base.obj L)),
    ringKrullDim ((normObj V₁ L).presheaf.stalk
        ((normMapGeom psi L (base.obj L) (emb L) (hemb L)).base y))
      ≤ ringKrullDim ((normObj V₂ (base.obj L)).presheaf.stalk y)
  /-- ★★原文の仮定 (a) —— `D_{K₁}` の引き戻しは `D_{K₂}` に入る。 -/
  hpull : ∀ (L : FinSub V₁.functionField Kbar₁)
    (w : PrimeDivisorPt (normObj V₂ (base.obj L))), w ∉ DLSet V₂ DK₂ (base.obj L) →
    ∀ hc : IsCodimOnePt (normObj V₁ L)
      ((normMapGeom psi L (base.obj L) (emb L) (hemb L)).base w.1),
      (⟨_, hc⟩ : PrimeDivisorPt (normObj V₁ L)) ∉ DLSet V₁ DK₁ L

namespace GeomPullDatum

variable {V₁ V₂ : Scheme.{u}} [IsIntegral V₁] [IsIntegral V₂]
  {Kbar₁ Kbar₂ : Type u} [Field Kbar₁] [Field Kbar₂]
  [Algebra V₁.functionField Kbar₁] [Algebra V₂.functionField Kbar₂]
  [∀ L : FinSub V₁.functionField Kbar₁, IsLocallyNoetherian (normObj V₁ L)]
  [∀ L : FinSub V₂.functionField Kbar₂, IsLocallyNoetherian (normObj V₂ L)]
  [∀ L : FinSub V₂.functionField Kbar₂, CompactSpace (normObj V₂ L)]
  {DK₁ : Set (PrimeDivisorPt V₁)} {DK₂ : Set (PrimeDivisorPt V₂)}
  {psi : V₂ ⟶ V₁} [IsDominant psi] (P : GeomPullDatum Kbar₁ Kbar₂ DK₁ DK₂ psi)

/-- ★★**`π_L : V₂[K₂(L)] ⟶ V₁[L]`** —— 底変換した正規化の間の射。 -/
noncomputable abbrev pi (L : FinSub V₁.functionField Kbar₁) :
    normObj V₂ (P.base.obj L) ⟶ normObj V₁ L :=
  normMapGeom psi L (P.base.obj L) (P.emb L) (P.hemb L)

/-- ★`π_L` は素因子の上で次元 `≤ 1`。 -/
theorem hdim (L : FinSub V₁.functionField Kbar₁)
    (w : PrimeDivisorPt (normObj V₂ (P.base.obj L))) :
    ringKrullDim ((normObj V₁ L).presheaf.stalk ((P.pi L).base w.1)) ≤ 1 :=
  w.2 ▸ P.hdimLe L w.1

/-- ★★**`π` の族は自然** —— `normMapGeom_naturality` に体の側の立方体を渡すだけ。 -/
theorem hsq {L M : FinSub V₁.functionField Kbar₁} (g : L ⟶ M) :
    normMap V₂ (P.base.map g) ≫ P.pi L = P.pi M ≫ normMap V₁ g :=
  normMapGeom_naturality psi g (P.base.map g) (P.emb L) (P.emb M) (P.hemb L) (P.hemb M)
    (by rw [specMapOf, specMapOf, ← Spec.map_comp, ← Spec.map_comp, P.hnat g])

/-- ★合成した射でも次元 `≤ 1`(`pullCoeff_square` に渡す形)。 -/
theorem hdimLM {L M : FinSub V₁.functionField Kbar₁} (g : L ⟶ M)
    (w : PrimeDivisorPt (normObj V₂ (P.base.obj M))) :
    ringKrullDim ((normObj V₁ L).presheaf.stalk
      ((normMap V₂ (P.base.map g) ≫ P.pi L).base w.1)) ≤ 1 :=
  le_trans (P.hdimLe L _)
    (le_trans (ringKrullDim_stalk_normMap_le V₂ (P.base.map g) w.1) (le_of_eq w.2))

/-! ## ★5. `ModelDataHomOver` と `Ψ : 𝒞₁ ⥤ 𝒞₂` -/

section Hom

variable [IsGalois V₁.functionField Kbar₁] [IsGalois V₂.functionField Kbar₂]
  [∀ L : FinSub V₁.functionField Kbar₁, CompactSpace (normObj V₁ L)]
  (hkq₁ : IsKQCartier V₁ DK₁
    (fun (L : FinSub V₁.functionField Kbar₁) _ => normObj_isNormalScheme V₁ L))
  (hkq₂ : IsKQCartier V₂ DK₂
    (fun (L : FinSub V₂.functionField Kbar₂) _ => normObj_isNormalScheme V₂ L))

/-- ★★★★★★★**[FrdI] `Theorem 6.2, (i)`** ——
幾何のデータの射は `Example 6.1` の model 同士の射を与える。

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism

★★5 フィールドの中身はすべて `Thm62Pull.lean` と本ファイルの ★1〜★3 で閉じている。 -/
noncomputable def toModelDataHomOver :
    ModelDataHomOver P.base.op (ex61ModelData V₁ DK₁ hkq₁)
      (ex61ModelData V₂ DK₂ hkq₂) where
  phiHom A := phiPullEff DK₁ DK₂ A.unop (P.base.obj A.unop) (P.pi A.unop)
    (P.hdim A.unop) (P.hpull A.unop)
  phiNat {A B} f x := phiPullEff_square DK₁ DK₂ f.unop (P.base.map f.unop)
    (P.pi A.unop) (P.pi B.unop) (P.hdim A.unop) (P.hpull A.unop)
    (P.hdim B.unop) (P.hpull B.unop) hkq₁ hkq₂ (P.hdimLM f.unop) (P.hsq f.unop) x
  bmonHom A := bPullBaseAdd DK₁ DK₂ A.unop (P.base.obj A.unop) (P.pi A.unop)
    (P.hdim A.unop) (P.hpull A.unop)
  bmonNat {A B} f x := bPullBaseAdd_square DK₁ DK₂ f.unop (P.base.map f.unop)
    (P.pi A.unop) (P.pi B.unop) (P.hdim A.unop) (P.hpull A.unop)
    (P.hdim B.unop) (P.hpull B.unop) (P.hsq f.unop) x
  divCompat d u := divBGeom_bPullBaseAdd DK₁ DK₂ (P.pi d.unop)
    (P.hdim d.unop) (P.hpull d.unop) hkq₁ hkq₂ u

/-- ★★★★★★★**`Ψ : 𝒞₁ ⥤ 𝒞₂`** —— `Theorem 6.2, (i)` の結論そのもの。

★底の圏の関手 `Ψ_𝒟 = base.op` の上にある。 -/
noncomputable def psiFunctor :
    ModelData.Obj (ex61ModelData V₁ DK₁ hkq₁) ⥤ ModelData.Obj (ex61ModelData V₂ DK₂ hkq₂) :=
  (P.toModelDataHomOver hkq₁ hkq₂).functor

@[simp] theorem psiFunctor_obj_base (A : ModelData.Obj (ex61ModelData V₁ DK₁ hkq₁)) :
    ((P.psiFunctor hkq₁ hkq₂).obj A).base = P.base.op.obj A.base := rfl

@[simp] theorem psiFunctor_map_base
    {A B : ModelData.Obj (ex61ModelData V₁ DK₁ hkq₁)} (φ : A ⟶ B) :
    ((P.psiFunctor hkq₁ hkq₂).map φ).base = P.base.op.map φ.base := rfl

@[simp] theorem psiFunctor_map_deg
    {A B : ModelData.Obj (ex61ModelData V₁ DK₁ hkq₁)} (φ : A ⟶ B) :
    ((P.psiFunctor hkq₁ hkq₂).map φ).deg = φ.deg := rfl

end Hom

end GeomPullDatum

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def phiPullEff_square.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — Φ₁ → Φ₂|𝒟₁ の自然性",
    sectionId := "frdi-thm-6-2" }

def phiPullEff_square.needs : List ProofObligation :=
  [ .citation "[ABC3]" "pullCoeff_square(引き戻しの関手性)"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullCoeff_square") 110,
    .citation "[ABC3]" "toWeilOnDL_phiPull / toWeilOnDL_phiPullBase(切り詰めと復元)"
      (.inProject "ABC3" "ABC3.Found.Divisor.toWeilOnDL_phiPullBase") 110 ]

def bPullBaseAdd_square.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — B₁ → B₂|𝒟₁ の自然性",
    sectionId := "frdi-thm-6-2" }

def bPullBaseAdd_square.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ffMap_square(可換な四角形では関数体の引き戻しが一致)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ffMap_square") 110 ]

def divBGeom_bPullBaseAdd.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — Div_B との両立",
    sectionId := "frdi-thm-6-2" }

def divBGeom_bPullBaseAdd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "divBHom_bPullBase(主因子の引き戻しは引き戻しの主因子)"
      (.inProject "ABC3" "ABC3.Found.Divisor.divBHom_bPullBase") 110,
    .citation "[ABC3]" "effSubGpHom_gpMap_restrictEff(Φ^gp ≃ Γ との両立)"
      (.inProject "ABC3" "ABC3.Found.Divisor.effSubGpHom_gpMap_restrictEff") 110 ]

def GeomPullDatum.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 幾何のデータの射の入力データ",
    sectionId := "frdi-thm-6-2" }

def GeomPullDatum.toModelDataHomOver.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 幾何のデータの射は model 同士の射を与える",
    sectionId := "frdi-thm-6-2" }

def GeomPullDatum.toModelDataHomOver.needs : List ProofObligation :=
  [ .citation "[ABC3]" "phiPullEff / phiPullEff_square(Φ の側)"
      (.inProject "ABC3" "ABC3.Found.Divisor.phiPullEff_square") 110,
    .citation "[ABC3]" "bPullBaseAdd / bPullBaseAdd_square(B の側)"
      (.inProject "ABC3" "ABC3.Found.Divisor.bPullBaseAdd_square") 110,
    .citation "[ABC3]" "divBGeom_bPullBaseAdd(Div_B との両立)"
      (.inProject "ABC3" "ABC3.Found.Divisor.divBGeom_bPullBaseAdd") 110,
    .citation "[ABC3]" "normMapGeom / normMapGeom_naturality(π_L の族)"
      (.inProject "ABC3" "ABC3.Found.Divisor.normMapGeom_naturality") 110 ]

def GeomPullDatum.psiFunctor.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — Ψ : 𝒞₁ ⥤ 𝒞₂",
    sectionId := "frdi-thm-6-2" }

def GeomPullDatum.psiFunctor.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ModelDataHomOver.functor(底の圏が違う版の関手化)"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelDataHomOver.functor") 110 ]

end ABC3.Found.Divisor
