/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop53Rlf

/-!
# [FrdI] Proposition 5.3 の 1-可換図式 —— 下の行の右の矢印 `(𝒞^un-tr)^pf ⟶ 𝒞^rlf`

原文 (FrdI p.103):
> if C is of Frobenius-isotropic type, then there is a natural 1-commutative

★`Proposition 5.3` 末尾の図式

```
𝒞       ⟶ 𝒞^istr       ⟶ 𝒞^pf
↓                          ↓
𝒞^un-tr ⟶ (𝒞^un-tr)^pf ⟶ 𝒞^rlf
```

の**下の行の右の矢印**を通す。★`Prop53Rlf.lean` で
`(𝒞^un-tr)^pf` の行き先も `𝒞^rlf` も**同じ形 `ScModelObj S`** で作ってあり、
違いは係数 `S` が `ℚ≥0` か `ℝ≥0` かだけである。したがって要るのは

  ★★**係数の拡大 `σ : S ⟶ S'` が誘導する model Frobenioid のあいだの関手**

だけであり、それが本ファイルの `scBaseFunctor` である。
`σ := NNRat.castHom NNReal` と取れば `(𝒞^un-tr)^pf ⟶ 𝒞^rlf` になる。

## ★下の行の 1-可換性

`untrToSc ℚ≥0 ⋙ scBaseFunctor = untrToSc ℝ≥0`(`untrToSc_comp_scBaseFunctor`)を
**関手の等式として**示す。★鍵は `scBase σ ∘ toSc = toSc`(`σ 1 = 1` だから)の 1 本で、
そこから `ModelDataHom` の等式(`untrToScHom_comp_scBaseHom`)が出る。

## ★★`ModelData` の射の合成(配管)

`ModelDataHom` には合成が無かったので、本ファイルで足す。
★★`(F ∘ G).functor = F.functor ⋙ G.functor` は **`rfl` ではない** ——
対象の側が `gpMap (g ∘ f) = gpMap g ∘ gpMap f`(`gpMap_comp`、群化の関手性)を
経由するからである。`eqToHom` を通すために
`ModelData.hom_eq_of_eqToHom`(4 成分が `HEq` で一致すれば
`eqToHom` を挟んだ形と等しい)を 1 本置く。

★★★これらは本来 `Thm52Change.lean` に置くべき在庫だが、
そこへ足すと下流の再ビルドが大きいので本ファイルに置く(**逸脱の記録**)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

/-! ## ★1. `ModelData` の射の合成 -/

section ModelDataComp

variable {D : Type u} [Category.{v} D] {D₂ : Type u} [Category.{v} D₂]

/-- ★★**4 成分が `HEq` で一致すれば `eqToHom` を挟んだ形と等しい**。

★`Functor.ext` の `h_map` を埋めるための道具である。
`ModelData.Obj` は `(base, cls)` の対で、`Hom` の 4 成分の**型は `cls` に依らない**ので、
`cls` だけが違う対象のあいだの `eqToHom` は 4 成分の上で恒等になる。 -/
theorem ModelData.hom_eq_of_eqToHom {M : ModelData.{v, u, w} D}
    {A B A' B' : ModelData.Obj M} (hA : A = A') (hB : B = B') (φ : A ⟶ B) (ψ : A' ⟶ B')
    (hbase : HEq φ.base ψ.base) (hdiv : HEq φ.div ψ.div) (hdeg : φ.deg = ψ.deg)
    (hu : HEq φ.u ψ.u) :
    φ = eqToHom hA ≫ ψ ≫ eqToHom hB.symm := by
  subst hA
  subst hB
  simp only [eqToHom_refl, Category.id_comp, Category.comp_id]
  exact ModelData.Hom.ext (eq_of_heq hbase) (eq_of_heq hdiv) hdeg (eq_of_heq hu)

/-- ★★`ModelDataHom` は**データ 2 本(`phiHom` / `bmonHom`)で決まる**。 -/
theorem ModelDataHom.ext' {M M' : ModelData.{v, u, w} D} {F G : ModelDataHom M M'}
    (hphi : ∀ d, F.phiHom d = G.phiHom d) (hbmon : ∀ d, F.bmonHom d = G.bmonHom d) : F = G := by
  cases F
  cases G
  congr 1
  · exact funext hphi
  · exact funext hbmon

/-- ★★★**`ModelData` の射の合成**。 -/
noncomputable def ModelDataHom.comp {M M' M'' : ModelData.{v, u, w} D}
    (F : ModelDataHom M M') (G : ModelDataHom M' M'') : ModelDataHom M M'' where
  phiHom d := (G.phiHom d).comp (F.phiHom d)
  phiNat f x := by
    show G.phiHom _ (F.phiHom _ (M.phi.map f x)) = _
    rw [F.phiNat, G.phiNat]
    rfl
  bmonHom d := (G.bmonHom d).comp (F.bmonHom d)
  bmonNat f x := by
    show G.bmonHom _ (F.bmonHom _ (M.bmon.map f x)) = _
    rw [F.bmonNat, G.bmonNat]
    rfl
  divCompat d u := by
    show M''.divB d (G.bmonHom d (F.bmonHom d u)) = _
    rw [G.divCompat, F.divCompat, ← AddMonoidHom.comp_apply, ← gpMap_comp]

@[simp] theorem ModelDataHom.comp_phiHom {M M' M'' : ModelData.{v, u, w} D}
    (F : ModelDataHom M M') (G : ModelDataHom M' M'') (d : D) :
    (F.comp G).phiHom d = (G.phiHom d).comp (F.phiHom d) := rfl

@[simp] theorem ModelDataHom.comp_bmonHom {M M' M'' : ModelData.{v, u, w} D}
    (F : ModelDataHom M M') (G : ModelDataHom M' M'') (d : D) :
    (F.comp G).bmonHom d = (G.bmonHom d).comp (F.bmonHom d) := rfl

/-- ★★対象の側 —— **`gpMap_comp`(群化の関手性)を経由する**ので `rfl` ではない。 -/
theorem ModelDataHom.comp_obj {M M' M'' : ModelData.{v, u, w} D}
    (F : ModelDataHom M M') (G : ModelDataHom M' M'') (A : ModelData.Obj M) :
    (F.comp G).obj A = G.obj (F.obj A) := by
  show (⟨A.base, gpMap _ ((G.phiHom A.base).comp (F.phiHom A.base)) A.cls⟩ : ModelData.Obj M'')
    = ⟨A.base, gpMap _ (G.phiHom A.base) (gpMap _ (F.phiHom A.base) A.cls)⟩
  rw [gpMap_comp]
  rfl

/-- ★★★**合成の関手は関手の合成**。 -/
theorem ModelDataHom.comp_functor {M M' M'' : ModelData.{v, u, w} D}
    (F : ModelDataHom M M') (G : ModelDataHom M' M'') :
    (F.comp G).functor = F.functor ⋙ G.functor :=
  CategoryTheory.Functor.ext (fun A => F.comp_obj G A) (fun A B _ =>
    ModelData.hom_eq_of_eqToHom (F.comp_obj G A) (F.comp_obj G B) _ _
      HEq.rfl HEq.rfl rfl HEq.rfl)

/-! ### ★底の圏が違う場合(`Corollary 5.4` の形)との合成 -/

theorem ModelData.hom_eq_of_eqToHom₂ {M : ModelData.{v, u, w} D₂}
    {A B A' B' : ModelData.Obj M} (hA : A = A') (hB : B = B') (φ : A ⟶ B) (ψ : A' ⟶ B')
    (hbase : HEq φ.base ψ.base) (hdiv : HEq φ.div ψ.div) (hdeg : φ.deg = ψ.deg)
    (hu : HEq φ.u ψ.u) :
    φ = eqToHom hA ≫ ψ ≫ eqToHom hB.symm := by
  subst hA
  subst hB
  simp only [eqToHom_refl, Category.id_comp, Category.comp_id]
  exact ModelData.Hom.ext (eq_of_heq hbase) (eq_of_heq hdiv) hdeg (eq_of_heq hu)

/-- ★★`ModelDataHomOver` も**データ 2 本で決まる**。 -/
theorem ModelDataHomOver.ext' {ΨB : D ⥤ D₂} {M : ModelData.{v, u, w} D}
    {M₂ : ModelData.{v, u, w} D₂} {F G : ModelDataHomOver ΨB M M₂}
    (hphi : ∀ d, F.phiHom d = G.phiHom d) (hbmon : ∀ d, F.bmonHom d = G.bmonHom d) : F = G := by
  cases F
  cases G
  congr 1
  · exact funext hphi
  · exact funext hbmon

/-- ★★**`M ⟶ M'` の後に `M' ⟶_{Ψ_𝒟} M₂`**。 -/
noncomputable def ModelDataHom.compOver {ΨB : D ⥤ D₂} {M M' : ModelData.{v, u, w} D}
    {M₂ : ModelData.{v, u, w} D₂}
    (F : ModelDataHom M M') (H : ModelDataHomOver ΨB M' M₂) : ModelDataHomOver ΨB M M₂ where
  phiHom d := (H.phiHom d).comp (F.phiHom d)
  phiNat f x := by
    show H.phiHom _ (F.phiHom _ (M.phi.map f x)) = _
    rw [F.phiNat, H.phiNat]
    rfl
  bmonHom d := (H.bmonHom d).comp (F.bmonHom d)
  bmonNat f x := by
    show H.bmonHom _ (F.bmonHom _ (M.bmon.map f x)) = _
    rw [F.bmonNat, H.bmonNat]
    rfl
  divCompat d u := by
    show M₂.divB (ΨB.obj d) (H.bmonHom d (F.bmonHom d u)) = _
    rw [H.divCompat, F.divCompat, ← AddMonoidHom.comp_apply, ← gpMap_comp]

/-- ★★**`M ⟶_{Ψ_𝒟} M₂` の後に `M₂ ⟶ M₂'`**。 -/
noncomputable def ModelDataHomOver.compHom {ΨB : D ⥤ D₂} {M : ModelData.{v, u, w} D}
    {M₂ M₂' : ModelData.{v, u, w} D₂}
    (H : ModelDataHomOver ΨB M M₂) (G : ModelDataHom M₂ M₂') : ModelDataHomOver ΨB M M₂' where
  phiHom d := (G.phiHom (ΨB.obj d)).comp (H.phiHom d)
  phiNat f x := by
    show G.phiHom _ (H.phiHom _ (M.phi.map f x)) = _
    rw [H.phiNat, G.phiNat]
    rfl
  bmonHom d := (G.bmonHom (ΨB.obj d)).comp (H.bmonHom d)
  bmonNat f x := by
    show G.bmonHom _ (H.bmonHom _ (M.bmon.map f x)) = _
    rw [H.bmonNat, G.bmonNat]
    rfl
  divCompat d u := by
    show M₂'.divB (ΨB.obj d) (G.bmonHom (ΨB.obj d) (H.bmonHom d u)) = _
    rw [G.divCompat, H.divCompat, ← AddMonoidHom.comp_apply, ← gpMap_comp]

@[simp] theorem ModelDataHom.compOver_phiHom {ΨB : D ⥤ D₂} {M M' : ModelData.{v, u, w} D}
    {M₂ : ModelData.{v, u, w} D₂}
    (F : ModelDataHom M M') (H : ModelDataHomOver ΨB M' M₂) (d : D) :
    (F.compOver H).phiHom d = (H.phiHom d).comp (F.phiHom d) := rfl

@[simp] theorem ModelDataHomOver.compHom_phiHom {ΨB : D ⥤ D₂} {M : ModelData.{v, u, w} D}
    {M₂ M₂' : ModelData.{v, u, w} D₂}
    (H : ModelDataHomOver ΨB M M₂) (G : ModelDataHom M₂ M₂') (d : D) :
    (H.compHom G).phiHom d = (G.phiHom (ΨB.obj d)).comp (H.phiHom d) := rfl

theorem ModelDataHom.compOver_obj {ΨB : D ⥤ D₂} {M M' : ModelData.{v, u, w} D}
    {M₂ : ModelData.{v, u, w} D₂}
    (F : ModelDataHom M M') (H : ModelDataHomOver ΨB M' M₂) (A : ModelData.Obj M) :
    (F.compOver H).obj A = H.obj (F.obj A) := by
  show (⟨ΨB.obj A.base, gpMap _ ((H.phiHom A.base).comp (F.phiHom A.base)) A.cls⟩ :
      ModelData.Obj M₂)
    = ⟨ΨB.obj A.base, gpMap _ (H.phiHom A.base) (gpMap _ (F.phiHom A.base) A.cls)⟩
  rw [gpMap_comp]
  rfl

theorem ModelDataHomOver.compHom_obj {ΨB : D ⥤ D₂} {M : ModelData.{v, u, w} D}
    {M₂ M₂' : ModelData.{v, u, w} D₂}
    (H : ModelDataHomOver ΨB M M₂) (G : ModelDataHom M₂ M₂') (A : ModelData.Obj M) :
    (H.compHom G).obj A = G.obj (H.obj A) := by
  show (⟨ΨB.obj A.base,
      gpMap _ ((G.phiHom (ΨB.obj A.base)).comp (H.phiHom A.base)) A.cls⟩ : ModelData.Obj M₂')
    = ⟨ΨB.obj A.base, gpMap _ (G.phiHom (ΨB.obj A.base)) (gpMap _ (H.phiHom A.base) A.cls)⟩
  rw [gpMap_comp]
  rfl

theorem ModelDataHom.compOver_functor {ΨB : D ⥤ D₂} {M M' : ModelData.{v, u, w} D}
    {M₂ : ModelData.{v, u, w} D₂}
    (F : ModelDataHom M M') (H : ModelDataHomOver ΨB M' M₂) :
    (F.compOver H).functor = F.functor ⋙ H.functor :=
  CategoryTheory.Functor.ext (fun A => F.compOver_obj H A) (fun A B _ =>
    ModelData.hom_eq_of_eqToHom₂ (F.compOver_obj H A) (F.compOver_obj H B) _ _
      HEq.rfl HEq.rfl rfl HEq.rfl)

theorem ModelDataHomOver.compHom_functor {ΨB : D ⥤ D₂} {M : ModelData.{v, u, w} D}
    {M₂ M₂' : ModelData.{v, u, w} D₂}
    (H : ModelDataHomOver ΨB M M₂) (G : ModelDataHom M₂ M₂') :
    (H.compHom G).functor = H.functor ⋙ G.functor :=
  CategoryTheory.Functor.ext (fun A => H.compHom_obj G A) (fun A B _ =>
    ModelData.hom_eq_of_eqToHom₂ (H.compHom_obj G A) (H.compHom_obj G B) _ _
      HEq.rfl HEq.rfl rfl HEq.rfl)

end ModelDataComp

/-! ## ★2. 係数の拡大 `S ⟶ S'` -/

section ScBase

variable {S S' : Type} [CommSemiring S] [CommSemiring S']

/-- ★★**係数の拡大** `S ⊗_ℕ M → S' ⊗_ℕ M`。 -/
noncomputable def scBase (σ : S →+* S') {M : Type w} [AddCommMonoid M] :
    ScT S M →+ ScT S' M :=
  (TensorProduct.map (σ : S →+ S').toNatLinearMap (LinearMap.id (R := ℕ) (M := M))).toAddMonoidHom

@[simp] theorem scBase_tmul (σ : S →+* S') {M : Type w} [AddCommMonoid M] (r : S) (m : M) :
    scBase σ (r ⊗ₜ m) = (σ r) ⊗ₜ[ℕ] m := rfl

/-- ★★**`σ 1 = 1` なので `1 ⊗ m` は動かない** —— 1-可換性の鍵。 -/
theorem scBase_toSc (σ : S →+* S') {M : Type w} [AddCommMonoid M] (m : M) :
    scBase σ (toSc (S := S) m) = toSc m := by
  show (σ 1) ⊗ₜ[ℕ] m = (1 : S') ⊗ₜ[ℕ] m
  rw [map_one]

/-- ★係数の拡大は `M` の側の射と可換。 -/
theorem scBase_scMap (σ : S →+* S') {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]
    (f : M →+ N) (x : ScT S M) :
    scBase σ (scMap f x) = scMap f (scBase σ x) := by
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul r m => rfl
  | add x y hx hy => simp [map_add, hx, hy]

/-- ★係数の拡大はスカラー倍を `σ` で送る。 -/
theorem scBase_smul (σ : S →+* S') {M : Type w} [AddCommMonoid M] (r : S) (x : ScT S M) :
    scBase σ (r • x) = σ r • scBase σ x := by
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul a m =>
      show (σ (r * a)) ⊗ₜ[ℕ] m = σ r • ((σ a) ⊗ₜ[ℕ] m)
      rw [map_mul, TensorProduct.smul_tmul', smul_eq_mul]
  | add x y hx hy => simp [smul_add, map_add, hx, hy]

/-! ### ★群化の側 -/

theorem gpMap_scBase_toScGp (σ : S →+* S') {M : Type w} [AddCommMonoid M] (x : Gp M) :
    gpMap _ (scBase σ) (toScGp (S := S) x) = toScGp (S := S') x := by
  have hc : (scBase σ (M := M)).comp (toSc : M →+ ScT S M) = (toSc : M →+ ScT S' M) := by
    ext y; exact scBase_toSc σ y
  show gpMap _ _ (gpMap _ _ x) = _
  rw [← AddMonoidHom.comp_apply, ← gpMap_comp, hc]
  rfl

theorem gpMap_scBase_sSmulGp (σ : S →+* S') {M : Type w} [AddCommMonoid M] (r : S)
    (x : Gp (ScT S M)) :
    gpMap _ (scBase σ) (sSmulGp r x) = sSmulGp (σ r) (gpMap _ (scBase σ) x) := by
  have hc : (scBase σ (M := M)).comp (DistribSMul.toAddMonoidHom (ScT S M) r)
      = (DistribSMul.toAddMonoidHom (ScT S' M) (σ r)).comp (scBase σ) := by
    ext y; exact scBase_smul σ r y
  rw [sSmulGp_apply, sSmulGp_apply, ← AddMonoidHom.comp_apply, ← gpMap_comp,
    ← AddMonoidHom.comp_apply, ← gpMap_comp, hc]

theorem gpMap_scBase_scMap (σ : S →+* S') {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]
    (f : M →+ N) (x : Gp (ScT S M)) :
    gpMap _ (scBase σ) (gpMap _ (scMap f) x)
      = gpMap _ (scMap f) (gpMap _ (scBase σ) x) := by
  have hc : (scBase σ (M := N)).comp (scMap (S := S) f)
      = (scMap (S := S') f).comp (scBase σ) := by
    ext y; exact scBase_scMap σ f y
  rw [← AddMonoidHom.comp_apply, ← gpMap_comp, ← AddMonoidHom.comp_apply, ← gpMap_comp, hc]

end ScBase

/-! ## ★3. 係数の拡大が誘導する model Frobenioid の射 -/

section ScBaseHom

variable {S S' : Type} [CommSemiring S] [CommSemiring S']
variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★**`S·Φ^birat` は係数の拡大で `S'·Φ^birat` に入る**。

★生成元 `r · (1 ⊗ y)` が `σ r · (1 ⊗ y)` に移るだけである
(`gpMap_scBase_sSmulGp` ＋ `gpMap_scBase_toScGp`)。 -/
theorem sPhiBiratOn_base (σ : S →+* S') (G : Frobenioid P) {d : D}
    {x : Gp (ScT S (Φ.val d))} (hx : x ∈ sPhiBiratOn S G d) :
    gpMap _ (scBase σ) x ∈ sPhiBiratOn S' G d := by
  have hle : sPhiBiratOn S G d
      ≤ AddSubgroup.comap (gpMap _ (scBase σ)) (sPhiBiratOn S' G d) := by
    refine AddSubgroup.closure_le _ |>.mpr ?_
    rintro _ ⟨r, y, hy, rfl⟩
    show gpMap _ (scBase σ) (sSmulGp r (toScGp y)) ∈ sPhiBiratOn S' G d
    rw [gpMap_scBase_sSmulGp, gpMap_scBase_toScGp]
    exact mem_sPhiBiratOn_gen G (σ r) hy
  exact hle hx

/-- ★★★★★**係数の拡大が誘導する `ModelData` の射**
`(Φ ⊗_S, S·Φ^birat) ⟶ (Φ ⊗_{S'}, S'·Φ^birat)`。 -/
noncomputable def scBaseHom (σ : S →+* S') (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (ScT S' (Φ.val A)))
    (hfsmD : IsOfFSMType D) :
    ModelDataHom (scModel S G hiso hfn hcharInj hint hfsmD)
      (scModel S' G hiso hfn hcharInj' hint' hfsmD) where
  phiHom _ := scBase σ
  phiNat f x := scBase_scMap σ (Φ.map f) x
  bmonHom d :=
    AddMonoidHom.codRestrict ((gpMap _ (scBase σ)).comp (sPhiBiratOn S G d).subtype) _
      (fun x => sPhiBiratOn_base σ G x.2)
  bmonNat := by
    intro A B f x
    exact Subtype.ext (gpMap_scBase_scMap σ (Φ.map f) _)
  divCompat _ _ := rfl

/-- ★★★★★**`Proposition 5.3` の図式の下の行の右の矢印**
`(Φ ⊗_S, S·Φ^birat) ⥤ (Φ ⊗_{S'}, S'·Φ^birat)`。

★`S = ℚ≥0`, `S' = ℝ≥0` なら **`(𝒞^un-tr)^pf ⟶ 𝒞^rlf`** である。 -/
noncomputable def scBaseFunctor (σ : S →+* S') (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (ScT S' (Φ.val A)))
    (hfsmD : IsOfFSMType D) :
    ScModelObj S G hiso hfn hcharInj hint hfsmD
      ⥤ ScModelObj S' G hiso hfn hcharInj' hint' hfsmD :=
  (scBaseHom σ G hiso hfn hcharInj hint hcharInj' hint' hfsmD).functor

end ScBaseHom

/-! ## ★4. 下の行の 1-可換性 -/

section BottomRow

variable {S S' : Type} [CommSemiring S] [CommSemiring S']
variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} [IsConnected D]

/-- ★★★★★**`ModelData` の射の段での 1-可換性** ——
`𝒞^un-tr → (𝒞^un-tr)^pf → 𝒞^rlf` は `𝒞^un-tr → 𝒞^rlf` に等しい。

★中身は `scBase σ (1 ⊗ m) = 1 ⊗ m`(`scBase_toSc`)だけである。 -/
theorem untrToScHom_comp_scBaseHom (σ : S →+* S') (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ.map α)))
    (hintS' : ∀ A : D, IsIntegralMonoid (ScT S' (Φ.val A)))
    (hfsmD : IsOfFSMType D) :
    (untrToScHom S Fc G hint hcharInj hintS hfsmD).comp
        (scBaseHom σ (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
          (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hcharInj' hintS' hfsmD)
      = untrToScHom S' Fc G hint hcharInj' hintS' hfsmD := by
  refine ModelDataHom.ext' (fun d => ?_) (fun d => ?_)
  · exact AddMonoidHom.ext (fun x => scBase_toSc σ x)
  · exact AddMonoidHom.ext (fun x => Subtype.ext (gpMap_scBase_toScGp σ _))

/-- ★★★★★★**`Proposition 5.3` の図式の下の行は 1-可換**(関手の等式として)。

```
𝒞^un-tr ⟶ (𝒞^un-tr)^pf ⟶ 𝒞^rlf   ＝   𝒞^un-tr ⟶ 𝒞^rlf
```
-/
theorem untrToSc_comp_scBaseFunctor (σ : S →+* S') (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ.map α)))
    (hintS' : ∀ A : D, IsIntegralMonoid (ScT S' (Φ.val A)))
    (hfsmD : IsOfFSMType D) :
    untrToSc S Fc G hint hcharInj hintS hfsmD
        ⋙ scBaseFunctor σ (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
          (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hcharInj' hintS' hfsmD
      = untrToSc S' Fc G hint hcharInj' hintS' hfsmD := by
  show (unTr_modelFrobenioid Fc G hint hfsmD).functor
      ⋙ ((untrToScHom S Fc G hint hcharInj hintS hfsmD).functor
        ⋙ (scBaseHom σ (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
          (fun Z => (unTr_isOfModelType Fc G).2 Z)
          hcharInj hintS hcharInj' hintS' hfsmD).functor) = _
  rw [← ModelDataHom.comp_functor, untrToScHom_comp_scBaseHom]
  rfl

end BottomRow

/-! ## ★5. `ℚ≥0 ⟶ ℝ≥0` —— 原文の `(𝒞^un-tr)^pf ⟶ 𝒞^rlf` -/

section PfToRlf

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★★★★**`Proposition 5.3` の図式の `(𝒞^un-tr)^pf ⟶ 𝒞^rlf`**。

★係数を `ℚ≥0` から `ℝ≥0` へ拡げる関手である。 -/
noncomputable def pfToRlfFunctor (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNRat) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (PfT (Φ.val A)))
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A)))
    (hfsmD : IsOfFSMType D) :
    CuntrPf G hiso hfn hcharInj hint hfsmD ⥤ Crlf G hiso hfn hcharInj' hint' hfsmD :=
  scBaseFunctor (NNRat.castHom NNReal) G hiso hfn hcharInj hint hcharInj' hint' hfsmD

end PfToRlf

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.3` の 1-可換図式の**下の行の右の矢印**
`(𝒞^un-tr)^pf ⟶ 𝒞^rlf` と、下の行の 1-可換性。 -/
def scBaseFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 1-可換図式の下の行 𝒞^un-tr ⥤ (𝒞^un-tr)^pf ⥤ 𝒞^rlf",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
