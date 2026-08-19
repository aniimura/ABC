/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def31Pf
import ABC3.Found.FrdI.Prop25iii
import ABC3.Found.FrdI.Prop32

/-!
# [FrdI] Proposition 5.5, (i) の中身 —— `𝒞^pf` の遷移写像は `α ↦ α^m`

原文 (FrdI p.104):
> (i) If A ∈Ob(Cistr) maps to an object Apf ∈Ob(Cpf), then the natural functor

★原文 `(i)` は `𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)` を主張し、証明を
「Frobenius-trivial な `A` については、base-identity な Frobenius 型自己射を考え、
`𝒞` が Frobenius-normalized 型であることを使えば **immediately**」で畳んでいる。

★★**その "immediately" の中身がこれ**である。

## ★なぜこれで済むのか

`𝒞^pf` の射は `Definition 3.1, (iii)` の**余極限**で定義され、
その遷移写像は `frobTransport`(`Proposition 1.10, (i)` の存在と**一意性**)で与えられる。
`A` が Frobenius-trivial なら、添字は base-identity な Frobenius 型自己射 `ζ_m`
(次数 `m`)で走らせられる。★そこで `𝒪^▷(A)` の元 `α` に対し遷移写像を計算すると、
**Frobenius-normalized**(`ζ ≫ α^m = α ≫ ζ`)と一意性から

  `frobTransport ζ ζ α = α ^ m`

★★したがって `𝒪^▷(A^pf)` は「`𝒪^▷(A)` を `α ↦ α^m` で繋いだ余極限」、
すなわち **perfection `𝒪^▷(A)^pf`** になる。

## ★★残り(記録)

余極限そのものの同定(`ℕ≥1` で添字づけた族が `IdxPf A A` の中で **cofinal** であること、
および `Pf (Additive (𝒪^▷(A)))` との突き合わせ)は別の葉である
—— 依存グラフの鎖 `prop55` の `p55i`。
★本ファイルが与えるのは**その計算の核**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★★★**`𝒞^pf` の遷移写像は `𝒪^▷(A)` の上で `α ↦ α^m`**。

★証明は `frobTransport` の**一意性**(`Proposition 1.10, (i)`)1 本:
Frobenius-normalized が `ζ ≫ α^m = α ≫ ζ` を与えるので、`α^m` が遷移写像に一致する。 -/
theorem frobTransport_otri_pow {A : C} (hfn : IsFrobeniusNormalized P A)
    (ζ : End A) (hζ : IsFrobeniusType P ((ζ : A ⟶ A))) (hζb : IsBaseIdentity P ζ)
    (α : End A) (hα : α ∈ OTri P A) :
    frobTransport (F := F) ((ζ : A ⟶ A)) hζ ((ζ : A ⟶ A)) hζ rfl ((α : A ⟶ A))
      = ((α ^ ((P.degFr ((ζ : A ⟶ A)) : ℕ+) : ℕ) : End A) : A ⟶ A) :=
  frobTransport_eq _ hζ _ hζ rfl _ _ (hfn ζ hζb α hα).symm

/-- ★遷移した先も `𝒪^▷(A)` に入る(部分単系だから)。 -/
theorem frobTransport_otri_mem {A : C} (hfn : IsFrobeniusNormalized P A)
    (ζ : End A) (hζ : IsFrobeniusType P ((ζ : A ⟶ A))) (hζb : IsBaseIdentity P ζ)
    (α : End A) (hα : α ∈ OTri P A) :
    (frobTransport (F := F) ((ζ : A ⟶ A)) hζ ((ζ : A ⟶ A)) hζ rfl ((α : A ⟶ A)) : End A)
      ∈ OTri P A := by
  rw [frobTransport_otri_pow hfn ζ hζ hζb α hα]
  exact (OTri P A).pow_mem hα _

/-- ★`𝒪^▷(A)` の可換性を `End A` の言葉で(`Commute` 版)。 -/
theorem otri_comm_end' {A : C} (hfn : IsFrobeniusNormalized P A) {x y : End A}
    (hx : x ∈ OTri P A) (hy : y ∈ OTri P A) : Commute x y :=
  congrArg Subtype.val (otri_mul_comm P hfn ⟨x, hx⟩ ⟨y, hy⟩)

/-- ★★**遷移写像は `𝒪^▷(A)` の単系準同型** —— `α ↦ α^m`。

★`𝒪^▷(A)` は Frobenius-normalized なら**可換**(`Proposition 2.5, (iii)`)なので、
`α ↦ α^m` は単系準同型になる。 -/
noncomputable def otriPowHom {A : C} (hfn : IsFrobeniusNormalized P A) (m : ℕ) :
    OTri P A →* OTri P A where
  toFun α := ⟨((α : End A) ^ m : End A), (OTri P A).pow_mem α.2 m⟩
  map_one' := Subtype.ext (one_pow m)
  map_mul' x y := Subtype.ext (by
    show ((x : End A) * (y : End A)) ^ m = ((x : End A) ^ m) * ((y : End A) ^ m)
    exact (otri_comm_end' hfn x.2 y.2).mul_pow m)

/-! ## ★1b. 同型による共役(一般形)

★`rtObj A 1` は `A` と同型だが等しくない。そこで `End` の共役を単系準同型として
一度だけ作り、必要な性質(base-identity・次数・`𝒪^▷`・Frobenius-normalized)が
すべて移ることを示しておく。 -/

/-- ★同型による共役 `End B →* End A`。 -/
noncomputable def conjEnd {A B : C} (u : A ⟶ B) [IsIso u] : End B →* End A where
  toFun x := u ≫ ((x : B ⟶ B)) ≫ inv u
  map_one' := by
    show u ≫ (𝟙 B) ≫ inv u = 𝟙 A
    rw [Category.id_comp, IsIso.hom_inv_id]
  map_mul' x y := by
    show u ≫ (((y : B ⟶ B)) ≫ ((x : B ⟶ B))) ≫ inv u
      = (u ≫ ((y : B ⟶ B)) ≫ inv u) ≫ (u ≫ ((x : B ⟶ B)) ≫ inv u)
    simp only [Category.assoc]
    rw [← Category.assoc (inv u) u, IsIso.inv_hom_id, Category.id_comp]

@[simp] theorem conjEnd_apply {A B : C} (u : A ⟶ B) [IsIso u] (x : End B) :
    ((conjEnd u x : End A) : A ⟶ A) = u ≫ ((x : B ⟶ B)) ≫ inv u := rfl

theorem conjEnd_baseIdentity {A B : C} (u : A ⟶ B) [IsIso u] {x : End B}
    (hx : IsBaseIdentity P x) : IsBaseIdentity P (conjEnd u x) := by
  haveI : IsIso (P.Base u) := isIso_Base_of_isIso u
  show P.Base (u ≫ ((x : B ⟶ B)) ≫ inv u) = P.Base (𝟙 A)
  rw [P.Base_comp, P.Base_comp, P.Base_id]
  have h1 : P.Base ((x : B ⟶ B)) = P.Base (𝟙 B) := hx
  rw [h1, P.Base_id, Category.id_comp, ← P.Base_comp, IsIso.hom_inv_id, P.Base_id]

theorem conjEnd_degFr {A B : C} (u : A ⟶ B) [IsIso u] (hu1 : P.degFr u = 1) (x : End B) :
    P.degFr (((conjEnd u x : End A)) : A ⟶ A) = P.degFr ((x : B ⟶ B)) := by
  show P.degFr (u ≫ ((x : B ⟶ B)) ≫ inv u) = _
  rw [P.degFr_comp, P.degFr_comp, hu1, degFr_inv_eq_one u hu1]
  simp

theorem conjEnd_otri {A B : C} (u : A ⟶ B) [IsIso u] (hu1 : P.degFr u = 1) {x : End B}
    (hx : x ∈ OTri P B) : (conjEnd u x) ∈ OTri P A :=
  ⟨conjEnd_baseIdentity u hx.1, (conjEnd_degFr u hu1 x).trans hx.2⟩

theorem conjEnd_inv_conjEnd {A B : C} (u : A ⟶ B) [IsIso u] (x : End B) :
    conjEnd (inv u) (conjEnd u x) = x := by
  show inv u ≫ (u ≫ ((x : B ⟶ B)) ≫ inv u) ≫ inv (inv u) = x
  rw [IsIso.inv_inv]
  simp only [Category.assoc]
  rw [← Category.assoc (inv u) u, IsIso.inv_hom_id, Category.id_comp]
  exact Category.comp_id _

/-- ★★**Frobenius-normalized は同型で移る**。

★共役 `conjEnd` が単系準同型であることと、`map_pow` だけで出る。 -/
theorem frobNormalized_of_iso {A B : C} (u : A ⟶ B) [IsIso u] (hu1 : P.degFr u = 1)
    (hfn : IsFrobeniusNormalized P A) : IsFrobeniusNormalized P B := by
  intro φ hφ β hβ
  have h := hfn (conjEnd u φ) (conjEnd_baseIdentity u hφ) (conjEnd u β)
    (conjEnd_otri u hu1 hβ)
  rw [conjEnd_degFr u hu1 φ] at h
  set d := ((P.degFr ((φ : B ⟶ B)) : ℕ+) : ℕ) with hd
  have h' : (conjEnd u β ^ d) * (conjEnd u φ) = (conjEnd u φ) * (conjEnd u β) := h
  have hback := congrArg (conjEnd (inv u)) h'
  rw [map_mul, map_mul, map_pow, conjEnd_inv_conjEnd, conjEnd_inv_conjEnd] at hback
  exact hback

/-! ## ★2. 添字を `ℕ≥1` に落とす —— cofinality の中身 -/

/-- ★★★**Frobenius-trivial 対象から出る Frobenius 型射の終域は、その対象自身と同型**。

★`Definition 1.3, (ii)` の**本質的一意性**(`frobDegUniq`)1 本で出る。 -/
theorem frobType_cod_iso_of_frobTrivial (F : FrobenioidCore P) {A A' : C}
    (hA : IsFrobeniusTrivial P A) (a : A ⟶ A') (ha : IsFrobeniusType P a) :
    ∃ (e : A' ⟶ A), IsIso e ∧ ∃ ζ : End A,
      IsBaseIdentity P ζ ∧ IsFrobeniusType P ((ζ : A ⟶ A)) ∧
      P.degFr ((ζ : A ⟶ A)) = P.degFr a ∧ a ≫ e = ((ζ : A ⟶ A)) := by
  obtain ⟨ζ, hdeg, hprop⟩ := hA
  obtain ⟨e, he, hcomp⟩ :=
    F.frobDegUniq A A' A a ((ζ (P.degFr a) : End A) : A ⟶ A) ha (hprop _).2 (hdeg _).symm
  exact ⟨e, he, ζ (P.degFr a), (hprop _).1, (hprop _).2, hdeg _, hcomp⟩

/-- ★★★★**`𝒞^pf` の添字はすべて「`ζ_m` の対」と同型** ——
`IdxPf A A` の任意の対象 `(A', B', a, b)` に対し、同型 `e : A' ≅ A`, `e' : B' ≅ A` が
あって `a ≫ e = b ≫ e' = ζ_m`(`m = deg_Fr a`)。

★★これが `Proposition 5.5, (i)` の余極限を **`ℕ≥1` へ落とす**ための cofinality の中身。
★`BiFr` の射は「等次数の Frobenius 型射の対」であり、同型は次数 1 の Frobenius 型射
なので、`(e, e')` はちょうど添字圏の射になる。 -/
theorem idxPf_iso_zeta (F : FrobenioidCore P) {A : C} (hA : IsFrobeniusTrivial P A)
    {A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : A ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    ∃ (e : A' ⟶ A) (e' : B' ⟶ A), IsIso e ∧ IsIso e' ∧ ∃ ζ : End A,
      IsBaseIdentity P ζ ∧ IsFrobeniusType P ((ζ : A ⟶ A)) ∧
      P.degFr ((ζ : A ⟶ A)) = P.degFr a ∧
      a ≫ e = ((ζ : A ⟶ A)) ∧ b ≫ e' = ((ζ : A ⟶ A)) := by
  obtain ⟨ζ, hdeg, hprop⟩ := hA
  obtain ⟨e, he, hcomp⟩ :=
    F.frobDegUniq A A' A a ((ζ (P.degFr a) : End A) : A ⟶ A) ha (hprop _).2 (hdeg _).symm
  obtain ⟨e', he', hcomp'⟩ :=
    F.frobDegUniq A B' A b ((ζ (P.degFr a) : End A) : A ⟶ A) hb (hprop _).2
      (by rw [hdeg, hd])
  exact ⟨e, e', he, he', ζ (P.degFr a), (hprop _).1, (hprop _).2, hdeg _, hcomp, hcomp'⟩

/-! ## ★3. 写像の側 —— `𝒪^▷(A) → 𝒪^▷(A^pf)` -/

/-- ★★★**自然な関手 `𝒞 → 𝒞^pf` は `𝒪^▷(A) → 𝒪^▷(A^pf)` を誘導する**。

★★これが `Proposition 5.5, (i)` の「the natural functor `𝒞 → 𝒞^pf` determines
a natural isomorphism `𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)`」の**写像の側**である。
★同型であることは余極限の同定を要する(鎖 `prop55` の `p55i-colim`)が、
その材料(遷移写像 `α ↦ α^m` と添字の cofinality)は本ファイルで揃っている。 -/
noncomputable def otriToPfRoot {A : C} :
    OTri P A →* OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F) where
  toFun α := ⟨(Functor.mapEnd A (toPfRoot (P := P) (F := F))) (α : End A), by
    refine ⟨?_, ?_⟩
    · show rootBase (toRootHom (F := F) (((α : End A) : A ⟶ A)))
        = rootBase (idRoot P F ⟨A, 1⟩)
      rw [rootBase_toRootHom, rootBase_id]
      have : P.Base (((α : End A) : A ⟶ A)) = P.Base (𝟙 A) := α.2.1
      rw [this, P.Base_id]
    · show rootDeg (toRootHom (F := F) (((α : End A) : A ⟶ A))) = 1
      rw [rootDeg_toRootHom]
      exact α.2.2⟩
  map_one' := Subtype.ext (map_one _)
  map_mul' x y := Subtype.ext (map_mul _ _ _)

/-! ## ★4. 同型による共役

★`HomRoot ⟨A,1⟩ ⟨A,1⟩` の添字は `rtObj A 1`(`A` と**同型**だが等しくはない)の上にある。
そこで `𝒪^▷` を同型 `rtExt A 1 : A ≅ rtObj A 1` で移す必要がある。 -/

/-- ★★同型による共役は `𝒪^▷` を移す。 -/
theorem otri_conj_iso {A B : C} (u : A ⟶ B) [IsIso u] (hu1 : P.degFr u = 1)
    {α : End A} (hα : α ∈ OTri P A) :
    ((inv u ≫ ((α : End A) : A ⟶ A) ≫ u : B ⟶ B) : End B) ∈ OTri P B := by
  have hiu1 : P.degFr (inv u) = 1 := degFr_inv_eq_one u hu1
  haveI hbu : IsIso (P.Base u) := isIso_Base_of_isIso u
  have hbinv : P.Base (inv u) = inv (P.Base u) := by
    refine IsIso.eq_inv_of_hom_inv_id ?_
    rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
  refine ⟨?_, ?_⟩
  · show P.Base (inv u ≫ ((α : End A) : A ⟶ A) ≫ u) = P.Base (𝟙 B)
    rw [P.Base_comp, P.Base_comp, P.Base_id, hbinv]
    have : P.Base (((α : End A) : A ⟶ A)) = P.Base (𝟙 A) := hα.1
    rw [this, P.Base_id, Category.id_comp, IsIso.inv_hom_id]
  · show P.degFr (inv u ≫ ((α : End A) : A ⟶ A) ≫ u) = 1
    rw [P.degFr_comp, P.degFr_comp, hiu1, hu1, hα.2]
    simp

/-! ## ★5. 逆向きの写像の材料 —— `α^{1/m}` を `𝒞^pf` の中に作る

★`HomRoot ⟨A,1⟩ ⟨A,1⟩` の添字は `IdxPf (rtObj A 1) (rtObj A 1)` にある。
`ζ`(次数 `m` の Frobenius 型自己射)を共役で `rtObj A 1` へ移して添字を作り、
そこに `α` を載せると、「`α` の `m` 乗根」にあたる `𝒞^pf` の自己射が得られる。 -/

/-- ★★同型による共役は Frobenius 型を保つ(isotropic 型のとき)。 -/
theorem frobType_conj_iso (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (u : A ⟶ B) [IsIso u] (hu0 : P.Div u = 0)
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) :
    IsFrobeniusType P (inv u ≫ ζ ≫ u) := by
  have hiu0 : P.Div (inv u) = 0 := Div_inv_eq_zero u hu0
  haveI hbu : IsIso (P.Base u) := isIso_Base_of_isIso u
  haveI hbiu : IsIso (P.Base (inv u)) := isIso_Base_of_isIso (inv u)
  haveI hbz : IsIso (P.Base ζ) := hζ.2
  refine ⟨⟨prop_1_4_i P _ (fun X _ => hiso X), ?_⟩, ?_⟩
  · have h1 : P.Div (ζ ≫ u) = 0 := by
      rw [P.Div_comp, hu0, hζ.1.2, map_zero, smul_zero, add_zero]
    show P.Div (inv u ≫ ζ ≫ u) = 0
    rw [P.Div_comp, h1, hiu0, map_zero, smul_zero, add_zero]
  · show IsIso (P.Base (inv u ≫ ζ ≫ u))
    rw [P.Base_comp, P.Base_comp]
    infer_instance

/-- ★`A` の自己射を `rtObj A 1` へ共役で移したもの。 -/
noncomputable def conjRt {A : C} (f : A ⟶ A) : rtObj P F A 1 ⟶ rtObj P F A 1 :=
  haveI := isIso_rtExt_one P F A
  inv (rtExt P F A 1) ≫ f ≫ rtExt P F A 1

theorem conjRt_frobType (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) :
    IsFrobeniusType P (conjRt (F := F) ζ) := by
  haveI := isIso_rtExt_one P F A
  exact frobType_conj_iso hiso _ ((rtExt_frobType P F A 1).1.2) hζ

theorem conjRt_base_id {A : C} {α : End A} (hα : α ∈ OTri P A) :
    P.Base (conjRt (F := F) ((α : A ⟶ A))) = 𝟙 _ := by
  haveI := isIso_rtExt_one P F A
  haveI : IsIso (P.Base (rtExt P F A 1)) := isIso_Base_of_isIso (rtExt P F A 1)
  show P.Base (inv (rtExt P F A 1) ≫ ((α : A ⟶ A)) ≫ rtExt P F A 1) = 𝟙 _
  rw [P.Base_comp, P.Base_comp]
  have h1 : P.Base ((α : A ⟶ A)) = P.Base (𝟙 A) := hα.1
  rw [h1, P.Base_id, Category.id_comp, ← P.Base_comp, IsIso.inv_hom_id, P.Base_id]

/-- ★`ζ` から作る `𝒞^pf` の添字。 -/
noncomputable def idxZeta (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) :
    IdxPf P F (rtObj P F A 1) (rtObj P F A 1) :=
  idxMk (conjRt (F := F) ζ) (conjRt (F := F) ζ)
    (conjRt_frobType hiso hζ) (conjRt_frobType hiso hζ) rfl

/-- ★★**「`α` の `m` 乗根」にあたる `𝒞^pf` の自己射**。 -/
noncomputable def otriPfMk (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) (α : End A) :
    HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨A, 1⟩ :=
  HomPf.mk (idxZeta hiso hζ) (conjRt (F := F) ((α : A ⟶ A)))

theorem otriPfMk_deg (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) (α : End A) (hα : α ∈ OTri P A) :
    rootDeg (otriPfMk (F := F) hiso hζ α) = 1 := by
  show pfDeg (HomPf.mk (idxZeta hiso hζ) (conjRt (F := F) ((α : A ⟶ A)))) = 1
  rw [pfDeg_mk]
  show P.degFr (conjRt (F := F) ((α : A ⟶ A))) = 1
  haveI := isIso_rtExt_one P F A
  show P.degFr (inv (rtExt P F A 1) ≫ ((α : A ⟶ A)) ≫ rtExt P F A 1) = 1
  rw [P.degFr_comp, P.degFr_comp, rtExt_degFr,
    degFr_inv_eq_one (rtExt P F A 1) (rtExt_degFr P F A 1), hα.2]
  simp

theorem otriPfMk_base (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) (α : End A) (hα : α ∈ OTri P A) :
    rootBase (otriPfMk (F := F) hiso hζ α) = 𝟙 _ := by
  haveI := isIso_rtExt_one P F A
  haveI hbr : IsIso (P.Base (rtExt P F A 1)) := isIso_Base_of_isIso (rtExt P F A 1)
  haveI hbz : IsIso (P.Base (conjRt (F := F) ζ)) := (conjRt_frobType hiso hζ).2
  have hrep : pfBase (otriPfMk (F := F) hiso hζ α) = 𝟙 _ := by
    show pfBase (HomPf.mk (idxZeta hiso hζ) (conjRt (F := F) ((α : A ⟶ A)))) = _
    rw [pfBase_mk]
    show P.Base (conjRt (F := F) ζ) ≫ P.Base (conjRt (F := F) ((α : A ⟶ A)))
        ≫ inv (P.Base (conjRt (F := F) ζ)) = 𝟙 _
    rw [conjRt_base_id hα, Category.id_comp, IsIso.hom_inv_id]
  show P.Base (rtExt P F A 1) ≫ pfBase (otriPfMk (F := F) hiso hζ α)
      ≫ inv (P.Base (rtExt P F A 1)) = 𝟙 _
  rw [hrep, Category.id_comp, IsIso.hom_inv_id]

/-- ★★★★**`α^{1/m}` は `𝒪^▷(A^pf)` の元である**。

★★これが `Proposition 5.5, (i)` の同型の**逆向きの写像の材料**である ——
`𝒪^▷(A)` の元 `α` と Frobenius 型自己射 `ζ`(次数 `m`)の対から、
`𝒞^pf` の base-identity かつ linear な自己射が作れる。 -/
theorem otriPfMk_mem (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) (α : End A) (hα : α ∈ OTri P A) :
    ((otriPfMk (F := F) hiso hζ α : End (⟨A, 1⟩ : PfRootObj P F)))
      ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F) := by
  refine ⟨?_, otriPfMk_deg hiso hζ α hα⟩
  show rootBase (otriPfMk (F := F) hiso hζ α) = rootBase (idRoot P F ⟨A, 1⟩)
  rw [otriPfMk_base hiso hζ α hα, rootBase_id]

/-! ## ★6. `Pf` の関係が `𝒞^pf` の中で成り立つこと

★`Pf (Additive (𝒪^▷(A)))` の元は「`α` の `m` 乗根」であり、
その同値関係は「`α^{c}` の `m·c` 乗根と同じ」である。
★★それが `𝒞^pf` の余極限の中で実際に成り立つ、というのが本節である。 -/

theorem conjRt_comp {A : C} (f g : A ⟶ A) :
    conjRt (F := F) (f ≫ g) = conjRt (F := F) f ≫ conjRt (F := F) g := by
  haveI := isIso_rtExt_one P F A
  show inv (rtExt P F A 1) ≫ (f ≫ g) ≫ rtExt P F A 1
    = (inv (rtExt P F A 1) ≫ f ≫ rtExt P F A 1) ≫ (inv (rtExt P F A 1) ≫ g ≫ rtExt P F A 1)
  simp only [Category.assoc]
  rw [← Category.assoc (rtExt P F A 1) (inv (rtExt P F A 1)), IsIso.hom_inv_id,
    Category.id_comp]

/-- ★★★**`rtObj A 1` の上でも遷移写像は `α ↦ α^m`** ——
共役で `A` の計算(`frobTransport_otri_pow`)に帰着する。

★`frobTransport` の**一意性**に、`hfn` を共役で移した四角形を渡すだけ。 -/
theorem frobTransport_conjRt_pow (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hfn : IsFrobeniusNormalized P A) {ζ : End A}
    (hζ : IsFrobeniusType P ((ζ : A ⟶ A))) (hζb : IsBaseIdentity P ζ)
    (α : End A) (hα : α ∈ OTri P A) :
    frobTransport (F := F) (conjRt (F := F) ((ζ : A ⟶ A))) (conjRt_frobType hiso hζ)
        (conjRt (F := F) ((ζ : A ⟶ A))) (conjRt_frobType hiso hζ) rfl
        (conjRt (F := F) ((α : A ⟶ A)))
      = conjRt (F := F) (((α ^ ((P.degFr ((ζ : A ⟶ A)) : ℕ+) : ℕ) : End A) : A ⟶ A)) := by
  refine frobTransport_eq _ _ _ _ rfl _ _ ?_
  have h := congrArg (conjRt (F := F)) (hfn ζ hζb α hα)
  rw [conjRt_comp, conjRt_comp] at h
  exact h.symm

/-- ★添字 `Z_ζ ⟶ Z_{ζ≫ξ}`(`ξ` を右から足す)。

★`BiFr` の射は「等次数の Frobenius 型射の対」なので `(conjRt ξ, conjRt ξ)` が使え、
`Under` の条件は `conjRt` の関手性(`conjRt_comp`)そのもの。 -/
noncomputable def idxZetaStep (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ ξ : A ⟶ A} (hζ : IsFrobeniusType P ζ) (hξ : IsFrobeniusType P ξ)
    (hcomp : IsFrobeniusType P (ζ ≫ ξ)) :
    (idxZeta (F := F) hiso hζ) ⟶ (idxZeta (F := F) hiso hcomp) :=
  Under.homMk (⟨(conjRt (F := F) ξ, conjRt (F := F) ξ),
      conjRt_frobType hiso hξ, conjRt_frobType hiso hξ, rfl⟩ :
      (⟨(rtObj P F A 1, rtObj P F A 1)⟩ : BiFr P F) ⟶ ⟨(rtObj P F A 1, rtObj P F A 1)⟩)
    (by
      refine WideSubcategory.hom_ext _ ?_
      exact Prod.ext (conjRt_comp ζ ξ).symm (conjRt_comp ζ ξ).symm)

/-- ★★★★**`Pf` の関係が `𝒞^pf` の中で成り立つ** ——
「`α` の `deg ζ` 乗根」と「`α^{deg ξ}` の `deg(ζ≫ξ)` 乗根」は**同じ元**である。

★★これが `Proposition 5.5, (i)` の同型の **well-defined 性**の中身である。 -/
theorem otriPfMk_step (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hfn : IsFrobeniusNormalized P A) {ζ ξ : A ⟶ A}
    (hζ : IsFrobeniusType P ζ) (hξ : IsFrobeniusType P ξ) (hξb : IsBaseIdentity P ξ)
    (hcomp : IsFrobeniusType P (ζ ≫ ξ))
    (α : End A) (hα : α ∈ OTri P A) :
    otriPfMk (F := F) hiso hcomp ((α ^ ((P.degFr ξ : ℕ+) : ℕ) : End A))
      = otriPfMk (F := F) hiso hζ α := by
  have h := HomPf.mk_map (idxZetaStep (F := F) hiso hζ hξ hcomp)
    (conjRt (F := F) ((α : A ⟶ A)))
  refine Eq.trans ?_ h
  show HomPf.mk (idxZeta hiso hcomp) (conjRt (F := F) ((α ^ ((P.degFr ξ : ℕ+) : ℕ) : End A)))
    = HomPf.mk (idxZeta hiso hcomp)
      (idxTransport P F (idxZetaStep (F := F) hiso hζ hξ hcomp)
        (conjRt (F := F) ((α : A ⟶ A))))
  congr 1
  show conjRt (F := F) ((α ^ ((P.degFr ξ : ℕ+) : ℕ) : End A))
    = frobTransport (F := F) (conjRt (F := F) ξ) _ (conjRt (F := F) ξ) _ _
        (conjRt (F := F) ((α : A ⟶ A)))
  exact (frobTransport_conjRt_pow hiso hfn hξ hξb α hα).symm

/-! ## ★7. 全射性の核 —— どの自己射も「`ζ` の添字」で代表される -/

/-- ★同型は Frobenius 型(isotropic 型のとき)。 -/
theorem frobType_of_isIso (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (u : A ⟶ B) [IsIso u] :
    IsFrobeniusType P u :=
  (prop_1_4_i_frobeniusType P u (fun X _ => hiso X)).mpr
    ⟨isIsometric_of_isIso P u, isIso_Base_of_isIso u⟩

/-- ★base-identity は共役で保たれる。 -/
theorem conjRt_baseIdentity {A : C} {f : End A} (hf : IsBaseIdentity P f) :
    IsBaseIdentity P (conjRt (F := F) ((f : A ⟶ A))) := by
  haveI := isIso_rtExt_one P F A
  haveI : IsIso (P.Base (rtExt P F A 1)) := isIso_Base_of_isIso (rtExt P F A 1)
  show P.Base (inv (rtExt P F A 1) ≫ ((f : A ⟶ A)) ≫ rtExt P F A 1) = P.Base (𝟙 _)
  rw [P.Base_comp, P.Base_comp, P.Base_id]
  have h1 : P.Base ((f : A ⟶ A)) = P.Base (𝟙 A) := hf
  rw [h1, P.Base_id, Category.id_comp, ← P.Base_comp, IsIso.inv_hom_id, P.Base_id]

/-- ★Frobenius 次数は共役で保たれる。 -/
theorem conjRt_degFr {A : C} (f : A ⟶ A) :
    P.degFr (conjRt (F := F) f) = P.degFr f := by
  haveI := isIso_rtExt_one P F A
  show P.degFr (inv (rtExt P F A 1) ≫ f ≫ rtExt P F A 1) = P.degFr f
  rw [P.degFr_comp, P.degFr_comp, rtExt_degFr,
    degFr_inv_eq_one (rtExt P F A 1) (rtExt_degFr P F A 1)]
  simp

/-- ★★`rtObj A 1` も Frobenius-trivial(`A` と同型だから)。 -/
theorem frobTrivial_rtObj (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A) : IsFrobeniusTrivial P (rtObj P F A 1) := by
  obtain ⟨ζ, hdeg, hprop⟩ := hA
  haveI := isIso_rtExt_one P F A
  refine ⟨{ toFun := fun n => conjRt (F := F) ((ζ n : End A) : A ⟶ A)
            map_one' := by
              show conjRt (F := F) ((ζ 1 : End A) : A ⟶ A) = 𝟙 _
              rw [map_one ζ]
              show inv (rtExt P F A 1) ≫ (𝟙 A) ≫ rtExt P F A 1 = 𝟙 _
              rw [Category.id_comp, IsIso.inv_hom_id]
            map_mul' := fun m n => by
              show conjRt (F := F) ((ζ (m * n) : End A) : A ⟶ A)
                = conjRt (F := F) ((ζ n : End A) : A ⟶ A)
                  ≫ conjRt (F := F) ((ζ m : End A) : A ⟶ A)
              rw [← conjRt_comp]
              exact congrArg (conjRt (F := F)) (map_mul ζ m n) },
    fun n => (conjRt_degFr _).trans (hdeg n),
    fun n => ⟨conjRt_baseIdentity (hprop n).1, conjRt_frobType hiso (hprop n).2⟩⟩

/-- ★`rtObj A 1` の Frobenius 型自己射から作る添字(`idxZeta` の一般形)。 -/
noncomputable def idxZeta' {A : C} {z : rtObj P F A 1 ⟶ rtObj P F A 1}
    (hz : IsFrobeniusType P z) : IdxPf P F (rtObj P F A 1) (rtObj P F A 1) :=
  idxMk z z hz hz rfl

/-- ★★★★**全射性の核** —— `𝒞^pf` の自己射はすべて「`ζ` の添字」で代表される。

★`HomPf.exists_rep` で任意の代表を取り、`idxPf_iso_zeta`(cofinality)で
その添字を `ζ` の添字へ移す。同型は次数 1 の Frobenius 型射なので、
`(e, e')` はちょうど添字圏の射になる。 -/
theorem homPf_exists_zeta_rep (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A) (f : HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨A, 1⟩) :
    ∃ (z : rtObj P F A 1 ⟶ rtObj P F A 1) (hz : IsFrobeniusType P z)
      (φ : rtObj P F A 1 ⟶ rtObj P F A 1), HomPf.mk (idxZeta' hz) φ = f := by
  obtain ⟨Z, φ₀, hZ⟩ := HomPf.exists_rep f
  obtain ⟨e, e', he, he', z, hzb, hzf, hzdeg, hae, hbe⟩ :=
    idxPf_iso_zeta F (frobTrivial_rtObj hiso hA) Z.hom.hom.1 Z.hom.property.1
      Z.hom.hom.2 Z.hom.property.2.1 Z.hom.property.2.2
  haveI := he
  haveI := he'
  have hu : Z ⟶ idxZeta' (F := F) hzf :=
    Under.homMk (⟨(e, e'), frobType_of_isIso hiso e, frobType_of_isIso hiso e',
        (degFr_of_isIso P e).trans (degFr_of_isIso P e').symm⟩ :
        Z.right ⟶ (⟨(rtObj P F A 1, rtObj P F A 1)⟩ : BiFr P F))
      (by
        refine WideSubcategory.hom_ext _ ?_
        exact Prod.ext hae hbe)
  refine ⟨z, hzf, idxTransport P F hu φ₀, ?_⟩
  rw [HomPf.mk_map hu φ₀]
  exact hZ

/-! ## ★8. 単射性の核 -/

/-- ★★`ζ` 添字どうしの射に沿う遷移は `φ ↦ φ^{deg}`。

★2 つの脚が一致することは `z₁` の epi 性(`𝒞` は totally epimorphic)から出る。
★そのあとは `frobTransport` の一意性に Frobenius-normalized の等式を渡すだけ。 -/
theorem idxTransport_zeta {A : C} (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    {z₁ z : rtObj P F A 1 ⟶ rtObj P F A 1}
    (hz₁ : IsFrobeniusType P z₁) (hz₁b : IsBaseIdentity P z₁)
    (hz : IsFrobeniusType P z) (hzb : IsBaseIdentity P z)
    (t : (idxZeta' (F := F) hz₁) ⟶ (idxZeta' (F := F) hz))
    {φ : End (rtObj P F A 1)} (hφ : φ ∈ OTri P (rtObj P F A 1)) :
    idxTransport P F t ((φ : rtObj P F A 1 ⟶ rtObj P F A 1))
      = ((φ ^ ((P.degFr (t.right.hom.1) : ℕ+) : ℕ) : End (rtObj P F A 1))
          : rtObj P F A 1 ⟶ rtObj P F A 1) := by
  have hw := Under.w t
  have hw1 : z₁ ≫ t.right.hom.1 = z :=
    congrArg (fun s : (idxZeta' (F := F) hz₁).right ⟶ _ => s.hom.1) hw
  have hw2 : z₁ ≫ t.right.hom.2 = z :=
    congrArg (fun s : (idxZeta' (F := F) hz₁).right ⟶ _ => s.hom.2) hw
  haveI : Epi z₁ := P.totEpiC _ _ _
  have heq12 : t.right.hom.1 = t.right.hom.2 := (cancel_epi z₁).mp (hw1.trans hw2.symm)
  have hkb : IsBaseIdentity P (t.right.hom.1) := by
    have h0 : P.Base z₁ ≫ P.Base (t.right.hom.1) = P.Base z := by
      rw [← P.Base_comp]; exact congrArg P.Base hw1
    have hb1 : P.Base z₁ = P.Base (𝟙 (rtObj P F A 1)) := hz₁b
    rw [hb1, P.Base_id, Category.id_comp] at h0
    exact h0.trans hzb
  refine frobTransport_eq _ _ _ _ _ _ _ ?_
  rw [← heq12]
  exact (hfn' (t.right.hom.1) hkb φ hφ).symm

/-- ★★★★**単射性の核** —— `ζ` 添字での 2 つの代表が `𝒞^pf` で等しければ、
`Pf` の同値関係(`k₁ · deg z₁ = k₂ · deg z₂` かつ `φ₁^{k₁} = φ₂^{k₂}`)が成り立つ。

★`HomPf.eq_iff` で共通の上界 `V` を取り、`idxPf_iso_zeta` でそれを `ζ` 添字へ移す。
そこへの 2 本の射に `idxTransport_zeta` を当てると、両辺が冪になる。 -/
theorem homPf_zeta_eq (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    {z₁ z₂ : rtObj P F A 1 ⟶ rtObj P F A 1}
    (hz₁ : IsFrobeniusType P z₁) (hz₁b : IsBaseIdentity P z₁)
    (hz₂ : IsFrobeniusType P z₂) (hz₂b : IsBaseIdentity P z₂)
    {φ₁ φ₂ : End (rtObj P F A 1)}
    (hφ₁ : φ₁ ∈ OTri P (rtObj P F A 1)) (hφ₂ : φ₂ ∈ OTri P (rtObj P F A 1))
    (heq : HomPf.mk (idxZeta' (F := F) hz₁) ((φ₁ : rtObj P F A 1 ⟶ rtObj P F A 1))
      = HomPf.mk (idxZeta' (F := F) hz₂) ((φ₂ : rtObj P F A 1 ⟶ rtObj P F A 1))) :
    ∃ k₁ k₂ : ℕ+, k₁ * P.degFr z₁ = k₂ * P.degFr z₂ ∧
      ((φ₁ ^ ((k₁ : ℕ+) : ℕ) : End (rtObj P F A 1)))
        = ((φ₂ ^ ((k₂ : ℕ+) : ℕ) : End (rtObj P F A 1))) := by
  obtain ⟨V, uu, vv, hUV⟩ := HomPf.eq_iff.mp heq
  obtain ⟨e, e', he, he', z, hzb, hzf, hzdeg, hae, hbe⟩ :=
    idxPf_iso_zeta F (frobTrivial_rtObj hiso hA) V.hom.hom.1 V.hom.property.1
      V.hom.hom.2 V.hom.property.2.1 V.hom.property.2.2
  haveI := he
  haveI := he'
  have ww : V ⟶ idxZeta' (F := F) hzf :=
    Under.homMk (⟨(e, e'), frobType_of_isIso hiso e, frobType_of_isIso hiso e',
        (degFr_of_isIso P e).trans (degFr_of_isIso P e').symm⟩ :
        V.right ⟶ (⟨(rtObj P F A 1, rtObj P F A 1)⟩ : BiFr P F))
      (by
        refine WideSubcategory.hom_ext _ ?_
        exact Prod.ext hae hbe)
  have h1 := idxTransport_zeta (F := F) hfn' hz₁ hz₁b hzf hzb (uu ≫ ww) hφ₁
  have h2 := idxTransport_zeta (F := F) hfn' hz₂ hz₂b hzf hzb (vv ≫ ww) hφ₂
  have hmid : idxTransport P F (uu ≫ ww) ((φ₁ : rtObj P F A 1 ⟶ rtObj P F A 1))
      = idxTransport P F (vv ≫ ww) ((φ₂ : rtObj P F A 1 ⟶ rtObj P F A 1)) := by
    rw [idxTransport_comp, idxTransport_comp, hUV]
  refine ⟨P.degFr ((uu ≫ ww).right.hom.1), P.degFr ((vv ≫ ww).right.hom.1), ?_, ?_⟩
  · have e1 : P.degFr (z₁ ≫ (uu ≫ ww).right.hom.1) = P.degFr z :=
      congrArg P.degFr (congrArg (fun s : (idxZeta' (F := F) hz₁).right ⟶ _ => s.hom.1)
        (Under.w (uu ≫ ww)))
    have e2 : P.degFr (z₂ ≫ (vv ≫ ww).right.hom.1) = P.degFr z :=
      congrArg P.degFr (congrArg (fun s : (idxZeta' (F := F) hz₂).right ⟶ _ => s.hom.1)
        (Under.w (vv ≫ ww)))
    rw [P.degFr_comp] at e1 e2
    exact e1.trans e2.symm
  · exact h1.symm.trans (hmid.trans h2)

/-! ## ★9. `𝒪^▷` の層での全射性

★★`Proposition 5.5, (i)` の同型の**全射性**を、`𝒪^▷` の元として述べる:
`𝒪^▷(A^pf)` の元はすべて「`𝒪^▷` の元の `m` 乗根」である。 -/

/-- ★`ζ` 添字での代表(型を `HomRoot` に固定したもの)。 -/
noncomputable def zetaMk {A : C} {z : rtObj P F A 1 ⟶ rtObj P F A 1}
    (hz : IsFrobeniusType P z) (φ : rtObj P F A 1 ⟶ rtObj P F A 1) :
    HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨A, 1⟩ :=
  HomPf.mk (idxZeta' (F := F) hz) φ

theorem rootDeg_zetaMk {A : C} {z : rtObj P F A 1 ⟶ rtObj P F A 1}
    (hz : IsFrobeniusType P z) (φ : rtObj P F A 1 ⟶ rtObj P F A 1) :
    rootDeg (zetaMk (F := F) hz φ) = P.degFr φ := by
  show pfDeg (HomPf.mk (idxZeta' (F := F) hz) φ) = _
  rw [pfDeg_mk]
  rfl

theorem baseIdentity_of_zetaMk {A : C} {z : rtObj P F A 1 ⟶ rtObj P F A 1}
    (hz : IsFrobeniusType P z) (hzb : IsBaseIdentity P z)
    (φ : rtObj P F A 1 ⟶ rtObj P F A 1)
    (h : rootBase (zetaMk (F := F) hz φ) = 𝟙 _) :
    P.Base φ = 𝟙 _ := by
  haveI := isIso_rtExt_one P F A
  haveI hbr : IsIso (P.Base (rtExt P F A 1)) := isIso_Base_of_isIso (rtExt P F A 1)
  haveI hbz : IsIso (P.Base z) := hz.2
  have hzid : P.Base z = 𝟙 _ := by
    rw [show P.Base z = P.Base (𝟙 (rtObj P F A 1)) from hzb, P.Base_id]
  have hinv : inv (P.Base z) = 𝟙 _ :=
    IsIso.inv_eq_of_hom_inv_id (by rw [hzid, Category.id_comp])
  have hrep : pfBase (HomPf.mk (idxZeta' (F := F) hz) φ) = P.Base φ := by
    rw [pfBase_mk]
    show P.Base z ≫ P.Base φ ≫ inv (P.Base z) = P.Base φ
    rw [hinv, Category.comp_id, hzid, Category.id_comp]
  have h2 : P.Base (rtExt P F A 1) ≫ pfBase (HomPf.mk (idxZeta' (F := F) hz) φ)
      ≫ inv (P.Base (rtExt P F A 1)) = 𝟙 _ := h
  rw [hrep] at h2
  have h3 := congrArg (fun t => t ≫ P.Base (rtExt P F A 1)) h2
  simp only [Category.assoc, IsIso.inv_hom_id, Category.comp_id, Category.id_comp] at h3
  have h4 : P.Base (rtExt P F A 1) ≫ P.Base φ
      = P.Base (rtExt P F A 1) ≫ 𝟙 (P.toElem.obj (rtObj P F A 1)).base := by
    rw [Category.comp_id]; exact h3
  exact (cancel_epi (P.Base (rtExt P F A 1))).mp h4

/-- ★添字の射が等しければ同じ元(束ね上げで使う)。 -/
theorem otriPfMk_congr (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ ζ' : A ⟶ A} (hζ : IsFrobeniusType P ζ) (hζ' : IsFrobeniusType P ζ')
    (h : ζ = ζ') {α α' : End A} (hα : α = α') :
    otriPfMk (F := F) hiso hζ α = otriPfMk (F := F) hiso hζ' α' := by
  subst h
  subst hα
  rfl

/-- ★★★★★**`𝒪^▷(A^pf)` の元はすべて「`𝒪^▷` の元の `m` 乗根」** ——
`Proposition 5.5, (i)` の同型の**全射性**(`𝒪^▷` の層で述べたもの)。

★`HomPf.exists_rep` ＋ `idxPf_iso_zeta`(cofinality)で代表を `ζ` 添字に移し、
`rootDeg` と `rootBase` の計算で代表が `𝒪^▷` に入ることを見る。 -/
theorem otri_pfRoot_exists_rep (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A) {f : End (⟨A, 1⟩ : PfRootObj P F)}
    (hf : f ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
    ∃ (z : rtObj P F A 1 ⟶ rtObj P F A 1) (hz : IsFrobeniusType P z),
      IsBaseIdentity P z ∧ ∃ φ : rtObj P F A 1 ⟶ rtObj P F A 1,
        (φ : End (rtObj P F A 1)) ∈ OTri P (rtObj P F A 1) ∧ zetaMk (F := F) hz φ = f := by
  obtain ⟨Z, φ₀, hZ⟩ := HomPf.exists_rep f
  obtain ⟨e, e', he, he', z, hzb, hzf, hzdeg, hae, hbe⟩ :=
    idxPf_iso_zeta F (frobTrivial_rtObj hiso hA) Z.hom.hom.1 Z.hom.property.1
      Z.hom.hom.2 Z.hom.property.2.1 Z.hom.property.2.2
  haveI := he
  haveI := he'
  have ww : Z ⟶ idxZeta' (F := F) hzf :=
    Under.homMk (⟨(e, e'), frobType_of_isIso hiso e, frobType_of_isIso hiso e',
        (degFr_of_isIso P e).trans (degFr_of_isIso P e').symm⟩ :
        Z.right ⟶ (⟨(rtObj P F A 1, rtObj P F A 1)⟩ : BiFr P F))
      (by
        refine WideSubcategory.hom_ext _ ?_
        exact Prod.ext hae hbe)
  have hmk : zetaMk (F := F) hzf (idxTransport P F ww φ₀) = f := by
    show HomPf.mk (idxZeta' (F := F) hzf) (idxTransport P F ww φ₀) = f
    rw [HomPf.mk_map ww φ₀]
    exact hZ
  refine ⟨z, hzf, hzb, idxTransport P F ww φ₀, ⟨?_, ?_⟩, hmk⟩
  · have hb : rootBase (zetaMk (F := F) hzf (idxTransport P F ww φ₀)) = 𝟙 _ := by
      rw [hmk]
      have : rootBase f = rootBase (idRoot P F (⟨A, 1⟩ : PfRootObj P F)) := hf.1
      rw [this, rootBase_id]
    exact (baseIdentity_of_zetaMk hzf hzb _ hb).trans (P.Base_id _).symm
  · have hd : rootDeg (zetaMk (F := F) hzf (idxTransport P F ww φ₀)) = 1 := by
      rw [hmk]; exact hf.2
    exact (rootDeg_zetaMk hzf _).symm.trans hd

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.5, (i)` の「immediately」の中身
(★**条つき**: 余極限の同定そのものは未実装)。 -/
def frobTransport_otri_pow.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 𝒞^pf の遷移写像は 𝒪^▷(A) の上で α ↦ α^m",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (i)` の cofinality の中身
(添字がすべて `ζ_m` の対と同型であること)。 -/
def idxPf_iso_zeta.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 𝒞^pf の添字はすべて ζ_m の対と同型",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (i)` の写像の側
(`𝒞 → 𝒞^pf` が `𝒪^▷(A) → 𝒪^▷(A^pf)` を誘導すること)。 -/
def otriToPfRoot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 𝒞 → 𝒞^pf が誘導する 𝒪^▷(A) → 𝒪^▷(A^pf)",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (i)` の逆向きの写像の材料
(`α` と次数 `m` の対から `𝒪^▷(A^pf)` の元を作る)。 -/
def otriPfMk_mem.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — α^{1/m} は 𝒪^▷(A^pf) の元",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (i)` の同型の well-defined 性
(`Pf` の関係が `𝒞^pf` の中で成り立つこと)。 -/
def otriPfMk_step.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — Pf の関係が 𝒞^pf の中で成り立つ",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (i)` の全射性の核。 -/
def homPf_exists_zeta_rep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 𝒞^pf の自己射はすべて ζ の添字で代表される",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (i)` の単射性の核。 -/
def homPf_zeta_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — ζ 添字の 2 代表が等しければ Pf の関係が成り立つ",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (i)` の全射性(`𝒪^▷` の層)。 -/
def otri_pfRoot_exists_rep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 𝒪^▷(A^pf) の元はすべて 𝒪^▷ の元の m 乗根",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
