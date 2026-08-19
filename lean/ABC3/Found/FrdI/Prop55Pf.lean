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

end ABC3.Found.FrdI
