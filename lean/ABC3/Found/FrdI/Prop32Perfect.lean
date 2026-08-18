/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Dict

/-!
# [FrdI] `Proposition 3.2, (iii)` の「perfect 型」

原文 (FrdI p.59):
> (iii) The category Cpf, equipped with the functor Cpf → FΦpf of the diagram of

原文 (FrdI p.59):
> (i), is a Frobenioid of perfect and isotropic type. Moreover, there is a natural

`Frobenioid` は `Prop32Equiv.lean` の `pfRoot_frobenioid` で、
`isotropic 型` は `Prop32Frob.lean` の `pfRoot_isOfIsotropicType` で済んでいる。
★このファイルは **perfect 型**(`Definition 1.2, (iv)`)を扱う。

## ★`IsPerfectObj` の 2 条(`MorphismTypes.lean`)

原文 (FrdI p.23):
> C such that for every n ∈N≥1, it holds that every B ∈Ob(C) base-isomorphic to

(a) 各次数 `n` について、`C` と底同型な対象 `B` は**次数 `n` の Frobenius 型射の
    終域として現れる**。
(b) 次数 `n` の Frobenius 型射 `φ₁ : B₁ ⟶ B₁'`・`φ₂ : B₂ ⟶ B₂'` と
    pre-step `ψ' : B₁' ⟶ B₂'` に対し、`φ₁ ≫ ψ' = ψ ≫ φ₂` を満たす
    **pre-step `ψ : B₁ ⟶ B₂` が一意に存在する**。

## ★★★測定(2026-08-18)—— (b) の一意性は `𝒞` では成り立たない

★`𝒞` では Frobenius 型射は **mono とは限らない**。
`χ₁ ≫ φ = χ₂ ≫ φ` からは(次数・底・零因子が一致するので)
**(vi) `faithfulUpToUnits` により単元のぶんだけのずれ**しか出ない。
★捻れ積 `𝔽_Φ ⋉ G` の模型では `u^n = 1` なる単元 `u` について
`φ ∘ u = φ` が実際に起き、mono は破れる。

★★★**ところが `𝒞^pf` では破れない。** 添字の遷移が Frobenius 型射なので、
**単元は遷移で `u ↦ u^d` と送られ、捻れのぶんが潰れる**。
これを補題の形にしたのが

  ★`homPf_cancel_frobType_idx` —— **Frobenius 型射を添字の遷移に吸収する**。

`ψ` を第 2 脚に持つ遷移 `u` を作ると `φ ≫ ψ = a ≫ idxTransport u φ` であり、
`a` は(pre-Frobenioid の全射性 `totEpiC` で)epi。よって
`φ ≫ ψ = φ' ≫ ψ` から遷移先で `idxTransport u φ = idxTransport u φ'` になる。

★これが `pfRoot_frobTypeMono`(**`𝒞^pf` の Frobenius 型射は mono**)を与え、
(b) の**一意性**が落ちる。

## ★(a) は「根を `n` 倍する」

  `(A, r·n) --Λ_{r·n}(rtExt A n)--> (A^{(n)}, r·n) ≅ (A, r)`

★`𝒞` の `frobDegSurj` は**始域**が自由だが、`𝒞^pf` では**終域**も自由になる
(`pfRoot_frobDegSurj_cod`)。これが perfection の第 1 の効き目である。

## ★★残り —— (b) の**存在**

★測定: `Hom_{𝒞^pf}((A,r),(B,s)) ≅ HomPf P F (rtObj A (e·s)) (rtObj B (e·r))` が
**任意の `e`** で成り立つ(`rtRootIso`)。したがって

  `Hom((A₁,r₁),(A₂,r₂))`  (`e := n` で取る)
  `Hom((A₁^{(n)},r₁),(A₂^{(n)},r₂))`  (`e := 1` で取り、`rtObj (rtObj A n) m ≅ rtObj A (n·m)`)

は**同じ `HomPf P F (rtObj A₁ (n·r₂)) (rtObj A₂ (n·r₁))` に写る**。
★★この同一視が (b) の存在を与えるはずである(未実装)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

section

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★1. 遷移が保つもの -/

/-- ★遷移は Frobenius 型射を保つ。 -/
theorem idxTransport_isFrobeniusType {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (h : IsFrobeniusType P φ) :
    IsFrobeniusType P (idxTransport P F u φ) :=
  prop_1_10_i_frobType_of P F u.right.property.1 u.right.property.2.1
    (idxTransport_spec (F := F) u φ) h

/-- ★遷移は等長性を保つ。 -/
theorem idxTransport_isIsometric {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (h : IsIsometric P φ) :
    IsIsometric P (idxTransport P F u φ) :=
  prop_1_10_i_isometric_of P u.right.property.1.1.2 u.right.property.2.1.1.2
    (idxTransport_spec (F := F) u φ) h

/-- ★遷移は底同型性を保つ。 -/
theorem idxTransport_isBaseIso {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (h : IsBaseIsomorphism P φ) :
    IsBaseIsomorphism P (idxTransport P F u φ) :=
  prop_1_10_i_baseIso_of P u.right.property.1.2 u.right.property.2.1.2
    (idxTransport_spec (F := F) u φ) h

/-! ## ★2. ★★★Frobenius 型射を添字の遷移に吸収する -/

/-- ★★★**`Hom^pf` では Frobenius 型射は右から消せる**。

★★手は「**Frobenius 型射を添字の遷移に吸収する**」こと ——
`ψ` を第 2 脚に持つ遷移 `u` を作ると `idxTransport u φ` は
`φ ≫ ψ = a ≫ idxTransport u φ` を満たす。`a` は(pre-Frobenioid の全射性で)epi
なので、`φ ≫ ψ = φ' ≫ ψ` から遷移先で 2 本が一致する。

★`𝒞` の中では Frobenius 型射は mono とは限らない(単元のぶんだけずれる)。
**`Hom^pf` の中では、そのずれが遷移で潰れる**——これが perfection の効き目である。 -/
theorem homPf_cancel_frobType_idx {A B : C} (Z : IdxPf P F A B)
    (φ φ' : Z.right.obj.1 ⟶ Z.right.obj.2) {E : C} (ψ : Z.right.obj.2 ⟶ E)
    (hψ : IsFrobeniusType P ψ) (h : φ ≫ ψ = φ' ≫ ψ) :
    HomPf.mk Z φ = HomPf.mk Z φ' := by
  obtain ⟨hz1, hz2, hzd⟩ := Z.hom.property
  obtain ⟨W₁, a, ha, had⟩ := F.frobDegSurj Z.right.obj.1 (P.degFr ψ)
  have hdeg : P.degFr (Z.hom.hom.1 ≫ a) = P.degFr (Z.hom.hom.2 ≫ ψ) := by
    rw [P.degFr_comp, P.degFr_comp, had, hzd]
  set u : Z ⟶ idxMk (P := P) (F := F) (Z.hom.hom.1 ≫ a) (Z.hom.hom.2 ≫ ψ)
      (IsFrobeniusType.comp P F hz1 ha) (IsFrobeniusType.comp P F hz2 hψ) hdeg :=
    Under.homMk (show Z.right ⟶ (⟨(W₁, E)⟩ : BiFr P F) from ⟨(a, ψ), ha, hψ, had⟩)
      (WideSubcategory.hom_ext _ rfl) with hu
  have h1 : φ ≫ ψ = a ≫ idxTransport P F u φ := idxTransport_spec (F := F) u φ
  have h2 : φ' ≫ ψ = a ≫ idxTransport P F u φ' := idxTransport_spec (F := F) u φ'
  haveI : Epi a := P.totEpiC _ _ _
  have hkey : idxTransport P F u φ = idxTransport P F u φ' :=
    (cancel_epi a).mp (h1.symm.trans (h.trans h2))
  rw [← HomPf.mk_map u φ, ← HomPf.mk_map u φ', hkey]

set_option maxHeartbeats 1000000 in
/-- ★★★**`Hom^pf` の中で Frobenius 型射は消せる**(3 脚の添字版)。

★`homPf_cancel_preStep` と同じ骨組みだが、最後の一手が違う ——
`𝒞` では Frobenius 型射は mono とは限らないので、
代わりに**添字を isotropic まで押し上げてから**
`homPf_cancel_frobType_idx` を当てる。 -/
theorem homPf_cancel_frobType (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (V : IdxPf3 P F A B E)
    (φ φ' : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2)
    (hψi : IsIsometric P ψ) (hψb : IsBaseIsomorphism P ψ)
    (h : HomPf.mk ((idx13 P F A B E).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F A B E).obj V) (φ' ≫ ψ)) :
    HomPf.mk ((idx12 P F A B E).obj V) φ = HomPf.mk ((idx12 P F A B E).obj V) φ' := by
  obtain ⟨V', t, t', ht⟩ := HomPf.eq_iff.mp h
  rw [idx_hom_ext t' t] at ht
  obtain ⟨V'', ⟨k⟩⟩ := exists_hom_of_final (idx13 P F A B E) V'
  set s : V ⟶ IsFiltered.max V V'' := IsFiltered.leftToMax V V'' with hs
  set r : V'' ⟶ IsFiltered.max V V'' := IsFiltered.rightToMax V V'' with hr
  have hm : t ≫ k ≫ (idx13 P F A B E).map r = (idx13 P F A B E).map s :=
    idx_hom_ext _ _
  have hA : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have hB : idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ' ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have key : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ) := by
    rw [hA, hB, ht]
  have hcan : idxTransport P F ((idx12 P F A B E).map s) φ
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ
      = idxTransport P F ((idx12 P F A B E).map s) φ'
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ := by
    rw [idxTransport_comp_pair (F := F) s φ ψ, idxTransport_comp_pair (F := F) s φ' ψ]
    exact key
  -- ★添字を「第 2 脚が isotropic」まで押し上げる
  obtain ⟨M, w, hiso⟩ := exists_idx3_isotropic2 (F := F) hfi (IsFiltered.max V V'')
  have p1 := idxTransport_comp_pair (F := F) w
    (idxTransport P F ((idx12 P F A B E).map s) φ)
    (idxTransport P F ((idx23 P F A B E).map s) ψ)
  have p2 := idxTransport_comp_pair (F := F) w
    (idxTransport P F ((idx12 P F A B E).map s) φ')
    (idxTransport P F ((idx23 P F A B E).map s) ψ)
  have hcanw := p1.trans ((congrArg (idxTransport P F ((idx13 P F A B E).map w)) hcan).trans
    p2.symm)
  have hfrob : IsFrobeniusType P
      (idxTransport P F ((idx23 P F A B E).map w)
        (idxTransport P F ((idx23 P F A B E).map s) ψ)) :=
    (isFrobeniusType_iff_of_isotropic F hiso _).mpr
      ⟨idxTransport_isIsometric _ _ (idxTransport_isIsometric _ _ hψi),
        idxTransport_isBaseIso _ _ (idxTransport_isBaseIso _ _ hψb)⟩
  have hend := homPf_cancel_frobType_idx ((idx12 P F A B E).obj M) _ _ _ hfrob hcanw
  rw [← HomPf.mk_map ((idx12 P F A B E).map s) φ, ← HomPf.mk_map ((idx12 P F A B E).map s) φ',
    ← HomPf.mk_map ((idx12 P F A B E).map w) (idxTransport P F ((idx12 P F A B E).map s) φ),
    ← HomPf.mk_map ((idx12 P F A B E).map w) (idxTransport P F ((idx12 P F A B E).map s) φ')]
  exact hend

set_option maxHeartbeats 1000000 in
/-- ★★★★**`𝒞^pf` の Frobenius 型射は monomorphism**。

★★これは `𝒞` では**成り立たない**——`𝒞` では単元のぶんだけずれる。
`𝒞^pf` では、そのずれが添字の遷移で潰れる(`homPf_cancel_frobType`)。
★これが `Definition 1.2, (iv)` の (b) の**一意性**を与える。 -/
theorem pfRoot_frobTypeMono (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) (hf : IsFrobeniusType (pfRootPre P F) f) : Mono f := by
  refine ⟨fun {W} u v huv => ?_⟩
  obtain ⟨V, φ, φ', ψ, hurep, hvrep, hfrep, hu, hv⟩ := compRoot_rep_pairL (F := F) u v f
  have hψ : IsIsometric P ψ ∧ IsBaseIsomorphism P ψ := by
    refine (isFrobeniusType_mk_iff (X := X) (Y := Y) hfi
      ((pushIdx (F := F)
        (rtLift P F X.obj (show Y.root * W.root = W.root * Y.root from mul_comm _ _))
        (rtLift_frobType P F X.obj _)
        (rtLift P F Y.obj (show X.root * W.root = W.root * X.root from mul_comm _ _))
        (rtLift_frobType P F Y.obj _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx23 P F _ _ _).obj V)) ψ).mp ?_
    rw [← rtRootIso_hom_mk (F := F) X.obj Y.obj _ _ ((idx23 P F _ _ _).obj V) ψ, ← hfrep]
    exact hf
  have h0 : (rtRootIso P F W.obj Y.obj
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
        (show X.root * W.root = X.root * W.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ))
      = (rtRootIso P F W.obj Y.obj
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
        (show X.root * W.root = X.root * W.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ)) := by
    rw [← hu, ← hv]; exact huv
  have h1 : HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ) :=
    ((Iso.hom_inv_id_apply _ _).symm.trans (congrArg _ h0)).trans (Iso.hom_inv_id_apply _ _)
  rw [hurep, hvrep, homPf_cancel_frobType (F := F) hfi V φ φ' ψ hψ.1 hψ.2 h1]

/-! ## ★3. `Definition 1.2, (iv)` の (a) —— 終域としての存在 -/

/-- ★★★**`𝒞^pf` では、どの対象も各次数の Frobenius 型射の「終域」として現れる**。

★★これが perfection の第 1 の効き目 —— `𝒞` では `frobDegSurj`(始域が自由)しか
無いが、`𝒞^pf` では**根を `n` 倍した対象から**同じ次数の Frobenius 型射が来る:

  `(A, r·n) --Λ_{r·n}(rtExt A n)--> (A^{(n)}, r·n) ≅ (A, r)`. -/
theorem pfRoot_frobDegSurj_cod (hfi : IsOfFrobeniusIsotropicType P) (B : PfRootObj P F)
    (n : ℕ+) :
    ∃ (B₀ : PfRootObj P F) (φ : B₀ ⟶ B),
      IsFrobeniusType (pfRootPre P F) φ ∧ (pfRootPre P F).degFr φ = n := by
  obtain ⟨e, he⟩ := pfRoot_exists_iso_root (F := F) B.obj B.root n (B.root * n) rfl
  haveI := he
  refine ⟨⟨B.obj, B.root * n⟩,
    lamHom (F := F) (B.root * n) (rtExt P F B.obj n) ≫ inv e, ?_, ?_⟩
  · exact IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi)
      (lamHom_isFrobeniusType (F := F) hfi (B.root * n) (rtExt P F B.obj n)
        (rtExt_frobType P F B.obj n))
      (isFrobeniusType_of_isIso (pfRootPre P F) (inv e))
  · rw [(pfRootPre P F).degFr_comp, lamHom_degFr, rtExt_degFr,
      show (pfRootPre P F).degFr (inv e) = 1 from isLinear_of_isIso (pfRootPre P F) (inv e),
      one_mul]

/-- ★★**(b) の一意性** —— pre-step は Frobenius 型射との合成で決まる。 -/
theorem pfRoot_preStep_uniq_of_comp_frobType (hfi : IsOfFrobeniusIsotropicType P)
    {B₁ B₂ B₂' : PfRootObj P F} (φ₂ : B₂ ⟶ B₂')
    (hφ₂ : IsFrobeniusType (pfRootPre P F) φ₂) (ψ ψ' : B₁ ⟶ B₂)
    (h : ψ ≫ φ₂ = ψ' ≫ φ₂) : ψ = ψ' :=
  haveI := pfRoot_frobTypeMono (F := F) hfi φ₂ hφ₂
  (cancel_mono φ₂).mp h

end

/-! ## ★出典の紐付け(条つき) -/

def pfRoot_frobDegSurj_cod.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59,
    item := "Proposition 3.2, (iii) — perfect 型の (a)",
    sectionId := "frdi-prop-3-2" }

def pfRoot_frobTypeMono.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59,
    item := "Proposition 3.2, (iii) — perfect 型の (b) の一意性",
    sectionId := "frdi-prop-3-2" }

end ABC3.Found.FrdI
